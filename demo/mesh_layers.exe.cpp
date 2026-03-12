#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_exchange.hpp"
#include "smesh_mesh.hpp"
#include "smesh_output.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "graph/layers.hpp"

#ifdef SMESH_ENABLE_MPI
#include "distributed/psdf.hpp"
#endif

// #error "Not implemented"

int main(int argc, char **argv) {
    auto ctx = smesh::initialize(argc, argv);
    SMESH_TRACE_SCOPE("mesh_layers.exe");

    {
        auto comm = ctx->communicator();

        if (argc != 4) {
            if (!comm->rank()) {
                std::cerr << "Usage: " << argv[0] << " <mesh> <surface> <output_folder>" << std::endl;
            }
            return 1;
        }

        auto mesh = smesh::Mesh::create_from_file(comm, smesh::Path(argv[1]));

        // Everyone reads the whole surface
        auto surface = smesh::Mesh::create_from_file(smesh::Communicator::self(), smesh::Path(argv[2]));

        const ptrdiff_t n_nodes = mesh->n_nodes();
        auto layer = smesh::create_host_buffer<uint8_t>(n_nodes);
        auto b_layer = layer->data();

        const uint8_t SSDF_MAX_LAYERS = smesh::Env::read("SSDF_MAX_LAYERS", 256);

        auto p = mesh->points()->data();
        auto selements = surface->elements(0)->data();
        auto sp = surface->points()->data();

        ssdf::points_inside_bvh(n_nodes,
                                p[0],
                                p[1],
                                p[2],
                                surface->n_elements(),
                                selements[0],
                                selements[1],
                                selements[2],
                                surface->n_nodes(),
                                sp[0],
                                sp[1],
                                sp[2],
                                b_layer);

#ifdef SMESH_ENABLE_MPI
        if (comm->size() > 1) {
            auto dist = mesh->distributed();
            auto exchange = smesh::Exchange::create_nodal(mesh, smesh::Exchange::ExchangeScope::GhostsAndAura);

            ptrdiff_t *n2e_ptr{nullptr};
            smesh::idx_t *n2e_idx{nullptr};

            ssdf::create_n2e(mesh->n_elements(0),
                             mesh->n_nodes(),
                             mesh->n_nodes_per_element(0),
                             mesh->elements(0)->data(),
                             &n2e_ptr,
                             &n2e_idx);

            auto queue = smesh::create_host_buffer<smesh::idx_t>(n_nodes);
            auto b_queue = queue->data();

            int queue_size = 0;
            for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                if (b_layer[i] == 1) {
                    b_layer[i] = 0;
                    b_queue[queue_size++] = static_cast<smesh::idx_t>(i);
                } else {
                    // Make sure that the point is not assigned to a layer
                    b_layer[i] = 255;
                }
            }

            const ptrdiff_t n_changing = n_nodes - dist->n_nodes_owned();
            auto changing = smesh::create_host_buffer<uint8_t>(n_changing);
            auto b_changing = changing->data();
            memcpy(b_changing, &b_layer[dist->n_nodes_owned()], sizeof(*b_layer) * n_changing);

            for (int i = 0; i < comm->size() * 1000; i++) {
                int cursor = 0;
                ssdf::layers_iterative(mesh->n_nodes_per_element(0),
                                       mesh->elements(0)->data(),
                                       b_layer,
                                       SSDF_MAX_LAYERS,
                                       //  Auxiliary arrays
                                       n2e_ptr,
                                       n2e_idx,
                                       &cursor,
                                       &queue_size,
                                       queue->size(),
                                       queue->data());

                exchange->gather(b_layer);

                const ptrdiff_t start = dist->n_nodes_owned();
                const ptrdiff_t end = n_nodes;
                ptrdiff_t n_changes = 0;
                for (ptrdiff_t i = start; i < end; i++) {
                    n_changes += b_layer[i] != b_changing[i - start];
                    b_changing[i - start] = b_layer[i];
                }

                n_changes = comm->sum(n_changes);
                if (!n_changes) {
                    if (!comm->rank()) {
                        printf("Finished at iteration %i\n", i);
                    }
                    break;
                }

                queue_size = 0;

                for (ptrdiff_t i = start; i < end; i++) {
                    b_queue[queue_size++] = i;
                }
            }

            free(n2e_ptr);
            free(n2e_idx);

            // Repeat layers and see if something changed, alltoall to synch
        } else
#endif
        {

            for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                b_layer[i] = b_layer[i] != 1;
            }

            ssdf::layers(mesh->n_elements(0),
                         mesh->n_nodes_per_element(0),
                         mesh->elements(0)->data(),
                         mesh->n_nodes(),
                         b_layer,
                         SSDF_MAX_LAYERS);
        }

        auto output = smesh::Output::create(mesh, smesh::Path(argv[3]));
        output->write_nodal("layer", smesh::astype<float>(layer));

        {
            auto sdf = smesh::create_host_buffer<smesh::geom_t>(n_nodes);
            auto b_sdf = sdf->data();
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < n_nodes; i++) {
                b_sdf[i] = std::numeric_limits<smesh::geom_t>::max();
            }

            {
                auto boundary = smesh::skin(mesh);
                auto be = boundary->elements(0)->data();
                auto bp = boundary->points()->data();

#ifdef SMESH_ENABLE_MPI
                auto dist = mesh->distributed();
                ssdf::pedf(comm->get(),
                           dist->n_nodes_owned(),
                           p[0],
                           p[1],
                           p[2],
                           boundary->n_elements(),
                           be[0],
                           be[1],
                           be[2],
                           boundary->n_nodes(),
                           bp[0],
                           bp[1],
                           bp[2],
                           b_sdf);
#else
                ssdf::edf_select(mesh->n_nodes(),
                                 p[0],
                                 p[1],
                                 p[2],
                                 boundary->n_elements(),
                                 be[0],
                                 be[1],
                                 be[2],
                                 boundary->n_nodes(),
                                 bp[0],
                                 bp[1],
                                 bp[2],
                                 b_sdf);
#endif
            }

            ssdf::edf_select(mesh->n_nodes(),
                             p[0],
                             p[1],
                             p[2],
                             surface->n_elements(),
                             selements[0],
                             selements[1],
                             selements[2],
                             surface->n_nodes(),
                             sp[0],
                             sp[1],
                             sp[2],
                             b_sdf);

#pragma omp parallel for
            for (ptrdiff_t i = 0; i < n_nodes; i++) {
                if (b_layer[i] == 0) {
                    b_sdf[i] = -b_sdf[i];
                }
            }

            output->write_nodal("sdf", sdf);
        }
    }

    return 0;
}
