#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_output.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"
#include "smesh_exchange.hpp"

#include "graph/layers.hpp"

#ifdef SMESH_ENABLE_MPI
#include "distributed/psdf.hpp"
#endif

// #error "Not implemented"

int main(int argc, char **argv) {
    auto ctx = smesh::initialize(argc, argv);
    SMESH_TRACE_SCOPE("mesh_distance_to_boundary.exe");

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
        auto out = smesh::create_host_buffer<uint8_t>(n_nodes);

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
                                out->data());

        auto b_out = out->data();
        for (ptrdiff_t i = 0; i < n_nodes; i++) {
            b_out[i] = !b_out[i];
        }

#ifdef SMESH_ENABLE_MPI
        {
            ssdf::layers(mesh->n_elements(0),
                         mesh->n_nodes_per_element(0),
                         mesh->elements(0)->data(),
                         mesh->n_nodes(),
                         b_out,
                         SSDF_MAX_LAYERS);

            auto exchange = smesh::Exchange::create_nodal(mesh, smesh::Exchange::ExchangeScope::GhostsOnly);
            exchange->gather(b_out);

            // Repeat layers and see if something changed, alltoall to synch
        }
#else

        {
            ssdf::layers(mesh->n_elements(0),
                         mesh->n_nodes_per_element(0),
                         mesh->elements(0)->data(),
                         mesh->n_nodes(),
                         out->data(),
                         SSDF_MAX_LAYERS);
        }
#endif

        auto output = smesh::Output::create(mesh, smesh::Path(argv[3]));
        output->write_nodal("layer", smesh::astype<float>(out));
    }

    return 0;
}
