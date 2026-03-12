#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_env.hpp"
#include "smesh_mesh.hpp"
#include "smesh_output.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "ssdf.hpp"

#ifdef SMESH_ENABLE_MPI
#include "distributed/psdf.hpp"
#endif

// #error "Not implemented"

int main(int argc, char **argv) {
    auto ctx = smesh::initialize(argc, argv);
    SMESH_TRACE_SCOPE("mesh_distance_to_boundary.exe");

    {
        auto comm = ctx->communicator();

        if (argc != 3) {
            if (!comm->rank()) {
                std::cerr << "Usage: " << argv[0] << " <mesh> <output_folder>" << std::endl;
            }
            return 1;
        }

        auto mesh = smesh::Mesh::create_from_file(comm, smesh::Path(argv[1]));
        auto surface = smesh::skin(mesh);
        auto out = smesh::create_host_buffer<smesh::geom_t>(mesh->n_nodes());

        {
            auto b_out = out->data();
            for (int i = 0; i < mesh->n_nodes(); i++) {
                b_out[i] = 1e8f;
            }
        }

        auto p = mesh->points()->data();
        auto sp = surface->points()->data();
        auto se = surface->elements(0)->data();

#ifdef SMESH_ENABLE_MPI

        bool update_enabled = smesh::Env::read("SSDF_ENABLE_UPDATE", false);

        if (update_enabled) {
            auto dist = mesh->distributed();
            ssdf::edf_select(dist->n_nodes_owned(),
                             p[0],
                             p[1],
                             p[2],
                             surface->n_elements(),
                             se[0],
                             se[1],
                             se[2],
                             surface->n_nodes(),
                             sp[0],
                             sp[1],
                             sp[2],
                             out->data());

            ssdf::pedf_update(comm->get(),
                              dist->n_nodes_owned(),
                              p[0],
                              p[1],
                              p[2],
                              surface->n_elements(),
                              se[0],
                              se[1],
                              se[2],
                              surface->n_nodes(),
                              sp[0],
                              sp[1],
                              sp[2],
                              out->data());
        } else {
            if (comm->size() > 1) {
                auto dist = mesh->distributed();
                ssdf::pedf(comm->get(),
                           dist->n_nodes_owned(),
                           p[0],
                           p[1],
                           p[2],
                           surface->n_elements(),
                           se[0],
                           se[1],
                           se[2],
                           surface->n_nodes(),
                           sp[0],
                           sp[1],
                           sp[2],
                           out->data());
            }
        }
#else
        {
            ssdf::edf_select(mesh->n_nodes(),
                             p[0],
                             p[1],
                             p[2],
                             surface->n_elements(),
                             se[0],
                             se[1],
                             se[2],
                             surface->n_nodes(),
                             sp[0],
                             sp[1],
                             sp[2],
                             out->data());
        }
#endif

        auto output = smesh::Output::create(mesh, smesh::Path(argv[2]));
        output->write_nodal("edf", out);
    }

    return 0;
}
