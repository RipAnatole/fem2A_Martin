#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose = false)
        {
            std::cout << "Simulation du probleme de Dirichlet pur\n";
            
            // Condition de Dirichlet u = x + y
            // Coefficient de diffusion k constant, k = 1
            
            // Initialisation des attributs
            Mesh mesh;
            mesh.load("data/square.mesh");
            SparseMatrix k_global(mesh.nb_vertices());
            //k_global.set_size(mesh.nb_vertices(), mesh.nb_vertices());
            std::vector <double> F_global(mesh.nb_vertices(), 0.);
            
            //On complete k_global, donc on a besoin que des triangles
            for (int i = 0; i < mesh.nb_triangles(); ++i) {
                ElementMapping my_el_map(mesh, false, i);
                ShapeFunctions my_shape(2, 1); 
                Quadrature quadrature;
                quadrature = quadrature.get_quadrature(2);
                DenseMatrix ke;
                
                assemble_elementary_matrix( my_el_map, my_shape, quadrature, FEM2A::Simu::unit_fct, ke ); //unit_fct car k est constant, k = 1, ce que renvoie unit_fct
                
                //On complete k_global
                local_to_global_matrix( mesh, i, ke, k_global );
            }
            
            //On met les conditions de dirichlet dans notre maillage
            std::vector< bool > att_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;
            mesh.set_attribute(unit_fct, 1, true);
            
            std::vector< double > imposed_values(mesh.nb_vertices());
            
            for (int i = 0; i < mesh.nb_vertices(); ++i) {
                imposed_values[i] = mesh.get_vertex(i).x + mesh.get_vertex(i).y; 
            } 
            
            apply_dirichlet_boundary_conditions(mesh, att_is_dirichlet, imposed_values, k_global, F_global);  
            
            //On lance la simulation
            std::vector< double > u(mesh.nb_vertices());
            solve(k_global, F_global, u);
            std::string export_name = "pure_dirichlet";
            mesh.save(export_name + ".mesh");
            save_solution(u, export_name + ".bb");
            
            
        }

    }

}
