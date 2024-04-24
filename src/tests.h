#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"
#include "simu.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        void affiche_vector(std::vector< double > int_vec) {
            for_each(int_vec.begin(), int_vec.end(),
                [](const int& n) { std::cout << n << "; "; });
            std::cout << std::endl;
        }


        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quad( int ordre, bool bord )
        {
            Quadrature my_quad; // il faut instancier avant d'appeler get_quadrature() 
            my_quad = Quadrature::get_quadrature( ordre, bord );
	    
	    std::cout << "nombre de points : " << my_quad.nb_points() << std::endl;
	    
	    double sum = 0; //les doubles sont bcp plus précis que les float
	    
	    for ( int i = 0; i < my_quad.nb_points(); ++i )
	    {
	        /*
	        std::cout << my_quad.point(i).x << " " << my_quad.point(i).y <<
	            "\npoids : " << my_quad.weight(i) << std::endl;
	        */ 
	        sum += my_quad.weight(i);
	    }
	    
	    std::cout << "somme des poids : " << sum << std::endl;
	    
            return true;
        }
        
        bool test_ElementMapping()
        {
            std::cout << "test coordonnées matrice\n";
            
            Mesh mesh;
            mesh.load("data/square.mesh");
            
            ElementMapping ElMapTr4(mesh, false, 4);
            ElementMapping ElMapBd4(mesh, true, 4);
            
            std::cout << "\ntest transformation, soit du mapping du point (xi=0.2 , eta=0.4)\n";
            vertex x_r; x_r.x = 0.2; x_r.y = 0.4;
            vertex r = ElMapTr4.transform(x_r);
            std::cout << r.x << " " << r.y << std::endl;
            
            std::cout << "\ntest Matrice Jacobienne du point (xi=0.2 , eta=0.4)\n";
            DenseMatrix J = ElMapTr4.jacobian_matrix(x_r);
            J.print();
            std::cout << "determinant :" << J.det_2x2() << std::endl;
            return true;
        }

        bool test_ShapeFunction()
        {
            std::cout << "test ShapeFunction\n";
            ShapeFunctions sf(1, 1);
            std::cout << "c'est teste\n";
            return true;
        }
        
        bool test_assemble_elementary_matrix()
        {
            std::cout << "\ntest assemblage elementaire Ke\n";
            Mesh mesh;
            mesh.load("data/square.mesh");
            ElementMapping my_el_map(mesh, false, 4);
            ShapeFunctions my_shape(2, 1);
        
            Quadrature quadrature0;
            quadrature0 = quadrature0.get_quadrature(0); //normalement, dans notre cas, on a le meme resultat pour une quadrature d'ordre 0 ou 2 car nos fonctions sont d'ordre 0
            std::cout << "test avec ordre de quadrature 0\n";
            DenseMatrix Ke;
            
            assemble_elementary_matrix( my_el_map, my_shape, quadrature0, FEM2A::Simu::unit_fct, Ke );
            Ke.print();
            
            Quadrature quadrature2;
            quadrature2 = quadrature2.get_quadrature(2);
            assemble_elementary_matrix( my_el_map, my_shape, quadrature2, FEM2A::Simu::unit_fct, Ke );
            
            std::cout << "test avec ordre de quadrature 2\n";
            Ke.print();
            
            //test assemblage global K
            std::cout << "\ntest assemblage global K\n";
            SparseMatrix K(mesh.nb_vertices()); //mettre autre valeur plus sure
            local_to_global_matrix(mesh, 4, Ke, K );
            K.print();
            return true;
        }
        
        bool test_apply_dirichlet_boundary_conditions() {
            Mesh mesh;
            mesh.load("data/square.mesh");
            /*
            const std::vector< bool > attributes_Dirichlet; // size: nb of attributes 
            const std::vector< double > Dirichlet_values; // size: nb of DOFs
            SparseMatrix K(mesh.nb_vertices());
            std::vector< double >& F(mesh.nb_vertices(), 0.);
            
            apply_dirichlet_boundary_conditions(mesh, attributes_Dirichlet, Dirichlet_values, K, F );
            */
            return true;
        }
            
        
        
        bool test_assemble_local_and_global_vector()
        {
            std::cout << "\ntest assemblage elementaire Fe\n";
            
            Mesh mesh;
            mesh.load("data/square.mesh");
            
            ElementMapping my_el_map(mesh, false, 4);

            ShapeFunctions my_shape(2, 1);
            Quadrature quadrature;
            quadrature = quadrature.get_quadrature(0);
            std::vector< double > Fe;
            
            assemble_elementary_vector(my_el_map, my_shape, quadrature, FEM2A::Simu::unit_fct, Fe);
            affiche_vector(Fe);
            
            //pour un bord
            std::cout << "pour un bord\n";
            ElementMapping my_el_map2(mesh, true, 4);
            assemble_elementary_vector(my_el_map2, my_shape, quadrature, FEM2A::Simu::unit_fct, Fe);
            affiche_vector(Fe);
              
            return true;      
        }
        
        bool test_solver()
        {
            std::cout << "le test de la solution est vide pour l'instant\n";
            return true;
        }

        
    }
}
