#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border ) //prend la valeur de border
    {
        /* 
        - border true : on est a un bord, donc un segment, donc 2 coordonnées, donc 2 vertices
        - border false : c'est un triangle, donc 3 coordonnées, donc 3 vertices 
        pour passer d'un indice en vecteur, il faut voir le COURS
        */
        // TODO, affecter les valeurs à vertice_, c'est un vecteur de vertex, de coordonnées x et y. si border = False : triangle, 3 vertexes = coordonnées des 3 vertices du système
        if ( border ) { //border
            for (int n = 0; n < 2; ++n) vertices_.push_back(M.get_edge_vertex(i,n));
            for (int n = 0; n < 2; ++n) std::cout << "x : " << vertices_[n].x << ", y : " << vertices_[n].y << std::endl; 
        }
        else { //triangle
            for (int n = 0; n < 3; ++n) vertices_.push_back(M.get_triangle_vertex(i, n));
            for (int n = 0; n < 3; ++n) std::cout << "x : " << vertices_[n].x << ", y :" << vertices_[n].y << std::endl;
        }
        
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        // TODO
        vertex r ;
        r.x = 0; r.y = 0;
        double phiZero; double phiUn; double phiDeux;
        if (border_) {
            phiZero = 1 - x_r.x; 
            phiUn = x_r.x;
            double phi[2] = {phiZero, phiUn};
            for (int i = 0; i<2; ++i) {
                r.x += phi[i]*vertices_[i].x;
                r.y += phi[i]*vertices_[i].y;
            } 
        } else {
            phiZero = 1 - x_r.x - x_r.y; 
            phiUn = x_r.x;
            phiDeux = x_r.y;
            double phi[3] = {phiZero, phiUn, phiDeux};
            for (int i = 0; i<3; ++i) {
                r.x += phi[i]*vertices_[i].x;
                r.y += phi[i]*vertices_[i].y;
            } 
        }
        return r ;
    }

    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        // TODO
        DenseMatrix J ;
        if ( border_ ) {
            J.set_size( 2, 1 );
            double a = -vertices_[0].x + vertices_[1].x;
            double b = -vertices_[0].y + vertices_[1].y;
            J.set(0, 0, a);
            J.set(1, 0, b);
        } else {
            J.set_size( 2, 2 );
            double a = -vertices_[0].x + vertices_[1].x;
            double b = -vertices_[0].x + vertices_[2].x;
            double c = -vertices_[0].y + vertices_[1].y;
            double d = -vertices_[0].y + vertices_[2].y;
            J.set(0, 0, a);
            J.set(0, 1, b);
            J.set(1, 0, c);
            J.set(1, 1, d);
        }
        return J ;
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        // TODO
        DenseMatrix J = jacobian_matrix(x_r);
        
        if ( border_ ) {
            double JTJ = J.get(0, 0)*J.get(0, 0) + J.get(1, 0)*J.get(1, 0);
            return std::sqrt(JTJ);
        } else {
            assert ( J.det_2x2() != 0 );
            return J.det_2x2();
            
        }
        
        return 0. ;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        /**
        * \brief Constructor of the ShapeFunctions
        * \param dim 1 for reference segment, 2 for reference triangle
        * \param order Should be 1 (linear functions only)
        */
        
        // TODO
        assert (order_ <= 1) ;
    }

    int ShapeFunctions::nb_functions() const
    {
        // TODO
        
        if ( dim_ == 2 ) return 3 ;
        return 2;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        // TODO
        
        if ( dim_ == 1 ) {
            double xi = x_r.x;
            switch(i) {
                case(0) : return 1-xi;
                case(1) : return xi;
            }
        } else {
            double xi = x_r.x; double eta = x_r.y;
            switch(i) {
                case(0) : return 1-xi-eta;
                case(1) : return xi;
                case(2) : return eta;
            }
        }
        
        /*
        vertex r ;
        r.x = 0; r.y = 0;
        double phiZero; double phiUn; double phiDeux;
        if (dim_ == 1) {
            phiZero = 1 - x_r.x; 
            phiUn = x_r.x;
            double phi[2] = {phiZero, phiUn};
            for (int i = 0; i<2; ++i) {
                r.x += phi[i]*vertices_[i].x;
                r.y += phi[i]*vertices_[i].y;
            } 
        } else {
            phiZero = 1 - x_r.x - x_r.y; 
            phiUn = x_r.x;
            phiDeux = x_r.y;
            double phi[3] = {phiZero, phiUn, phiDeux};
            for (int i = 0; i<3; ++i) {
                r.x += phi[i]*vertices_[i].x;
                r.y += phi[i]*vertices_[i].y;
            } 
        }
        */
        return 0. ; // should not be reached
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        // TODO
        vec2 g ;
        
        if ( dim_ == 1 ) { //border
            switch(i) {
                case(0) : g.x = -1.; break;
                case(1) : g.x = 1.; break; 
            }
            g.y = 0.;
        } else { //triangle
            switch(i) {
                case(0) : g.x = -1.; g.y = -1.; break;
                case(1) : g.x = 1.; g.y = 0.; break;
                case(2) : g.x = 0.; g.y = 1.; break;
            }
        }
        
        /*
        // border
        if(dim_==1){
            if(i ==1){
                g.x = -1;
                g.y = 0;
                }
            else{
                g.x = 1;
                g.y = 0;
                }
        }// Pour un triangle
        else {
            switch (i){
            case 0:
                g.x = -1;
                g.y = -1;
                break;
            case 1:
                g.x = 1;
                g.y = 0;
                break;
            case 2:
                g.x = 0;
                g.y = 1;
            }
        }
        */
        /*
        assert ( J.det_2x2() != 0 );
    
        double det = J.det_2x2()
        
        DenseMatrix Jt = J.invert_2x2()
        
    
        double a = J.get(0, 0)/det;
        double b = J.get(0, 1)/det;
        double c = J.get(1, 0)/det;
        double d = J.get(1, 1)/det;
    
        J.set(0, 0, d);
        J.set(0, 1, -b);
        J.set(1, 0, -c);
        J.set(1, 1, a);
    
        DenseMatrix Jf = invert_2x2()
        */
        
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),
        DenseMatrix& Ke )
    {
        // TODO
        //jacobian_matrix( vertex x_r )
    
        Ke.set_size(reference_functions.nb_functions(), reference_functions.nb_functions());
    
        for (int i = 0; i < reference_functions.nb_functions() ; ++i) { 
            for (int j = 0; j < reference_functions.nb_functions(); ++j) {
                
            //Ke ne se fait que dans les triangles et pas sur les bords, pas besoin de faire un cas bords_
                double valeur = 0;

                //On récupère la matrice jacobienne
                
                for (int q = 0; q < quadrature.nb_points(); ++q) {
                
                    double wq = quadrature.weight(q);
                    
                    double ksiq = quadrature.point(q).x; 
                    double etaq = quadrature.point(q).y;
                    //vertex point = quadrature.point(q);
                    
                    //On recupere l'element Me(ksiq, etaq)
                    vertex Me = elt_mapping.transform(quadrature.point(q));
                    double kq = coefficient(Me);
                    
                    //On recupere la jacobienne
                    DenseMatrix J = elt_mapping.jacobian_matrix( quadrature.point(q) );
                    
                    //On recupere la transposee de l'inverse de la jacobienne
                    DenseMatrix JinvT = J.invert_2x2().transpose();
                    
                    //On met la jacobienne en valeurs positives
                    /*for (int height = 0; height < J.height(); ++height) {
                        for (int width = 0; width < J.width(); ++width) {
                            if (J.get(height, width) < 0) {
                                J.set(height, width, -J.get(height, width));
                            }
                        }    
                    }
                    */
                    double detJ = J.det_2x2();
                
                    //On recupere les gradiants de phi
                    vec2 phi_grad_i = reference_functions.evaluate_grad( i, quadrature.point(q) );
                    
                    vec2 phi_grad_j = reference_functions.evaluate_grad( j, quadrature.point(q) );   

                    //On complete l'expression de Ke(i, j)
                    valeur += wq * kq * dot ( JinvT.mult_2x2_2( phi_grad_i ), JinvT.mult_2x2_2( phi_grad_j ) ) * detJ;   
                    
                }
            Ke.set(i, j, valeur);             
            }
        }
    }



    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        std::cout << "Ke -> K" << '\n';
        // TODO
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (source term)" << '\n';
        // TODO
        
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        // TODO
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {
        std::cout << "Fe -> F" << '\n';
        // TODO
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        std::cout << "apply dirichlet boundary conditions" << '\n';
        // TODO
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
