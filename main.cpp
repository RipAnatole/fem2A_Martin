#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quadrature = false;
    const bool t_ElementMapping = false;
    const bool t_test_assemble_elementary_matrix = false;
    const bool t_test_assemble_vector_matrix = true;
    const bool t_solver = false;

    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if( t_quadrature ){
        Tests::test_quad(0, false);
        Tests::test_quad(2, false);
        Tests::test_quad(4, false);
        Tests::test_quad(6, false);
        }
    if( t_ElementMapping ) Tests::test_ElementMapping();
    if ( t_test_assemble_elementary_matrix ) Tests::test_assemble_elementary_matrix();
    if ( t_test_assemble_vector_matrix ) Tests::test_assemble_local_and_global_vector();
    if( t_solver ) Tests::test_solver();

}

void run_simu()
{

    const bool simu_pure_dirichlet = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) {
        Simu::pure_dirichlet_pb("data/square.mesh", verbose);
    }
}

int main( int argc, const char * argv[] )
/* 	argc = nombre d'arguments
		ex : ./build/fem2a => argc = 0
	argv = liste des arguments
		argc[0] = le 1er argument
 		argc[1] = le 2e argument
*/ 
{
    /* Command line parsing */
    /* Il boucle sur le nombre d'arguments. 
    Il remplit arguments avec les arguments donnés */
    
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    /* S'il n'y a pas d'arguments ou (||) on a utilisé l'argument -h,
    il effectue runtest*/
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked, cad -t */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked, cad -s */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
