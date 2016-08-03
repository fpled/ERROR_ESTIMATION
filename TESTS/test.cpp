//
// C++ Implementation: test_cpp
//
// Description: Test
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "../build/problem_error_estimation/all_in_one.h"
#include "../Mesh.h"
#include "../Material_properties.h"
#include "../Boundary_conditions.h"
#include "../Display.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"
#include "../VERIFICATION/Verification.h"
#include "../DISCRETIZATION_ERROR/Calcul_discretization_error.h"
#include "../Calcul_global_error_estimation.h"
#include "../Calcul_goal_oriented_error_estimation.h"

#include "LMT/include/containers/gnuplot.h"
#include "LMT/include/containers/matlabplot.h"
#include "LMT/include/containers/matcholamd.h"
#include "LMT/include/containers/conjugate_gradient.h"
#include "LMT/include/containers/MatWithTinyBlocks.h"

#include "LMT/include/util/MKL_direct_solver.h"
#include "LMT/include/util/MKL_iterative_solver.h"
#include "LMT/include/util/MUMPS_solver.h"

#include "LMT/include/mesh/meshcaracstd.h"

using namespace LMT;
using namespace std;

/*
int main( int argc, char **argv ) {
    TicToc t_total;
    t_total.start();
    static const unsigned dim = 3;
    static const bool wont_add_nz = true;
    typedef Mesh<Mesh_carac_error_estimation<double,dim> > TM;
    typedef Formulation<TM,FormulationElasticity,DefaultBehavior,double,wont_add_nz> TF;
    typedef TM::Pvec Pvec;
    typedef TM::TNode::T T;
    static const string structure = "hashin_32";
    static const string mesh_size = "fine";
    static const string loading = "Step-2";
    static const unsigned deg_p = 1;
    static const string boundary_condition_D = "penalty";

    static const string method = "EET";
    static const unsigned cost_function = 0;
    static const T penalty_val_N = 1e6;
    static const string solver = "LDL";
    static const string solver_minimisation = "UMFPACK";

    static const bool enhancement_with_geometric_criterium = 0;
    static const string geometric_criterium = "radius_ratio";
    static const T val_geometric_criterium = 0.34;
    static const bool enhancement_with_estimator_criterium = 0;
    static const T val_estimator_criterium = 0.8;

    static const bool verif_eq = 1;
    static const bool verif_compatibility_conditions = 1;
    static const T tol_compatibility_conditions = 1e-6;
    static const bool verif_eq_force_fluxes = 1;
    static const T tol_eq_force_fluxes = 1e-6;

    static const bool verif_solver = 1;
    static const T tol_solver = 1e-6;
    static const bool verif_solver_enhancement = 1;
    static const T tol_solver_enhancement = 1e-6;
    static const bool verif_solver_minimisation = 1;
    static const T tol_solver_minimisation = 1e-6;
    static const bool verif_solver_minimisation_enhancement = 1;
    static const T tol_solver_minimisation_enhancement = 1e-6;

    /// ------------------------------------------------------- ///
    /// Construction de la solution elements finis du pb direct ///
    /// ------------------------------------------------------- ///

    display_pb( dim, structure, deg_p  );

    /// Maillage du pb direct
    /// ---------------------
    TM m;
    set_mesh( m, structure, mesh_size, loading, deg_p );

    /// Formulation du pb direct
    /// ------------------------
    TF f( m );

    /// Proprietes materiaux du pb direct
    /// ---------------------------------
    set_material_properties( f, m, structure );

    /// Conditions aux limites du pb direct
    /// -----------------------------------
    set_constraints( f, m, boundary_condition_D, "direct", structure, loading );
    set_load_conditions( m, structure, loading, mesh_size );

    /// Resolution du pb direct
    /// -----------------------
    cout << "Resolution du pb direct" << endl << endl;
    TicToc t;
    t.start();
    f.solve();
    t.stop();
    cout << "temps de calcul de la resolution du pb direct = " << t.res << endl << endl;

//    f.allocate_matrices();
//    f.shift();
//    f.assemble();
////    f.solve_system();
//    f.get_initial_conditions();
//    f.update_variables();
//    f.call_after_solve();

    /// Verification de l'equilibre du pb direct
    /// ----------------------------------------
    if ( verif_eq )
        check_equilibrium( f, "direct" );

    /// Calcul de la norme du champ de deplacement approche du pb direct
    /// ----------------------------------------------------------------
    calcul_norm_dep( m, f, "direct" );

    /// Affichage
    /// ---------
    display( m, "test" );

//    cout << m.type_elements() << endl;
}
/*

/// ---------- ///
/// refinement ///
/// ---------- ///
/*
int main( int argc, char **argv ) {
    static const unsigned dim = 2;
    static const bool wont_add_nz = true;
//    typedef Mesh<MeshCaracStd<dim,2> > TM;
    typedef Mesh<Mesh_carac_error_estimation<double,dim> > TM;
//    typedef Formulation<TM,elasticity_iso> TF;
    typedef Formulation<TM,FormulationElasticity,DefaultBehavior,double,wont_add_nz> TF;
    typedef TM::Pvec Pvec;
    typedef TM::TNode::T T;

    TM m;
    T lx = 1., ly = 1., lz = 1.;
//    make_rect( m, Triangle(), Pvec( 0., 0. ), Pvec( lx, ly ), Pvec( 5, 5 ) );
//    make_rect( m, Quad(), Pvec( 0., 0. ), Pvec( lx, ly ), Pvec( 5, 5 ) );
//    make_rect( m, Tetra(), Pvec( 0., 0., 0. ), Pvec( lx, ly, lz ), Pvec( 5, 5, 5 ) );
//    make_rect( m, Hexa(), Pvec( 0., 0., 0. ), Pvec( lx, ly, lz ), Pvec( 5, 5, 5 ) ); replace_Hexa_by_Tetra( m );
//    make_rect( m, Triangle(), Pvec( 0., 0. ), Pvec( lx, ly ), Pvec( 21, 5 ) );
//    make_rect( m, Tetra(), Pvec( 0., 0., 0. ), Pvec( lx, ly, lz ), Pvec( 21, 5, 5 ) );
   display_mesh_carac( m );

    TF f( m );

    for(unsigned i = 0; i < m.node_list.size(); ++i ) {
        if ( ( m.node_list[ i ].pos[ 0 ] < 1e-6 ) or ( m.node_list[ i ].pos[ 0 ] > lx - 1e-6 ) ) {
            for( int d = 0; d < dim; ++d )
                f.add_constraint( "sin( node[" + to_string( i ) + "].dep[" + to_string( d ) + "] ) - " + to_string( 0.1 * m.node_list[ i ].pos[ 0 ] * ( d == 0 ) ), 1e5 );
        }
    }

    f.solve();

//    display_mesh( m );

//    PRINT( generate( m.node_list, ExtractDM<pos_DM>() ) );
//    PRINT( generate( m.node_list, ExtractDMi<pos_DM>( 0 ) ) );
//    PRINT( generate( m.node_list, ExtractDMi<pos_DM>( 1 ) ) );

    bool spread_cut = true;

    refinement_if_constraints( m, f, spread_cut );

//    while ( refinement_if_length_sup( m, 0.05, spread_cut ) );

//    while ( refinement_point( m, 0.01, 0.2, Pvec( 0.2, 0.5 ), spread_cut ) );

//    while ( refinement_circle( m, 0.01, 0.2, Pvec( 0.3, 0.5 ), 0.1, spread_cut ) );

//    for( unsigned i = 0 ; i < m.node_list.size(); ++i )
//        m.node_list[i].phi_domain = sin( std::sqrt( i ) * 5. );
//    level_set_refinement( m, ExtractDM< phi_domain_DM >(), spread_cut );

//    for( unsigned i = 0 ; i < m.node_list.size(); ++i )
//        m.node_list[i].phi_domain = sin( std::sqrt( i ) * 5. );
//    level_set_cut( m, ExtractDM< phi_domain_DM >(), spread_cut );

//    for( unsigned i = 0 ; i < m.node_list.size(); ++i )
//        m.node_list[i].theta_nodal = sin( std::sqrt( i ) * 5. );
//    refinement_if_nodal_field_sup( m, ExtractDM< theta_nodal_DM >(), 0.75, spread_cut );

//    for( unsigned i = 0 ; i < m.node_list.size(); ++i )
//        m.node_list[i].theta_nodal = sin( std::sqrt( i ) * 5. );
//    refinement_if_constraints_or_nodal_field_sup( m, f, ExtractDM< theta_nodal_DM >(), 0.75, spread_cut );

//    for( unsigned n = 0 ; n < m.elem_list.size(); ++n )
//        m.elem_list[n]->set_field( "theta_elem_EET", sin( std::sqrt( n ) * 5. ) );
//    refinement_if_elem_field_sup( m, theta_elem_EET_DM(), 0.75, spread_cut );

//    for( unsigned n = 0 ; n < m.elem_list.size(); ++n )
//        m.elem_list[n]->set_field( "theta_elem_EET", sin( std::sqrt( n ) * 5. ) );
//    refinement_if_constraints_or_elem_field_sup( m, f, theta_elem_EET_DM(), 0.75, spread_cut );

//    for( unsigned n = 0 ; n < m.elem_list.size(); ++n )
//        m.elem_list[n]->set_field( "theta_elem_EET", sin( std::sqrt( n ) * 5. ) );
//    smoothing( m, ExtractDM< theta_nodal_DM >(), ExtractDM< theta_elem_EET_DM >() );
//    refinement_if_nodal_field_sup( m, ExtractDM< theta_nodal_DM >(), 0.75, spread_cut );

//    divide( m );

//    divide_element( m, Vec<unsigned>( 5, 50, 80, 100, 200 ) );

//    replace_Quad_by_Triangle( m );
//    replace_Hexa_by_Tetra( m );
//    replace_Wedge_by_Tetra( m );

    display_mesh_carac( m );
    display_mesh( m );

    return 0;
}
*/

/// ------ ///
/// HDF 2D ///
/// ------ ///
/*
int main() {
    Hdf hdf( "/Users/Op/Documents/Recherche/LMT/Code_Cpp_FEM_Error_Estimation/ERROR_ESTIMATION/DATA/square-32x32.hdf5" );

    Vec<int> s;
    Tens3<double> u, tau;

//    hdf.read_size( "/u", s );
//    u.resize( s );
//    hdf.read_data( "/u", u.ptr(), s, s );
//    PRINT( s );
//    PRINT( u( 0, 0, 0 ) );
//    PRINT( u( 1, 5, 10 ) );

    hdf.read( "/u", u );
    PRINT( u( 0, 0, 0 ) );
    PRINT( u( 1, 5, 10 ) );

//    hdf.read_size( "/tau", s );
//    tau.resize( s );
//    hdf.read_data( "/tau", tau.ptr(), s, s );
//    PRINT( s );
//    PRINT( tau( 0, 0, 0 ) );
//    PRINT( tau( 1, 5, 10 ) );

    hdf.read( "/tau", tau );
    PRINT( tau( 0, 0, 0 ) );
    PRINT( tau( 1, 5, 10 ) );
}
*/

/// ------ ///
/// HDF 3D ///
/// ------ ///

int main() {
    Hdf hdf( "/Users/Op/Documents/Recherche/LMT/Code_Cpp_FEM_Error_Estimation/ERROR_ESTIMATION/DATA/hashin-32x32x32.hdf5" );

    Vec<int> s;
    unsigned grid_size;
    double k0, k1, k2, k3, nu;
    Tens3<double> f1, f2, f3;
    Tens4<double> u, tau;

    hdf.read_tag( "/", "grid_size", grid_size );
    PRINT( grid_size );
    hdf.read_tag( "/", "k0", k0 );
    PRINT( k0 );
    hdf.read_tag( "/", "k1", k1 );
    PRINT( k1 );
    hdf.read_tag( "/", "k2", k2 );
    PRINT( k2 );
    hdf.read_tag( "/", "k3", k3 );
    PRINT( k3 );
    hdf.read_tag( "/", "nu", nu );
    PRINT( nu );

//    hdf.read_size( "/f1", s );
//    f1.resize( s );
//    hdf.read_data( "/f1", f1.ptr(), s, s );
//    PRINT( s );
//    PRINT( f1( 0, 0, 0 ) );
//    PRINT( f1( 1, 5, 10 ) );

    hdf.read( "/f1", f1 );
    PRINT( f1( 0, 0, 0 ) );
    PRINT( f1( 1, 5, 10 ) );

//    hdf.read_size( "/f2", s );
//    f2.resize( s );
//    hdf.read_data( "/f2", f2.ptr(), s, s );
//    PRINT( s );
//    PRINT( f2( 0, 0, 0 ) );
//    PRINT( f2( 1, 5, 10 ) );

    hdf.read( "/f2", f2 );
    PRINT( f2( 0, 0, 0 ) );
    PRINT( f2( 1, 5, 10 ) );

    f3.resize(grid_size);
    f3.set(1.);
    f3 -= f1 + f2;
    PRINT( f3( 12, 8, 25 ) );
    PRINT( 1 - f1( 12, 8, 25 ) - f2( 12, 8, 25 ) );

//    hdf.read_size( "/filtered/u", s );
//    u.resize( s );
//    hdf.read_data( "/filtered/u", u.ptr(), s, s );
//    PRINT( s );
//    PRINT( u( 0, 0, 0, 0 ) );
//    PRINT( u( 1, 5, 10, 15 ) );

    hdf.read( "/filtered/u", u );
    PRINT( u( 0, 0, 0, 0 ) );
    PRINT( u( 1, 5, 10, 15 ) );

//    hdf.read_size( "/filtered/tau", s );
//    tau.resize( s );
//    hdf.read_data( "/filtered/tau", tau.ptr(), s, s );
//    PRINT( s );
//    PRINT( tau( 0, 0, 0, 0 ) );
//    PRINT( tau( 1, 5, 10, 15 ) );

    hdf.read( "/filtered/tau", tau );
    PRINT( tau( 0, 0, 0, 0 ) );
    PRINT( tau( 1, 5, 10, 15 ) );

}


/// ---------- ///
/// eig_lapack ///
/// ---------- ///
/*
#include "LMT/include/containers/eig_lapack.h"

int main ( int argc, char **argv ) {
    typedef double T;

    Mat<T, Sym<> > Ms( 3, 3, 0. );
    Ms.diag() = 5.;
    Ms( 1, 0 ) = 2.;
    Ms( 2, 0 ) = 1.;
    Ms( 2, 1 ) = -0.3;

    Mat<T> A( Ms );
    PRINTN( A );

    Vec<T> eig_val;
    get_eig_val_sym( A, eig_val );
    PRINT( eig_val );

    Mat<T> eig_vec;
    get_eig_sym( A, eig_val, eig_vec );
    PRINT( eig_val );
    PRINTN( eig_vec );

    PRINTN( eig_vec.row( 0 ) );
    PRINTN( A * eig_vec.row( 0 ) - eig_val[ 0 ] * eig_vec.row( 0 ) );
    PRINTN( A * trans( eig_vec ) - trans( eig_vec ) * diag( eig_val ) );
    PRINTN( A - trans( eig_vec ) * diag( eig_val ) * eig_vec );

    Ms.set( 0. );
    Ms.diag() = 1.;
    Ms( 1, 1 ) = 2.;
    Mat<T> B( Ms );
    PRINTN( B );

    get_eig_val_gen( A, B, eig_val );
    PRINT( eig_val );

    get_eig_gen( A, B, eig_val, eig_vec );
    PRINT( eig_val );
    PRINTN( eig_vec );

    PRINTN( eig_vec.row( 0 ) );
    PRINTN( A * eig_vec.row( 0 ) - eig_val[ 0 ] * B * eig_vec.row( 0 ) );
    PRINTN( A * trans( eig_vec ) - B * trans( eig_vec ) * diag( eig_val ) );

    return 0.;
}
*/

/// -------------------------- ///
/// eigen_problem_using_lapack ///
/// -------------------------- ///
/*
#include "LMT/include/containers/eigen_problem_using_lapack.h"

int main ( int argc, char **argv ) {
    typedef double T;

    Mat<T, Sym<> > Ms( 3, 3, 0. );
    Ms.diag() = 5.;
    Ms( 1, 0 ) = 2.;
    Ms( 2, 0 ) = 1.;
    Ms( 2, 1 ) = -0.3;

    Mat<T> A( Ms );
    PRINTN( A );

    Vec<T> eig_val;
    eigen_values_using_lapack( A, eig_val );
    PRINT( eig_val );

    Mat<T> eig_vec;
    eigen_problem_using_lapack( A, eig_val, eig_vec );
    PRINT( eig_val );
    PRINTN( eig_vec );

    PRINTN( eig_vec.col( 0 ) );
    PRINTN( A * eig_vec.col( 0 ) - eig_val[ 0 ] * eig_vec.col( 0 ) );
    PRINTN( A * eig_vec - trans( eig_vec ) * diag( eig_val ) );
    PRINTN( A - eig_vec * diag( eig_val ) * trans( eig_vec ) );

    Ms.set( 0. );
    Ms.diag() = 1.;
    Ms( 1, 1 ) = 2.;
    Mat<T> B( Ms );
    PRINTN( B );

    generalized_eigen_values_using_lapack( A, B, eig_val );
    PRINT( eig_val );

    generalized_eigen_problem_using_lapack( A, B, eig_val, eig_vec );
    PRINT( eig_val );
    PRINTN( eig_vec );

    PRINTN( eig_vec.col( 0 ) );
    PRINTN( A * eig_vec.col( 0 ) - eig_val[ 0 ] * B * eig_vec.col( 0 ) );
    PRINTN( A * eig_vec - B * eig_vec * diag( eig_val ) );

    return 0.;
}
*/

/// --------------- ///
/// rand / RAND_MAX ///
/// --------------- ///
/*
int main ( int argc, char **argv ) {
    double random = rand();
    cout << random << endl;
    double random_max = RAND_MAX;
    cout << random_max << endl;
    cout << random/random_max << endl;
    return 0;
}
*/

/// --------- ///
/// print env ///
/// --------- ///
/*
int main( int argc, char **argv, char** env ) {
    /// Print all environment variables
    while (*env)
        printf("%s\n", *env++);
    return 0;
}
*/
