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

//int main( int argc, char **argv ) {
//    TicToc t_total;
//    t_total.start();
//    static const unsigned dim = 2;
//    static const bool wont_add_nz = true;
//    typedef Mesh<Mesh_carac_error_estimation<double,dim> > TM;
//    typedef Formulation<TM,FormulationElasticity,DefaultBehavior,double,wont_add_nz> TF;
//    typedef TM::Pvec Pvec;
//    typedef TM::TNode::T T;
//    static const string structure = "structure_crack";
//    static const string mesh_size = "coarse";
//    static const string loading = "pull";
//    static const unsigned deg_p = 1;
//    static const unsigned deg_k = 3;
//    static const string boundary_condition_D = "penalty";

//    display_pb( dim, structure, deg_p  );

//    /// Maillage du pb direct
//    /// ---------------------
//    TM m; // declaration d'un maillage de type TM
//    set_mesh( m, structure, mesh_size, loading, deg_p );

//    /// Formulation du pb direct
//    /// ------------------------
//    TF f( m ); // creation d'une formulation du type TF avec le maillage m

//    /// Proprietes materiaux du pb direct
//    /// ---------------------------------
//    set_material_properties( f, m, structure );

//    /// Conditions aux limites du pb direct
//    /// -----------------------------------
//    set_constraints( f, m, boundary_condition_D, "direct", structure, loading );
////    check_constraints( f );
//    set_load_conditions( m, structure, loading, mesh_size );

//    /// Resolution du pb direct
//    /// -----------------------
//    cout << "Resolution du pb direct" << endl << endl;
//    TicToc t;
//    t.start();
//    f.solve();
//    t.stop();
//    cout << "Temps de calcul de la resolution du pb direct = " << t.res << endl << endl;

////    f.allocate_matrices();
////    f.shift();
////    f.assemble();
//////    f.solve_system();
////    f.get_initial_conditions();
////    f.update_variables();
////    f.call_after_solve();

//    /// Verification de l'equilibre du pb direct
//    /// ----------------------------------------
//    check_equilibrium( f, "direct" );

////    cout << m.type_elements() << endl;

//    TM m_ref = m;

//    divide( m_ref );

////    TF f_ref( m_ref );
////    set_material_properties( f_ref, m_ref, structure );
////    set_constraints( f_ref, m_ref, boundary_condition_D, "direct", structure, loading );
////    set_load_conditions( m_ref, structure, loading, mesh_size );
////    f_ref.solve();

//    display( m_ref );

//}

/// ---------- ///
/// refinement ///
/// ---------- ///
int main( int argc, char **argv ) {
    static const unsigned dim = 2;
//    typedef Mesh<MeshCaracStd<dim,2> > TM;
    typedef Mesh<Mesh_carac_error_estimation<double,dim> > TM;
    typedef TM::Pvec Pvec;
    typedef TM::TNode::T T;

    TM m;
//    m.add_node( TM::Pvec( 0, 1 ) );
//    PRINT( generate( m.node_list, ExtractDM<pos_DM>() ) );
//    PRINT( generate( m.node_list, ExtractDMi<pos_DM>( 0 ) ) );
//    PRINT( generate( m.node_list, ExtractDMi<pos_DM>( 1 ) ) );
    make_rect( m, Quad(), Pvec( 0, 0 ), Pvec( 1., 1. ), Pvec( 20, 20 ) );

    bool spread_cut = false;

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

    for( unsigned n = 0 ; n < m.elem_list.size(); ++n )
        m.elem_list[n]->set_field( "theta_elem_EET", sin( std::sqrt( n ) * 5. ) );
    refinement_if_elem_field_sup( m, ExtractDM< theta_elem_EET_DM >(), 0.75, spread_cut );

//    for( unsigned n = 0 ; n < m.elem_list.size(); ++n )
//        m.elem_list[n]->set_field( "theta_elem_EET", sin( std::sqrt( n ) * 5. ) );
//    smoothing( m, ExtractDM< theta_nodal_DM >(), ExtractDM< theta_elem_EET_DM >() );
//    refinement_if_nodal_field_sup( m, ExtractDM< theta_nodal_DM >(), 0.75, spread_cut );

//    divide( m );

//    divide_element( m, Vec<unsigned>( 5, 50, 80, 100, 200 ) );

    display_mesh( m );

    return 0;
}

/// ---------- ///
/// eig_lapack ///
/// ---------- ///
//#include "LMT/include/containers/eig_lapack.h"

//int main ( int argc, char **argv ) {
//    typedef double T;

//    Mat<T, Sym<> > Ms( 3, 3, 0. );
//    Ms.diag() = 5.;
//    Ms( 1, 0 ) = 2.;
//    Ms( 2, 0 ) = 1.;
//    Ms( 2, 1 ) = -0.3;

//    Mat<T> A( Ms );
//    PRINTN( A );

//    Vec<T> eig_val;
//    get_eig_val_sym( A, eig_val );
//    PRINT( eig_val );

//    Mat<T> eig_vec;
//    get_eig_sym( A, eig_val, eig_vec );
//    PRINT( eig_val );
//    PRINTN( eig_vec );

//    PRINTN( eig_vec.row( 0 ) );
//    PRINTN( A * eig_vec.row( 0 ) - eig_val[ 0 ] * eig_vec.row( 0 ) );
//    PRINTN( A * trans( eig_vec ) - trans( eig_vec ) * diag( eig_val ) );
//    PRINTN( A - trans( eig_vec ) * diag( eig_val ) * eig_vec );

//    Ms.set( 0. );
//    Ms.diag() = 1.;
//    Ms( 1, 1 ) = 2.;
//    Mat<T> B( Ms );
//    PRINTN( B );

//    get_eig_val_gen( A, B, eig_val );
//    PRINT( eig_val );

//    get_eig_gen( A, B, eig_val, eig_vec );
//    PRINT( eig_val );
//    PRINTN( eig_vec );

//    PRINTN( eig_vec.row( 0 ) );
//    PRINTN( A * eig_vec.row( 0 ) - eig_val[ 0 ] * B * eig_vec.row( 0 ) );
//    PRINTN( A * trans( eig_vec ) - B * trans( eig_vec ) * diag( eig_val ) );

//    return 0.;
//}

/// -------------------------- ///
/// eigen_problem_using_lapack ///
/// -------------------------- ///
//#include "LMT/include/containers/eigen_problem_using_lapack.h"

//int main ( int argc, char **argv ) {
//    typedef double T;

//    Mat<T, Sym<> > Ms( 3, 3, 0. );
//    Ms.diag() = 5.;
//    Ms( 1, 0 ) = 2.;
//    Ms( 2, 0 ) = 1.;
//    Ms( 2, 1 ) = -0.3;

//    Mat<T> A( Ms );
//    PRINTN( A );

//    Vec<T> eig_val;
//    eigen_values_using_lapack( A, eig_val );
//    PRINT( eig_val );

//    Mat<T> eig_vec;
//    eigen_problem_using_lapack( A, eig_val, eig_vec );
//    PRINT( eig_val );
//    PRINTN( eig_vec );

//    PRINTN( eig_vec.col( 0 ) );
//    PRINTN( A * eig_vec.col( 0 ) - eig_val[ 0 ] * eig_vec.col( 0 ) );
//    PRINTN( A * eig_vec - trans( eig_vec ) * diag( eig_val ) );
//    PRINTN( A - eig_vec * diag( eig_val ) * trans( eig_vec ) );

//    Ms.set( 0. );
//    Ms.diag() = 1.;
//    Ms( 1, 1 ) = 2.;
//    Mat<T> B( Ms );
//    PRINTN( B );

//    generalized_eigen_values_using_lapack( A, B, eig_val );
//    PRINT( eig_val );

//    generalized_eigen_problem_using_lapack( A, B, eig_val, eig_vec );
//    PRINT( eig_val );
//    PRINTN( eig_vec );

//    PRINTN( eig_vec.col( 0 ) );
//    PRINTN( A * eig_vec.col( 0 ) - eig_val[ 0 ] * B * eig_vec.col( 0 ) );
//    PRINTN( A * eig_vec - B * eig_vec * diag( eig_val ) );

//    return 0.;
//}

/// --------------- ///
/// rand / RAND_MAX ///
/// --------------- ///
//int main ( int argc, char **argv ) {
//    double random = rand();
//    cout << random << endl;
//    double random_max = RAND_MAX;
//    cout << random_max << endl;
//    cout << random/random_max << endl;
//    return 0;
//}

/// --------- ///
/// print env ///
/// --------- ///
//int main( int argc, char **argv, char** env ) {
//    /// Print all environment variables
//    while (*env)
//        printf("%s\n", *env++);
//    return 0;
//}
