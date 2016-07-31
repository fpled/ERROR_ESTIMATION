//
// C++ Interface: Calcul_goal_oriented_estimation
//
// Description: calcul des bornes d'erreur sur une quantite d'interet locale
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_goal_oriented_error_estimation_h
#define Calcul_goal_oriented_error_estimation_h

#include "LMT/include/mesh/make_mesh_from_mask.h"
#include "LMT/include/containers/eig_lapack.h"
#include "LMT/include/containers/QR_lapack.h"
#include "LMT/include/mesh/mesh_transformation.h"
#include "INTEREST_QUANTITY/Calcul_interest_quantity.h"

using namespace LMT;
using namespace std;

/// Calcul standard des bornes d'erreur sur la quantite d'interet locale (avec ou sans introduction de sigma_hat_m)
/// ---------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T>
void calcul_standard_local_error_bounds( TM &m, TM &m_adjoint, const TF &f, const TF &f_adjoint, const string &method, const T &theta, const T &theta_adjoint, const Vec<T> &theta_adjoint_elem, const Vec<unsigned> &correspondance_elem_m_adjoint_to_elem_m, const Vec< Vec<T> > &dep_hat, const T &I_h, const T &I_hh, const bool want_introduction_sigma_hat_m = true ) {

    cout << "Bornes d'erreur standard sur la quantite d'interet locale" << endl;
    cout << "---------------------------------------------------------" << endl << endl;

    if ( not want_introduction_sigma_hat_m ) {
        cout << "xi = theta * theta_adjoint" << endl;
        cout << "   = " << theta * theta_adjoint << endl << endl;
    }
    else {
        cout << "xi = 1/2 * ( theta * theta_adjoint )" << endl;
        cout << "   = " << theta * theta_adjoint / 2 << endl << endl;
    }
    
    Vec<T> theta_elem_proj_on_adjoint;
    theta_elem_proj_on_adjoint.resize( m_adjoint.elem_list.size() );
    theta_elem_proj_on_adjoint.set( 0. );
    
    Calcul_Error_Estimate_Proj_On_Adjoint<TM, TF, T> calcul_error_estimate_proj_on_adjoint;
    calcul_error_estimate_proj_on_adjoint.m = &m;
    calcul_error_estimate_proj_on_adjoint.m_adjoint = &m_adjoint;
    calcul_error_estimate_proj_on_adjoint.f = &f;
    calcul_error_estimate_proj_on_adjoint.f_adjoint = &f_adjoint;
    calcul_error_estimate_proj_on_adjoint.method = &method;
    calcul_error_estimate_proj_on_adjoint.correspondance_elem_m_adjoint_to_elem_m = &correspondance_elem_m_adjoint_to_elem_m;
    calcul_error_estimate_proj_on_adjoint.dep_hat = &dep_hat;
    calcul_error_estimate_proj_on_adjoint.theta_elem_proj_on_adjoint = &theta_elem_proj_on_adjoint;
    
    apply( m_adjoint.elem_list, calcul_error_estimate_proj_on_adjoint );

    T sum_theta_direct_adjoint = 0.;

    for (unsigned n=0;n<m_adjoint.elem_list.size();++n)
        sum_theta_direct_adjoint += sqrt( theta_elem_proj_on_adjoint[ n ] * theta_adjoint_elem[ n ] );
    
    T xi_inf = 0.;
    T xi_sup = 0.;
    
    if ( not want_introduction_sigma_hat_m ) {
        cout << "somme des contributions elementaires aux estimateurs d'erreur globale :" << endl;
        cout << "sum( sqrt( theta_elem_proj_on_adjoint^2 * theta_adjoint_elem^2 ) ) = " << sum_theta_direct_adjoint << endl << endl;
        
        xi_inf = I_h + I_hh - theta * theta_adjoint;
        xi_sup = I_h + I_hh + theta * theta_adjoint;
        
        cout << "borne inferieure :" << endl;
        cout << "xi_inf = I_h + I_hh - ( theta * theta_adjoint )" << endl;
        cout << "       = " << xi_inf << endl << endl;
        
        cout << "borne superieure :" << endl;
        cout << "xi_sup = I_h + I_hh + ( theta * theta_adjoint )" << endl;
        cout << "       = " << xi_sup << endl << endl;
    }
    else {
        cout << "somme des contributions elementaires aux estimateurs d'erreur globale :" << endl;
        cout << "1/2 * sum( sqrt( theta_elem_proj_on_adjoint^2 * theta_adjoint_elem^2 ) ) = " << sum_theta_direct_adjoint / 2 << endl << endl;
        
        xi_inf = I_h + I_hh - theta * theta_adjoint / 2;
        xi_sup = I_h + I_hh + theta * theta_adjoint / 2;
        
        cout << "borne inferieure :" << endl;
        cout << "xi_inf = I_h + I_hh - 1/2 * ( theta * theta_adjoint )" << endl;
        cout << "       = " << xi_inf << endl << endl;
        
        cout << "borne superieure :" << endl;
        cout << "xi_sup = I_h + I_hh + 1/2 * ( theta * theta_adjoint )" << endl;
        cout << "       = " << xi_sup << endl << endl;
    }
}

/// Calcul ameliore des bornes d'erreur sur la quantite d'interet locale (avec ou sans introduction de sigma_hat_m)
/// ---------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class Pvec>
void calcul_enhanced_local_error_bounds( TM &m, TM &m_adjoint, const TF &f, const TF &f_adjoint, const string &structure, const unsigned &deg_p, const string &method, const string &method_adjoint, const string &local_improvement, const string &shape, const T &k_min, const T &k_max, const T &k_opt, const string &interest_quantity, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list_interest_quantity, const unsigned &node_interest_quantity, const Pvec &pos_interest_quantity, const Pvec &pos_crack_tip, const T &radius_Ri, const T &radius_Re, const bool &spread_cut, const T &theta, const T &theta_adjoint, const Vec< Vec<T> > &dep_hat, const Vec< Vec<T> > &dep_adjoint_hat, const T &I_h, const T &I_hh, const string &integration_k, const unsigned &integration_nb_points, const bool want_introduction_sigma_hat_m = true, const bool want_solve_eig_local_improvement = false, const bool use_mask_eig_local_improvement = false, const bool display_vtu_lambda = false, const bool display_vtu_adjoint_lambda = false, const string &prefix = "paraview_direct", const string &prefix_adjoint = "paraview_adjoint", const bool disp = false ) {

    static const unsigned dim = TM::dim;

    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    typedef Mat<T, Sym<> > TMatSym;
    typedef Mat<T, Gen<>, SparseLine<> > TMatGenSparse;
    typedef Mat<T, Gen<> > TMatGen;
    typedef Mat<T, Diag<> > TMatDiag;
    
    cout << "Amelioration des bornes d'erreur sur la quantite d'interet locale" << endl;
    cout << "-----------------------------------------------------------------" << endl << endl;

    cout << "technique d'amelioration de l'erreur locale : " << local_improvement << endl;
    cout << "forme geometrique des domaines homothetiques = " << shape << endl;
    if ( local_improvement == "steklov" ) {
        cout << "parametres des domaines homothetiques : k_min = " << k_min << endl << " et k_max = " << k_max << endl << endl;
    }
    else if ( local_improvement == "rayleigh" )
        cout << "parametre du domaine homothetique : k_opt = " << k_opt << endl << endl;
    
    /// Construction du centre et de la taille du domaine homothetique 
    /// --------------------------------------------------------------
    Vec<T> domain_length;
    Pvec domain_center( 0. );
    construct_center_length_domain( m, deg_p, shape, interest_quantity, pointwise_interest_quantity, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_Re, domain_center, domain_length );

    /// Calcul de la constante dans l'amelioration
    /// ------------------------------------------
    
    T h = 0.;
    if ( not want_solve_eig_local_improvement ) {
        if ( disp )
            cout << "Calcul analytique de la constante h" << endl << endl;
        
        if ( local_improvement == "steklov" ) {
            calcul_h_alternativeontype( m, Number< AreSameType< typename ExtractDM<la_DM>::ReturnType<TM>::T, void >::res >(), Number< AreSameType< typename ExtractDM<mu_DM>::ReturnType<TM>::T, void >::res >(), h, shape );
            cout << "constante de Steklov : h = " << h << endl;
        }
        else if ( local_improvement == "rayleigh" ) {
            if ( dim == 2 )
                h = 2;
            else if ( dim == 3 )
                h = 3;
            else
                cerr << "forme " << shape << " non implementee pour amelioration locale basee sur rayleigh..." << endl << endl;
            cout << "constante de Rayleigh : h = " << h << endl;
        }
        else
            cerr << "amelioration locale basee sur " << local_improvement << " non implementee..." << endl << endl;
        cout << "constante h dans l'amelioration : h = " << h << endl << endl;
    }
    else {
        if ( disp )
            cout << "Calcul numerique de la constante h" << endl << endl;
        
        if ( use_mask_eig_local_improvement ) {
            
            string filename_mask;
            
            if ( dim ==2 ) {
                if ( shape == "circle" )
                    filename_mask = "./MASK/circle.png";
                else if ( shape == "half_circle" )
                    filename_mask = "./MASK/hal_circle.png";
                else if ( shape == "circle_crack" )
                    filename_mask = "./MASK/circle_crack.png";
                else if ( shape == "rectangle" )
                    filename_mask = "./MASK/rectangle.png";
                else
                    cerr << "forme " << shape << " non implementee pour le calcul de la constante dans l'amelioration..." << endl << endl;
            }
            
            Pvec translation( -1000. );
            
            Pvec scale_factor( 1./1000 );
            if ( local_improvement == "rayleigh" )
                scale_factor *= k_opt;
            
            T criterium_eq_gen_eig_pb = 1e-7;
            T residual = 1.;
            T elem_size_mask = 200;
            unsigned iter = 0;
            
            Vec<T> generalized_eig_val;
            Mat<T> generalized_eig_vec;

            while( residual > criterium_eq_gen_eig_pb ) {

                ++iter;
                
                TM m_mask;
                make_mesh_from_mask( filename_mask, m_mask, elem_size_mask );
                translate_mesh( m_mask, translation );
                scale_mesh( m_mask, scale_factor );
//                display( m_mask, shape );
                
                if ( disp )
                    cout << "Construction de la solution numerique du pb aux valeurs propres generalisees a l'iteration " << iter << endl << endl;

                if ( m_mask.node_list.size() ) {
                    if ( remove_lonely_nodes( m_mask ) )
                        cerr << "Des noeuds seuls ont ete retires du maillage du pb aux valeurs propres generalisees..." << endl << endl;
                    if ( not m_mask.elem_list.size() ) {
                        cerr << "Arret brutal, car il n'y a aucun element dans le maillage du pb aux valeurs propres generalisees..." << endl << endl;
                        throw "Baleinou sous caillou...";
                    }
                    if ( disp )
                        display_mesh_carac( m_mask );
                }
                else {
                    cerr << "Arret brutal, car il n'y a aucun noeud dans le maillage du pb aux valeurs propres generalisees..." << endl << endl;
                    throw "Ourson sous gravillon...";
                }
                
                TF f_mask( m_mask );

                set_material_properties( f_mask, m_mask, structure );
                
                if ( local_improvement == "steklov" )
                    m_mask.phi_local_improvement_steklov = 1;
                else if ( local_improvement == "rayleigh" )
                    m_mask.phi_local_improvement_rayleigh = 1;
                else
                    cerr << "amelioration locale basee sur " << local_improvement << " non implementee..." << endl << endl;
                
                f_mask.allocate_matrices();
                
                unsigned n = f_mask.vectors[0].size();
                unsigned r = unsigned(dim*(dim+1)/2);
                
                TMatSymSparse A;
                TMatSymSparse B;
                TMatGenSparse C;
                Vec<Vec<T> > vectors;
                A.resize( n );
                B.resize( n );
                C.resize( r, n );
                
                f_mask.assemble( A, B, vectors );
                f_mask.assemble( C, vectors );
                
                Vec<double> V;
                V.resize( n, 0. );
                
                for (unsigned i=0;i<m_mask.node_list.size();++i) {
                    for (unsigned d=0;d<dim;++d)
                        V[ i * dim + d ] = m_mask.node_list[ i ].pos[ d ];
                }
                
                h = dot( V, A * V ) / dot( V, B * V );
                
                cout << "constante h = " << h << endl;
                
                if ( disp )
                    cout << "residu associe aux conditions d'equilibre = " << C * V << endl;
                
                B.diag() += 1e-12;
                
                TicToc t_eig_gen;
                t_eig_gen.start();
                
                get_eig_gen( A, B, generalized_eig_val, generalized_eig_vec );
                generalized_eig_val.resize( n - r );
                if ( disp ) {
                    cout << "nb valeurs propres generalisees = " << generalized_eig_val.size() << endl;
                    cout << "valeurs propres generalisees : " << endl;
                    cout << generalized_eig_val << endl << endl;
                }
                
                t_eig_gen.stop();
                if ( disp )
                    cout << "temps de calcul de la resolution du pb aux valeurs propres generalisees = " << t_eig_gen.res << endl << endl;
                
                T generalized_eig_val_eq = generalized_eig_val[ 0 ];
                Vec<T> residual_vec = C * generalized_eig_vec.row( 0 );
                residual = norm_2( residual_vec );
                for (unsigned i=0;i<generalized_eig_val.size();++i) {
                    if ( norm_2( C * generalized_eig_vec.row( i ) ) < residual and generalized_eig_val[ i ] > 1e-10 ) {
                        residual_vec = C * generalized_eig_vec.row( i );
                        residual = norm_2( residual_vec );
                        generalized_eig_val_eq = generalized_eig_val[ i ];
                    }
                }
                if ( disp ) {
                    cout << "constante h = " << generalized_eig_val_eq << endl;
                    cout << "residu associe aux conditions d'equilibre = " << residual_vec << endl;
                    cout << "norme du residu sur les conditions d'equilibre = " << residual << endl << endl;
                }
                elem_size_mask /= 2.;
            }
        
//            if ( local_improvement == "steklov" )
//                cout << "constante de Steklov : h = " << max( generalized_eig_val ) << endl;
//            else if ( local_improvement == "rayleigh" )
//                cout << "constante de Rayleigh : h = " << min( generalized_eig_val ) << endl;
        }
        else {
            
            TM m_mask;
            if ( dim == 2 ) {
                if ( shape == "circle" ) {
//                    read_msh_2( m_mask, "MESH/CIRCLE_2D/circle_very_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_2D/circle_coarse_Triangle.msh" );
                    read_msh_2( m_mask, "MESH/CIRCLE_2D/circle_fine_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_2D/circle_very_fine_Triangle.msh" );
                }
                else if ( shape == "half_circle" ) {
//                    read_msh_2( m_mask, "MESH/HALF_CIRCLE_2D/half_circle_very_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/HALF_CIRCLE_2D/half_circle_coarse_Triangle.msh" );
                    read_msh_2( m_mask, "MESH/HALF_CIRCLE_2D/half_circle_fine_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/HALF_CIRCLE_2D/half_circle_very_fine_Triangle.msh" );
                }
                else if ( shape == "circle_crack" ) {
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_very_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_coarse_Triangle.msh" );
                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_fine_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_very_fine_Triangle.msh" );

//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_15_very_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_15_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_15_fine_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_15_very_fine_Triangle.msh" );

//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_30_very_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_30_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_30_fine_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_30_very_fine_Triangle.msh" );

//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_45_very_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_45_coarse_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_45_fine_Triangle.msh" );
//                    read_msh_2( m_mask, "MESH/CIRCLE_CRACK_2D/circle_crack_45_very_fine_Triangle.msh" );
                }
                else if ( shape == "square" or shape == "rectangle" )
                    make_rect( m_mask, Triangle(), Pvec( -1., -1. ), Pvec( 1., 1. ), Pvec( 101, 101 ) );
                else
                    cerr << "forme " << shape << " non implementee pour le calcul de la constante dans l'amelioration..." << endl << endl;
            }
            else if ( dim ==3 ) {
                if ( shape == "sphere" ) {
//                    read_msh_2( m_mask, "MESH/SPHERE_3D/sphere_very_coarse_Tetra.msh" );
//                    read_msh_2( m_mask, "MESH/SPHERE_3D/sphere_coarse_Tetra.msh" );
                    read_msh_2( m_mask, "MESH/SPHERE_3D/sphere_fine_Tetra.msh" );
//                    read_msh_2( m_mask, "MESH/SPHERE_3D/sphere_very_fine_Tetra.msh" );
                }
                else if ( shape == "half_sphere" ) {
//                    read_msh_2( m_mask, "MESH/HALF_SPHERE_3D/half_sphere_very_coarse_Tetra.msh" );
//                    read_msh_2( m_mask, "MESH/HALF_SPHERE_3D/half_sphere_coarse_Tetra.msh" );
                    read_msh_2( m_mask, "MESH/HALF_SPHERE_3D/half_sphere_fine_Tetra.msh" );
//                    read_msh_2( m_mask, "MESH/HALF_SPHERE_3D/half_sphere_very_fine_Tetra.msh" );
                }
                else if ( shape == "cube" or shape == "cuboid" )
                    make_rect( m_mask, Tetra(), Pvec( -1., -1., -1. ), Pvec( 1., 1., 1. ), Pvec( 11, 11, 11 ) );
                else
                    cerr << "forme " << shape << " non implementee pour le calcul de la constante dans l'amelioration..." << endl << endl;
            }
            else
                cerr << "dim non implementee..." << endl << endl;
            
            Pvec scale_factor( 1. );
//            if ( local_improvement == "rayleigh" )
//                scale_factor *= k_opt;
            scale_mesh( m_mask, scale_factor );
//            display( m_mask, shape );
            
            if ( disp )
                cout << "Construction de la solution numerique du pb aux valeurs propres generalisees" << endl << endl;

            if ( m_mask.node_list.size() ) {
                if ( remove_lonely_nodes( m_mask ) )
                    cerr << "Des noeuds seuls ont ete retires du maillage du pb aux valeurs propres generalisees..." << endl << endl;
                if ( not m_mask.elem_list.size() ) {
                    cerr << "Arret brutal, car il n'y a aucun element dans le maillage du pb aux valeurs propres generalisees..." << endl << endl;
                    throw "Baleinou sous caillou...";
                }
                if ( disp )
                    display_mesh_carac( m );
            }
            else {
                cerr << "Arret brutal, car il n'y a aucun noeud dans le maillage du pb aux valeurs propres generalisees..." << endl << endl;
                throw "Ourson sous gravillon...";
            }
            
            TF f_mask( m_mask );
            
            if ( local_improvement == "steklov" )
                m_mask.phi_local_improvement_steklov = 1;
            else if ( local_improvement == "rayleigh" )
                m_mask.phi_local_improvement_rayleigh = 1;
            else
                cerr << "amelioration locale basee sur " << local_improvement << " non implementee..." << endl << endl;
            
            f_mask.allocate_matrices();
            
            unsigned n = f_mask.vectors[0].size();
            unsigned r = unsigned(dim*(dim+1)/2);
            
            TMatSymSparse A;
            TMatSymSparse B;
            TMatGenSparse C;
            Vec<Vec<T> > vectors;
            A.resize( n );
            B.resize( n );
            C.resize( r, n );
            
            f_mask.assemble( A, B, vectors );
            f_mask.assemble( C, vectors );
            
            Vec<double> V;
            V.resize( n, 0. );
            
            for (unsigned i=0;i<m_mask.node_list.size();++i) {
                for (unsigned d=0;d<dim;++d)
                    V[ i * dim + d ] = m_mask.node_list[ i ].pos[ d ];
            }
            
            h = dot( V, A * V ) / dot( V, B * V );
            
            cout << "constante h :" << endl;
            cout << "h = " << h << endl << endl;
            
            if ( disp )
                cout << "residu associe aux conditions d'equilibre = " << C * V << endl;
            
            B.diag() += 1e-12;
            
            TicToc t_eig_gen;
            t_eig_gen.start();
            
            Vec<T> generalized_eig_val;
            Mat<T> generalized_eig_vec;
            get_eig_gen( A, B, generalized_eig_val, generalized_eig_vec );
            generalized_eig_val.resize( n - r );
            if ( disp ) {
                cout << "nb valeurs propres generalisees = " << generalized_eig_val.size() << endl;
                cout << "valeurs propres generalisees : " << endl;
                cout << generalized_eig_val << endl << endl;
            }
            
            t_eig_gen.stop();
            if ( disp )
                cout << "temps de calcul de la resolution du pb aux valeurs propres generalisees = " << t_eig_gen.res << endl << endl;
            
            T generalized_eig_val_eq = generalized_eig_val[ 0 ];
            Vec<T> residual_vec = C * generalized_eig_vec.row( 0 );
            T residual = norm_2( residual_vec );
            for (unsigned i=0;i<generalized_eig_val.size();++i) {
                if ( norm_2( C * generalized_eig_vec.row( i ) ) < residual and generalized_eig_val[ i ] > 1e-10 ) {
                    residual_vec = C * generalized_eig_vec.row( i );
                    residual = norm_2( residual_vec );
                    generalized_eig_val_eq = generalized_eig_val[ i ];
                }
            }
            cout << "constante h = " << generalized_eig_val_eq << endl;
            cout << "residu associe aux conditions d'equilibre = " << residual_vec << endl;
            cout << "norme du residu sur les conditions d'equilibre = " << residual << endl << endl;
    
//            if ( local_improvement == "steklov" )
//                cout << "constante de Steklov : h = " << max( generalized_eig_val ) << endl;
//            else if ( local_improvement == "rayleigh" )
//                cout << "constante de Rayleigh : h = " << min( generalized_eig_val ) << endl;
        }
    }
    
    /// Maillages extraits du maillage direct autour de la quantite d'interet
    /// ---------------------------------------------------------------------
    TM m_lambda_min, m_lambda_max, m_lambda_opt;
    if ( local_improvement == "steklov" ) {
        set_mesh_domain( m_lambda_min, m, shape, k_min, domain_length, domain_center, spread_cut );
        set_mesh_domain( m_lambda_max, m, shape, k_max, domain_length, domain_center, spread_cut );
    }
    else if ( local_improvement == "rayleigh" )
        set_mesh_domain( m_lambda_opt, m, shape, k_opt, domain_length, domain_center, spread_cut );
    
    /// Maillages extraits du maillage adjoint autour de la quantite d'interet
    /// ----------------------------------------------------------------------
    TM m_adjoint_lambda_min, m_adjoint_lambda_opt;
    if ( local_improvement == "steklov" )
        set_mesh_domain( m_adjoint_lambda_min, m_adjoint, shape, k_min, domain_length, domain_center, spread_cut );
    else if ( local_improvement == "rayleigh" )
        set_mesh_domain( m_adjoint_lambda_opt, m_adjoint, shape, k_opt, domain_length, domain_center, spread_cut );

    /// Formulation des pbs extraits associes aux pbs direct/adjoint
    /// ------------------------------------------------------------
    TF f_lambda_min( m_lambda_min );
    TF f_lambda_max( m_lambda_max );
    TF f_adjoint_lambda_min( m_adjoint_lambda_min );
    
    TF f_lambda_opt( m_lambda_opt );
    TF f_adjoint_lambda_opt( m_adjoint_lambda_opt );
    
    /// --------------------------------------------------------------------------------------------------------------------------------------------------///
    /// Construction d'un champ de contrainte admissible et Calcul d'un estimateur d'erreur globale associe aux maillages extraits des pbs direct/adjoint ///
    /// --------------------------------------------------------------------------------------------------------------------------------------------------///
    
    if ( disp ) {
        cout << "Construction d'un champ de contrainte admissible par element" << endl;
        if ( local_improvement == "steklov" )
            cout << "Calcul d'un estimateur d'erreur globale sur les structures extraites associees aux pbs direct/adjoint de type " << shape << " pour le parametre min lambda_min = " << k_min << " et le parametre max lambda_max = " << k_max << endl;
        else if ( local_improvement == "rayleigh" )
            cout << "Calcul d'un estimateur d'erreur globale sur la structure extraite associee aux pbs direct/adjoint de type " << shape << " pour le parametre lambda = " << k_opt << endl;
        cout << endl;
    }
    
    T theta_lambda_min, theta_lambda_max, theta_adjoint_lambda_min;
    T theta_lambda_opt, theta_adjoint_lambda_opt;
    
    Vec< Vec<T> > dep_hat_lambda_min, dep_hat_lambda_max, dep_hat_lambda_opt;
    
    if ( local_improvement == "steklov" ) {
        dep_hat_lambda_min.resize( m_lambda_min.elem_list.size() );
        dep_hat_lambda_max.resize( m_lambda_max.elem_list.size() );
        calcul_error_estimate_lambda( m, m_lambda_min, f, f_lambda_min, "direct", method, shape, k_min, theta_lambda_min, dep_hat, dep_hat_lambda_min );
        calcul_error_estimate_lambda( m, m_lambda_max, f, f_lambda_max, "direct", method, shape, k_max, theta_lambda_max, dep_hat, dep_hat_lambda_max );

        const string prefix_lambda_min = prefix + "_lambda_min_" + to_string( k_min );
        const string prefix_lambda_max = prefix + "_lambda_max_" + to_string( k_max );
        if ( display_vtu_lambda ) {
            display( m_lambda_min, prefix_lambda_min );
            display( m_lambda_max, prefix_lambda_max );
        }
        else {
            save( m_lambda_min, prefix_lambda_min );
            save( m_lambda_max, prefix_lambda_max );
        }
    }
    else if ( local_improvement == "rayleigh" ) {
        dep_hat_lambda_opt.resize( m_lambda_opt.elem_list.size() );
        calcul_error_estimate_lambda( m, m_lambda_opt, f, f_lambda_opt, "direct", method, shape, k_opt, theta_lambda_opt, dep_hat, dep_hat_lambda_opt );

        const string prefix_lambda_opt = prefix + "_lambda_opt_" + to_string( k_opt );
        if ( display_vtu_lambda )
            display( m_lambda_opt, prefix_lambda_opt );
        else
            save( m_lambda_opt, prefix_lambda_opt );
    }
    
    Vec< Vec<T> > dep_adjoint_hat_lambda_min, dep_adjoint_hat_lambda_opt;

    if ( local_improvement == "steklov" ) {
        dep_adjoint_hat_lambda_min.resize( m_adjoint_lambda_min.elem_list.size() );
        calcul_error_estimate_lambda( m_adjoint, m_adjoint_lambda_min, f_adjoint, f_adjoint_lambda_min, "adjoint", method_adjoint, shape, k_min, theta_adjoint_lambda_min, dep_adjoint_hat, dep_adjoint_hat_lambda_min );

        const string prefix_adjoint_lambda_min = prefix_adjoint + "_lambda_min_" + to_string( k_min );
        if ( display_vtu_adjoint_lambda )
            display( m_adjoint_lambda_min, prefix_adjoint_lambda_min );
        else
            save( m_adjoint_lambda_min, prefix_adjoint_lambda_min );
    }
    else if ( local_improvement == "rayleigh" ) {
        dep_adjoint_hat_lambda_opt.resize( m_adjoint_lambda_opt.elem_list.size() );
        calcul_error_estimate_lambda( m_adjoint, m_adjoint_lambda_opt, f_adjoint, f_adjoint_lambda_opt, "adjoint", method, shape, k_opt, theta_adjoint_lambda_opt, dep_adjoint_hat, dep_adjoint_hat_lambda_opt );

        string prefix_adjoint_lambda_opt = prefix_adjoint + "_lambda_opt_" + to_string( k_opt );
        if ( display_vtu_adjoint_lambda )
            display( m_adjoint_lambda_opt, prefix_adjoint_lambda_opt );
        else
            save( m_adjoint_lambda_opt, prefix_adjoint_lambda_opt );
    }
    
    T gamma;
    calcul_gamma( m, m_adjoint, m_adjoint_lambda_opt, f, f_adjoint, f_adjoint_lambda_opt, method, local_improvement, shape, k_min, k_max, k_opt, theta_lambda_min, theta_lambda_max, h, domain_center, domain_length, spread_cut, dep_hat, dep_adjoint_hat, integration_k, integration_nb_points, gamma );
    
    T I_hhh;
    if ( local_improvement == "steklov" )
        calcul_correction_interest_quantity_lambda( m_lambda_min, m_adjoint_lambda_min, f_lambda_min, f_adjoint_lambda_min, interest_quantity, method, method_adjoint, dep_hat_lambda_min, dep_adjoint_hat_lambda_min, I_hhh, want_introduction_sigma_hat_m );
    
    cout << "Bornes d'erreur ameliorees sur la quantite d'interet locale" << endl;
    cout << "-----------------------------------------------------------" << endl << endl;
    
    if ( local_improvement == "steklov" ) {
        T chi = 0.;
        if ( not want_introduction_sigma_hat_m ) {
            chi = theta_adjoint_lambda_min * sqrt( pow( k_min / k_max, 1/h ) * pow( theta, 2 ) + gamma ) + theta * sqrt( pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_min, 2 ) );
            cout << "chi = " << theta_adjoint_lambda_min << " * ( " << pow( k_min / k_max, 1/h ) * pow( theta, 2 ) << " + " << gamma << " )^1/2 + " << theta << " * ( " << sqrt( pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_min, 2 ) ) << " )" << endl;
            cout << "    = " << theta_adjoint_lambda_min * sqrt( pow( k_min / k_max, 1/h ) * pow( theta, 2 ) + gamma ) << " + " << theta * sqrt( pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_min, 2 ) ) << endl;
            cout << "    = " << chi << endl << endl;
        }
        else {
            chi = theta_adjoint_lambda_min * sqrt( pow( k_min / k_max, 1/h ) * pow( theta + theta_lambda_max, 2 ) / 4 + gamma ) + theta * sqrt( pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_min, 2 ) ) / 2;
            cout << "chi = " << theta_adjoint_lambda_min << " * ( " << pow( k_min / k_max, 1/h ) * pow( theta + theta_lambda_max, 2 ) / 4 << " + " << gamma << " )^1/2 + 1/2 * " << theta << " * ( " << sqrt( pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_min, 2 ) ) << " )" << endl;
            cout << "    = " << theta_adjoint_lambda_min * sqrt( pow( k_min / k_max, 1/h ) * pow( theta + theta_lambda_max, 2 ) / 4 + gamma ) << " + " << theta * sqrt( pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_min, 2 ) ) / 2 << endl;
            cout << "    = " << chi << endl << endl;
        }
        T chi_inf = I_h + I_hh + I_hhh - chi;
        T chi_sup = I_h + I_hh + I_hhh + chi;
        
        cout << "borne inferieure :" << endl;
        cout << "chi_inf = I_h + I_hh + I_hhh - chi" << endl;
        cout << "        = " << chi_inf << endl << endl;

        cout << "borne superieure :" << endl;
        cout << "chi_sup = I_h + I_hh + I_hhh + chi" << endl;
        cout << "        = " << chi_sup << endl << endl;
    }
    else if ( local_improvement == "rayleigh" ) {
        T zeta = 0.;
        if ( want_introduction_sigma_hat_m ) {
            zeta = theta / 2 * sqrt( pow( gamma, 2 ) + pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_opt, 2 ) ) + theta_lambda_opt / 2 * ( gamma + theta_adjoint_lambda_opt );
            cout << "zeta = " << theta / 2 << " * ( " << pow( gamma, 2 ) << " + " << pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_opt, 2 ) << " )^1/2 + " << theta_lambda_opt / 2 << " * ( " << gamma << " + " << theta_adjoint_lambda_opt << " )" << endl;
            cout << "     = " << theta / 2 * sqrt( pow( gamma, 2 ) + pow( theta_adjoint, 2 ) - pow( theta_adjoint_lambda_opt, 2 ) ) << " + " << theta_lambda_opt / 2 * ( gamma + theta_adjoint_lambda_opt ) << endl;
            cout << "     = "<< zeta << endl << endl;
        }
        else {
            cerr << "TODO..." << endl << endl;
        }
        T zeta_inf = I_h + I_hh - zeta;
        T zeta_sup = I_h + I_hh + zeta;
        
        cout << "borne inferieure :" << endl;
        cout << "zeta_inf = I_h + I_hh - zeta" << endl;
        cout << "         = " << zeta_inf << endl << endl;
        
        cout << "borne superieure :" << endl;
        cout << "zeta_sup = I_h + I_hh + zeta" << endl;
        cout << "         = " << zeta_sup << endl << endl;
    }
    
}

#endif // Calcul_goal_oriented_error_estimation_h
