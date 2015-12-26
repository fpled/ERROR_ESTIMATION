//
// C++ Interface: Calcul_global_error_estimation
//
// Description: construction d'un champ de contrainte admissible et calcul d'un estimateur d'erreur globale
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2015
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_global_error_estimation_homog_h
#define Calcul_global_error_estimation_homog_h

#include "EET/Construct_standard_force_fluxes_EET.h"
#include "EESPT/Construct_standard_force_fluxes_EESPT.h"
#include "ENHANCEMENT/Construct_enhanced_force_fluxes_EET_EESPT.h"
#include "ECRE/Construct_K_hat.h"
#include "ECRE/Construct_F_hat.h"
#include "ECRE/Construct_dep_hat.h"
#include "ECRE/Calcul_error_estimate_prolongation_condition.h"
#include "ECRE/Calcul_error_estimate_prolongation_condition_homog.h"
#include "SPET/Calcul_error_estimate_partition_unity.h"
#include "SPET/Calcul_error_estimate_partition_unity_homog.h"
#include "VERIFICATION/Verification.h"
#include "CRITERIUM_ENHANCEMENT/Construct_criterium_enhancement.h"
#include "Display.h"

using namespace LMT;
using namespace std;

/// Calcul d'un estimateur d'erreur globale
/// ---------------------------------------
template<class TF, class TM, class T>
void calcul_global_error_estimation( TF &f, TM &m, const string &pb, const string &method, const unsigned &cost_function, const T &penalty_val_N, const string &solver, const string &solver_minimisation, const bool &enhancement_with_geometric_criterium, const bool &enhancement_with_estimator_criterium, const string &geometric_criterium, const T &val_geometric_criterium, const T &val_estimator_criterium, T &theta, T &theta_init, Vec<T> &theta_elem, Vec<T> &theta_elem_init, Vec< Vec<T> > &dep_hat, const bool verif_compatibility_conditions = false, const T tol_compatibility_conditions = 1e-6, const bool verif_eq_force_fluxes = false, const T tol_eq_force_fluxes = 1e-6, const bool verif_solver = false, const T tol_solver = 1e-6, const bool verif_solver_enhancement = false, const T tol_solver_enhancement = 1e-6, const bool verif_solver_minimisation = false, const T tol_solver_minimisation = 1e-6, const bool verif_solver_minimisation_enhancement = false, const T tol_solver_minimisation_enhancement = 1e-6, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool want_local_enrichment = false, const bool debug_geometry = false, const bool debug_force_fluxes = false, const bool debug_force_fluxes_enhancement = false, const bool debug_criterium_enhancement = false, const bool debug_error_estimate = false, const bool debug_local_effectivity_index = false, const bool debug_method = false, const bool debug_method_enhancement = false ) {
    
    /// ----------- ///
    /// Methode EET ///
    /// ----------- ///
    
    if ( method.find("EET") != string::npos ) {
        
        display_method( pb, "EET", cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, solver, solver_minimisation );
        
        TicToc t_EET;
        t_EET.start();
        
        Vec< Vec< Vec<T> > > vec_force_fluxes;
        Vec< Mat<T, Sym<> > > K_hat;
        Vec< Vec<T> > F_hat;
        
        theta = 0.;
        theta_init = 0.;
        
        Vec<bool> elem_flag_enh;
        Vec<bool> face_flag_enh;
        Vec<bool> elem_flag_bal;
        Vec<unsigned> elem_list_enh;
        Vec<unsigned> face_list_enh;
        Vec<unsigned> elem_list_bal;
        
        bool enhancement = 0;
        bool balancing = 0;
        
        /// --------------------------------------------------------------------------- ///
        /// Critere d'amelioration geometrique de la construction des densites d'effort ///
        /// --------------------------------------------------------------------------- ///
        
        Vec<T> geometric_ratio;
        
        if ( enhancement_with_geometric_criterium ) {
            enhancement = 1;
            construct_geometric_criterium( m, geometric_criterium, geometric_ratio, debug_criterium_enhancement );
        }
        
        /// -------------------------------------------------------------------------------- ///
        /// Critere d'amelioration sur l'estimateur de la construction des densites d'effort ///
        /// -------------------------------------------------------------------------------- ///
        
        Vec<T> estimator_ratio;
        
        if ( enhancement == 0 or enhancement_with_estimator_criterium ) {
           
            enhancement = 0;
           
            /// Construction des densites d'effort standard
            /// -------------------------------------------
            
            Vec< Vec< Vec<T> > > vec_force_fluxes_standard;
            construct_standard_force_fluxes_EET( m, f, pb, cost_function, enhancement, face_flag_enh, solver_minimisation, vec_force_fluxes_standard, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation, verif_compatibility_conditions, tol_compatibility_conditions, debug_geometry, debug_force_fluxes, debug_method );
            
            /// Verification de l'equilibre des densites d'effort standard
            /// ----------------------------------------------------------
            
            if ( verif_eq_force_fluxes )
                check_equilibrium_force_fluxes( m, f, pb, vec_force_fluxes_standard, tol_eq_force_fluxes, want_local_enrichment, debug_geometry );
            
            /// Resolution des problemes locaux associes aux densites d'effort standard
            /// -----------------------------------------------------------------------
            
            construct_K_hat( m, f, "EET", K_hat, debug_method );
            
            balancing = 0;
            construct_F_hat( m, f, pb, balancing, elem_flag_bal, elem_flag_enh, vec_force_fluxes_standard, F_hat, want_local_enrichment, debug_method, debug_geometry );
            
            construct_dep_hat( m, f, "EET", solver, K_hat, F_hat, dep_hat, verif_solver, tol_solver, debug_method );
            
            /// Construction d'un champ de contrainte admissible par element, Calcul d'un estimateur d'erreur globale
            /// -----------------------------------------------------------------------------------------------------
            
            calcul_error_estimate_prolongation_condition( m, f, pb, "EET", theta, theta_init, theta_elem, theta_elem_init, dep_hat, want_global_discretization_error, want_local_discretization_error, debug_error_estimate, debug_local_effectivity_index, debug_method, debug_method_enhancement );
            
            if ( enhancement_with_estimator_criterium ) {
           
                enhancement = 1;
               
                /// Construction du critere d'amelioration sur l'estimateur d'erreur
                /// ----------------------------------------------------------------
                
                construct_estimator_criterium( m, estimator_ratio, theta_elem, debug_criterium_enhancement );
               
            }
            
        }
        
        if ( enhancement ) {
            
            /// Application du critere d'amelioration
            /// -------------------------------------
            
            apply_criterium_enhancement( m, "EET", enhancement_with_estimator_criterium, enhancement_with_geometric_criterium, estimator_ratio, geometric_ratio, val_estimator_criterium, val_geometric_criterium, elem_flag_enh, face_flag_enh, elem_flag_bal, elem_list_enh, face_list_enh, elem_list_bal, debug_criterium_enhancement );
            
            geometric_ratio.free();
            estimator_ratio.free();
            
            /// Construction de la partie standard des densites d'effort
            /// --------------------------------------------------------
            
            construct_standard_force_fluxes_EET( m, f, pb, cost_function, enhancement, face_flag_enh, solver_minimisation, vec_force_fluxes, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation, verif_compatibility_conditions, tol_compatibility_conditions, debug_geometry, debug_force_fluxes, debug_method );
            
            /// Resolution des problemes locaux associes a la partie standard des densites d'effort avec procedure d'equilibrage
            /// ----------------------------------------------------------------------------------------------------------------
            
            if ( enhancement_with_estimator_criterium == 0 )
                construct_K_hat( m, f, "EET", K_hat, debug_method );
            
            balancing = 1;
            construct_F_hat( m, f, pb, balancing, elem_flag_bal, elem_flag_enh, vec_force_fluxes, F_hat, want_local_enrichment, debug_method, debug_geometry );
            
            construct_dep_hat( m, f, "EET", solver, K_hat, F_hat, dep_hat, verif_solver, tol_solver, debug_method );
            
            /// Construction de la partie amelioree des densites d'effort
            /// ---------------------------------------------------------
            
            construct_enhanced_force_fluxes_EET_EESPT( m, f, "EET", elem_flag_enh, face_flag_enh, elem_flag_bal, elem_list_enh, face_list_enh, elem_list_bal, K_hat, dep_hat, vec_force_fluxes, solver, solver_minimisation, verif_solver_enhancement, tol_solver_enhancement, verif_solver_minimisation_enhancement, tol_solver_minimisation_enhancement, debug_force_fluxes_enhancement, debug_method_enhancement, debug_geometry );
            
            /// Verification de l'equilibre des densites d'effort standard + ameliorees
            /// ------------------------------------------------------------------------
            if ( verif_eq_force_fluxes )
                check_equilibrium_force_fluxes( m, f, pb, vec_force_fluxes, tol_eq_force_fluxes, want_local_enrichment, debug_geometry );
            
            /// Resolution des problemes locaux associes aux densites d'effort standard + ameliorees sans procedure d'equilibrage
            /// -----------------------------------------------------------------------------------------------------------------
            
            balancing = 0;
            construct_F_hat( m, f, pb, balancing, elem_flag_bal, elem_flag_enh, vec_force_fluxes, F_hat, want_local_enrichment, debug_method_enhancement, debug_geometry );
            
            construct_dep_hat( m, f, "EET", solver, K_hat, F_hat, dep_hat, verif_solver_enhancement, tol_solver_enhancement, debug_method_enhancement );
            
            /// Construction d'un champ de contrainte admissible par element, Calcul d'un estimateur d'erreur globale
            /// -----------------------------------------------------------------------------------------------------
            
            calcul_error_estimate_prolongation_condition( m, f, pb, "EET", theta, theta_init, theta_elem, theta_elem_init, dep_hat, want_global_discretization_error, want_local_discretization_error, debug_error_estimate, debug_local_effectivity_index, debug_method, debug_method_enhancement );
            
        }
        
        t_EET.stop();
        cout << "Temps de calcul total pour la technique EET : " << t_EET.res << endl << endl;
        
    }

    /// ------------ ///
    /// Methode SPET ///
    /// ------------ ///
    
    if ( method.find("SPET") != string::npos ) {
        
        display_method( pb, "SPET", cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, solver, solver_minimisation );
        
        TicToc t_SPET;
        t_SPET.start();
        
        /// Calcul d'un estimateur d'erreur globale
        /// ---------------------------------------
        
        theta = 0.;
        calcul_error_estimate_partition_unity( m, f, pb, solver, "SPET", theta, theta_init, theta_elem, theta_elem_init, dep_hat, verif_solver, tol_solver, want_global_discretization_error, want_local_discretization_error, want_local_enrichment, debug_geometry, debug_error_estimate, debug_local_effectivity_index, debug_method );
        
        t_SPET.stop();
        cout << "Temps de calcul total pour la technique SPET : " << t_SPET.res << endl << endl;
        
    }
    
    /// ------------- ///
    /// Methode EESPT ///
    /// ------------- ///
    
    if ( method.find("EESPT") != string::npos ) {
        
        display_method( pb, "EESPT", cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, solver, solver_minimisation );
        
        TicToc t_EESPT;
        t_EESPT.start();
        
        Vec< Vec< Vec<T> > > vec_force_fluxes;
        Vec< Mat<T, Sym<> > > K_hat;
        Vec< Vec<T> > F_hat;
        
        theta = 0.;
        
        Vec<bool> elem_flag_enh;
        Vec<bool> face_flag_enh;
        Vec<bool> elem_flag_bal;
        Vec<unsigned> elem_list_enh;
        Vec<unsigned> face_list_enh;
        Vec<unsigned> elem_list_bal;
        
        bool enhancement = 0;
        bool balancing = 0;
        
        /// --------------------------------------------------------------------------- ///
        /// Critere d'amelioration geometrique de la construction des densites d'effort ///
        /// --------------------------------------------------------------------------- ///
        
        Vec<T> geometric_ratio;
        
        if ( enhancement_with_geometric_criterium ) {
            enhancement = 1;
            construct_geometric_criterium( m, geometric_criterium, geometric_ratio, debug_criterium_enhancement );
        }
        
        /// -------------------------------------------------------------------------------- ///
        /// Critere d'amelioration sur l'estimateur de la construction des densites d'effort ///
        /// -------------------------------------------------------------------------------- ///
        
        Vec<T> estimator_ratio;
        
        if ( enhancement == 0 or enhancement_with_estimator_criterium ) {
           
            enhancement = 0;
            
            /// Construction des densites d'effort standard
            /// -------------------------------------------
            
            Vec< Vec< Vec<T> > > vec_force_fluxes_standard;
            construct_standard_force_fluxes_EESPT( m, f, cost_function, enhancement, face_flag_enh, solver_minimisation, penalty_val_N, pb, vec_force_fluxes_standard, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation,  debug_geometry, debug_force_fluxes, debug_method );
            
            /// Verification de l'equilibre des densites d'effort standard
            /// ----------------------------------------------------------
            if ( verif_eq_force_fluxes )
                check_equilibrium_force_fluxes( m, f, pb, vec_force_fluxes_standard, tol_eq_force_fluxes, want_local_enrichment, debug_geometry );
            
            /// Resolution des problemes locaux associes aux densites d'effort standard
            /// -----------------------------------------------------------------------
            
            construct_K_hat( m, f, "EESPT", K_hat, debug_method );

            balancing = 0;
            construct_F_hat( m, f, pb, balancing, elem_flag_bal, elem_flag_enh, vec_force_fluxes_standard, F_hat, want_local_enrichment, debug_method, debug_geometry );
            
            construct_dep_hat( m, f, "EESPT", solver, K_hat, F_hat, dep_hat, verif_solver, tol_solver, debug_method );
            
            /// Construction d'un champ de contrainte admissible par element, Calcul d'un estimateur d'erreur globale
            /// -----------------------------------------------------------------------------------------------------
            
            calcul_error_estimate_prolongation_condition( m, f, pb, "EESPT", theta, theta_init, theta_elem, theta_elem_init, dep_hat, want_global_discretization_error, want_local_discretization_error, debug_error_estimate, debug_local_effectivity_index, debug_method, debug_method_enhancement );
            
            if ( enhancement_with_estimator_criterium ) {
                
                enhancement = 1;
               
                /// Construction du critere d'amelioration sur l'estimateur d'erreur globale
                /// ------------------------------------------------------------------------
                
                construct_estimator_criterium( m, estimator_ratio, theta_elem, debug_criterium_enhancement );
               
            }
            
        }
        
        if ( enhancement ) {
            
            /// Application du critere d'amelioration
            /// -------------------------------------
            
            apply_criterium_enhancement( m, "EESPT", enhancement_with_estimator_criterium, enhancement_with_geometric_criterium, estimator_ratio, geometric_ratio, val_estimator_criterium, val_geometric_criterium, elem_flag_enh, face_flag_enh, elem_flag_bal, elem_list_enh, face_list_enh, elem_list_bal, debug_criterium_enhancement );
            
            geometric_ratio.free();
            estimator_ratio.free();
            
            /// Construction de la partie standard des densites d'effort
            /// --------------------------------------------------------
            
            construct_standard_force_fluxes_EESPT( m, f, cost_function, enhancement, face_flag_enh, solver_minimisation, penalty_val_N, pb, vec_force_fluxes, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation, debug_geometry, debug_force_fluxes, debug_method );
            
            /// Resolution des problemes locaux associes a la partie standard des densites d'effort avec procedure d'equilibrage
            /// ----------------------------------------------------------------------------------------------------------------
            
            if ( enhancement_with_estimator_criterium == 0 )
                construct_K_hat( m, f, "EESPT", K_hat, debug_method );
            
            balancing = 1;
            construct_F_hat( m, f, pb, balancing, elem_flag_bal, elem_flag_enh, vec_force_fluxes, F_hat, want_local_enrichment, debug_method, debug_geometry );
            
            construct_dep_hat( m, f, "EESPT", solver, K_hat, F_hat, dep_hat, verif_solver, tol_solver, debug_method );
            
            /// Construction de la partie amelioree des densites d'effort
            /// ---------------------------------------------------------
            
            construct_enhanced_force_fluxes_EET_EESPT( m, f, "EESPT", elem_flag_enh, face_flag_enh, elem_flag_bal, elem_list_enh, face_list_enh, elem_list_bal, K_hat, dep_hat, vec_force_fluxes, solver, solver_minimisation, verif_solver_enhancement, tol_solver_enhancement, verif_solver_minimisation_enhancement, tol_solver_minimisation_enhancement, debug_force_fluxes_enhancement, debug_method_enhancement, debug_geometry );
            
            /// Verification de l'equilibre des densites d'effort standard + ameliorees
            /// ------------------------------------------------------------------------
            if ( verif_eq_force_fluxes )
                check_equilibrium_force_fluxes( m, f, pb, vec_force_fluxes, tol_eq_force_fluxes, want_local_enrichment, debug_geometry );
            
            /// Resolution des problemes locaux associes aux densites d'effort standard + ameliorees sans procedure d'equilibrage
            /// -----------------------------------------------------------------------------------------------------------------
            
            balancing = 0;
            construct_F_hat( m, f, pb, balancing, elem_flag_bal, elem_flag_enh, vec_force_fluxes, F_hat, want_local_enrichment, debug_method_enhancement, debug_geometry );
            
            construct_dep_hat( m, f, "EESPT", solver, K_hat, F_hat, dep_hat, verif_solver_enhancement, tol_solver_enhancement, debug_method_enhancement );
            
            /// Construction d'un champ de contrainte admissible par element, Calcul d'un estimateur d'erreur globale
            /// -----------------------------------------------------------------------------------------------------
            
            calcul_error_estimate_prolongation_condition( m, f, pb, "EESPT", theta, theta_init, theta_elem, theta_elem_init, dep_hat, want_global_discretization_error, want_local_discretization_error, debug_error_estimate, debug_local_effectivity_index, debug_method, debug_method_enhancement );
            
        }
        
        t_EESPT.stop();
        cout << "Temps de calcul total pour la technique EESPT : " << t_EESPT.res << endl << endl;
        
    }
}

#endif // Calcul_global_error_estimation_h
