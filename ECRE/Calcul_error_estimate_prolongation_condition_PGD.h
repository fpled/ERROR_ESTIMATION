//
// C++ Interface: Calcul_error_estimate_prolongation_condition
//
// Description: calcul d'un champ de contrainte admissible et d'un estimateur theta de l'erreur globale pour les methodes basees sur la condition de prolongement (EET,EESPT)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_error_estimate_prolongation_condition_h
#define Calcul_error_estimate_prolongation_condition_h

#include "ECRE.h"

using namespace LMT;
using namespace std;
/*
/// Construction d'un champ de contrainte admissible et Calcul d'un estimateur d'erreur globale pour les methodes basees sur la condition de prolongement (EET,EESPT)
/// -----------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T>
void calcul_error_estimate_prolongation_condition( TM &m, const TF &f, const string &pb, const string &method, T &theta, Vec<T> &theta_elem, const Vec< Vec<T> > &dep_hat, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool debug_error_estimate = false, const bool debug_local_effectivity_index = false, const bool debug_method = false, const bool debug_method_enhancement = false ) {
    
    const Vec< Vec<T> > &dep_psi, const Vec< Vec<T> > &dep_lambda, const unsigned &nb_modes, const TMAT &K_s, const TMAT &K_unk_s, const TMAT &K_k_s, const TMAT &K_unk_p, const TMAT &K_k_p, const Vec<T> &F_s, const Vec<T> &F_p, const Vec<unsigned> &elem_list_PGD_unknown_param,
	
	m, f, pb, "EET", theta, theta_elem, dep_part_hat, dep_psi_hat, debug_method, debug_method_enhancement, debug_ecre_theta, debug_local_effectivity_index, want_global_discretization_error, want_local_discretization_error

    /// ------------------------------------------------------------------------------------------------------- ///
    /// Construction d'un champ de contrainte admissible par element et Calcul d'un estimateur d'erreur globale ///
    /// ------------------------------------------------------------------------------------------------------- ///
    
    cout << "------------------------------------------------------------" << endl;
    cout << "Construction d'un champ de contrainte admissible par element" << endl;
    cout << "Calcul d'un estimateur d'erreur globale" << endl;
    cout << "------------------------------------------------------------" << endl << endl;

    theta = 0.;
    theta_elem.resize( m.elem_list.size() );
    theta_elem.set( 0. );

    Calc_Elem_Error_Estimate_EET_EESPT<T,TMAT> calc_elem_error_estimate_EET_EESPT;
    calc_elem_error_estimate_EET_EESPT.dep_hat = &dep_hat;
    calc_elem_error_estimate_EET_EESPT.method = &method;
    calc_elem_error_estimate_EET_EESPT.theta_elem = &theta_elem;
    calc_elem_error_estimate_EET_EESPT.dep_psi = &dep_psi;
    calc_elem_error_estimate_EET_EESPT.dep_lambda = &dep_lambda;
    calc_elem_error_estimate_EET_EESPT.dep_part = &dep_part;
    calc_elem_error_estimate_EET_EESPT.kappa = &kappa;
    calc_elem_error_estimate_EET_EESPT.K_k_p = &K_k_p;
    calc_elem_error_estimate_EET_EESPT.K_unk_p = K_unk_p;
    calc_elem_error_estimate_EET_EESPT.elem_list_PGD_unknown_param = &elem_list_PGD_unknown_param;
    calc_elem_error_estimate_EET_EESPT.nb_modes = &nb_modes;

    apply( m.elem_list, calc_elem_error_estimate_EET_EESPT, m, f, theta );

    if ( debug_error_estimate or debug_method or debug_method_enhancement ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            T ecre_elem = theta_elem[ n ] / 2;
            cout << "contribution a la mesure globale de l'erreur en relation de comportement au carre de l'element " << n << " :" << endl;
            cout << "ecre_elem^2 = " << ecre_elem << endl;
            cout << "contribution a l'estimateur d'erreur globale au carre de l'element " << n << " :" << endl;
            cout << "theta_elem^2 = " << theta_elem[ n ] << endl;
        }
    }

    T ecre = theta / 2.;
    cout << "mesure globale de l'erreur en relation de comportement au carre :" << endl;
    cout << "ecre^2 = " << ecre << endl;

    theta = sqrt( theta );
    ecre = sqrt( ecre );
    if ( method == "EET" ) {
        m.ecre_EET = ecre;
        m.theta_EET = theta;
    }
    if ( method == "EESPT" ) {
        m.ecre_EESPT = ecre;
        m.theta_EESPT = theta;
    }
    cout << "mesure globale de l'erreur en relation de comportement :" << endl;
    cout << "ecre = " << ecre << endl;
    cout << "estimateur d'erreur globale :" << endl;
    cout << "theta = " << theta << endl << endl;

    if ( pb == "direct" and want_global_discretization_error ) {
        T eff_index = theta / m.discretization_error;
        if ( method == "EET" )
            m.eff_index_EET = eff_index;
        if ( method == "EESPT" )
            m.eff_index_EESPT = eff_index;
        cout << "indice d'efficacite global :" << endl;
        cout << "eta = " << eff_index << endl << endl;
    }

    if ( pb == "direct" and want_local_discretization_error ) {
        Vec<T> eff_index_elem;
        eff_index_elem.resize( m.elem_list.size(), 0. );
        eff_index_elem.set( 0. );

        apply( m.elem_list, Calcul_Elem_Effectivity_Index(), method, eff_index_elem );

        if ( debug_local_effectivity_index or debug_method or debug_method_enhancement ) {
            for (unsigned n=0;n<m.elem_list.size();++n) {
                cout << "indice d'efficacite local de l'element " << n << " :" << endl;
                cout << "eta_elem = theta_elem / e_elem" << endl;
                cout << "         = " << eff_index_elem[ n ] << endl;
            }
        }
    }
}

#endif // Calcul_error_estimate_prolongation_condition_h
