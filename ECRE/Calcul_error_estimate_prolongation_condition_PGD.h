//
// C++ Interface: Calcul_error_estimate_prolongation_condition_PGD
//
// Description: calcul d'un champ de contrainte admissible et d'un estimateur theta de l'erreur globale pour les methodes basees sur la condition de prolongement (EET,EESPT) dans le cadre des methodes PGD
//
//
// Author: Pled Florent <florent.pled@univ-paris-est.fr>, (C) 2017
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_error_estimate_prolongation_condition_PGD_h
#define Calcul_error_estimate_prolongation_condition_PGD_h

#include "ECRE.h"
#include "../DISCRETIZATION_ERROR/Discretization_error.h"

using namespace LMT;
using namespace std;

/// Construction d'un champ de contrainte admissible & Calcul d'un estimateur d'erreur globale pour les methodes basees sur la condition de prolongement (EET,EESPT)
/// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class TV, class TVV, class TTVV, class TTVVV, class TTTVVV>
void calcul_error_estimate_prolongation_condition_PGD( TM &m, TF &f, const string &pb, T &theta, T &theta_PGD, T &theta_dis, TV &theta_elem, TV &theta_elem_PGD, TV &theta_elem_dis, const TTVV &dep_space, const TTVVV &dep_param, const TV &dep_space_FE_part, const TTTVVV &dep_space_FE, const TVV &dep_space_hat_part, const TTTVVV &dep_space_hat, const Vec< Vec<unsigned> > &elem_group, const unsigned &mode, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool disp = false ) {
    
    /// ------------------------------------------------------------------------------------------------------ ///
    /// Construction d'un champ de contrainte admissible par element & Calcul d'un estimateur d'erreur globale ///
    /// ------------------------------------------------------------------------------------------------------ ///
    
    if ( disp ) {
        cout << "Construction d'un champ de contrainte admissible & Calcul d'un estimateur d'erreur globale" << endl;
        cout << "------------------------------------------------------------------------------------------" << endl << endl;
    }
    
    theta = 0.;
    theta_PGD = 0.;
    theta_dis = 0.;
    
    theta_elem.resize( m.elem_list.size() );
    theta_elem_PGD.resize( m.elem_list.size() );
    theta_elem_dis.resize( m.elem_list.size() );
    theta_elem.set( 0. );
    theta_elem_PGD.set( 0. );
    theta_elem_dis.set( 0. );
    
    // PGD approximation
    f.vectors[0].set( 0. );
    for (unsigned n=0;n<mode+1;++n) {
        TV dep_mode = dep_space[ n ];
        for (unsigned p=0;p<elem_group.size()-1;++p)
            dep_mode *= dep_param[ p ][ n ][ 0 ];
        f.vectors[0] +=  dep_mode;
    }
    
    // FE approximation
    TVV dep_FE;
    dep_FE.resize( elem_group.size() );
    for (unsigned g=0;g<elem_group.size();++g) {
        dep_FE[ g ] = dep_space_FE_part;
        for (unsigned n=0;n<mode+1;++n) {
            TV dep_FE_mode = dep_space_FE[ n ][ g ];
            for (unsigned p=0;p<elem_group.size()-1;++p)
                dep_FE_mode *= dep_param[ p ][ n ][ 0 ];
            dep_FE[ g ] += dep_FE_mode;
        }
    }
    
    Calc_Elem_Error_Estimate_PGD<TV,TVV> calc_elem_error_estimate_PGD;
    calc_elem_error_estimate_PGD.dep = &dep_FE;
    calc_elem_error_estimate_PGD.elem_group = &elem_group;
    calc_elem_error_estimate_PGD.theta_elem_PGD = &theta_elem_PGD;
    
    apply( m.elem_list, calc_elem_error_estimate_PGD, m, f, theta_PGD );
    
    // Admissible approximation
    TVV dep_hat = dep_space_hat_part;
    for (unsigned n=0;n<mode+1;++n) {
        TVV dep_hat_mode = dep_space_hat[ n ];
        for (unsigned p=0;p<elem_group.size()-1;++p)
            dep_hat_mode *= dep_param[ p ][ n ][ 0 ];
        dep_hat += dep_hat_mode;
    }
    
    Calc_Elem_Error_Estimate<TV,TVV> calc_elem_error_estimate;
    calc_elem_error_estimate.dep_hat = &dep_hat;
    calc_elem_error_estimate.theta_elem = &theta_elem;
    
    apply( m.elem_list, calc_elem_error_estimate, m, f, theta );
    
    apply( m.elem_list, Calc_Elem_Error_Estimate_Dis(), theta_elem_dis );
    theta_dis = theta - theta_PGD;
        
    if ( disp ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            T ecre_elem = theta_elem[ n ] / 2.;
            T ecre_elem_PGD = theta_elem_PGD[ n ] / 2.;
            T ecre_elem_dis = theta_elem_dis[ n ] / 2.;
//            cout << "contribution a la mesure globale de l'erreur en relation de comportement au carre de l'element " << n << " :" << endl;
//            cout << "ecre_elem^2 = " << ecre_elem << endl;
//            cout << "ecre_elem_PGD^2 = " << ecre_elem << endl;
//            cout << "ecre_elem_dis^2 = " << ecre_elem << endl << endl;
            
            cout << "contribution a l'estimateur d'erreur globale au carre de l'element " << n << " :" << endl;
            cout << "theta_elem^2 = " << theta_elem[ n ] << endl;
            cout << "theta_elem_PGD^2 = " << theta_elem_PGD[ n ] << endl;
            cout << "theta_elem_dis^2 = " << theta_elem_dis[ n ] << endl << endl;
        }
    }
    
    T ecre = theta / 2.;
    T ecre_PGD = theta_PGD / 2.;
    T ecre_dis = theta_dis / 2.;
//    cout << "mesure globale de l'erreur en relation de comportement au carre :" << endl;
//    cout << "ecre^2 = " << ecre << endl;
//    cout << "ecre_PGD^2 = " << ecre_PGD << endl;
//    cout << "ecre_dis^2 = " << ecre_dis << endl << endl;
    
    cout << "estimateur d'erreur globale au carre :" << endl;
    cout << "theta^2     = " << theta << endl;
    cout << "theta_PGD^2 = " << theta_PGD << endl;
    cout << "theta_dis^2 = " << theta_dis << endl << endl;
    
    ecre = sqrt( ecre );
    ecre_PGD = sqrt( ecre_PGD );
    ecre_dis = sqrt( ecre_dis );
    theta = sqrt( theta );
    theta_PGD = sqrt( theta_PGD );
    theta_dis = sqrt( theta_dis );
    
    m.ecre = ecre;
    m.ecre_PGD = ecre_PGD;
    m.ecre_dis = ecre_dis;
    m.error_estimate = theta;
    m.error_estimate_PGD = theta_PGD;
    m.error_estimate_dis = theta_dis;
    
//    cout << "mesure globale de l'erreur en relation de comportement :" << endl;
//    cout << "ecre     = " << ecre << endl;
//    cout << "ecre_PGD = " << ecre_PGD << endl;
//    cout << "ecre_dis = " << ecre_dis << endl << endl;
    
    cout << "estimateur d'erreur globale :" << endl;
    cout << "theta     = " << theta << endl;
    cout << "theta_PGD = " << theta_PGD << endl;
    cout << "theta_dis = " << theta_dis << endl;
    cout << "theta     / norm(u_h) = " << theta / m.norm_dep * 100. << " %" << endl;
    cout << "theta_PGD / norm(u_h) = " << theta_PGD / m.norm_dep * 100. << " %" << endl;
    cout << "theta_dis / norm(u_h) = " << theta_dis / m.norm_dep * 100. << " %" << endl << endl;
    
    if ( pb == "direct" and want_global_discretization_error ) {
        T eff_index = theta / m.discretization_error;
        m.eff_index = eff_index;
        cout << "indice d'efficacite global :" << endl;
        cout << "eta = theta / e" << endl;
        cout << "    = " << eff_index << endl << endl;
    }
    
    if ( pb == "direct" and want_local_discretization_error ) {
        TV eff_index_elem;
        eff_index_elem.resize( m.elem_list.size(), 0. );
        eff_index_elem.set( 0. );
        
        apply( m.elem_list, Calcul_Elem_Effectivity_Index(), eff_index_elem );
        
        if ( disp ) {
            for (unsigned n=0;n<m.elem_list.size();++n) {
                cout << "indice d'efficacite local de l'element " << n << " :" << endl;
                cout << "eta_elem = theta_elem / e_elem" << endl;
                cout << "         = " << eff_index_elem[ n ] << endl;
            }
            cout << endl;
        }
    }
}

#endif // Calcul_error_estimate_prolongation_condition_PGD_h
