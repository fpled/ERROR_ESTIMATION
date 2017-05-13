//
// C++ Interface: Calcul_error_estimate_partition_unity_PGD
//
// Description: calcul d'un champ admissible, calcul d'un estimateur theta de l'erreur globale pour la methode basee sur la partition de l'unite (SPET) dans le cadre des methodes PGD
//
//
// Author: Pled Florent <florent.pled@univ-paris-est.fr>, (C) 2017
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_error_estimate_partition_unity_PGD_h
#define Calcul_error_estimate_partition_unity_PGD_h

#include "SPET.h"
#include "../ECRE/ECRE.h"
#include "../DISCRETIZATION_ERROR/Discretization_error.h"

using namespace LMT;
using namespace std;

/// Calcul d'un champ de contrainte admissible & Calcul d'un estimateur theta de l'erreur globale pour la methode basee sur la partition de l'unite (SPET)
/// ------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class TV, class TVV, class TTVV, class TTVVV, class TTTVVV>
void calcul_error_estimate_partition_unity_PGD( TM &m, TF &f, const string &pb, const string &method, T &theta, T &theta_PGD, T &theta_dis, TV &theta_elem, TV &theta_elem_PGD, TV &theta_elem_dis, const TTVV &dep_space, const TTVVV &dep_param, const TV &dep_space_part, const TTTVVV &dep_space_FE, const TTTVVV &dep_space_hat, const TVV &dep_space_part_hat, const Vec< Vec<unsigned> > &elem_group, const unsigned &nb_modes, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool disp = false ) {
    
    /// ------------------------------------------------------------------------------------------------------ ///
    /// Construction d'un champ de contrainte admissible par element & Calcul d'un estimateur d'erreur globale ///
    /// ------------------------------------------------------------------------------------------------------ ///
    
    if ( disp ) {
        cout << "Construction d'un champ de contrainte admissible par element & Calcul d'un estimateur d'erreur globale" << endl;
        cout << "------------------------------------------------------------------------------------------------------" << endl << endl;
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
    
    f.vectors[0].set( 0. );
    for (unsigned n=0;n<nb_modes;++n) {
        TV dep_mode = dep_space[ n ];
        for (unsigned p=0;p<elem_group.size()-1;++p)
            dep_mode *= dep_param[ p ][ n ][ 0 ];
        f.vectors[0] +=  dep_mode;
    }
    
    TVV dep;
    dep.resize( elem_group.size() );
    for (unsigned g=0;g<elem_group.size();++g) {
        dep[ g ] = dep_space_part;
        for (unsigned n=0;n<nb_modes;++n) {
            TV dep_mode = dep_space_FE[ n ][ g ];
            for (unsigned p=0;p<elem_group.size()-1;++p)
                dep_mode *= dep_param[ p ][ n ][ 0 ];
            dep[ g ] += dep_mode;
        }
    }
    
    TVV dep_hat = dep_space_part_hat;
    for (unsigned n=0;n<nb_modes;++n) {
        TVV dep_hat_mode = dep_space_hat[ n ];
        for (unsigned p=0;p<elem_group.size()-1;++p)
            dep_hat_mode *= dep_param[ p ][ n ][ 0 ];
        dep_hat += dep_hat_mode;
    }
    
    Calc_Elem_Error_Estimate_FE<TV,TVV> calc_elem_error_estimate_FE;
    calc_elem_error_estimate_FE.dep = &dep;
    calc_elem_error_estimate_FE.elem_group = &elem_group;
    calc_elem_error_estimate_FE.method = &method;
    calc_elem_error_estimate_FE.theta_elem = &theta_elem_PGD;
    
    apply( m.elem_list, calc_elem_error_estimate_FE, m, f, theta_PGD );
    
    Calcul_Elem_Error_Estimate_SPET<TV,TVV> calcul_elem_error_estimate_SPET;
    calcul_elem_error_estimate_SPET.dep_hat = &dep_hat;
    calcul_elem_error_estimate_SPET.theta_elem = &theta_elem;
    
    apply( m.elem_list, calcul_elem_error_estimate_SPET, m, f, theta );
    
    theta_dis = theta - theta_PGD;
    for (unsigned n=0;n<m.elem_list.size();++n)
        theta_elem_dis[ n ] = theta_elem[ n ] - theta_elem_PGD[ n ];
    
    if ( disp ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "contribution a l'estimateur d'erreur globale au carre de l'element " << n << " :" << endl;
            cout << "theta_elem^2 = " << theta_elem[ n ] << endl;
            cout << "theta_elem_PGD^2 = " << theta_elem_PGD[ n ] << endl;
            cout << "theta_elem_dis^2 = " << theta_elem_dis[ n ] << endl << endl;
        }
        cout << endl;
    }
    
    theta = sqrt( theta );
    theta_PGD = sqrt( theta_PGD );
    theta_dis = sqrt( theta_dis );
    m.theta_SPET = theta;
    cout << "estimateur d'erreur globale :" << endl;
    cout << "theta = " << theta << endl;
    cout << "theta_PGD = " << theta_PGD << endl;
    cout << "theta_dis = " << theta_dis << endl;
    cout << "theta     / norm(u_h) = " << theta / m.norm_dep * 100. << " %" << endl;
    cout << "theta_PGD / norm(u_h) = " << theta_PGD / m.norm_dep * 100. << " %" << endl;
    cout << "theta_dis / norm(u_h) = " << theta_dis / m.norm_dep * 100. << " %" << endl << endl;
    
    if ( pb == "direct" and want_global_discretization_error ) {
        m.eff_index_SPET = theta / m.discretization_error;
        cout << "indice d'efficacite global :" << endl;
        cout << "eta = theta / e" << endl;
        cout << "    = " << m.eff_index_SPET << endl << endl;
    }
    
    if ( pb == "direct" and want_local_discretization_error ) {
        TV eff_index_elem;
        eff_index_elem.resize( m.elem_list.size() );
        eff_index_elem.set( 0. );
        
        apply( m.elem_list, Calcul_Elem_Effectivity_Index(), method, eff_index_elem );
        
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

#endif // Calcul_error_estimate_partition_unity_PGD_h
