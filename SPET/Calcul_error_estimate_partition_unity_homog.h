//
// C++ Interface: Calcul_error_estimate_partition_unity_homog
//
// Description: calcul d'un champ admissible, calcul d'un estimateur theta de l'erreur globale pour la methode basee sur la partition de l'unite (SPET)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_error_estimate_partition_unity_homog_h
#define Calcul_error_estimate_partition_unity_homog_h

#include "SPET.h"
#include "../DISCRETIZATION_ERROR/Discretization_error.h"

using namespace LMT;
using namespace std;

/// Calcul d'un champ de contrainte admissible & Calcul d'un estimateur theta de l'erreur globale pour la methode basee sur la partition de l'unite (SPET)
/// ------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class TV, class TVV>
void calcul_error_estimate_partition_unity( TM &m, const TF &f, const string &pb, const string &method, T &theta, T &theta_init, TV &theta_elem, TV &theta_elem_init, TVV &dep_hat, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool disp = false ) {
    
    /// ------------------------------------------------------------------------------------------------------ ///
    /// Construction d'un champ de contrainte admissible par element & Calcul d'un estimateur d'erreur globale ///
    /// ------------------------------------------------------------------------------------------------------ ///
    
    if ( disp ) {
        cout << "Construction d'un champ de contrainte admissible par element & Calcul d'un estimateur d'erreur globale" << endl;
        cout << "------------------------------------------------------------------------------------------------------" << endl << endl;
    }
    
    theta = 0.;
    
    theta_elem.resize( m.elem_list.size() );
    theta_elem.set( 0. );
    
    Calcul_Elem_Error_Estimate_Init_SPET<T,TV,TVV> calcul_elem_error_estimate_init_SPET;
    calcul_elem_error_estimate_init_SPET.dep_hat = &dep_hat;
    calcul_elem_error_estimate_init_SPET.theta_elem = &theta_elem;
    calcul_elem_error_estimate_init_SPET.theta_elem_init = &theta_elem_init;
    calcul_elem_error_estimate_init_SPET.theta_init = &theta_init;
    
    apply( m.elem_list, calcul_elem_error_estimate_init_SPET, m, f, theta );
    
    if ( disp ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "contribution a l'estimateur d'erreur globale au carre de l'element " << n << " :" << endl;
            cout << "theta_elem^2 = " << theta_elem[ n ] << endl;
            cout << "theta_elem_init^2 = " << theta_elem_init[ n ] << endl;
        }
        cout << endl;
    }
    
    theta = sqrt( theta );
    m.theta_SPET = theta;
    cout << "estimateur d'erreur globale :" << endl;
    cout << "theta = " << theta << endl;
    cout << "theta_init = " << theta_init << endl;
    cout << "theta / norm(u_h) = " << theta / m.norm_dep * 100. << " %" << endl;
    cout << "theta_init / norm(u_h_init) = " << theta_init / m.norm_dep_init * 100. << " %" << endl << endl;
    
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

#endif // Calcul_error_estimate_partition_unity_homog_h
