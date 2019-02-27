//
// C++ Interface: Calcul_discretization_error
//
// Description: calcul de la mesure de l'erreur de discretisation globale et locale
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_discretization_error_h
#define Calcul_discretization_error_h

#include "Discretization_error.h"
#include "../LMT/include/containers/apply_ij.h"

using namespace LMT;
using namespace std;

/// Calcul de la norme du champ de deplacement
/// ------------------------------------------
template<class TM, class TF>
void calcul_norm_dep( TM &m, const TF &f, const string &pb, const bool disp = false ) {
    
    typedef typename TM::Pvec Pvec;
    typedef typename TM::TNode::T T;
    typedef Vec<T> TV;
    
    if ( disp ) {
        if ( pb == "direct" ) {
            cout << "Calcul de la norme du champ de deplacement approche" << endl << endl;
        }
        else if ( pb == "reference" ) {
            cout << "Calcul de la norme du champ de deplacement de reference" << endl << endl;
        }
    }
    
    T norm_dep = 0.;
    TV norm_dep_elem;
    norm_dep_elem.resize( m.elem_list.size() );
    norm_dep_elem.set( 0. );
    
    Add_Elem_Norm_Dep<T,TV> add_elem_norm_dep;
    add_elem_norm_dep.norm_dep_elem = &norm_dep_elem;
    
    apply( m.elem_list, add_elem_norm_dep, m, f, norm_dep );
    
    m.norm_dep = sqrt( norm_dep );
    
    if ( pb == "direct" ) {
//        cout << "norme du champ de deplacement approche au carre :" << endl;
//        cout << "norm(u_h)^2 = " << norm_dep << endl << endl;
        
        cout << "norme du champ de deplacement approche :" << endl;
        cout << "norm(u_h) = " << m.norm_dep << endl << endl;
    }
    else if ( pb == "reference" ) {
//        cout << "norme du champ de deplacement de reference au carre :" << endl;
//        cout << "norm(u_ex)^2 = " << norm_dep << endl << endl;
        
        cout << "norme du champ de deplacement de reference :" << endl;
        cout << "norm(u_ex) = " << m.norm_dep << endl << endl;
    }
}

/// Calcul de la norme du champ de deplacement initial
/// --------------------------------------------------
template<class TM, class TF>
void calcul_norm_dep_init( TM &m, const TF &f, const string &pb, const bool disp = false ) {
    
    typedef typename TM::Pvec Pvec;
    typedef typename TM::TNode::T T;
    typedef Vec<T> TV;
    
    if ( disp ) {
        if ( pb == "direct" )
            cout << "Calcul de la norme du champ de deplacement initial approche" << endl << endl;
        else if ( pb == "reference" )
            cout << "Calcul de la norme du champ de deplacement initial de reference" << endl << endl;
    }
    
    T norm_dep_init = 0.;
    TV norm_dep_elem_init;
    norm_dep_elem_init.resize( m.elem_list.size() );
    norm_dep_elem_init.set( 0. );
    
    Add_Elem_Norm_Dep_Init<T,TV> add_elem_norm_dep_init;
    add_elem_norm_dep_init.norm_dep_elem_init = &norm_dep_elem_init;
    
    apply( m.elem_list, add_elem_norm_dep_init, m, f, norm_dep_init );
    
    m.norm_dep_init = sqrt( norm_dep_init );
    
    if ( pb == "direct" ) {
//        cout << "norme du champ de deplacement approche initial au carre :" << endl;
//        cout << "norm(u_h)_init^2 = " << norm_dep_init << endl << endl;
        
        cout << "norme du champ de deplacement initial approche :" << endl;
        cout << "norm(u_h)_init = " << m.norm_dep_init << endl << endl;
    }
    else if ( pb == "reference" ) {
//        cout << "norme du champ de deplacement de reference initial au carre :" << endl;
//        cout << "norm(u_ex)_init^2 = " << norm_dep_init << endl << endl;
        
        cout << "norme du champ de deplacement de reference initial :" << endl;
        cout << "norm(u_ex)_init = " << m.norm_dep_init << endl << endl;
    }
}

/// Calcul de la mesure de l'erreur de discretisation globale et locale
/// -------------------------------------------------------------------
template<class TM, class TF>
void calcul_discretization_error( TM &m, TM &m_ref, const TF &f, const TF &f_ref, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool want_solve_ref = false, const bool disp = false ) {
    
    typedef typename TM::Pvec Pvec;
    typedef typename TM::TNode::T T;
    typedef Vec<T> TV;
    
    T discretization_error = 0.;
    TV discretization_error_elem;
    discretization_error_elem.resize( m.elem_list.size() );
    discretization_error_elem.set( 0. );
    
    if ( want_global_discretization_error or want_local_discretization_error ) {
        
        cout << "Mesure de l'erreur de discretisation" << endl;
        cout << "------------------------------------" << endl << endl;
        
        if ( want_solve_ref ) {
            
            bool homogeneous_constraints = 1;
            
            for (unsigned nc=0;nc<f.constraints.size();++nc) {
                if ( not ( f.constraints[nc].res == 0 ) ) {
                    homogeneous_constraints = 0;
                    cerr << "Conditions aux limites de Dirichlet non homogènes..." << endl << endl;
                    break;
                }
            }
            
            if ( want_global_discretization_error and homogeneous_constraints ) {
                
                TicToc t;
                t.start();
                
                /// Calcul de la norme du champ de deplacement de reference
                /// -------------------------------------------------------
                
                calcul_norm_dep( m_ref, f_ref, "reference" );
                
                /// Calcul exact de la mesure de l'erreur de discretisation globale
                /// ---------------------------------------------------------------
                
                if ( disp )
                    cout << "Calcul exact de la mesure de l'erreur de discretisation globale" << endl << endl;
                
                T discretization_error = pow( m_ref.norm_dep, 2 ) - pow( m.norm_dep, 2 );
//                cout << "mesure de l'erreur de discretisation globale au carre :" << endl;
//                cout << "e^2 = norm(u_ex - u_h)^2 = norm(u_ex)^2 - norm(u_h)^2 = " << discretization_error << endl << endl;
                
                discretization_error = sqrt( discretization_error );
                m.discretization_error = discretization_error;
                cout << "mesure de l'erreur de discretisation globale :" << endl;
                cout << "e = norm(u_ex - u_h) = sqrt( norm(u_ex)^2 - norm(u_h)^2 ) = " << m.discretization_error << endl;
                cout << "e / norm(u_h) = " << m.discretization_error / m.norm_dep * 100. << " %" << endl << endl;
                
                t.stop();
                cout << "temps de calcul de la mesure de l'erreur de discretisation globale = " << t.res << " s" << endl << endl;
            }
            
            if ( want_local_discretization_error ) {
                
                /// Calcul approche des mesures (elementaires) de l'erreur de discretisation locale
                /// -------------------------------------------------------------------------------
                
                TicToc t;
                t.start();
                
                if ( disp )
                    cout << "Calcul approche des mesures (elementaires) de l'erreur de discretisation locale" << endl << endl;
                
                Calcul_Elem_Discretization_Error<TM, TF, TV> calcul_elem_discretization_error;
                calcul_elem_discretization_error.m = &m;
                calcul_elem_discretization_error.m_ref = &m_ref;
                calcul_elem_discretization_error.f = &f;
                calcul_elem_discretization_error.f_ref = &f_ref;
                calcul_elem_discretization_error.discretization_error_elem = &discretization_error_elem;
                
                apply_ij( m.elem_list, m_ref.elem_list, calcul_elem_discretization_error );
                
                apply( m.elem_list, Set_Elem_Discretization_Error(), discretization_error_elem );
                
//                for (unsigned n=0;n<m.elem_list.size();++n) {
//                    cout << "mesure de l'erreur de discretisation locale au carre sur l'element " << n << " :" << endl;
//                    cout << "e_elem^2 = norm(u_ex - u_h)_elem^2 = " << discretization_error_elem[ n ] << endl;
//                }
                
                T sum_discretization_error_elem = 0.;
                for (unsigned n=0;n<m.elem_list.size();++n)
                    sum_discretization_error_elem += discretization_error_elem[ n ];
                    
//                cout << "indicateur d'erreur globale au carre :" << endl;
//                cout << "e^2 = norm(u_ex - u_h)^2 = " << sum_discretization_error_elem << endl << endl;
                
                cout << "indicateur d'erreur globale :" << endl;
                cout << "e = norm(u_ex - u_h) ≈ " << sqrt( sum_discretization_error_elem ) << endl;
                cout << "e / norm(u_h) ≈ " << sqrt( sum_discretization_error_elem ) / m.norm_dep * 100. << " %" << endl << endl;
                
                t.stop();
                cout << "temps de calcul des mesures (elementaires) de l'erreur de discretisation locale = " << t.res << " s" << endl << endl;
                
                /// Calcul approche de la mesure de l'erreur de discretisation globale comme somme de contributions locales
                /// -------------------------------------------------------------------------------------------------------
                
                if ( want_global_discretization_error and not( homogeneous_constraints ) )
                    m.discretization_error = sqrt( sum_discretization_error_elem );
            }
        }
        else {
            discretization_error = m.discretization_error;
            
//            cout << "mesure de l'erreur de discretisation globale au carre :" << endl;
//            cout << "e^2 = norm(u_ex - u_h)^2 = " << pow( m.discretization_error, 2 ) << endl << endl;
            
            cout << "mesure de l'erreur de discretisation globale :" << endl;
            cout << "e = norm(u_ex - u_h) = " << m.discretization_error << endl;
            cout << "e / norm(u_h) = " << m.discretization_error / m.norm_dep * 100. << " %" << endl << endl;
            
            T sum_discretization_error_elem = 0.;
            for (unsigned n=0;n<m.elem_list.size();++n) {
                discretization_error_elem[ n ] = m.elem_list[ n ]->get_field( "discretization_error_elem", StructForType<T>() );
                sum_discretization_error_elem += discretization_error_elem[ n ];
            }
//            cout << "indicateur d'erreur globale au carre :" << endl;
//            cout << "e^2 = norm(u_ex - u_h)^2 ≈ " << sum_discretization_error_elem << endl << endl;
            
            cout << "indicateur d'erreur globale :" << endl;
            cout << "e = norm(u_ex - u_h) ≈ " << sqrt( sum_discretization_error_elem ) << endl;
            cout << "e / norm(u_h) ≈ " << sqrt( sum_discretization_error_elem ) / m.norm_dep * 100. << " %" << endl << endl;
        }
    }
}

#endif // Calcul_discretization_error_h
