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

/// Calcul de la mesure de l'erreur de discretisation globale et locale
/// -------------------------------------------------------------------
template<class TM, class TF>
void calcul_discretization_error( TM &m, const TM &m_ref, const TF &f, const TF &f_ref, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool want_solve_ref = false, const bool debug_discretization_error = false ) {
    
    typedef typename TM::Pvec Pvec;
    typedef typename TM::TNode::T T;
    
    T discretization_error = 0.;
    Vec<T> discretization_error_elem;
    discretization_error_elem.resize( m.elem_list.size() );
    discretization_error_elem.set( 0. );
    
    if ( want_global_discretization_error or want_local_discretization_error ) {

        cout << "------------------------------------" << endl;
        cout << "Mesure de l'erreur de discretisation" << endl;
        cout << "------------------------------------" << endl << endl;

        if ( want_solve_ref ) {
            
            bool all_Dirichlet_boundary_conditions_equal_to_zero = 1;

            for (unsigned nc=0;nc<f.constraints.size();++nc) {
                if ( f.constraints[nc].res == 0 ) {
                }
                else {
                    all_Dirichlet_boundary_conditions_equal_to_zero = 0;
                    cerr << "Presence de conditions aux limites de type Dirichlet non imposees a 0..." << endl << endl;
                    break;
                }
            }

            if ( want_global_discretization_error and all_Dirichlet_boundary_conditions_equal_to_zero ) {

                TicToc t_global_discretization_error;
                t_global_discretization_error.start();

                /// Calcul de la norme du champ de deplacement approche
                /// ---------------------------------------------------

                cout << "Calcul de la norme du champ de deplacement approche" << endl << endl;

                T norm_dep = 0.;

                apply( m.elem_list, Add_Elem_Norm_Dep(), m, f, norm_dep );

                cout << "norme du champ de deplacement approche au carre :" << endl;
                cout << "||u_h||^2 = " << norm_dep << endl << endl;

                /// Calcul de la norme du champ de deplacement de reference
                /// -------------------------------------------------------

                cout << "Calcul de la norme du champ de deplacement de reference" << endl << endl;

                T norm_dep_ref = 0.;

                apply( m_ref.elem_list, Add_Elem_Norm_Dep(), m_ref, f_ref, norm_dep_ref );

                cout << "norme du champ de deplacement de reference au carre :" << endl;
                cout << "||u_ex||^2 = " << norm_dep_ref << endl << endl;

                /// Calcul exact de la mesure de l'erreur de discretisation globale
                /// ---------------------------------------------------------------

                cout << "Calcul exact de la mesure de l'erreur de discretisation globale" << endl << endl;

                T discretization_error = norm_dep_ref - norm_dep;
                m.discretization_error = sqrt(discretization_error);

                cout << "mesure de l'erreur de discretisation globale au carre :" << endl;
                cout << "e^2 = " << discretization_error << endl << endl;

                cout << "mesure de l'erreur de discretisation globale :" << endl;
                cout << "e = " << m.discretization_error << endl << endl;

                t_global_discretization_error.stop();
                cout << "Temps de calcul de la mesure de l'erreur de discretisation globale : " << t_global_discretization_error.res << endl << endl;
            }

            if ( want_local_discretization_error ) {

                /// Calcul approche des mesures (elementaires) de l'erreur de discretisation locale
                /// -------------------------------------------------------------------------------

                TicToc t_local_discretization_error;
                t_local_discretization_error.start();

                cout << "Calcul approche des mesures (elementaires) de l'erreur de discretisation locale" << endl << endl;

                Calcul_Elem_Discretization_Error<TM, TF, T> calcul_elem_discretization_error;
                calcul_elem_discretization_error.m = &m;
                calcul_elem_discretization_error.m_ref = &m_ref;
                calcul_elem_discretization_error.f = &f;
                calcul_elem_discretization_error.f_ref = &f_ref;
                calcul_elem_discretization_error.discretization_error_elem = &discretization_error_elem;

                apply_ij( m.elem_list, m_ref.elem_list, calcul_elem_discretization_error );

                apply( m.elem_list, Set_Elem_Discretization_Error(), discretization_error_elem );

                if ( debug_discretization_error ) {
                    for (unsigned n=0;n<m.elem_list.size();++n) {
                        cout << "mesure de l'erreur de discretisation locale au carre sur l'element " << n << " :" << endl;
                        cout << "e_elem^2 = " << discretization_error_elem[ n ] << endl;
                    }
                }

                T sum_discretization_error_elem = 0.;
                for (unsigned n=0;n<m.elem_list.size();++n) {
                    sum_discretization_error_elem += discretization_error_elem[ n ];
                }
                cout << "somme des mesures (elementaires) de l'erreur de discretisation locale au carre :" << endl;
                cout << "sum( e_elem^2 ) = " << sum_discretization_error_elem << endl << endl;

                t_local_discretization_error.stop();
                cout << "Temps de calcul des mesures (elementaires) de l'erreur de discretisation locale : " << t_local_discretization_error.res << endl << endl;

                /// Calcul approche de la mesure de l'erreur de discretisation globale comme somme de contributions locales
                /// -------------------------------------------------------------------------------------------------------

                if ( want_global_discretization_error and all_Dirichlet_boundary_conditions_equal_to_zero == 0 ) {

                    cout << "Calcul approche de la mesure de l'erreur de discretisation globale" << endl << endl;

                    cout << "mesure de l'erreur de discretisation globale au carre :" << endl;
                    cout << "e^2 = " << sum_discretization_error_elem << endl << endl;

                    m.discretization_error = sqrt( sum_discretization_error_elem );

                    cout << "mesure de l'erreur de discretisation globale :" << endl;
                    cout << "e = " << discretization_error << endl << endl;
                }
            }
        }
        else {
            discretization_error = m.discretization_error;

            cout << "mesure de l'erreur de discretisation globale au carre :" << endl;
            cout << "e^2 = " << pow( discretization_error, 2 ) << endl << endl;

            cout << "mesure de l'erreur de discretisation globale :" << endl;
            cout << "e = " << discretization_error << endl << endl;
            
            T sum_discretization_error_elem = 0.;
            for (unsigned n=0;n<m.elem_list.size();++n) {
                discretization_error_elem[ n ] = m.elem_list[ n ]->get_field( "discretization_error_elem", StructForType<T>() );
                sum_discretization_error_elem += discretization_error_elem[ n ];
            }
            cout << "somme des mesures (elementaires) de l'erreur de discretisation locale au carre :" << endl;
            cout << "sum( e_elem^2 ) = " << sum_discretization_error_elem << endl << endl;
        }
    }
}

#endif // Calcul_discretization_error_h
