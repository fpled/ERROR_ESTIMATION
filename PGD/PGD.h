//
// C++ Interface: PGD
//
// Description: fonctions relatives a la decomposition PGD
//
//
// Author: Pled Florent These 2011 <pled@lmt.ens-cachan.fr>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef PGD_h
#define PGD_h

using namespace LMT;
using namespace std;


/// Creation des proprietes materiaux
///----------------------------------
template<class TF, class TM, class TV>
void define_unknown_parameter_zone( TF &f, TM &m, const string &structure, TV &list_elems_PGD_unknown_parameter ) {

    static const unsigned dim = TM::dim;
    typedef typename TM::TNode::T T;

    if ( m.node_list.size() ) {
        /// Dimension 2
        ///------------
        if ( dim == 2 ) {
            /// Plaque rectangulaire 2D en flexion
            ///-----------------------------------
            if ( structure == "plate_flexion" ) {
                for (unsigned i=0;i<m.elem_list.size();++i) {
                    if ( center( *m.elem_list[i] )[1] < 0.5 )
                        list_elems_PGD_unknown_parameter.push_back( m.elem_list[i]->number );
                }
            }
            /// Inclusions circulaires 2D
            ///--------------------------
            if (structure == "circular_inclusions") {
                for (unsigned i=0;i<m.elem_list.size();++i) {
                    if ( pow(center( *m.elem_list[i] )[0] - 0.2, 2) + pow(center( *m.elem_list[i] )[1] - 0.2, 2) < pow(0.1 + 1e-6, 2) or pow(center( *m.elem_list[i] )[0] - 0.6, 2) + pow(center( *m.elem_list[i] )[1] - 0.3, 2) + pow(center( *m.elem_list[i] )[2] - 0.5, 2) < pow(0.1 + 1e-6, 2) or pow(center( *m.elem_list[i] )[0] - 0.4, 2) + pow(center( *m.elem_list[i] )[1] - 0.7, 2) + pow(center( *m.elem_list[i] )[2] - 0.8, 2) < pow(0.1 + 1e-6, 2) ) // ( x - 0.2 )^2 + ( y - 0.2 )^2 + ( z - 0.2 )^2 = (0.1)^2 or ( x - 0.6 )^2 + ( y - 0.3 )^2 + ( z - 0.5 )^2 = (0.1)^2 or ( x - 0.4 )^2 + ( y - 0.7 )^2 + ( z - 0.8 )^2 = (0.1)^2
                        list_elems_PGD_unknown_parameter.push_back( m.elem_list[i]->number );
                }
            }
        }
        /// Dimension 3
        ///------------
        else if ( dim == 3 ) {
            /// Inclusions spheriques 3D
            ///-------------------------
            if (structure == "spherical_inclusions") {
                for (unsigned i=0;i<m.elem_list.size();++i) {
                    if ( pow(center( *m.elem_list[i] )[0] - 0.2, 2) + pow(center( *m.elem_list[i] )[1] - 0.2, 2) + pow(center( *m.elem_list[i] )[2] - 0.2, 2) < pow(0.1 + 1e-6, 2) or pow(center( *m.elem_list[i] )[0] - 0.6, 2) + pow(center( *m.elem_list[i] )[1] - 0.3, 2) + pow(center( *m.elem_list[i] )[2] - 0.5, 2) < pow(0.1 + 1e-6, 2) or pow(center( *m.elem_list[i] )[0] - 0.4, 2) + pow(center( *m.elem_list[i] )[1] - 0.7, 2) + pow(center( *m.elem_list[i] )[2] - 0.8, 2) < pow(0.1 + 1e-6, 2) ) // ( x - 0.2 )^2 + ( y - 0.2 )^2 + ( z - 0.2 )^2 = (0.1)^2 or ( x - 0.6 )^2 + ( y - 0.3 )^2 + ( z - 0.5 )^2 = (0.1)^2 or ( x - 0.4 )^2 + ( y - 0.7 )^2 + ( z - 0.8 )^2 = (0.1)^2
                        list_elems_PGD_unknown_parameter.push_back( m.elem_list[i]->number );
                }
            }
        }
    }
}

/// Construction et resolution du pb en espace
///-------------------------------------------
template<class TM, class TF, class T, class TV, class TMAT, class TVVV>
void solve_PGD_space( TM &m, TF &f, const unsigned &n, const unsigned &k, const Vec<unsigned> &nb_iterations, const TV &F_s, const TV &F_p, const TMAT &K_unk_p, const TMAT &K_k_p, const bool &want_iterative_solver, const T &iterative_criterium, const Vec<unsigned> &list_elems_PGD_unknown_parameter, const TVVV &lambda, TVVV &psi ) {
    /// Construction du pb en espace
    ///-----------------------------
    T gamma_s = dot( F_p, lambda[ n ][ k ] );
    f.sollicitation = F_s * gamma_s;
    for (unsigned i=0;i<n;++i) {
        T alpha_s_i_unk = dot( lambda[ i ][ nb_iterations[ i ] ], K_unk_p * lambda[ n ][ k ] );
        T alpha_s_i_k = dot( lambda[ i ][ nb_iterations[ i ] ], K_k_p * lambda[ n ][ k ] );
        for (unsigned j=0;j<m.elem_list.size();++j) {
            if ( find( list_elems_PGD_unknown_parameter, _1 == j ) )
                m.elem_list[j]->set_field( "phi_elem_PGD_unknown_parameter", alpha_s_i_unk );
            else
                m.elem_list[j]->set_field( "phi_elem_PGD_unknown_parameter", alpha_s_i_k );
        }
        f.assemble( true, false );
        f.sollicitation -= f.matrices(Number<0>()) * psi[ i ][ nb_iterations[ i ] ];
    }
    T alpha_s_unk = dot( lambda[ n ][ k ], K_unk_p * lambda[ n ][ k ] );
    T alpha_s_k = dot( lambda[ n ][ k ], K_k_p * lambda[ n ][ k ] );
    for (unsigned i=0;i<m.elem_list.size();++i) {
        if ( find( list_elems_PGD_unknown_parameter, _1 == i ) )
            m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", alpha_s_unk );
        else
            m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", alpha_s_k );
    }
    f.assemble( true, false );
    
    /// Resolution du pb en espace
    ///---------------------------
    if ( want_iterative_solver == 0 )
        f.solve_system();
    else
        f.solve_system( iterative_criterium );
    f.update_variables();
    f.call_after_solve();
    
    /// Fonction psi en espace
    ///-----------------------
    psi[ n ][ k ] = f.vectors[0];
}

/// Construction et resolution du pb en parametre
///----------------------------------------------
template<class TM_param, class TF_unknown_param, class TV, class TMAT, class TVVV>
void solve_PGD_param( TM_param &m_param, TF_unknown_param &f_unknown_param, const unsigned &n, const unsigned &k, const Vec<unsigned> &nb_iterations, const TV &F_p, const TMAT &K_unk_p, const TMAT &K_k_p, const TV &F_s, const TMAT &K_unk_s, const TMAT &K_k_s, const TVVV &psi, TVVV &lambda ) {
    typedef typename TM_param::TNode::T T;
    
    T gamma_p = dot( F_s, psi[ n ][ k ] );
    T alpha_p_unk = dot( psi[ n ][ k ], K_unk_s * psi[ n ][ k ] );
    T alpha_p_k = dot( psi[ n ][ k ], K_k_s * psi[ n ][ k ] );
    /// Construction du pb en parametre
    ///--------------------------------
    f_unknown_param.sollicitation = F_p * gamma_p;
    for (unsigned i=0;i<n;++i) {
        T alpha_p_i_unk = dot( psi[ i ][ nb_iterations[ i ] ], K_unk_s * psi[ n ][ k ] );
        T alpha_p_i_k = dot( psi[ i ][ nb_iterations[ i ] ], K_k_s * psi[ n ][ k ] );
        f_unknown_param.sollicitation -= ( alpha_p_i_unk * K_unk_p + alpha_p_i_k * K_k_p ) * lambda[ i ][ nb_iterations[ i ] ];
    }
    f_unknown_param.matrices(Number<0>()) = TMAT( alpha_p_unk * K_unk_p + alpha_p_k * K_k_p );
    
    /// Resolution du pb en parametre
    ///------------------------------
    f_unknown_param.solve_system();
    f_unknown_param.update_variables();
    f_unknown_param.call_after_solve();
    
    /// Fonction lambda en parametre
    ///-----------------------------
    lambda[ n ][ k ] = f_unknown_param.vectors[0];
    
//     /// Resolution explicite du pb en parametre
//     ///----------------------------------------
//     for (unsigned j=0;j<m_param.node_list.size();++j) {
//         T function = 1 + m_param.node_list[ j ].pos[ 0 ];
//         lambda[ n ][ k ][ j ] = gamma_p;
//         for (unsigned i=0;i<n;++i) {
//             T alpha_p_i_unk = dot( psi[ i ][ nb_iterations[ i ] ], K_unk_s * psi[ n ][ k ] );
//             T alpha_p_i_k = dot( psi[ i ][ nb_iterations[ i ] ], K_k_s * psi[ n ][ k ] );
//             lambda[ n ][ k ][ j ] -= ( alpha_p_i_unk * function + alpha_p_i_k ) * lambda[ i ][ nb_iterations[ i ] ][ j ];
//         }
//         lambda[ n ][ k ][ j ] /= alpha_p_unk * function + alpha_p_k;
//     }
}

/// Construction des coefficients alpha, gamma associes au pb spatial
///------------------------------------------------------------------
template<class TE_param, class TM_param, class T> 
void construct_space_pb( const TE_param &elem_param, const TM_param &m_param, T &alpha_s, T &gamma_s ) {}

struct Construct_Space_Pb {
    template<class TE_param, class TM_param, class T> void operator()( const TE_param &elem_param, const TM_param &m_param, T &alpha_s, T &gamma_s ) const {
        construct_space_pb( elem_param, m_param, alpha_s, gamma_s );
    }
};

#endif // PGD_h
