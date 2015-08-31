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
/// ---------------------------------
template<class TM>
void partition_elem_list( TM &m, const string &structure, Vec< Vec<unsigned> > &elem_group ) {

    static const unsigned dim = TM::dim;

    if ( m.node_list.size() ) {
        /// Dimension 2
        /// -----------
        if ( dim == 2 ) {
            /// Plaque rectangulaire 2D en flexion
            /// ----------------------------------
            if ( structure == "plate_flexion" ) {
                elem_group.resize(2);
                for (unsigned n=0;n<m.elem_list.size();++n) {
                    if ( center( *m.elem_list[n] )[1] < 0.5 )
                        elem_group[0].push_back( m.elem_list[n]->number );
                    else
                        elem_group[1].push_back( m.elem_list[n]->number );
                }
            }
            /// Inclusions circulaires 2D
            /// -------------------------
            else if (structure == "circular_inclusions") {
                elem_group.resize(4);
                for (unsigned n=0;n<m.elem_list.size();++n) {
                    if ( pow(center( *m.elem_list[n] )[0] - 0.2, 2) + pow(center( *m.elem_list[n] )[1] - 0.2, 2) < pow(0.1 + 1e-6, 2) ) // ( x - 0.2 )^2 + ( y - 0.2 )^2 = (0.1)^2
                        elem_group[0].push_back( m.elem_list[n]->number );
                    else if ( pow(center( *m.elem_list[n] )[0] - 0.6, 2) + pow(center( *m.elem_list[n] )[1] - 0.3, 2) < pow(0.1 + 1e-6, 2) ) // ( x - 0.6 )^2 + ( y - 0.3 )^2 = (0.1)^2
                        elem_group[1].push_back( m.elem_list[n]->number );
                    else if ( pow(center( *m.elem_list[n] )[0] - 0.4, 2) + pow(center( *m.elem_list[n] )[1] - 0.7, 2) < pow(0.1 + 1e-6, 2) ) // ( x - 0.4 )^2 + ( y - 0.7 )^2 = (0.1)^2
                        elem_group[2].push_back( m.elem_list[n]->number );
                    else
                        elem_group[3].push_back( m.elem_list[n]->number );
                }
            }
        }
        /// Dimension 3
        /// -----------
        else if ( dim == 3 ) {
            /// Inclusions spheriques 3D
            /// ------------------------
            if (structure == "spherical_inclusions") {
                elem_group.resize(4);
                for (unsigned n=0;n<m.elem_list.size();++n) {
                    if ( pow(center( *m.elem_list[n] )[0] - 0.2, 2) + pow(center( *m.elem_list[n] )[1] - 0.2, 2) + pow(center( *m.elem_list[n] )[2] - 0.2, 2) < pow(0.1 + 1e-6, 2) ) // ( x - 0.2 )^2 + ( y - 0.2 )^2 + ( z - 0.2 )^2 = (0.1)^2
                        elem_group[0].push_back( m.elem_list[n]->number );
                    else if ( pow(center( *m.elem_list[n] )[0] - 0.6, 2) + pow(center( *m.elem_list[n] )[1] - 0.3, 2) + pow(center( *m.elem_list[n] )[2] - 0.5, 2) < pow(0.1 + 1e-6, 2) ) // ( x - 0.6 )^2 + ( y - 0.3 )^2 + ( z - 0.5 )^2 = (0.1)^2
                        elem_group[1].push_back( m.elem_list[n]->number );
                    else if ( pow(center( *m.elem_list[n] )[0] - 0.4, 2) + pow(center( *m.elem_list[n] )[1] - 0.7, 2) + pow(center( *m.elem_list[n] )[2] - 0.8, 2) < pow(0.1 + 1e-6, 2) ) // ( x - 0.4 )^2 + ( y - 0.7 )^2 + ( z - 0.8 )^2 = (0.1)^2
                        elem_group[2].push_back( m.elem_list[n]->number );
                    else
                        elem_group[3].push_back( m.elem_list[n]->number );
                }
            }
        }
    }
}

/// Construction et resolution du pb en espace
/// ------------------------------------------
template<class TM, class TF, class T, class TV, class TVV, class TMATVV, class TTVVV, class TTVV>
void solve_space( TM &m, TF &f, const unsigned &n, const TV &F_space, const TVV &F_param, const TMATVV &K_param, const Vec< Vec<unsigned> > &elem_group, const TTVVV &dep_param, TTVV &dep_space, const bool want_iterative_solver = false, const T iterative_criterium = 1e-3, const bool want_normalization = false ) {

    /// Construction du pb en espace
    /// ----------------------------
    f.sollicitation = F_space;
    for (unsigned p=0;p<elem_group.size()-1;++p)
        f.sollicitation *= dot( F_param[p], dep_param[ p ][ n ] );
    for (unsigned i=0;i<n;++i) {
        for (unsigned g=0;g<elem_group.size();++g) {
            T alpha = 1.;
            for (unsigned p=0;p<elem_group.size()-1;++p) {
                if ( g == p )
                    alpha *= dot( dep_param[ p ][ i ], K_param[p][1] * dep_param[ p ][ n ] );
                else
                    alpha *= dot( dep_param[ p ][ i ], K_param[p][0] * dep_param[ p ][ n ] );
            }
            for (unsigned j=0;j<elem_group[g].size();++j)
                m.elem_list[ elem_group[g][j] ]->set_field( "alpha", alpha );
        }
        f.assemble( true, false );
        f.sollicitation -= f.matrices(Number<0>()) * dep_space[ i ];
    }

    for (unsigned g=0;g<elem_group.size();++g) {
        T alpha = 1.;
        for (unsigned p=0;p<elem_group.size()-1;++p) {
            if ( g == p )
                alpha *= dot( dep_param[ p ][ n ], K_param[p][1] * dep_param[ p ][ n ] );
            else
                alpha *= dot( dep_param[ p ][ n ], K_param[p][0] * dep_param[ p ][ n ] );
        }
        for (unsigned j=0;j<elem_group[g].size();++j)
            m.elem_list[ elem_group[g][j] ]->set_field( "alpha", alpha );
    }
    f.assemble( true, false );
    
    /// Resolution du pb en espace
    /// --------------------------
    if ( want_iterative_solver == 0 )
        f.solve_system();
    else
        f.solve_system( iterative_criterium );
    f.update_variables();
    f.call_after_solve();
    
    /// Fonction en espace
    /// ------------------
    dep_space[ n ] = f.vectors[0];

    /// Normalisation
    /// -------------
    if ( want_normalization ) {
        for (unsigned j=0;j<m.elem_list.size();++j)
            m.elem_list[j]->set_field( "alpha", 1. );
        f.assemble( true, false );
        dep_space[ n ] /= sqrt( dot( dep_space[ n ], f.matrices(Number<0>()) * dep_space[ n ] ) );
    }
}

/// Construction et resolution du pb en parametre
/// ---------------------------------------------
template<class TM_param, class TF_param, class TV, class TVV, class TMATV, class TMATVV, class TTVV, class TTVVV>
void solve_param( TM_param &m_param, TF_param &f_param, const unsigned &p, const unsigned &n, const TV &F_space, const TVV &F_param, const TMATV &K_space, const TMATVV &K_param, const Vec< Vec<unsigned> > &elem_group, const TTVV &dep_space, TTVVV &dep_param, const bool want_normalization = false ) {

    typedef typename TM_param::TNode::T T;
    typedef Mat<T, Sym<>, SparseLine<> > TMAT;
    
    /// Construction du pb en parametre
    /// -------------------------------
    T gamma = dot( F_space, dep_space[ n ] );
    for (unsigned q=0;q<elem_group.size()-1;++q) {
        if ( q != p )
            gamma *= dot( F_param[q], dep_param[ q ][ n ] );
    }
    f_param.sollicitation = F_param[p] * gamma;
    for (unsigned i=0;i<n;++i) {
        for (unsigned g=0;g<elem_group.size();++g) {
            T alpha = dot( dep_space[ i ], K_space[g] * dep_space[ n ] );
            for (unsigned q=0;q<elem_group.size()-1;++q) {
                if ( q != p ) {
                    if ( g == q )
                        alpha *= dot( dep_param[ q ][ i ], K_param[q][1] * dep_param[ q ][ n ] );
                    else
                        alpha *= dot( dep_param[ q ][ i ], K_param[q][0] * dep_param[ q ][ n ] );
                }
            }
            if ( g == p )
                f_param.sollicitation -= alpha * K_param[p][1] * dep_param[ p ][ i ];
            else
                f_param.sollicitation -= alpha * K_param[p][0] * dep_param[ p ][ i ];
        }
    }
    TMAT K = f_param.matrices(Number<0>());
    K.clear();
    for (unsigned g=0;g<elem_group.size();++g) {
        T alpha = dot( dep_space[ n ], K_space[g] * dep_space[ n ] );
        for (unsigned q=0;q<elem_group.size()-1;++q) {
            if ( q != p ) {
                if ( g == q )
                    alpha *= dot( dep_param[ q ][ n ], K_param[q][1] * dep_param[ q ][ n ] );
                else
                    alpha *= dot( dep_param[ q ][ n ], K_param[q][0] * dep_param[ q ][ n ] );
            }
        }
        if ( g == p )
            K += alpha * K_param[p][1];
        else
            K += alpha * K_param[p][0];
    }
    f_param.matrices(Number<0>()) = K;
    
    /// Resolution du pb en parametre
    /// -----------------------------
    f_param.solve_system();
    f_param.update_variables();
    f_param.call_after_solve();
    
    /// Resolution explicite du pb en parametre
    /// ---------------------------------------
//    for (unsigned j=0;j<m_param.node_list.size();++j) {
//        m_param.node_list[ j ].fun = m_param.node_list[ j ].pos[ 0 ];
//        dep_param[ p ][ n ][ j ] = gamma;
//        for (unsigned i=0;i<n;++i) {
//            for (unsigned g=0;g<elem_group.size();++g) {
//                T alpha = dot( dep_space[ i ], K_space[g] * dep_space[ n ] );
//                for (unsigned q=0;q<elem_group.size()-1;++q) {
//                    if ( q != p ) {
//                        if ( g == q )
//                            alpha *= dot( dep_param[ q ][ i ], K_param[q][1] * dep_param[ q ][ n ] );
//                        else
//                            alpha *= dot( dep_param[ q ][ i ], K_param[q][0] * dep_param[ q ][ n ] );
//                    }
//                }
//                if ( g == p )
//                    dep_param[ p ][ n ][ j ] -= alpha * m_param.node_list[ j ].fun * dep_param[ p ][ i ][ j ];
//                else
//                    dep_param[ p ][ n ][ j ] -= alpha * dep_param[ p ][ i ][ j ];
//            }
//        }
//        T delta = 0.;
//        for (unsigned g=0;g<elem_group.size();++g) {
//            T alpha = dot( dep_space[ n ], K_space[g] * dep_space[ n ] );
//            for (unsigned q=0;q<elem_group.size()-1;++q) {
//                if ( q != p ) {
//                    if ( g == q )
//                        alpha *= dot( dep_param[ q ][ n ], K_param[q][1] * dep_param[ q ][ n ] );
//                    else
//                        alpha *= dot( dep_param[ q ][ n ], K_param[q][0] * dep_param[ q ][ n ] );
//                }
//            }
//            if ( g == p )
//                delta += alpha * m_param.node_list[ j ].fun;
//            else
//                delta += alpha;
//        }
//        dep_param[ p ][ n ][ j ] /= delta;
//    }

    /// Fonction en parametre
    /// ---------------------
    dep_param[ p ][ n ] = f_param.vectors[0];

    /// Normalisation
    /// -------------
    if ( want_normalization )
        dep_param[ p ][ n ] /= sqrt( dot( dep_param[ p ][ n ], K_param[p][1] * dep_param[ p ][ n ] ) );
}

/// Calcul de l'indicateur de stagnation
/// ------------------------------------
template<class TM, class TF, class T, class TV, class TVV, class TMATVV, class TTVV, class TTVVV>
void calc_stagnation_indicator( TM &m, TF &f, const unsigned &n, const TMATVV &K_param, const Vec< Vec<unsigned> > &elem_group, const TTVVV &dep_param, const TVV &dep_param_old, const TTVV &dep_space, const TV &dep_space_old, T &stagnation_indicator ) {
    stagnation_indicator = 0.;
    for (unsigned g=0;g<elem_group.size();++g) {
        T alpha = 1.;
        for (unsigned p=0;p<elem_group.size()-1;++p) {
            if ( g == p )
                alpha *= dot( dep_param[ p ][ n ], K_param[p][1] * dep_param[ p ][ n ] );
            else
                alpha *= dot( dep_param[ p ][ n ], K_param[p][0] * dep_param[ p ][ n ] );
        }
        for (unsigned j=0;j<elem_group[g].size();++j)
            m.elem_list[ elem_group[g][j] ]->set_field( "alpha", alpha );
    }
    f.assemble( true, false );
    stagnation_indicator += dot( dep_space[ n ], f.matrices(Number<0>()) * dep_space[ n ] );
    for (unsigned g=0;g<elem_group.size();++g) {
        T alpha = 1.;
        for (unsigned p=0;p<elem_group.size()-1;++p) {
            if ( g == p )
                alpha *= dot( dep_param_old[ p ], K_param[p][1] * dep_param_old[ p ] );
            else
                alpha *= dot( dep_param_old[ p ], K_param[p][0] * dep_param_old[ p ] );
        }
        for (unsigned j=0;j<elem_group[g].size();++j)
            m.elem_list[ elem_group[g][j] ]->set_field( "alpha", alpha );
    }
    f.assemble( true, false );
    stagnation_indicator += dot( dep_space_old, f.matrices(Number<0>()) * dep_space_old );
    for (unsigned g=0;g<elem_group.size();++g) {
        T alpha = 1.;
        for (unsigned p=0;p<elem_group.size()-1;++p) {
            if ( g == p )
                alpha *= dot( dep_param[ p ][ n ], K_param[p][1] * dep_param_old[ p ] );
            else
                alpha *= dot( dep_param[ p ][ n ], K_param[p][0] * dep_param_old[ p ] );
        }
        for (unsigned j=0;j<elem_group[g].size();++j)
            m.elem_list[ elem_group[g][j] ]->set_field( "alpha", alpha );
    }
    f.assemble( true, false );
    stagnation_indicator -= 2 * dot( dep_space[ n ], f.matrices(Number<0>()) * dep_space_old );
    stagnation_indicator = sqrt( fabs( stagnation_indicator ) );
}

/// Calcul de l'indicateur d'erreur
/// -------------------------------
template<class TM, class TF, class T, class TV, class TVV, class TMATVV, class TTVV, class TTVVV>
void calc_error_indicator( TM &m, TF &f, const unsigned &n, const TV &F_space, const TVV &F_param, const TMATVV &K_param, const Vec< Vec<unsigned> > &elem_group, const TTVVV &dep_param, const TTVV &dep_space, T &residual, T &error_indicator ) {
    residual = 0.;
    T sollicitation = 0.;
    for (unsigned i=0;i<n+1;++i) {
        T sollicitation_mode = dot( F_space, dep_space[ i ] );
        for (unsigned p=0;p<elem_group.size()-1;++p)
            sollicitation_mode *= dot( F_param[p], dep_param[ p ][ i ] );
        sollicitation += sollicitation_mode;
        residual -= sollicitation_mode;
        for (unsigned j=0;j<n+1;++j) {
            for (unsigned g=0;g<elem_group.size();++g) {
                T alpha = 1.;
                for (unsigned p=0;p<elem_group.size()-1;++p) {
                    if ( g == p )
                        alpha *= dot( dep_param[ p ][ i ], K_param[p][1] * dep_param[ p ][ j ] );
                    else
                        alpha *= dot( dep_param[ p ][ i ], K_param[p][0] * dep_param[ p ][ j ] );
                }
                for (unsigned l=0;l<elem_group[g].size();++l)
                    m.elem_list[ elem_group[g][l] ]->set_field( "alpha", alpha );
            }
            f.assemble( true, false );
            residual += dot( dep_space[ i ], f.matrices(Number<0>()) * dep_space[ j ] );
        }
    }
    error_indicator = fabs( residual ) / fabs( sollicitation );
}

/// Construction des coefficients alpha, gamma associes au pb spatial
/// -----------------------------------------------------------------
template<class TE_param, class TM_param, class T> 
void construct_space_pb( const TE_param &elem_param, const TM_param &m_param, T &alpha_s, T &gamma_s ) {}

struct Construct_Space_Pb {
    template<class TE_param, class TM_param, class T> void operator()( const TE_param &elem_param, const TM_param &m_param, T &alpha_s, T &gamma_s ) const {
        construct_space_pb( elem_param, m_param, alpha_s, gamma_s );
    }
};

#endif // PGD_h
