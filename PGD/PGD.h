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

#include "../Material_properties.h"

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
            else if ( structure == "circular_inclusions" ) {
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
            if ( structure == "spherical_inclusions" ) {
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

/// Construction des opérateurs et du second membre en espace
/// ---------------------------------------------------------
template<class TM, class TF, class TV, class TMatV>
void assemble_space( TM &m, TF &f, TV &F_space, TMatV &K_space, const Vec< Vec<unsigned> > &elem_group ) {
    f.allocate_matrices();
    f.shift();
    f.assemble( false, true );
    F_space = f.sollicitation;
    K_space.resize( elem_group.size() );
    for (unsigned g=0;g<elem_group.size();++g) {
        for (unsigned j=0;j<m.elem_list.size();++j) {
            if ( find( elem_group[g], _1 == j ) )
                m.elem_list[j]->set_field( "alpha", 1. );
            else
                m.elem_list[j]->set_field( "alpha", 0. );
        }
        f.assemble( true, false );
        K_space[g] = f.matrices(Number<0>());
    }
    set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
}

/// Construction des opérateurs et seconds membres en parametre
/// -----------------------------------------------------------
template<class TM_param, class TF_param, class TVV, class TMatVV>
void assemble_param( TM_param &m_param, TF_param &f_param, TVV &F_param, TMatVV &K_param, const Vec< Vec<unsigned> > &elem_group ) {
    F_param.resize( elem_group.size()-1 );
    K_param.resize( elem_group.size()-1 );
    for (unsigned p=0;p<elem_group.size()-1;++p) {
        f_param[p].allocate_matrices();
        f_param[p].shift();
        f_param[p].assemble( false, true );
        F_param[p] = f_param[p].sollicitation;
        m_param[p].phi = 0;
        f_param[p].assemble( true, false );
        K_param[p][0] = f_param[p].matrices(Number<0>());
        m_param[p].phi = 1;
        f_param[p].assemble( true, false );
        K_param[p][1] = f_param[p].matrices(Number<0>());
    }
}

/// Construction et resolution du pb en espace
/// ------------------------------------------
template<class TM, class TF, class TV, class TVV, class TMatVV, class TTVVV, class TTVV>
void solve_space( TM &m, TF &f, const unsigned &mode, const TV &F_space, const TVV &F_param, const TMatVV &K_param, const Vec< Vec<unsigned> > &elem_group, const TTVVV &dep_param, TTVV &dep_space, const bool want_normalization = false ) {
    
    typedef typename TM::TNode::T T;
    
    /// Construction du pb en espace
    /// ----------------------------
    f.sollicitation = F_space;
    for (unsigned p=0;p<elem_group.size()-1;++p)
        f.sollicitation *= dot( F_param[p], dep_param[ p ][ mode ] );
    for (unsigned i=0;i<mode;++i) {
        for (unsigned g=0;g<elem_group.size();++g) {
            T alpha = 1.;
            for (unsigned p=0;p<elem_group.size()-1;++p) {
                if ( p == g )
                    alpha *= dot( dep_param[ p ][ i ], K_param[p][1] * dep_param[ p ][ mode ] );
                else
                    alpha *= dot( dep_param[ p ][ i ], K_param[p][0] * dep_param[ p ][ mode ] );
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
            if ( p == g )
                alpha *= dot( dep_param[ p ][ mode ], K_param[p][1] * dep_param[ p ][ mode ] );
            else
                alpha *= dot( dep_param[ p ][ mode ], K_param[p][0] * dep_param[ p ][ mode ] );
        }
        for (unsigned j=0;j<elem_group[g].size();++j)
            m.elem_list[ elem_group[g][j] ]->set_field( "alpha", alpha );
    }
    f.assemble( true, false );
    set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
    
    /// Resolution du pb en espace
    /// --------------------------
    f.solve_system();
    f.update_variables();
    
    /// Fonction en espace
    /// ------------------
    dep_space[ mode ] = f.vectors[0];
    
    /// Normalisation
    /// -------------
    if ( want_normalization ) {
        // set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
        f.assemble( true, false );
        dep_space[ mode ] /= sqrt( dot( dep_space[ mode ], f.matrices(Number<0>()) * dep_space[ mode ] ) );
    }
}

/// Construction et resolution du pb en parametre
/// ---------------------------------------------
template<class TM_param, class TF_param, class TV, class TVV, class TMatV, class TMatVV, class TTVV, class TTVVV>
void solve_param( TM_param &m_param, TF_param &f_param, const unsigned &p, const unsigned &mode, const TV &F_space, const TVV &F_param, const TMatV &K_space, const TMatVV &K_param, const Vec< Vec<unsigned> > &elem_group, const TTVV &dep_space, TTVVV &dep_param, const bool want_normalization = false ) {
    
    typedef typename TM_param::TNode::T T;
    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    
    /// Construction du pb en parametre
    /// -------------------------------
    T gamma = dot( F_space, dep_space[ mode ] );
    for (unsigned q=0;q<elem_group.size()-1;++q) {
        if ( q != p )
            gamma *= dot( F_param[q], dep_param[ q ][ mode ] );
    }
    f_param.sollicitation = F_param[p] * gamma;
    for (unsigned i=0;i<mode;++i) {
        for (unsigned g=0;g<elem_group.size();++g) {
            T alpha = dot( dep_space[ i ], K_space[g] * dep_space[ mode ] );
            for (unsigned q=0;q<elem_group.size()-1;++q) {
                if ( q != p ) {
                    if ( g == q )
                        alpha *= dot( dep_param[ q ][ i ], K_param[q][1] * dep_param[ q ][ mode ] );
                    else
                        alpha *= dot( dep_param[ q ][ i ], K_param[q][0] * dep_param[ q ][ mode ] );
                }
            }
            if ( g == p )
                f_param.sollicitation -= alpha * K_param[p][1] * dep_param[ p ][ i ];
            else
                f_param.sollicitation -= alpha * K_param[p][0] * dep_param[ p ][ i ];
        }
    }
    TMatSymSparse K = f_param.matrices(Number<0>());
    K.clear();
    for (unsigned g=0;g<elem_group.size();++g) {
        T alpha = dot( dep_space[ mode ], K_space[g] * dep_space[ mode ] );
        for (unsigned q=0;q<elem_group.size()-1;++q) {
            if ( q != p ) {
                if ( g == q )
                    alpha *= dot( dep_param[ q ][ mode ], K_param[q][1] * dep_param[ q ][ mode ] );
                else
                    alpha *= dot( dep_param[ q ][ mode ], K_param[q][0] * dep_param[ q ][ mode ] );
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
    
    /// Resolution explicite du pb en parametre
    /// ---------------------------------------
//    for (unsigned j=0;j<m_param.node_list.size();++j) {
//        m_param.node_list[ j ].fun = m_param.node_list[ j ].pos[ 0 ];
//        dep_param[ p ][ mode ][ j ] = gamma;
//        for (unsigned i=0;i<mode;++i) {
//            for (unsigned g=0;g<elem_group.size();++g) {
//                T alpha = dot( dep_space[ i ], K_space[g] * dep_space[ mode ] );
//                for (unsigned q=0;q<elem_group.size()-1;++q) {
//                    if ( q != p ) {
//                        if ( q == g )
//                            alpha *= dot( dep_param[ q ][ i ], K_param[q][1] * dep_param[ q ][ mode ] );
//                        else
//                            alpha *= dot( dep_param[ q ][ i ], K_param[q][0] * dep_param[ q ][ mode ] );
//                    }
//                }
//                if ( g == p )
//                    dep_param[ p ][ mode ][ j ] -= alpha * m_param.node_list[ j ].fun * dep_param[ p ][ i ][ j ];
//                else
//                    dep_param[ p ][ mode ][ j ] -= alpha * dep_param[ p ][ i ][ j ];
//            }
//        }
//        T delta = 0.;
//        for (unsigned g=0;g<elem_group.size();++g) {
//            T alpha = dot( dep_space[ mode ], K_space[g] * dep_space[ mode ] );
//            for (unsigned q=0;q<elem_group.size()-1;++q) {
//                if ( q != p ) {
//                    if ( q == g )
//                        alpha *= dot( dep_param[ q ][ mode ], K_param[q][1] * dep_param[ q ][ mode ] );
//                    else
//                        alpha *= dot( dep_param[ q ][ mode ], K_param[q][0] * dep_param[ q ][ mode ] );
//                }
//            }
//            if ( g == p )
//                delta += alpha * m_param.node_list[ j ].fun;
//            else
//                delta += alpha;
//        }
//        dep_param[ p ][ mode ][ j ] /= delta;
//    }
    
    /// Fonction en parametre
    /// ---------------------
    dep_param[ p ][ mode ] = f_param.vectors[0];
    
    /// Normalisation
    /// -------------
    if ( want_normalization )
        dep_param[ p ][ mode ] /= sqrt( dot( dep_param[ p ][ mode ], K_param[p][1] * dep_param[ p ][ mode ] ) );
}

/// Calcul de l'indicateur de stagnation
/// ------------------------------------
template<class TM, class TF, class T, class TV, class TVV, class TMatVV, class TTVV, class TTVVV>
void calc_stagnation_indicator( TM &m, TF &f, const unsigned &mode, const TMatVV &K_param, const Vec< Vec<unsigned> > &elem_group, const TTVV &dep_space, const TTVVV &dep_param, const TV &dep_space_old, const TVV &dep_param_old, T &stagnation_indicator ) {
    stagnation_indicator = 0.;
    for (unsigned g=0;g<elem_group.size();++g) {
        T alpha = 1.;
        for (unsigned p=0;p<elem_group.size()-1;++p) {
            if ( g == p )
                alpha *= dot( dep_param[ p ][ mode ], K_param[p][1] * dep_param[ p ][ mode ] );
            else
                alpha *= dot( dep_param[ p ][ mode ], K_param[p][0] * dep_param[ p ][ mode ] );
        }
        for (unsigned j=0;j<elem_group[g].size();++j)
            m.elem_list[ elem_group[g][j] ]->set_field( "alpha", alpha );
    }
    f.assemble( true, false );
    stagnation_indicator += dot( dep_space[ mode ], f.matrices(Number<0>()) * dep_space[ mode ] );
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
                alpha *= dot( dep_param[ p ][ mode ], K_param[p][1] * dep_param_old[ p ] );
            else
                alpha *= dot( dep_param[ p ][ mode ], K_param[p][0] * dep_param_old[ p ] );
        }
        for (unsigned j=0;j<elem_group[g].size();++j)
            m.elem_list[ elem_group[g][j] ]->set_field( "alpha", alpha );
    }
    f.assemble( true, false );
    stagnation_indicator -= 2 * dot( dep_space[ mode ], f.matrices(Number<0>()) * dep_space_old );
    stagnation_indicator = sqrt( fabs( stagnation_indicator ) );
    
    set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
}

/// Calcul de l'indicateur d'erreur
/// -------------------------------
template<class TM, class TF, class T, class TV, class TVV, class TMatVV, class TTVV, class TTVVV>
void calc_error_indicator( TM &m, TF &f, const unsigned &mode, const TV &F_space, const TVV &F_param, const TMatVV &K_param, const Vec< Vec<unsigned> > &elem_group, const TTVV &dep_space, const TTVVV &dep_param, T &residual, T &error_indicator ) {
    residual = 0.;
    T sollicitation = 0.;
    for (unsigned i=0;i<mode+1;++i) {
        T sollicitation_mode = dot( F_space, dep_space[ i ] );
        for (unsigned p=0;p<elem_group.size()-1;++p)
            sollicitation_mode *= dot( F_param[p], dep_param[ p ][ i ] );
        sollicitation += sollicitation_mode;
        residual -= sollicitation_mode;
        for (unsigned j=0;j<mode+1;++j) {
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
    
    set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
}

/// Construction du champ de deplacement EF en espace
/// -------------------------------------------------
template<class TM, class TVV, class TMatVV, class TTVVV, class TTVV, class TV, class TTTVV>
void construct_dep_space_FE( const TM &m, const TVV &F_param, const TMatVV &K_param, const Vec< Vec<unsigned> > &elem_group, unsigned &mode, const TTVVV &dep_param, const TTVV &dep_space, const TV &dep_space_FE_part, TTTVV &dep_space_FE ) {
    
    typedef typename TM::TNode::T T;
    
    dep_space_FE.resize( elem_group.size() );
    for (unsigned g=0;g<elem_group.size();++g) {
        dep_space_FE[ g ] = - dep_space_FE_part;
        for (unsigned p=0;p<elem_group.size()-1;++p)
            dep_space_FE[ g ] *= dot( F_param[p], dep_param[ p ][ mode ] );
        for (unsigned i=0;i<mode+1;++i) {
            T alpha = 1.;
            for (unsigned p=0;p<elem_group.size()-1;++p) {
                if ( p == g )
                    alpha *= dot( dep_param[ p ][ mode ], K_param[p][1] * dep_param[ p ][ i ] );
                else
                    alpha *= dot( dep_param[ p ][ mode ], K_param[p][0] * dep_param[ p ][ i ] );
            }
            dep_space_FE[ g ] += alpha * dep_space[ i ];
        }
    }
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


/// Evaluation de la solution PGD pour un jeu connu de parametres
/// -------------------------------------------------------------
template<class TM_param, class TM, class TF, class TVV, class TTVV, class TTVVV>
void eval_PGD( TM_param &m_param, TM &m, TF &f, const string &pb, const string &structure, const string &boundary_condition_D, const string &loading, const string &mesh_size, const Vec< Vec<unsigned> > &elem_group, const unsigned &nb_vals, const TVV &vals_param, const unsigned &mode, const TTVV &dep_space, const TTVVV &dep_param, const string &filename = "paraview", const bool display_pvd = false ) {
    
    typedef typename TM::TNode::T T;
    
    Vec< DisplayParaview > dp;
    dp.resize( nb_vals );
    Vec<string> lp("all");
//    lp.push_back( "dep" );
//    lp.push_back( "young_eff" );
    
    for (unsigned i=0;i<nb_vals;++i) {
        Vec<unsigned> ind;
        ind.resize( elem_group.size()-1 );
        // generate pseudo-random integral number ind[p] in the range between 0 and m_param[p].node_list.size()-1
        for (unsigned p=0;p<elem_group.size()-1;++p)
            ind[p] = rand() % m_param[p].node_list.size();
        
        /// Proprietes materiaux du pb direct
        /// ---------------------------------
        set_material_properties( f, m, structure );
        
        /// Conditions aux limites du pb direct
        /// -----------------------------------
        set_constraints( f, m, boundary_condition_D, "direct", structure, loading );
        set_load_conditions( m, structure, loading, mesh_size );
        
        for (unsigned p=0;p<elem_group.size()-1;++p) {
            for (unsigned j=0;j<elem_group[p].size();++j)
                m.elem_list[ elem_group[p][j] ]->set_field( "alpha", vals_param[p][ ind[p] ] );
        }
        for (unsigned j=0;j<elem_group.back().size();++j)
            m.elem_list[ elem_group.back()[j] ]->set_field( "alpha", 1. );
        
        // Reference FE solution
        f.solve();
        
        string filename_ = filename + "_space_REF_params";
        for (unsigned p=0;p<elem_group.size()-1;++p)
            filename_ += '_' + to_string( vals_param[p][ ind[p] ] );
        dp[ i ].add_mesh_iter( m, filename_, lp, 0 );
        
        // PGD solution
        f.vectors[0].set( 0. );
        for (unsigned n=0;n<mode+1;++n) {
            Vec<T> dep_mode = dep_space[ n ];
            for (unsigned p=0;p<elem_group.size()-1;++p)
                dep_mode *= dep_param[ p ][ n ][ ind[p] ];
            f.vectors[0] +=  dep_mode;
        }
        f.update_variables();
        f.call_after_solve();
        
        filename_ = filename + "_space_PGD_params";
        for (unsigned p=0;p<elem_group.size()-1;++p)
            filename_ += '_' + to_string( vals_param[p][ ind[p] ] );
        dp[ i ].add_mesh_iter( m, filename_, lp, 1 );
        
        filename_ = filename + "_space_params";
        for (unsigned p=0;p<elem_group.size()-1;++p)
            filename_ += '_' + to_string( vals_param[p][ ind[p] ] );
        if ( display_pvd )
            dp[ i ].exec( filename_ );
        else
            dp[ i ].make_pvd_file( filename_ );
    }
    
    set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
}

#endif // PGD_h
