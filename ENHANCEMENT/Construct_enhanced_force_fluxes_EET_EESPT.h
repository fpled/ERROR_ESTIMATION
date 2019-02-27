//
// C++ Interface: Construct_enhanced_force_fluxes_EET_EESPT
//
// Description: construction amelioree des densites d'effort par les methodes EET et EESPT
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_enhanced_force_fluxes_EET_EESPT_h
#define Construct_enhanced_force_fluxes_EET_EESPT_h

#include "Enhancement_EET_EESPT.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"
#include "../TOOLS/Algebre.h"

using namespace LMT;
using namespace std;

/// Construction amelioree des densites d'effort par les methodes EET et EESPT
/// --------------------------------------------------------------------------
template<class TM, class TF, class T, class TVV, class TVVV, class TMatV>
void construct_enhanced_force_fluxes_EET_EESPT( TM &m, const TF &f, const string &method, const Vec<bool> &elem_flag_enh, const Vec<bool> &face_flag_enh, const Vec<bool> &elem_flag_bal, const Vec<unsigned> &elem_list_enh, const Vec<unsigned> &face_list_enh, const Vec<unsigned> &elem_list_bal, TMatV &K_hat, const TVV &dep_hat, TVVV &force_fluxes, const string &solver, const string &solver_minimisation, const bool verif_solver_enhancement = false, const T tol_solver_enhancement = 1e-6, const bool verif_solver_minimisation_enhancement = false, const T tol_solver_minimisation_enhancement = 1e-6, const bool disp = false ) {
    
    static const unsigned dim = TM::dim;
    
    typedef Mat<T, Gen<>, SparseLine<> > TMatGenSparse;
    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    typedef Mat<T, Sym<>, SparseCholMod > TMatSymSparseCholMod;
    typedef Mat<T, Gen<>, SparseUMFPACK > TMatSymSparseUMFPACK;
    typedef Mat<T, Gen<>, SparseLU > TMatGenSparseLU;
    
    TicToc t_force_fluxes_enh;
    t_force_fluxes_enh.start();
    
    Vec<unsigned> node_cpt_face;
    Vec< Vec<unsigned> > node_list_face;
    construct_nodes_connected_to_face( m, node_cpt_face, node_list_face );
    
    Vec<unsigned> elem_cpt_node;
    Vec< Vec<unsigned> > elem_list_node;
    construct_elems_connected_to_node( m, elem_cpt_node, elem_list_node );
    
    elem_list_node.free();
    
    Vec<bool> correspondance_node_to_vertex_node;
    Vec<unsigned> connect_node_to_vertex_node;
    unsigned nb_vertex_nodes = match_node_to_vertex_node( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node );
    
    connect_node_to_vertex_node.free();
    
    Vec< Vec<unsigned> > face_type;
    construct_face_type( m, f, face_type );
    
    if ( disp ) {
        cout << "Construction amelioree des densites d'effort" << endl;
        cout << "--------------------------------------------" << endl << endl;
    }
    
    unsigned nb_unk_enh = 0;
    
    for (unsigned k=0;k<face_list_enh.size();++k) {
        nb_unk_enh += m.sub_mesh(Number<1>()).elem_list[ face_list_enh[ k ] ]->nb_nodes_virtual() * dim;
    }
    
    /// Construction des vecteurs F_hat_enh[ n ] pour chaque element ameliore n du maillage
    /// -----------------------------------------------------------------------------------
    
    if ( disp ) {
        cout << "Calcul des vecteurs nb_unk_local_enh" << endl << endl;
        cout << "Construction des vecteurs F_hat_enh" << endl << endl;
    }
    
    Vec<unsigned> nb_unk_local_enh;
    nb_unk_local_enh.resize( elem_list_enh.size() );
    
    Vec< Vec< Vec<T> > > F_hat_enh;
    F_hat_enh.resize( elem_list_enh.size() );
    
    Calcul_Elem_Vector_F_hat_enh calc_elem_vector_F_hat_enh;
    calc_elem_vector_F_hat_enh.elem_flag_enh = &elem_flag_enh;
    calc_elem_vector_F_hat_enh.elem_list_enh = &elem_list_enh;
    calc_elem_vector_F_hat_enh.nb_unk_local_enh = &nb_unk_local_enh;
    
    apply( m.elem_list, calc_elem_vector_F_hat_enh, m, f, F_hat_enh );
    
    if ( disp ) {
        for (unsigned n=0;n<elem_list_enh.size();++n) {
            cout << "vecteur F_hat_enh de l'element " << elem_list_enh[ n ] << " =" << endl;
            for (unsigned i=0;i<nb_unk_local_enh[ n ];++i)
                cout << F_hat_enh[ n ][ i ] << endl;
        }
    }
    
    /// Construction des vecteurs dep_hat_enh[ n ] pour chaque element ameliore n du maillage
    /// -------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs dep_hat_enh" << endl << endl;
    
    Vec< Vec< Vec<T> > > dep_hat_enh;
    dep_hat_enh.resize( elem_list_enh.size() );
    
    for (unsigned n=0;n<elem_list_enh.size();++n) {
        dep_hat_enh[ n ].resize( nb_unk_local_enh[ n ] );
        for (unsigned i=0;i<nb_unk_local_enh[ n ];++i) {
            if ( solver == "CholMod" ) {
                #ifdef WITH_CHOLMOD
                TMatSymSparse K_hat_sym = K_hat[ elem_list_enh[ n ] ];
                TMatSymSparseCholMod K_hat_CholMod = K_hat_sym;
                K_hat_CholMod.get_factorization();
                dep_hat_enh[ n ][ i ] = K_hat_CholMod.solve( F_hat_enh[ n ][ i ] );
                #endif
            }
            else if ( solver == "LDL" ) {
                #ifdef WITH_LDL
                TMatSymSparse K_hat_LDL = K_hat[ elem_list_enh[ n ] ];
                dep_hat_enh[ n ][ i ] = F_hat_enh[ n ][ i ];
                LDL_solver ls;
                Vec< Vec<T> > Ker;
                Vec<int> Pivots;
//                bool want_free_matrix = false;
                ls.get_factorization( K_hat_LDL, Ker, Pivots/*, want_free_matrix*/ );
                ls.solve( dep_hat_enh[ n ][ i ] );
                #endif
            }
            else if ( solver == "UMFPACK" ) {
                #ifdef WITH_UMFPACK
                TMatSymSparseUMFPACK K_hat_UMFPACK = K_hat[ elem_list_enh[ n ] ];
                K_hat_UMFPACK.get_factorization();
                dep_hat_enh[ n ][ i ] = K_hat_UMFPACK.solve( F_hat_enh[ n ][ i ] );
                #endif
            }
            else if ( solver == "CholFactorize" ) {
                TMatSymSparse K_hat_Chol = K_hat[ elem_list_enh[ n ] ];
                chol_factorize( K_hat_Chol );
                solve_using_chol_factorize( K_hat_Chol, F_hat_enh[ n ][ i ], dep_hat_enh[ n ][ i ] );
            }
            else if ( solver == "LUFactorize" ) {
                TMatSymSparse K_hat_sym = K_hat[ elem_list_enh[ n ] ];
                TMatGenSparseLU K_hat_LU = K_hat_sym;
//                Vec<int> vector_permutation;
                lu_factorize( K_hat_LU/*, vector_permutation*/ );
                solve_using_lu_factorize( K_hat_LU, /*vector_permutation,*/ F_hat_enh[ n ][ i ], dep_hat_enh[ n ][ i ] );
            }
            else if ( solver == "Inv" ) {
                TMatSymSparse K_hat_Inv = K_hat[ elem_list_enh[ n ] ];
                dep_hat_enh[ n ][ i ] = inv( K_hat_Inv ) * F_hat_enh[ n ][ i ];
            }
            else {
                cerr << "Bing. Error : solveur " << solver << " pour la resolution des pbs locaux non implemente" << endl << endl;
            }
        }
    }
    
    if ( disp ) {
        for (unsigned n=0;n<elem_list_enh.size();++n) {
            cout << "vecteur dep_hat_enh de l'element " << elem_list_enh[ n ] << " =" << endl;
            for (unsigned i=0;i<nb_unk_local_enh[ n ];++i)
                cout << dep_hat_enh[ n ][ i ] << endl;
        }
    }
    if ( verif_solver_enhancement ) {
        if ( disp )
            cout << "Verification de la resolution des pbs locaux par element : tolerance = " << tol_solver_enhancement << endl << endl;
        for (unsigned n=0;n<elem_list_enh.size();++n) {
            for (unsigned i=0;i<nb_unk_local_enh[ n ];++i) {
                T residual = norm_2( K_hat[ elem_list_enh[ n ] ] * dep_hat_enh[ n ][ i ] - F_hat_enh[ n ][ i ] );
                T b = norm_2( F_hat_enh[ n ][ i ] );
                if ( residual / b > tol_solver_enhancement ) {
                    cout << "residu pour l'element " << elem_list_enh[ n ] << " avec cas de charge " << i << " :" << endl;
//                    cout << "K_hat * dep_hat_enh - F_hat_enh =" << endl;
//                    cout << K_hat[ elem_list_enh[ n ] ] * dep_hat_enh[ n ][ i ] - F_hat_enh[ n ][ i ] << endl;
                    cout << "norme du residu = " << residual << endl;
                    cout << "norme du residu relatif = " << residual / b << endl << endl;
                }
            }
        }
    }
    
    F_hat_enh.free();
    
    /// Construction des matrices A_local_enh[ n ] pour chaque element ameliore n du maillage
    /// -------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices A_local_enh" << endl << endl;
    
    Vec< TMatGenSparse > A_local_enh;
    A_local_enh.resize( elem_list_enh.size() );
    
    for (unsigned n=0;n<elem_list_enh.size();++n) {
        A_local_enh[ n ].resize( nb_unk_local_enh[ n ] );
    }
    
    Calcul_Elem_Matrix_A_enh<T> calcul_elem_matrix_A_enh;
    calcul_elem_matrix_A_enh.elem_flag_enh = &elem_flag_enh;
    calcul_elem_matrix_A_enh.elem_list_enh = &elem_list_enh;
    calcul_elem_matrix_A_enh.dep_hat_enh = &dep_hat_enh;
    
    apply( m.elem_list, calcul_elem_matrix_A_enh, m, f, A_local_enh );
    
    if ( disp ) {
        for (unsigned n=0;n<elem_list_enh.size();++n) {
            cout << "dimension de la matrice A_local_enh de l'element " << elem_list_enh[ n ] << " = ( " << nb_unk_local_enh[ n ] << ", " << nb_unk_local_enh[ n ] << " )" << endl;
            cout << "matrice A_local_enh de l'element " << elem_list_enh[ n ] << " =" << endl;
            cout << A_local_enh[ n ] << endl << endl;
        }
    }
    
    /// Construction des vecteurs d_local_enh[ n ] pour chaque element ameliore n du maillage
    /// -------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs d_local_enh" << endl << endl;
    
    Vec< Vec<T> > d_local_enh;
    d_local_enh.resize( elem_list_enh.size() );
    
    for (unsigned n=0;n<elem_list_enh.size();++n) {
        d_local_enh[ n ].resize( nb_unk_local_enh[ n ] );
        d_local_enh[ n ].set( 0. );
    }
    
    Calcul_Elem_Vector_d_enh<T> calcul_elem_vector_d_enh;
    calcul_elem_vector_d_enh.elem_flag_enh = &elem_flag_enh;
    calcul_elem_vector_d_enh.elem_list_enh = &elem_list_enh;
    calcul_elem_vector_d_enh.dep_hat = &dep_hat;
    calcul_elem_vector_d_enh.dep_hat_enh = &dep_hat_enh;
    
    apply( m.elem_list, calcul_elem_vector_d_enh, m, f, d_local_enh );
    
    if ( disp ) {
        for (unsigned n=0;n<elem_list_enh.size();++n) {
            cout << "dimension du vecteur d_local_enh de l'element " << elem_list_enh[ n ] << " = " << nb_unk_local_enh[ n ] << endl;
            cout << "vecteur d_local_enh de l'element " << elem_list_enh[ n ] << " =" << endl;
            cout << d_local_enh[ n ] << endl << endl;
        }
    }
    
    dep_hat_enh.free();
    
    /// Construction de la matrice C_enh
    /// --------------------------------
    
    if ( disp )
        cout << "Construction de la matrice C_enh" << endl << endl;
    
    unsigned nb_eq_f_surf_enh = 0;
    
    for (unsigned k=0;k<face_list_enh.size();++k) {
        for (unsigned d=0;d<dim;d++) {
            if ( face_type[ face_list_enh[ k ] ][ d ] == 2 ) {
                nb_eq_f_surf_enh += node_cpt_face[ face_list_enh[ k ] ];
            }
        }
    }
    
    unsigned cpt_eq_f_surf_enh = 0;
    
    TMatGenSparse C_enh;
    C_enh.resize( nb_eq_f_surf_enh, nb_unk_enh );
    
    Calcul_Global_Matrix_C_enh calcul_global_matrix_C_enh;
    calcul_global_matrix_C_enh.face_type = &face_type;
    calcul_global_matrix_C_enh.face_flag_enh = &face_flag_enh;
    calcul_global_matrix_C_enh.face_list_enh = &face_list_enh;
    
    apply( m.sub_mesh(Number<1>()).elem_list, calcul_global_matrix_C_enh, cpt_eq_f_surf_enh, C_enh );
    
    if ( disp ) {
        cout << "dimension de la matrice C_enh = ( " << nb_eq_f_surf_enh << ", " << nb_unk_enh << " )" << endl;
        cout << "matrice C_enh =" << endl;
        cout << C_enh << endl << endl;
    }
    
    /// Construction du vecteur global q_enh
    /// ------------------------------------
    
    if ( disp )
        cout << "Construction du vecteur q_enh" << endl << endl;
    
    Vec<T> q_enh;
    q_enh.resize( nb_eq_f_surf_enh );
    q_enh.set( 0. );
    
    cpt_eq_f_surf_enh = 0;
    
    Calcul_Global_Vector_q_enh calcul_global_vector_q_enh;
    calcul_global_vector_q_enh.face_type = &face_type;
    calcul_global_vector_q_enh.elem_flag_enh = &elem_flag_enh;
    calcul_global_vector_q_enh.cpt_eq_f_surf_enh = &cpt_eq_f_surf_enh;
    
    apply( m.elem_list, calcul_global_vector_q_enh, m, f, q_enh );
    
    if ( disp ) {
        cout << "dimension du vecteur q_enh = " << nb_eq_f_surf_enh << endl;
        cout << "vecteur q_enh =" << endl;
        cout << q_enh << endl << endl;
    }
    
    face_type.free();
    
    /// Construction des matrices L_local_enh[ n ] pour chaque element ameliore n du maillage
    /// -------------------------------------------------------------------------------------
    
    if ( disp ) {
        cout << "Calcul des vecteurs nb_eq_f_vol_local_enh" << endl << endl;
        cout << "Construction des matrices L_local_enh" << endl << endl;
    }
    
    Vec<unsigned> nb_unk_local_bal;
    nb_unk_local_bal.resize( elem_list_bal.size() );
    
    Vec<unsigned> nb_eq_f_vol_local_enh;
    nb_eq_f_vol_local_enh.resize( elem_list_bal.size() );
    
    Vec< TMatGenSparse > L_local_enh;
    L_local_enh.resize( elem_list_bal.size() );
    
    Calcul_Elem_Matrix_L_enh calcul_elem_matrix_L_enh;
    calcul_elem_matrix_L_enh.face_flag_enh = &face_flag_enh;
    calcul_elem_matrix_L_enh.elem_flag_bal = &elem_flag_bal;
    calcul_elem_matrix_L_enh.elem_list_bal = &elem_list_bal;
    calcul_elem_matrix_L_enh.nb_unk_local_bal = &nb_unk_local_bal;
    calcul_elem_matrix_L_enh.nb_eq_f_vol_local_enh = &nb_eq_f_vol_local_enh;
    
    apply( m.elem_list, calcul_elem_matrix_L_enh, m, f, L_local_enh );
    
    if ( disp ) {
        for (unsigned n=0;n<elem_list_bal.size();++n) {
            cout << "dimension de la matrice L_local_enh de l'element " << elem_list_bal[ n ] << " = ( " << nb_eq_f_vol_local_enh[ n ] << ", " << nb_unk_local_bal[ n ] << " )" << endl;
            cout << "matrice L_local_enh de l'element " << elem_list_bal[ n ] << " =" << endl;
            cout << L_local_enh[ n ] << endl << endl;
        }
    }
    
    /// Construction des vecteurs b_local_enh[ n ] pour chaque element ameliore n du maillage
    /// -------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs b_local_enh" << endl << endl;
    
    Vec< Vec<T> > b_local_enh;
    b_local_enh.resize( elem_list_bal.size() );
    
    Calcul_Elem_Vector_b_enh<T> calcul_elem_vector_b_enh;
    calcul_elem_vector_b_enh.node_list_face = &node_list_face;
    calcul_elem_vector_b_enh.elem_cpt_node = &elem_cpt_node;
    calcul_elem_vector_b_enh.force_fluxes = &force_fluxes;
    calcul_elem_vector_b_enh.elem_flag_bal = &elem_flag_bal;
    calcul_elem_vector_b_enh.elem_list_bal = &elem_list_bal;
    
    apply( m.elem_list, calcul_elem_vector_b_enh, m, f, b_local_enh );
    
    if ( disp ) {
        for (unsigned n=0;n<elem_list_bal.size();++n) {
            cout << "dimension du vecteur b_local_enh de l'element " << elem_list_bal[ n ] << " = " << nb_eq_f_vol_local_enh[ n ] << endl;
            cout << "vecteur b_local_enh de l'element " << elem_list_bal[ n ] << " =" << endl;
            cout << b_local_enh[ n ] << endl << endl;
        }
    }
    
    elem_cpt_node.free();
    
    /// Construction de la matrice globale A_enh
    /// ----------------------------------------
    
    if ( disp )
        cout << "Construction de la matrice A_enh" << endl << endl;
    
    TMatGenSparse A_enh;
    A_enh.resize( nb_unk_enh );
    
    Calcul_Global_Matrix_A_enh calcul_global_matrix_A_enh;
    calcul_global_matrix_A_enh.node_list_face = &node_list_face;
    calcul_global_matrix_A_enh.elem_flag_enh = &elem_flag_enh;
    calcul_global_matrix_A_enh.elem_list_enh = &elem_list_enh;
    calcul_global_matrix_A_enh.face_list_enh = &face_list_enh;
    
    apply( m.elem_list, calcul_global_matrix_A_enh, m, A_local_enh, A_enh );
    
    if ( disp ) {
        cout << "dimension de la matrice A_enh = ( " << nb_unk_enh << ", " << nb_unk_enh << " )" << endl;
        cout << "matrice A_enh =" << endl;
        cout << A_enh << endl << endl;
    }
    
    A_local_enh.free();
    
    /// Construction du vecteur global d_enh
    /// ------------------------------------
    
    if ( disp )
        cout << "Construction du vecteur d_enh" << endl << endl;
    
    Vec<T> d_enh;
    d_enh.resize( nb_unk_enh );
    d_enh.set( 0. );
    
    Calcul_Global_Vector_d_enh calcul_global_vector_d_enh;
    calcul_global_vector_d_enh.node_list_face = &node_list_face;
    calcul_global_vector_d_enh.elem_flag_enh = &elem_flag_enh;
    calcul_global_vector_d_enh.elem_list_enh = &elem_list_enh;
    calcul_global_vector_d_enh.face_list_enh = &face_list_enh;
    
    apply( m.elem_list, calcul_global_vector_d_enh, m, d_local_enh, d_enh );
    
    if ( disp ) {
        cout << "dimension du vecteur d_enh = " << nb_unk_enh << endl;
        cout << "vecteur d_enh =" << endl;
        cout << d_enh << endl << endl;
    }
    
    d_local_enh.free();
    
    /// Construction de la matrice globale L_enh
    /// ----------------------------------------
    
    if ( disp )
        cout << "Construction de la matrice L_enh" << endl << endl;
    
    unsigned nb_eq_f_vol_enh = 0;
    
    for (unsigned n=0;n<elem_list_bal.size();++n) {
        nb_eq_f_vol_enh += nb_eq_f_vol_local_enh[ n ];
    }
    
    TMatGenSparse L_enh;
    L_enh.resize( nb_eq_f_vol_enh, nb_unk_enh );
    
    Calcul_Global_Matrix_L_enh calcul_global_matrix_L_enh;
    calcul_global_matrix_L_enh.node_list_face = &node_list_face;
    calcul_global_matrix_L_enh.elem_flag_bal = &elem_flag_bal;
    calcul_global_matrix_L_enh.elem_list_bal = &elem_list_bal;
    calcul_global_matrix_L_enh.face_flag_enh = &face_flag_enh;
    calcul_global_matrix_L_enh.face_list_enh = &face_list_enh;
    calcul_global_matrix_L_enh.nb_eq_f_vol_local_enh = &nb_eq_f_vol_local_enh;
    
    apply( m.elem_list, calcul_global_matrix_L_enh, m, L_local_enh, L_enh );
    
    if ( disp ) {
        cout << "dimension de la matrice L_enh = ( " << nb_eq_f_vol_enh << ", " << nb_unk_enh << " )" << endl;
        cout << "matrice L_enh =" << endl;
        cout << L_enh << endl << endl;
    }
    
    L_local_enh.free();
    
    /// Construction du vecteur global b_enh
    /// ------------------------------------
    
    if ( disp )
        cout << "Construction du vecteur b_enh" << endl << endl;
    
    Vec<T> b_enh;
    b_enh.resize( nb_eq_f_vol_enh );
    b_enh.set( 0. );
    
    Calcul_Global_Vector_b_enh calcul_global_vector_b_enh;
    calcul_global_vector_b_enh.elem_flag_bal = &elem_flag_bal;
    calcul_global_vector_b_enh.elem_list_bal = &elem_list_bal;
    calcul_global_vector_b_enh.nb_eq_f_vol_local_enh = &nb_eq_f_vol_local_enh;
    
    apply( m.elem_list, calcul_global_vector_b_enh, m, b_local_enh, b_enh );
    
    if ( disp ) {
        cout << "dimension du vecteur b_enh = " << nb_eq_f_vol_enh << endl;
        cout << "vecteur b_enh =" << endl;
        cout << b_enh << endl << endl;
    }
    
    b_local_enh.free();
    nb_eq_f_vol_local_enh.free();
    
    /// Construction de la matrice globale P_enh
    /// ----------------------------------------
    
    if ( disp )
        cout << "Construction de la matrice P_enh" << endl << endl;
    
    unsigned nb_eq_proj_f_surf_enh = 0;
    
    for (unsigned k=0;k<face_list_enh.size();++k) {
        for (unsigned i=0;i<node_cpt_face[ face_list_enh[ k ] ];++i) {
            if ( not correspondance_node_to_vertex_node[ node_list_face[ face_list_enh[ k ] ][ i ] ] ) { // si le ieme noeud de la face k est non sommet
                nb_eq_proj_f_surf_enh += dim;
            }
        }
    }
    
    unsigned cpt_eq_proj_f_surf_enh = 0;
    
    TMatGenSparse P_enh;
    P_enh.resize( nb_eq_proj_f_surf_enh, nb_unk_enh );
    
    Calcul_Global_Matrix_P_enh calcul_global_matrix_P_enh;
    calcul_global_matrix_P_enh.face_flag_enh = &face_flag_enh;
    calcul_global_matrix_P_enh.face_list_enh = &face_list_enh;
    
    apply( m.sub_mesh(Number<1>()).elem_list, calcul_global_matrix_P_enh, cpt_eq_proj_f_surf_enh, P_enh );
    
    if ( disp ) {
        cout << "dimension de la matrice P_enh = ( " << nb_eq_proj_f_surf_enh << ", " << nb_unk_enh << " )" << endl;
        cout << "matrice P_enh =" << endl;
        cout << P_enh << endl << endl;
    }
    
    correspondance_node_to_vertex_node.free();
    node_cpt_face.free();
    node_list_face.free();
    
    /// Construction de la matrice globale C_tot_enh et du vecteur global q_tot_enh
    /// ---------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction de la matrice C_tot_enh et du vecteur q_tot_enh" << endl << endl;
    
    TMatGenSparse C_tot_enh;
    Vec<T> q_tot_enh;
    
    C_tot_enh.resize( nb_eq_f_surf_enh + nb_eq_proj_f_surf_enh + nb_eq_f_vol_enh, nb_unk_enh );
    q_tot_enh.resize( nb_eq_f_surf_enh + nb_eq_proj_f_surf_enh + nb_eq_f_vol_enh );
    q_tot_enh.set( 0. );
    
    Vec<unsigned> vec_unk = range( nb_unk_enh );
    Vec<unsigned> vec_eq_f_surf = range( nb_eq_f_surf_enh );
    Vec<unsigned> vec_eq_proj_f_surf = range( nb_eq_proj_f_surf_enh );
    Vec<unsigned> vec_eq_f_vol = range( nb_eq_f_vol_enh );
    
    Vec<unsigned> vec_eq_f_surf_to_eq_f_vol_plus_eq_proj_f_surf = range( nb_eq_f_surf_enh, nb_eq_f_surf_enh + nb_eq_proj_f_surf_enh );
    Vec<unsigned> vec_eq_f_surf_plus_eq_proj_f_surf_to_eq_f_surf_plus_eq_proj_f_surf_plus_eq_f_vol = range( nb_eq_f_surf_enh + nb_eq_proj_f_surf_enh, nb_eq_f_surf_enh + nb_eq_proj_f_surf_enh + nb_eq_f_vol_enh );
    
    C_tot_enh( vec_eq_f_surf, vec_unk ) = C_enh( vec_eq_f_surf, vec_unk ) * 1.;
    C_tot_enh( vec_eq_f_surf_to_eq_f_vol_plus_eq_proj_f_surf, vec_unk ) = P_enh( vec_eq_proj_f_surf, vec_unk ) * 1.;
    C_tot_enh( vec_eq_f_surf_plus_eq_proj_f_surf_to_eq_f_surf_plus_eq_proj_f_surf_plus_eq_f_vol, vec_unk ) = L_enh( vec_eq_f_vol, vec_unk ) * 1.;
    q_tot_enh[ vec_eq_f_surf ] = q_enh[ vec_eq_f_surf ] * 1.;
    q_tot_enh[ vec_eq_f_surf_plus_eq_proj_f_surf_to_eq_f_surf_plus_eq_proj_f_surf_plus_eq_f_vol ] = b_enh[ vec_eq_f_vol ] * 1.;
    
    if ( disp ) {
        cout << "dimension de la matrice C_tot_enh = ( " << C_tot_enh.nb_rows() << ", " << C_tot_enh.nb_cols() << " )" << endl;
        cout << "matrice C_tot_enh =" << endl;
        cout << C_tot_enh << endl << endl;
        cout << "dimension du vecteur q_tot_enh = " << q_tot_enh.size() << endl;
        cout << "vecteur q_tot_enh =" << endl;
        cout << q_tot_enh << endl << endl;
    }
    
    L_enh.clear();
    b_enh.free();
    C_enh.clear();
    q_enh.free();
    P_enh.clear();
    vec_eq_f_surf.free();
    vec_eq_proj_f_surf.free();
    vec_eq_f_vol.free();
    vec_eq_f_surf_to_eq_f_vol_plus_eq_proj_f_surf.free();
    vec_eq_f_surf_plus_eq_proj_f_surf_to_eq_f_surf_plus_eq_proj_f_surf_plus_eq_f_vol.free();
    
    /// Construction de la matrice globale C_tot_tilde_enh et du vecteur global q_tot_tilde_enh
    /// ---------------------------------------------------------------------------------------
    
    if ( disp ) {
        cout << "Suppression du noyau de la matrice C_tot_enh" << endl << endl;
        cout << "Calcul du vecteur nb_eq_indep_enh" << endl << endl;
    }
    
    Vec<unsigned> eq_indep_enh;
    
    /// Orthonormalisation de Schmidt
    /// -----------------------------
//    Mat<T> C_tot_orth_enh = C_tot_enh;
    
//    TicToc t_orth_schmidt;
//    t_orth_schmidt.start();
    
//    orthonormalisation_schmidt_row( C_tot_orth_enh, eq_indep_enh, 1.e-10 );
    
//    t_orth_schmidt.stop();
//    if ( disp )
//        cout << "temps de calcul de l'orthonormalisation de schmidt de la matrice C_tot_enh associee au pb global de minimisation = " << t_orth_schmidt.res << " s" << endl << endl;
    
//    C_tot_orth_enh.clear();
    
    /// Factorisation LDL avec detection des pivots nuls
    /// ------------------------------------------------
    eq_indep_enh = range( nb_eq_f_surf_enh + nb_eq_proj_f_surf_enh + nb_eq_f_vol_enh );
    
    TicToc t_detect_pivots;
    t_detect_pivots.start();
    
    TMatSymSparse C_tot_2_enhancement = C_tot_enh * trans( C_tot_enh );
    
    LDL_solver ls_enh;
    Vec< Vec<T> > Ker_enh;
    Vec<int> Pivots_enh;
    ls_enh.get_factorization( C_tot_2_enhancement, Ker_enh, Pivots_enh );
    
    t_detect_pivots.stop();
    if ( disp )
        cout << "temps de calcul de la detection des pivots nuls de la matrice C_tot_enh associee au pb global de minimisation = " << t_detect_pivots.res << " s" << endl << endl;
    
    if ( Pivots_enh.size() ) { // si le nombre de pivots nuls est non nul
        for (unsigned n=0;n<Pivots_enh.size();++n) {
            remove_if( eq_indep_enh, _1 == _2, Pivots_enh[ n ] );
        }
    }
    
    if ( disp )
        cout << "indices des pivots nuls dans la matrice C_tot_enh = " << Pivots_enh << endl;
    
    C_tot_2_enhancement.clear();
    Ker_enh.free();
    Pivots_enh.free();
    
    unsigned nb_eq_indep_enh = eq_indep_enh.size();
    
    if ( disp ) {
        cout << "nombre de lignes independantes dans la matrice C_tot_enh = " << nb_eq_indep_enh << endl;
        cout << "indices des lignes independantes dans la matrice C_tot_enh = " << eq_indep_enh << endl << endl;
    }
    
    if ( disp )
        cout << "Construction de la matrice C_tot_tilde_enh et du vecteur q_tot_tilde_enh" << endl << endl;
    
    Mat<T, Gen<>, SparseLine<Row> > C_tot_tilde_enh;
    Vec<T> q_tot_tilde_enh;
    
    C_tot_tilde_enh.resize( nb_eq_indep_enh, nb_unk_enh );
    q_tot_tilde_enh.resize( nb_eq_indep_enh );
    q_tot_tilde_enh.set( 0. );
    
    for (unsigned n=0;n<nb_eq_indep_enh;++n) {
        C_tot_tilde_enh.row( n ) = C_tot_enh.row( eq_indep_enh[ n ] );
        q_tot_tilde_enh[ n ] += q_tot_enh[ eq_indep_enh[ n ] ];
    }
    
    if ( disp ) {
        cout << "dimension de la matrice C_tot_tilde_enh = ( " << C_tot_tilde_enh.nb_rows() << ", " << C_tot_tilde_enh.nb_cols() << " )" << endl;
        cout << "matrice C_tot_tilde_enh =" << endl;
        cout << C_tot_tilde_enh << endl << endl;
        cout << "dimension du vecteur q_tot_tilde_enh = " << q_tot_tilde_enh.size() << endl;
        cout << "vecteur q_tot_tilde_enh =" << endl;
        cout << q_tot_tilde_enh << endl << endl;
    }
    
    C_tot_enh.clear();
    q_tot_enh.free();
    eq_indep_enh.free();
    
    /// Construction de la matrice globale K_enh et du vecteur global F_enh
    /// -------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction de la matrice K_enh et du vecteur F_enh" << endl << endl;
    
    TMatGenSparse K_enh;
    Vec<T> F_enh;
    Vec<T> U_enh;
    
    TicToc t_construct_K_F_enh;
    t_construct_K_F_enh.start();
    
    K_enh.resize( nb_unk_enh + nb_eq_indep_enh );
//    K_enh.set( 0. );
    F_enh.resize( nb_unk_enh + nb_eq_indep_enh );
    F_enh.set( 0. );
    U_enh.resize( nb_unk_enh + nb_eq_indep_enh );
    U_enh.set( 0. );
    
    Vec<unsigned> vec_eq_indep = range( nb_eq_indep_enh );
    Vec<unsigned> vec_unk_to_unk_plus_eq_indep = range( nb_unk_enh, nb_unk_enh + nb_eq_indep_enh );
    TMatGenSparse trans_C_tot_tilde_enh = trans( C_tot_tilde_enh );
    K_enh( vec_unk, vec_unk ) = A_enh( vec_unk, vec_unk ) * 1.;
    K_enh( vec_unk_to_unk_plus_eq_indep, vec_unk ) = C_tot_tilde_enh( vec_eq_indep, vec_unk ) * 1.;
    K_enh( vec_unk, vec_unk_to_unk_plus_eq_indep ) = trans_C_tot_tilde_enh( vec_unk, vec_eq_indep ) * 1.;
    F_enh[ vec_unk ] = d_enh[ vec_unk ] * 1.;
    F_enh[ vec_unk_to_unk_plus_eq_indep ] = q_tot_tilde_enh[ vec_eq_indep ] * 1.;
    
    t_construct_K_F_enh.stop();
    if ( disp )
        cout << "temps de calcul du remplissage de la matrice K_enh et du vecteur F_enh associes au pb global de minimisation = " << t_construct_K_F_enh.res << " s" << endl << endl;
    
    if ( disp ) {
        cout << "dimension de la matrice K_enh = ( " << K_enh.nb_rows() << ", " << K_enh.nb_cols() << " ) "<< endl;
        cout << "matrice K_enh =" << endl;
        cout << K_enh << endl << endl;
        cout << "dimension du vecteur F_enh = " << F_enh.size() << endl;
        cout << "vecteur F_enh =" << endl;
        cout << F_enh << endl << endl;
    }
    
    A_enh.clear();
    d_enh.free();
    C_tot_tilde_enh.clear();
    trans_C_tot_tilde_enh.clear();
    q_tot_tilde_enh.free();
    vec_eq_indep.free();
    vec_unk_to_unk_plus_eq_indep.free();
    
    /// Resolution du probleme global de minimsation K_enh * U_enh = F_enh
    /// Construction du vecteur U_enh
    /// ------------------------------------------------------------------
    
    if ( disp ) {
        cout << "Resolution du pb global de minimsation K_enh * U_enh = F_enh a " << nb_unk_enh << " inconnues et " << nb_eq_indep_enh << " equations independantes imposees" << endl;
        cout << "Construction du vecteur U_enh" << endl << endl;
    }
    
    TicToc t_solve_minimization_enh;
    t_solve_minimization_enh.start();
    
    if ( solver_minimisation == "LDL" ) {
        #ifdef WITH_LDL
        TMatSymSparse K_LDL = K_enh;
        U_enh = F_enh;
        LDL_solver ls;
        Vec< Vec<T> > Ker;
        Vec<int> Pivots;
        ls.get_factorization( K_LDL, Ker, Pivots );
        ls.solve( U_enh );
        #endif
    }
    else if ( solver_minimisation == "UMFPACK" ) {
        #ifdef WITH_UMFPACK
        TMatSymSparseUMFPACK K_UMFPACK = K_enh;
        K_UMFPACK.get_factorization();
        U_enh = K_UMFPACK.solve( F_enh );
        #endif
    }
    else if ( solver_minimisation == "Inv" ) {
        TMatSymSparse K_Inv = K_enh;
        U_enh = inv( K_Inv ) * F_enh;
    }
    else if ( solver_minimisation == "LUFactorize" ) {
        TMatSymSparse K_sym = K_enh;
        TMatGenSparseLU K_LU = K_sym;
//        Vec<int> vector_permutation;
        lu_factorize( K_LU/*, vector_permutation*/ );
        solve_using_lu_factorize( K_LU, /*vector_permutation,*/ F_enh, U_enh );
    }
    else {
        cerr << "Bing. Error : solveur " << solver_minimisation << " pour la minimisation non implemente" << endl << endl;
    }
    
    t_solve_minimization_enh.stop();
    if ( disp )
        cout << "temps de calcul de la resolution du pb global de minimisation = " << t_solve_minimization_enh.res << " s" << endl << endl;
    
    if ( disp ) {
        cout << "dimension du vecteur U_enh = " << U_enh.size() << endl;
        cout << "vecteur U_enh =" << endl;
        cout << U_enh << endl << endl;
    }
    if ( verif_solver_minimisation_enhancement ) {
        if ( disp )
            cout << "Verification de la resolution du pb global de minimsation : tolerance = " << tol_solver_minimisation_enhancement << endl << endl;
        T residual = norm_2( K_enh * U_enh - F_enh );
        T b = norm_2( F_enh );
        if ( residual / b > tol_solver_minimisation_enhancement ) {
//            cout << "residu K_enh * U_enh - F_enh :" << endl;
//            cout << K_enh * U_enh - F_enh << endl;
            cout << "norme du residu = " << residual << endl;
            cout << "norme du residu relatif = " << residual / b << endl << endl;
        }
    }
    
    K_enh.clear();
    F_enh.free();
    
    /// Construction du vecteur force_fluxes_enhancement_global
    /// -------------------------------------------------------
    
    if ( disp )
        cout << "Construction du vecteur force_fluxes_enhancement_global" << endl << endl;
    
    Vec<T> force_fluxes_enhancement_global;
    
    force_fluxes_enhancement_global.resize( nb_unk_enh );
    force_fluxes_enhancement_global.set( 0. );
    
    force_fluxes_enhancement_global[ vec_unk ] = U_enh[ vec_unk ] * 1.;
    
    if ( disp ) {
        cout << "dimension du vecteur force_fluxes_enhancement_global = " << nb_unk_enh << endl;
        cout << "vecteur force_fluxes_enhancement_global =" << endl;
        cout << force_fluxes_enhancement_global << endl << endl;
    }
    
    vec_unk.free();
    U_enh.free();
    
    /// Construction des vecteurs force_fluxes_enhancement[ k ][ d ] pour chaque face k du maillage et chaque direction d
    /// -----------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs force_fluxes_enhancement" << endl << endl;
    
    Vec< Vec< Vec<T> > > force_fluxes_enhancement;
    
    force_fluxes_enhancement.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        force_fluxes_enhancement[ k ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            force_fluxes_enhancement[ k ][ d ].resize( m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
            force_fluxes_enhancement[ k ][ d ].set( 0. );
        }
    }
    
    for (unsigned k=0;k<face_list_enh.size();++k) {
        for (unsigned d=0;d<dim;++d) {
            for (unsigned j=0;j<m.sub_mesh(Number<1>()).elem_list[face_list_enh[ k ]]->nb_nodes_virtual();++j) {
                force_fluxes_enhancement[ face_list_enh[ k ] ][ d ][ j ] += force_fluxes_enhancement_global[ k * m.sub_mesh(Number<1>()).elem_list[face_list_enh[ k ]]->nb_nodes_virtual() * dim + j * dim + d ];
            }
        }
    }
    
    if ( disp ) {
        for (unsigned k=0;k<face_list_enh.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur des densites d'effort ameliorees associe a la face " << face_list_enh[ k ] << " dans la direction " << d << " = " << m.sub_mesh(Number<1>()).elem_list[face_list_enh[ k ]]->nb_nodes_virtual() << endl;
                cout << "vecteur des densites d'effort ameliorees associe a la face " << face_list_enh[ k ] << " dans la direction " << d << " =" << endl;
                cout << force_fluxes_enhancement[ face_list_enh[ k ] ][ d ] << endl << endl;
            }
        }
    }
    
    /// Construction des vecteurs force_fluxes_enhancement[ k ][ d ] pour chaque face k du maillage et chaque direction d
    /// -----------------------------------------------------------------------------------------------------------------
    
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        for (unsigned d=0;d<dim;++d) {
            for (unsigned j=0;j<m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual();++j) {
                force_fluxes[ k ][ d ][ j ] += force_fluxes_enhancement[ k ][ d ][ j ];
            }
        }
    }
    
    if ( disp ) {
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur des densites d'effort standard + ameliorees associe a la face " << k << " dans la direction " << d << " = " << m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << endl;
                cout << "vecteur des densites d'effort standard + ameliorees associe a la face " << k << " dans la direction " << d << " =" << endl;
                cout << force_fluxes[ k ][ d ] << endl << endl;
            }
        }
    }
    
    t_force_fluxes_enh.stop();
    cout << "temps de calcul de la construction amelioree des densites d'effort = " << t_force_fluxes_enh.res << " s" << endl << endl;
}

#endif // Construct_enhanced_force_fluxes_EET_EESPT_h
