//
// C++ Interface: Construct_enhanced_force_fluxes_EET_EESPT
//
// Description: construction des densites d'effort ameliorees par les methodes EET et EESPT
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
#include "../GEOMETRY/Calcul_geometry.h"
#include "../TOOLS/Algebre.h"

using namespace LMT;
using namespace std;

/// Construction amelioree des densites d'effort par les methodes EET et EESPT
///---------------------------------------------------------------------------
template<class TM, class TF, class T>
void construct_enhanced_force_fluxes_EET_EESPT( TM &m, const TF &f, const string &method, const Vec<bool> &flag_elem_enh, const Vec<bool> &flag_face_enh, const Vec<bool> &flag_elem_bal, const Vec<unsigned> &list_elems_enh, const Vec<unsigned> &list_faces_enh, const Vec<unsigned> &list_elems_bal, Vec< Mat<T, Sym<> > > &K_hat, const Vec< Vec<T> > &dep_hat, Vec< Vec< Vec<T> > > &vec_force_fluxes, const string &solver, const string &solver_minimisation, const bool &debug_method_enhancement, const bool &verif_solver_enhancement, const T &tol_solver_enhancement, const bool &verif_solver_minimisation_enhancement, const T &tol_solver_minimisation_enhancement, const bool &debug_geometry, const bool &debug_force_fluxes_enhancement ) {

    static const unsigned dim = TM::dim;

    TicToc t_construct_force_fluxes_enh;
    t_construct_force_fluxes_enh.start();

    Vec<unsigned> cpt_nodes_face;
    Vec< Vec<unsigned> > list_nodes_face;
    construct_nodes_connected_to_face( m, cpt_nodes_face, list_nodes_face, debug_geometry );

    Vec<unsigned> cpt_elems_node;
    Vec< Vec<unsigned> > list_elems_node;
    construct_elems_connected_to_node( m, cpt_elems_node, list_elems_node, debug_geometry );

    list_elems_node.free();

    Vec<bool> correspondance_node_to_vertex_node;
    Vec<unsigned> connect_node_to_vertex_node;
    unsigned nb_vertex_nodes = match_node_to_vertex_node( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, debug_geometry );

    connect_node_to_vertex_node.free();

    Vec< Vec<unsigned> > type_face;
    construct_type_face( m, f, type_face, debug_geometry );

    cout << "--------------------------------------------" << endl;
    cout << "Construction amelioree des densites d'effort" << endl;
    cout << "--------------------------------------------" << endl << endl;

    unsigned nb_unk_enh = 0;

    for (unsigned k=0;k<list_faces_enh.size();++k) {
        nb_unk_enh += m.sub_mesh(Number<1>()).elem_list[ list_faces_enh[ k ] ]->nb_nodes_virtual() * dim;
    }

    /// Construction des vecteurs F_hat_enh[ n ] pour chaque element ameliore n du maillage
    ///------------------------------------------------------------------------------------

    cout << "Calcul des vecteurs nb_unk_local_enh" << endl << endl;

    cout << "Construction des vecteurs F_hat_enh" << endl << endl;

    Vec<unsigned> nb_unk_local_enh;
    nb_unk_local_enh.resize( list_elems_enh.size() );

    Vec< Vec< Vec<T> > > F_hat_enh;
    F_hat_enh.resize( list_elems_enh.size() );

    Calcul_Elem_Vector_F_hat_enh calc_elem_vector_F_hat_enh;
    calc_elem_vector_F_hat_enh.flag_elem_enh = &flag_elem_enh;
    calc_elem_vector_F_hat_enh.list_elems_enh = &list_elems_enh;
    calc_elem_vector_F_hat_enh.nb_unk_local_enh = &nb_unk_local_enh;

    apply( m.elem_list, calc_elem_vector_F_hat_enh, m, f, F_hat_enh );

    if ( debug_method_enhancement ) {
        for (unsigned n=0;n<list_elems_enh.size();++n) {
            cout << "vecteur F_hat_enh de l'element " << list_elems_enh[ n ] << " :" << endl << endl;
            for (unsigned i=0;i<nb_unk_local_enh[ n ];++i) {
                cout << "avec cas de charge " << i << " :" << endl;
                cout << F_hat_enh[ n ][ i ] << endl << endl;
            }
        }
    }

    /// Construction des vecteurs dep_hat_enh[ n ] pour chaque element ameliore n du maillage
    ///--------------------------------------------------------------------------------------

    cout << "Construction des vecteurs dep_hat_enh" << endl << endl;

    Vec< Vec< Vec<T> > > dep_hat_enh;
    dep_hat_enh.resize( list_elems_enh.size() );

    for (unsigned n=0;n<list_elems_enh.size();++n) {
        dep_hat_enh[ n ].resize( nb_unk_local_enh[ n ] );
        for (unsigned i=0;i<nb_unk_local_enh[ n ];++i) {
            if ( solver == "CholMod" ) {
#ifdef WITH_CHOLMOD
                Mat<T, Sym<>, SparseLine<> > K_hat_sym = K_hat[ list_elems_enh[ n ] ];
                Mat<T, Sym<>, SparseCholMod > K_hat_CholMod = K_hat_sym;
                K_hat_CholMod.get_factorization();
                dep_hat_enh[ n ][ i ] = K_hat_CholMod.solve( F_hat_enh[ n ][ i ] );
#endif
            }
            else if ( solver == "LDL" ) {
#ifdef WITH_LDL
                Mat<T, Sym<>, SparseLine<> > K_hat_LDL = K_hat[ list_elems_enh[ n ] ];
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
                Mat<T, Gen<>, SparseUMFPACK > K_hat_UMFPACK = K_hat[ list_elems_enh[ n ] ];
                K_hat_UMFPACK.get_factorization();
                dep_hat_enh[ n ][ i ] = K_hat_UMFPACK.solve( F_hat_enh[ n ][ i ] );
#endif
            }
            else if ( solver == "CholFactorize" ) {
                Mat<T, Sym<>, SparseLine<> > K_hat_Chol = K_hat[ list_elems_enh[ n ] ];
                chol_factorize( K_hat_Chol );
                solve_using_chol_factorize( K_hat_Chol, F_hat_enh[ n ][ i ], dep_hat_enh[ n ][ i ] );
            }
            else if ( solver == "LUFactorize" ) {
                Mat<T, Sym<>, SparseLine<> > K_hat_sym = K_hat[ list_elems_enh[ n ] ];
                Mat<T, Gen<>, SparseLU > K_hat_LU = K_hat_sym;
//                Vec<int> vector_permutation;
                lu_factorize( K_hat_LU/*, vector_permutation*/ );
                solve_using_lu_factorize( K_hat_LU, /*vector_permutation,*/ F_hat_enh[ n ][ i ], dep_hat_enh[ n ][ i ] );
            }
            else if ( solver == "Inv" ) {
                Mat<T, Sym<>, SparseLine<> > K_hat_Inv = K_hat[ list_elems_enh[ n ] ];
                dep_hat_enh[ n ][ i ] = inv( K_hat_Inv ) * F_hat_enh[ n ][ i ];
            }
            else {
                cerr << "Bing. Error : type de solveur pour la resolution locaux non implemente" << endl << endl;
            }
        }
    }

    if ( debug_method_enhancement ) {
        for (unsigned n=0;n<list_elems_enh.size();++n) {
            cout << "vecteur dep_hat_enh de l'element " << list_elems_enh[ n ] << " :" << endl << endl;
            for (unsigned i=0;i<nb_unk_local_enh[ n ];++i) {
                cout << "avec cas de charge " << i << " :" << endl;
                cout << dep_hat_enh[ n ][ i ] << endl << endl;
            }
        }
    }
    if ( verif_solver_enhancement ) {
        cout << "Verification de la resolution des problemes locaux associes a la partie amelioree des densites d'effort pour la technique " << method << " : tolerance = " << tol_solver_enhancement << endl << endl;
        for (unsigned n=0;n<list_elems_enh.size();++n) {
            for (unsigned i=0;i<nb_unk_local_enh[ n ];++i) {
                if ( norm_2( K_hat[ list_elems_enh[ n ] ] * dep_hat_enh[ n ][ i ] - F_hat_enh[ n ][ i ] ) / norm_2( F_hat_enh[ n ][ i ] ) > tol_solver_enhancement ) {
                    cout << "residu pour l'element " << list_elems_enh[ n ] << " avec cas de charge " << i << " :" << endl << endl;
                    cout << "K_hat * dep_hat_enh - F_hat_enh =" << endl;
                    cout << K_hat[ list_elems_enh[ n ] ] * dep_hat_enh[ n ][ i ] - F_hat_enh[ n ][ i ] << endl << endl;
                    cout << "norme 2 du residu = " << norm_2( K_hat[ list_elems_enh[ n ] ] * dep_hat_enh[ n ][ i ] - F_hat_enh[ n ][ i ] ) << endl << endl;
                    cout << "norme 2 du residu relatif = " << norm_2( K_hat[ list_elems_enh[ n ] ] * dep_hat_enh[ n ][ i ] - F_hat_enh[ n ][ i ] ) / norm_2( F_hat_enh[ n ][ i ] ) << endl << endl;
                }
            }
        }
    }

    F_hat_enh.free();

    /// Construction des matrices A_local_enh[ n ] pour chaque element ameliore n du maillage
    ///--------------------------------------------------------------------------------------

    cout << "Construction des matrices A_local_enh" << endl << endl;

    Vec< Mat<T, Gen<>, SparseLine<> > > A_local_enh;
    A_local_enh.resize( list_elems_enh.size() );

    for (unsigned n=0;n<list_elems_enh.size();++n) {
        A_local_enh[ n ].resize( nb_unk_local_enh[ n ] );
    }

    Calcul_Elem_Matrix_A_enh<T> calcul_elem_matrix_A_enh;
    calcul_elem_matrix_A_enh.flag_elem_enh = &flag_elem_enh;
    calcul_elem_matrix_A_enh.list_elems_enh = &list_elems_enh;
    calcul_elem_matrix_A_enh.dep_hat_enh = &dep_hat_enh;

    apply( m.elem_list, calcul_elem_matrix_A_enh, m, f, A_local_enh );

    if ( debug_method_enhancement ) {
        for (unsigned n=0;n<list_elems_enh.size();++n) {
            cout << "dimension de la matrice A_local_enh de l'element " << list_elems_enh[ n ] << " : ( " << nb_unk_local_enh[ n ] << ", " << nb_unk_local_enh[ n ] << " )" << endl;
            cout << "matrice A_local_enh de l'element " << list_elems_enh[ n ] << " :" << endl;
            cout << A_local_enh[ n ] << endl << endl;
        }
    }

    /// Construction des vecteurs d_local_enh[ n ] pour chaque element ameliore n du maillage
    ///--------------------------------------------------------------------------------------

    cout << "Construction des vecteurs d_local_enh" << endl << endl;

    Vec< Vec<T> > d_local_enh;
    d_local_enh.resize( list_elems_enh.size() );

    for (unsigned n=0;n<list_elems_enh.size();++n) {
        d_local_enh[ n ].resize( nb_unk_local_enh[ n ] );
        d_local_enh[ n ].set( 0. );
    }

    Calcul_Elem_Vector_d_enh<T> calcul_elem_vector_d_enh;
    calcul_elem_vector_d_enh.flag_elem_enh = &flag_elem_enh;
    calcul_elem_vector_d_enh.list_elems_enh = &list_elems_enh;
    calcul_elem_vector_d_enh.dep_hat = &dep_hat;
    calcul_elem_vector_d_enh.dep_hat_enh = &dep_hat_enh;

    apply( m.elem_list, calcul_elem_vector_d_enh, m, f, d_local_enh );

    if ( debug_method_enhancement ) {
        for (unsigned n=0;n<list_elems_enh.size();++n) {
            cout << "dimension du vecteur d_local_enh de l'element " << list_elems_enh[ n ] << " : " << nb_unk_local_enh[ n ] << endl;
            cout << "vecteur d_local_enh de l'element " << list_elems_enh[ n ] << " :" << endl;
            cout << d_local_enh[ n ] << endl << endl;
        }
    }

    dep_hat_enh.free();

    /// Construction de la matrice C_enh
    ///---------------------------------

    cout << "Construction de la matrice C_enh" << endl << endl;

    unsigned nb_eq_f_surf_enh = 0;

    for (unsigned k=0;k<list_faces_enh.size();++k) {
        for (unsigned d=0;d<dim;d++) {
            if ( type_face[ list_faces_enh[ k ] ][ d ] == 2 ) {
                nb_eq_f_surf_enh += cpt_nodes_face[ list_faces_enh[ k ] ];
            }
        }
    }

    unsigned cpt_eq_f_surf_enh = 0;

    Mat<T, Gen<>, SparseLine<> > C_enh;
    C_enh.resize( nb_eq_f_surf_enh, nb_unk_enh );

    Calcul_Global_Matrix_C_enh calcul_global_matrix_C_enh;
    calcul_global_matrix_C_enh.type_face = &type_face;
    calcul_global_matrix_C_enh.flag_face_enh = &flag_face_enh;
    calcul_global_matrix_C_enh.list_faces_enh = &list_faces_enh;

    apply( m.sub_mesh(Number<1>()).elem_list, calcul_global_matrix_C_enh, cpt_eq_f_surf_enh, C_enh );

    if ( debug_method_enhancement ) {
        cout << "dimension de la matrice C_enh : ( " << nb_eq_f_surf_enh << ", " << nb_unk_enh << " )" << endl;
        cout << "matrice C_enh :" << endl;
        cout << C_enh << endl << endl;
    }

    /// Construction du vecteur global q_enh
    ///-------------------------------------

    cout << "Construction du vecteur q_enh" << endl << endl;

    Vec<T> q_enh;
    q_enh.resize( nb_eq_f_surf_enh );
    q_enh.set( 0. );

    cpt_eq_f_surf_enh = 0;

    Calcul_Global_Vector_q_enh calcul_global_vector_q_enh;
    calcul_global_vector_q_enh.type_face = &type_face;
    calcul_global_vector_q_enh.flag_elem_enh = &flag_elem_enh;
    calcul_global_vector_q_enh.cpt_eq_f_surf_enh = &cpt_eq_f_surf_enh;

    apply( m.elem_list, calcul_global_vector_q_enh, m, f, q_enh );

    if ( debug_method_enhancement ) {
        cout << "dimension du vecteur q_enh : " << nb_eq_f_surf_enh << endl;
        cout << "vecteur q_enh :" << endl;
        cout << q_enh << endl << endl;
    }

    type_face.free();
    
    /// Construction des matrices L_local_enh[ n ] pour chaque element ameliore n du maillage
    ///--------------------------------------------------------------------------------------

    cout << "Calcul des vecteurs nb_eq_f_vol_local_enh" << endl << endl;

    cout << "Construction des matrices L_local_enh" << endl << endl;

    Vec<unsigned> nb_unk_local_bal;
    nb_unk_local_bal.resize( list_elems_bal.size() );

    Vec<unsigned> nb_eq_f_vol_local_enh;
    nb_eq_f_vol_local_enh.resize( list_elems_bal.size() );

    Vec< Mat<T, Gen<>, SparseLine<> > > L_local_enh;
    L_local_enh.resize( list_elems_bal.size() );

    Calcul_Elem_Matrix_L_enh calcul_elem_matrix_L_enh;
    calcul_elem_matrix_L_enh.flag_face_enh = &flag_face_enh;
    calcul_elem_matrix_L_enh.flag_elem_bal = &flag_elem_bal;
    calcul_elem_matrix_L_enh.list_elems_bal = &list_elems_bal;
    calcul_elem_matrix_L_enh.nb_unk_local_bal = &nb_unk_local_bal;
    calcul_elem_matrix_L_enh.nb_eq_f_vol_local_enh = &nb_eq_f_vol_local_enh;

    apply( m.elem_list, calcul_elem_matrix_L_enh, m, f, L_local_enh );

    if ( debug_method_enhancement ) {
        for (unsigned n=0;n<list_elems_bal.size();++n) {
            cout << "dimension de la matrice L_local_enh de l'element " << list_elems_bal[ n ] << " : ( " << nb_eq_f_vol_local_enh[ n ] << ", " << nb_unk_local_bal[ n ] << " )" << endl;
            cout << "matrice L_local_enh de l'element " << list_elems_bal[ n ] << " :" << endl;
            cout << L_local_enh[ n ] << endl << endl;
        }
    }

    /// Construction des vecteurs b_local_enh[ n ] pour chaque element ameliore n du maillage
    ///--------------------------------------------------------------------------------------

    cout << "Construction des vecteurs b_local_enh" << endl << endl;

    Vec< Vec<T> > b_local_enh;
    b_local_enh.resize( list_elems_bal.size() );

    Calcul_Elem_Vector_b_enh<T> calcul_elem_vector_b_enh;
    calcul_elem_vector_b_enh.list_nodes_face = &list_nodes_face;
    calcul_elem_vector_b_enh.cpt_elems_node = &cpt_elems_node;
    calcul_elem_vector_b_enh.vec_force_fluxes = &vec_force_fluxes;
    calcul_elem_vector_b_enh.flag_elem_bal = &flag_elem_bal;
    calcul_elem_vector_b_enh.list_elems_bal = &list_elems_bal;

    apply( m.elem_list, calcul_elem_vector_b_enh, m, f, b_local_enh );

    if ( debug_method_enhancement ) {
        for (unsigned n=0;n<list_elems_bal.size();++n) {
            cout << "dimension du vecteur b_local_enh de l'element " << list_elems_bal[ n ] << " : " << nb_eq_f_vol_local_enh[ n ] << endl;
            cout << "vecteur b_local_enh de l'element " << list_elems_bal[ n ] << " :" << endl;
            cout << b_local_enh[ n ] << endl << endl;
        }
    }

    cpt_elems_node.free();

    /// Construction de la matrice globale A_enh
    ///-----------------------------------------

    cout << "Construction de la matrice A_enh" << endl << endl;

    Mat<T, Gen<>, SparseLine<> > A_enh;
    A_enh.resize( nb_unk_enh );

    Calcul_Global_Matrix_A_enh calcul_global_matrix_A_enh;
    calcul_global_matrix_A_enh.list_nodes_face = &list_nodes_face;
    calcul_global_matrix_A_enh.flag_elem_enh = &flag_elem_enh;
    calcul_global_matrix_A_enh.list_elems_enh = &list_elems_enh;
    calcul_global_matrix_A_enh.list_faces_enh = &list_faces_enh;

    apply( m.elem_list, calcul_global_matrix_A_enh, m, A_local_enh, A_enh );

    if ( debug_method_enhancement ) {
        cout << "dimension de la matrice A_enh : ( " << nb_unk_enh << ", " << nb_unk_enh << " )" << endl;
        cout << "matrice A_enh" << " :" << endl;
        cout << A_enh << endl << endl;
    }

    A_local_enh.free();

    /// Construction du vecteur global d_enh
    ///-------------------------------------

    cout << "Construction du vecteur d_enh" << endl << endl;

    Vec<T> d_enh;
    d_enh.resize( nb_unk_enh );
    d_enh.set( 0. );

    Calcul_Global_Vector_d_enh calcul_global_vector_d_enh;
    calcul_global_vector_d_enh.list_nodes_face = &list_nodes_face;
    calcul_global_vector_d_enh.flag_elem_enh = &flag_elem_enh;
    calcul_global_vector_d_enh.list_elems_enh = &list_elems_enh;
    calcul_global_vector_d_enh.list_faces_enh = &list_faces_enh;

    apply( m.elem_list, calcul_global_vector_d_enh, m, d_local_enh, d_enh );

    if ( debug_method_enhancement ) {
        cout << "dimension du vecteur d_enh : " << nb_unk_enh << endl;
        cout << "vecteur d_enh" << " :" << endl;
        cout << d_enh << endl << endl;
    }

    d_local_enh.free();

    /// Construction de la matrice globale L_enh
    ///-----------------------------------------

    cout << "Construction de la matrice L_enh" << endl << endl;

    unsigned nb_eq_f_vol_enh = 0;

    for (unsigned n=0;n<list_elems_bal.size();++n) {
        nb_eq_f_vol_enh += nb_eq_f_vol_local_enh[ n ];
    }

    Mat<T, Gen<>, SparseLine<> > L_enh;
    L_enh.resize( nb_eq_f_vol_enh, nb_unk_enh );

    Calcul_Global_Matrix_L_enh calcul_global_matrix_L_enh;
    calcul_global_matrix_L_enh.list_nodes_face = &list_nodes_face;
    calcul_global_matrix_L_enh.flag_elem_bal = &flag_elem_bal;
    calcul_global_matrix_L_enh.list_elems_bal = &list_elems_bal;
    calcul_global_matrix_L_enh.flag_face_enh = &flag_face_enh;
    calcul_global_matrix_L_enh.list_faces_enh = &list_faces_enh;
    calcul_global_matrix_L_enh.nb_eq_f_vol_local_enh = &nb_eq_f_vol_local_enh;

    apply( m.elem_list, calcul_global_matrix_L_enh, m, L_local_enh, L_enh );

    if ( debug_method_enhancement ) {
        cout << "dimension de la matrice L_enh : ( " << nb_eq_f_vol_enh << ", " << nb_unk_enh << " )" << endl;
        cout << "matrice L_enh :" << endl;
        cout << L_enh << endl << endl;
    }

    L_local_enh.free();

    /// Construction du vecteur global b_enh
    ///-------------------------------------

    cout << "Construction du vecteur b_enh" << endl << endl;

    Vec<T> b_enh;
    b_enh.resize( nb_eq_f_vol_enh );
    b_enh.set( 0. );

    Calcul_Global_Vector_b_enh calcul_global_vector_b_enh;
    calcul_global_vector_b_enh.flag_elem_bal = &flag_elem_bal;
    calcul_global_vector_b_enh.list_elems_bal = &list_elems_bal;
    calcul_global_vector_b_enh.nb_eq_f_vol_local_enh = &nb_eq_f_vol_local_enh;

    apply( m.elem_list, calcul_global_vector_b_enh, m, b_local_enh, b_enh );

    if ( debug_method_enhancement ) {
        cout << "dimension du vecteur b_enh : " << nb_eq_f_vol_enh << endl;
        cout << "vecteur b_enh :" << endl;
        cout << b_enh << endl << endl;
    }

    b_local_enh.free();
    nb_eq_f_vol_local_enh.free();

    /// Construction de la matrice globale P_enh
    ///-----------------------------------------

    cout << "Construction de la matrice P_enh" << endl << endl;

    unsigned nb_eq_proj_f_surf_enh = 0;

    for (unsigned k=0;k<list_faces_enh.size();++k) {
        for (unsigned i=0;i<cpt_nodes_face[ list_faces_enh[ k ] ];++i) {
            if ( correspondance_node_to_vertex_node[ list_nodes_face[ list_faces_enh[ k ] ][ i ] ] == 0 ) { // si le ieme noeud de la face k est non sommet
                nb_eq_proj_f_surf_enh += dim;
            }
        }
    }

    unsigned cpt_eq_proj_f_surf_enh = 0;

    Mat<T, Gen<>, SparseLine<> > P_enh;
    P_enh.resize( nb_eq_proj_f_surf_enh, nb_unk_enh );

    Calcul_Global_Matrix_P_enh calcul_global_matrix_P_enh;
    calcul_global_matrix_P_enh.flag_face_enh = &flag_face_enh;
    calcul_global_matrix_P_enh.list_faces_enh = &list_faces_enh;

    apply( m.sub_mesh(Number<1>()).elem_list, calcul_global_matrix_P_enh, cpt_eq_proj_f_surf_enh, P_enh );

    if ( debug_method_enhancement ) {
        cout << "dimension de la matrice P_enh : ( " << nb_eq_proj_f_surf_enh << ", " << nb_unk_enh << " )" << endl;
        cout << "matrice P_enh :" << endl;
        cout << P_enh << endl << endl;
    }

    correspondance_node_to_vertex_node.free();
    cpt_nodes_face.free();
    list_nodes_face.free();

    /// Construction de la matrice globale C_tot_enh et du vecteur global q_tot_enh
    ///----------------------------------------------------------------------------

    cout << "Construction de la matrice C_tot_enh et du vecteur q_tot_enh" << endl << endl;

    Mat<T, Gen<>, SparseLine<> > C_tot_enh;
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

    if ( debug_method_enhancement ) {
        cout << "dimension de la matrice C_tot_enh : ( " << C_tot_enh.nb_rows() << ", " << C_tot_enh.nb_cols() << " )" << endl;
        cout << "matrice C_tot_enh :" << endl;
        cout << C_tot_enh << endl << endl;
        cout << "dimension du vecteur q_tot_enh : " << q_tot_enh.size() << endl;
        cout << "vecteur q_tot_enh :" << endl;
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
    ///----------------------------------------------------------------------------------------

    cout << "Suppression du noyau de la matrice C_tot_enh" << endl;
    cout << "Calcul du vecteur nb_eq_indep_enh" << endl << endl;

    Vec<unsigned> eq_indep_enh;

    /// Orthonormalisation de Schmidt
    ///------------------------------
    //     Mat<T> C_tot_orth_enh = C_tot_enh;
    //
    //     TicToc t_orth_schmidt;
    //     t_orth_schmidt.start();
    //
    //     orthonormalisation_schmidt_row( C_tot_orth_enh, eq_indep_enh, 1.e-10 );
    //
    //     t_orth_schmidt.stop();
    //     cout << "Temps de calcul de l'orthonormalisation de schmidt de la matrice C_tot_enh associee au probleme global de minimisation pour la technique " << method << " : " << t_orth_schmidt.res << endl << endl;
    //
    //     C_tot_orth_enh.clear();

    /// Factorisation LDL avec detection des pivots nuls
    ///-------------------------------------------------
    eq_indep_enh = range( nb_eq_f_surf_enh + nb_eq_proj_f_surf_enh + nb_eq_f_vol_enh );

    TicToc t_detect_pivots;
    t_detect_pivots.start();

    Mat<T, Sym<>, SparseLine<> > C_tot_2_enhancement = C_tot_enh * trans( C_tot_enh );

    LDL_solver ls_enh;
    Vec< Vec<T> > Ker_enh;
    Vec<int> Pivots_enh;
    ls_enh.get_factorization( C_tot_2_enhancement, Ker_enh, Pivots_enh );

    t_detect_pivots.stop();
    cout << "Temps de calcul de la detection des pivots nuls de la matrice C_tot_enh associee au probleme global de minimisation pour la technique " << method << " : " << t_detect_pivots.res << endl << endl;

    if ( Pivots_enh.size() ) { // si le nombre de pivots nuls est non nul
        for (unsigned n=0;n<Pivots_enh.size();++n) {
            remove_if( eq_indep_enh, _1 == _2, Pivots_enh[ n ] );
        }
    }

    if ( debug_method_enhancement ) {
        cout << "indices des pivots nuls dans la matrice C_tot_enh " << Pivots_enh << endl;
    }

    C_tot_2_enhancement.clear();
    Ker_enh.free();
    Pivots_enh.free();

    unsigned nb_eq_indep_enh = eq_indep_enh.size();

    if ( debug_method_enhancement ) {
        cout << "nombre de lignes independantes dans la matrice C_tot_enh : " << nb_eq_indep_enh << endl;
        cout << "indices des lignes independantes dans la matrice C_tot_enh : " << eq_indep_enh << endl << endl;
    }

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

    if ( debug_method_enhancement ) {
        cout << "dimension de la matrice C_tot_tilde_enh : ( " << C_tot_tilde_enh.nb_rows() << ", " << C_tot_tilde_enh.nb_cols() << " )" << endl;
        cout << "matrice C_tot_tilde_enh :" << endl;
        cout << C_tot_tilde_enh << endl << endl;
        cout << "dimension du vecteur q_tot_tilde_enh : " << q_tot_tilde_enh.size() << endl;
        cout << "vecteur q_tot_tilde_enh :" << endl;
        cout << q_tot_tilde_enh << endl << endl;
    }

    C_tot_enh.clear();
    q_tot_enh.free();
    eq_indep_enh.free();

    /// Construction de la matrice globale K_enh et du vecteur global F_enh
    ///--------------------------------------------------------------------

    cout << "Construction de la matrice K_enh et du vecteur F_enh" << endl << endl;

    Mat<T, Gen<>, SparseLine<> > K_enh;
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
    Mat<T, Gen<>, SparseLine<> > trans_C_tot_tilde_enh = trans( C_tot_tilde_enh );
    K_enh( vec_unk, vec_unk ) = A_enh( vec_unk, vec_unk ) * 1.;
    K_enh( vec_unk_to_unk_plus_eq_indep, vec_unk ) = C_tot_tilde_enh( vec_eq_indep, vec_unk ) * 1.;
    K_enh( vec_unk, vec_unk_to_unk_plus_eq_indep ) = trans_C_tot_tilde_enh( vec_unk, vec_eq_indep ) * 1.;
    F_enh[ vec_unk ] = d_enh[ vec_unk ] * 1.;
    F_enh[ vec_unk_to_unk_plus_eq_indep ] = q_tot_tilde_enh[ vec_eq_indep ] * 1.;

    t_construct_K_F_enh.stop();
    cout << "Temps de calcul de remplissage de la matrice K_enh et du vecteur F_enh associes au probleme global de minimisation pour la technique " << method << " : " << t_construct_K_F_enh.res << endl << endl;

    if ( debug_method_enhancement ) {
        cout << "dimension de la matrice K_enh : ( " << K_enh.nb_rows() << ", " << K_enh.nb_cols() << " ) "<< endl;
        cout << "matrice K_enh :" << endl;
        cout << K_enh << endl << endl;
        cout << "dimension du vecteur F_enh : " << F_enh.size() << endl;
        cout << "vecteur F_enh :" << endl;
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
    ///-------------------------------------------------------------------

    cout << "Resolution du probleme global de minimsation K_enh * U_enh = F_enh a " << nb_unk_enh << " inconnues et " << nb_eq_indep_enh << " equations independantes imposees" << endl;
    cout << "Construction du vecteur U_enh" << endl << endl;

    TicToc t_solve_minimization_enhancement;
    t_solve_minimization_enhancement.start();

    if ( solver_minimisation == "LDL" ) {
#ifdef WITH_LDL
        Mat<T, Sym<>, SparseLine<> > K_LDL = K_enh;
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
        Mat<T, Gen<>, SparseUMFPACK > K_UMFPACK = K_enh;
        K_UMFPACK.get_factorization();
        U_enh = K_UMFPACK.solve( F_enh );
#endif
    }
    else if ( solver_minimisation == "Inv" ) {
        Mat<T, Sym<>, SparseLine<> > K_Inv = K_enh;
        U_enh = inv( K_Inv ) * F_enh;
    }
    else if ( solver_minimisation == "LUFactorize" ) {
        Mat<T, Sym<>, SparseLine<> > K_sym = K_enh;
        Mat<T, Gen<>, SparseLU > K_LU = K_sym;
//        Vec<int> vector_permutation;
        lu_factorize( K_LU/*, vector_permutation*/ );
        solve_using_lu_factorize( K_LU, /*vector_permutation,*/ F_enh, U_enh );
    }
    else {
        cerr << "Bing. Error : type de solveur pour la minimisation globale non implemente" << endl << endl;
    }

    t_solve_minimization_enhancement.stop();
    cout << "Temps de calcul de la resolution du probleme global de minimisation pour la technique " << method << " amelioree : " << t_solve_minimization_enhancement.res << endl << endl;

    if ( debug_method_enhancement ) {
        cout << "dimension du vecteur U_enh : " << U_enh.size() << endl;
        cout << "vecteur U_enh :" << endl;
        cout << U_enh << endl << endl;
    }
    if ( verif_solver_minimisation_enhancement ) {
        cout << "Verification de la resolution du probleme global de minimsation pour la technique " << method << " amelioree" << endl << endl;
        if ( norm_2( K_enh * U_enh - F_enh ) / norm_2( F_enh ) > tol_solver_minimisation_enhancement ) {
            cout << "residu K_enh * U_enh - F_enh =" << endl;
            cout << K_enh * U_enh - F_enh << endl << endl;
            cout << "norme 2 du residu =" << endl;
            cout << norm_2( K_enh * U_enh - F_enh ) / norm_2( F_enh ) << endl << endl;
        }
    }

    K_enh.clear();
    F_enh.free();

    /// Construction du vecteur vec_force_fluxes_enhancement_global
    ///------------------------------------------------------------

    cout << "Construction du vecteur vec_force_fluxes_enhancement_global" << endl << endl;

    Vec<T> vec_force_fluxes_enhancement_global;

    vec_force_fluxes_enhancement_global.resize( nb_unk_enh );
    vec_force_fluxes_enhancement_global.set( 0. );

    vec_force_fluxes_enhancement_global[ vec_unk ] = U_enh[ vec_unk ] * 1.;

    if ( debug_method_enhancement ) {
        cout << "dimension du vecteur vec_force_fluxes_enhancement_global : " << nb_unk_enh << endl;
        cout << "vecteur vec_force_fluxes_enhancement_global :" << endl;
        cout << vec_force_fluxes_enhancement_global << endl << endl;
    }

    vec_unk.free();
    U_enh.free();

    /// Construction des vecteurs vec_force_fluxes_enhancement[ k ][ d ] pour chaque face k du maillage et chaque direction d
    /// Construction des matrices mat_force_fluxes_enhancement[ k ] pour chaque face k du maillage
    ///----------------------------------------------------------------------------------------------------------------------

    cout << "Construction des vecteurs vec_force_fluxes_enhancement et des matrices force_fluxes_enhancement" << endl << endl;

    Vec< Vec< Vec<T> > > vec_force_fluxes_enhancement;

    vec_force_fluxes_enhancement.resize( m.sub_mesh(Number<1>()).elem_list.size() );

    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        vec_force_fluxes_enhancement[ k ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            vec_force_fluxes_enhancement[ k ][ d ].resize( m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
            vec_force_fluxes_enhancement[ k ][ d ].set( 0. );
        }
    }

    Vec< Mat<T> > mat_force_fluxes_enhancement;
    mat_force_fluxes_enhancement.resize( m.sub_mesh(Number<1>()).elem_list.size() );

    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        mat_force_fluxes_enhancement[ k ].resize( m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual(), dim );
        mat_force_fluxes_enhancement[ k ].set( 0. );
    }

    for (unsigned k=0;k<list_faces_enh.size();++k) {
        for (unsigned d=0;d<dim;++d) {
            for (unsigned j=0;j<m.sub_mesh(Number<1>()).elem_list[list_faces_enh[ k ]]->nb_nodes_virtual();++j) {
                vec_force_fluxes_enhancement[ list_faces_enh[ k ] ][ d ][ j ] += vec_force_fluxes_enhancement_global[ k * m.sub_mesh(Number<1>()).elem_list[list_faces_enh[ k ]]->nb_nodes_virtual() * dim + j * dim + d ];
            }
        }
    }

    for (unsigned k=0;k<list_faces_enh.size();++k) {
        for (unsigned d=0;d<dim;++d) {;
            for (unsigned j=0;j<m.sub_mesh(Number<1>()).elem_list[list_faces_enh[ k ]]->nb_nodes_virtual();++j) {
                mat_force_fluxes_enhancement[ list_faces_enh[ k ] ]( j, d ) += vec_force_fluxes_enhancement[ list_faces_enh[ k ] ][ d ][ j ];
            }
        }
    }

    if ( debug_force_fluxes_enhancement or debug_method_enhancement ) {
        for (unsigned k=0;k<list_faces_enh.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur des densites d'effort ameliorees associe a la face " << list_faces_enh[ k ] << " dans la direction " << d << " : " << m.sub_mesh(Number<1>()).elem_list[list_faces_enh[ k ]]->nb_nodes_virtual() << endl;
                cout << "vecteur des densites d'effort ameliorees associe a la face " << list_faces_enh[ k ] << " dans la direction " << d << " :" << endl;
                cout << vec_force_fluxes_enhancement[ list_faces_enh[ k ] ][ d ] << endl << endl;
            }
        }
    }
    if ( debug_force_fluxes_enhancement or debug_method_enhancement ) {
        for (unsigned k=0;k<list_faces_enh.size();++k) {
            cout << "dimension de la matrice des densites d'effort ameliorees associe a la face " << list_faces_enh[ k ] << " : ( " << m.sub_mesh(Number<1>()).elem_list[list_faces_enh[ k ]]->nb_nodes_virtual() << ", " << dim << " )" << endl;
            cout << "matrice des densites d'effort ameliorees associe a la face " << list_faces_enh[ k ] << " :" << endl;
            cout << mat_force_fluxes_enhancement[ list_faces_enh[ k ] ] << endl << endl;
        }
    }

    /// Construction des vecteurs vec_force_fluxes_enhancement[ k ][ d ] pour chaque face k du maillage et chaque direction d
    /// Construction des matrices mat_force_fluxes_enhancement[ k ] pour chaque face k du maillage
    ///----------------------------------------------------------------------------------------------------------------------

    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        for (unsigned d=0;d<dim;++d) {
            for (unsigned j=0;j<m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual();++j) {
                vec_force_fluxes[ k ][ d ][ j ] += vec_force_fluxes_enhancement[ k ][ d ][ j ];
            }
        }
    }

    if ( debug_force_fluxes_enhancement or debug_method_enhancement ) {
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur des densites d'effort standard + ameliorees associe a la face " << k << " dans la direction " << d << " : " << m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << endl;
                cout << "vecteur des densites d'effort standard + ameliorees associe a la face " << k << " dans la direction " << d << " :" << endl;
                cout << vec_force_fluxes[ k ][ d ] << endl << endl;
            }
        }
    }

    t_construct_force_fluxes_enh.stop();
    cout << "Temps de calcul de la construction amelioree des densites d'effort pour la technique " << method << " : " << t_construct_force_fluxes_enh.res << endl << endl;
}

#endif // Construct_enhanced_force_fluxes_EET_EESPT_h
