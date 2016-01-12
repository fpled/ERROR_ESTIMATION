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
#include "Construct_connectivity_patch.h"
#include "../GEOMETRY/Calcul_connectivity.h"
#include "../DISCRETIZATION_ERROR/Discretization_error.h"

using namespace LMT;
using namespace std;

/// Calcul d'un champ de contrainte admissible, calcul d'un estimateur theta de l'erreur globale pour la methode basee sur la partition de l'unite (SPET)
/// -----------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T>
void calcul_error_estimate_partition_unity( TM &m, const TF &f, const string &pb, const string &solver, const string &method, T &theta, T &theta_init, Vec<T> &theta_elem, Vec<T> &theta_elem_init, Vec< Vec<T> > &E, const bool verif_solver = false, const T tol_solver = 1e-6, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool want_local_enrichment = false, const bool debug_mesh = false, const bool debug_error_estimate = false, const bool debug_local_effectivity_index = false, const bool debug_method = false ) {
    
    static const unsigned dim = TM::dim;
    
    Vec<unsigned> child_cpt;
    Vec< Vec<unsigned> > child_list;
    construct_child( m, child_cpt, child_list, debug_mesh );
    
    Vec<bool> correspondance_node_to_vertex_node;
    Vec<unsigned> connect_node_to_vertex_node;
    unsigned nb_vertex_nodes = match_node_to_vertex_node( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, debug_mesh );
    
    Vec<unsigned> elem_cpt_vertex_node;
    Vec< Vec<unsigned> > elem_list_vertex_node;
    construct_elems_connected_to_vertex_node( m, nb_vertex_nodes, correspondance_node_to_vertex_node, connect_node_to_vertex_node, elem_cpt_vertex_node, elem_list_vertex_node, debug_mesh );
    
    correspondance_node_to_vertex_node.free();
    
    Vec< Vec<unsigned> > face_type;
    construct_face_type( m, f, face_type, debug_mesh );
    
    /// Construction de la table de connectivite de chaque patch
    /// --------------------------------------------------------
    
    Vec< Vec<unsigned> > face_list_patch;
    Vec<unsigned> nb_points_face;
    Vec<unsigned> nb_points_elem;
    Vec<unsigned> nb_points_patch;
    Vec< Vec< Vec<unsigned> > > patch_face;
    Vec< Vec< Vec<unsigned> > > patch_elem;
    
    construct_connectivity_patch( m, nb_vertex_nodes, face_list_patch, child_cpt, child_list, elem_cpt_vertex_node, elem_list_vertex_node, nb_points_face, nb_points_elem, nb_points_patch, patch_face, patch_elem, debug_method );
    
    child_cpt.free();
    child_list.free();
    elem_cpt_vertex_node.free();
    
    /// Construction des contraintes cinematiques pour chaque noeud sommet j du maillage et chaque direction d
    /// ------------------------------------------------------------------------------------------------------
    
    if ( debug_method )
        cout << "Construction des contraintes cinematiques" << endl << endl;
    
    Vec< Vec< Vec<unsigned> > > constrained_points_list_patch;
    constrained_points_list_patch.resize( nb_vertex_nodes );
    
    Vec<unsigned> nb_constraints_patch;
    nb_constraints_patch.resize( nb_vertex_nodes );
    nb_constraints_patch.set( 0 );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        constrained_points_list_patch[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            for (unsigned k=0;k<face_list_patch[ j ].size();++k) {
                if ( face_type[ face_list_patch[ j ][ k ] ][ d ] == 1 ) { // si face_list_patch[ j ][ k ] (k eme face connectee au patch associe au noeud sommet j) est une face de type Dirichlet (en deplacement) dans la direction d
                    for (unsigned p=0;p<nb_points_face[ face_list_patch[ j ][ k ] ];++p) {
                        constrained_points_list_patch[ j ][ d ].push_back( patch_face[ j ][ k ][ p ] );
                    }
                }
            }
            remove_doubles( constrained_points_list_patch[ j ][ d ] );
            nb_constraints_patch[ j ] += constrained_points_list_patch[ j ][ d ].size();
        }
    }
    
    if ( debug_method ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "liste des points bloques dans le patch associe au noeud sommet " << j << " dans la direction " << d << " = " << constrained_points_list_patch[ j ][ d ] << endl;
            }
            cout << "nb de contraintes sur le patch associe au noeud sommet " << j << " = " << nb_constraints_patch[ j ] << endl << endl;
        }
    }
    
    face_list_patch.free();
    nb_points_face.free();
    patch_face.free();
    nb_constraints_patch.free();
    
    /// Construction des matrices K[ j ] pour chaque noeud sommet j du maillage
    /// -----------------------------------------------------------------------
    
    if ( debug_method )
        cout << "Construction des matrices K" << endl << endl;
    
    Vec< Mat<T, Sym<> > > K; 
    K.resize( nb_vertex_nodes );
    
    TicToc t_construct_K;
    t_construct_K.start();
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        K[ j ].resize( nb_points_patch[ j ] * dim );
        K[ j ].set( 0. );
    }
    
    Calcul_Vertex_Nodal_Matrix_K calcul_vertex_nodal_matrix_K;
    calcul_vertex_nodal_matrix_K.elem_list_vertex_node = &elem_list_vertex_node;
    calcul_vertex_nodal_matrix_K.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_vertex_nodal_matrix_K.patch_elem = &patch_elem;
    
    apply( m.elem_list, calcul_vertex_nodal_matrix_K, m, f, K );
    
    /// Prise en compte des conditions aux limites dans les matrices K[ j ] pour chaque noeud sommet j du maillage
    /// ----------------------------------------------------------------------------------------------------------
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        Vec<unsigned> vec_unk_patch = range( nb_points_patch[ j ] * dim );
        Vec<unsigned> vec_null;
        vec_null.resize( nb_points_patch[ j ] * dim );
        vec_null.set( 0 );
        for (unsigned d=0;d<dim;++d) {
            for (unsigned p=0;p<constrained_points_list_patch[ j ][ d ].size();++p) {
                K[ j ].row( constrained_points_list_patch[ j ][ d ][ p ] * dim + d ) = vec_null;
                K[ j ].col( constrained_points_list_patch[ j ][ d ][ p ] * dim + d ) = vec_null;
                K[ j ]( constrained_points_list_patch[ j ][ d ][ p ] * dim + d, constrained_points_list_patch[ j ][ d ][ p ] * dim + d ) = 1;
            }
        }
    }
    
    t_construct_K.stop();
    if ( debug_method )
        cout << "Temps de calcul du remplissage des matrices associees aux pbs locaux auto-equilibres par patch = " << t_construct_K.res << endl << endl;
    
    if ( debug_method ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "dimension de la matrice K associee au noeud sommet " << j << " = ( " << nb_points_patch[ j ] * dim << ", " << nb_points_patch[ j ] * dim << " )" << endl;
            cout << "matrice K associee au noeud sommet " << j << " =" << endl;
            cout << K[ j ] << endl << endl;
        }
    }
    
    /// Construction des vecteurs F[ j ] pour chaque noeud sommet j du maillage
    /// -----------------------------------------------------------------------
    
    if ( debug_method )
        cout << "Construction des vecteurs F" << endl << endl;
    
    Vec< Vec<T> > F;
    F.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        F[ j ].resize( nb_points_patch[ j ] * dim );
        F[ j ].set( 0. );
    }
    
    Calcul_Vertex_Nodal_Vector_F calcul_vertex_nodal_vector_F;
    calcul_vertex_nodal_vector_F.elem_list_vertex_node = &elem_list_vertex_node;
    calcul_vertex_nodal_vector_F.face_type = &face_type;
    calcul_vertex_nodal_vector_F.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_vertex_nodal_vector_F.patch_elem = &patch_elem;
    calcul_vertex_nodal_vector_F.pb = &pb;
    calcul_vertex_nodal_vector_F.want_local_enrichment = &want_local_enrichment;
    
    apply( m.elem_list, calcul_vertex_nodal_vector_F, m, f, F );
    
    face_type.free();
    
    /// Prise en compte des conditions aux limites dans les vecteurs F[ j ] pour chaque noeud sommet j du maillage
    /// ----------------------------------------------------------------------------------------------------------
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            for (unsigned p=0;p<constrained_points_list_patch[ j ][ d ].size();++p) {
                F[ j ][ constrained_points_list_patch[ j ][ d ][ p ] * dim + d ] = 0;
            }
        }
    }
    
    if ( debug_method ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "dimension du vecteur F associe au noeud sommet " << j << " = " << nb_points_patch[ j ] * dim << endl;
            cout << "vecteur F associe au noeud sommet " << j << " =" << endl;
            cout << F[ j ] << endl << endl;
        }
    }
    
    constrained_points_list_patch.free();
    
    cout << "Resolution des pbs locaux auto-equilibres par patch" << endl;
    cout << "---------------------------------------------------" << endl << endl;
    
    /// Resolution des problemes locaux K[ j ] * U[ j ] = F[ j ] pour chaque noeud sommet j du maillage
    /// Construction des vecteurs U[ j ] pour chaque noeud sommet j du maillage
    /// -------------------------------------------------------------------------------------------------
    
    if ( debug_method ) {
        cout << "Resolution des pbs locaux K * U = F" << endl;
        cout << "Construction des vecteurs U" << endl << endl;
    }
    
    Vec< Vec<T> > U;
    U.resize( nb_vertex_nodes );
    
    TicToc t_solve_local_SPET;
    t_solve_local_SPET.start();
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        U[ j ].resize( nb_points_patch[ j ] * dim );
        U[ j ].set( 0. );
    }
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        if ( solver == "CholMod" ) {
            #ifdef WITH_CHOLMOD
            Mat<T, Sym<>, SparseLine<> > K_sym = K[ j ];
            Mat<T, Sym<>, SparseCholMod > K_Cholmod = K_sym;
            U[ j ] = K_Cholmod.solve( F[ j ] );
            K_sym.clear();
            K_Cholmod.clear();
            #endif
        }
        else if ( solver == "LDL" ) {
            #ifdef WITH_LDL
            Mat<T, Sym<>, SparseLine<> > K_LDL = K[ j ];
            U[ j ] = F[ j ];
            LDL_solver ls;
            Vec< Vec<T> > Ker;
            Vec<int> Pivots;
            ls.get_factorization( K_LDL, Ker, Pivots );
            ls.solve( U[ j ] );
            Ker.free();
            Pivots.free();
            K_LDL.clear();
            #endif
        }
        else if ( solver == "UMFPACK" ) {
            #ifdef WITH_UMFPACK
            Mat<T, Gen<>, SparseUMFPACK > K_UMFPACK = K[ j ];
            K_UMFPACK.get_factorization();
            U[ j ] = K_UMFPACK.solve( F[ j ] );
            K_UMFPACK.clear();
            #endif
        }
        else if ( solver == "CholFactorize" ) {
            Mat<T, Sym<>, SparseLine<> > K_Chol = K[ j ];
            chol_factorize( K_Chol );
            solve_using_chol_factorize( K_Chol, F[ j ], U[ j ] );
            K_Chol.clear();
        }
        else if ( solver == "LUFactorize" ) {
            Mat<T, Sym<>, SparseLine<> > K_sym = K[ j ];
            Mat<T, Gen<>, SparseLU > K_LU = K_sym;
//             Vec<int> vector_permutation;
            lu_factorize( K_LU/*, vector_permutation*/ );
            solve_using_lu_factorize( K_LU, /*vector_permutation,*/ F[ j ], U[ j ] );
            K_sym.clear();
            K_LU.clear();
        }
        else if ( solver == "Inv" ) {
            Mat<T, Sym<>, SparseLine<> > K_sym = K[ j ];
            U[ j ] = inv( K_sym ) * F[ j ];
            K_sym.clear();
        }
        else {
            cerr << "Bing. Error : solveur " << solver << " pour la resolution des pbs locaux non implemente" << endl << endl;
        }
        if ( debug_method == 0 and verif_solver == 0 ) {
            K[ j ].clear();
            F[ j ].free();
        }
    }
    
    t_solve_local_SPET.stop();
    cout << "Temps de calcul de la resolution des pbs locaux auto-equilibres par patch = " << t_solve_local_SPET.res << endl << endl;
    
    if ( debug_method ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "dimension du vecteur U associee au noeud sommet " << j << " = " << nb_points_patch[ j ] * dim << endl;
            cout << "vecteur U associe au noeud sommet " << j << " =" << endl;
            cout << U[ j ] << endl << endl;
        }
    }
     if ( verif_solver ) {
         if ( debug_method )
             cout << "Verification de la resolution des problemes locaux auto-equilibres par patch : tolerance = " << tol_solver << endl << endl;
         for (unsigned j=0;j<nb_vertex_nodes;++j) {
             T residual = norm_2( K[ j ] * U[ j ] - F[ j ] );
             T b = norm_2( F[ j ] );
            if ( residual / b > tol_solver ) {
                cout << "residu associe au noeud sommet " << j << " :" << endl;
//                cout << "K * U - F =" << endl;
//                cout << K[ j ] * U[ j ] - F[ j ] << endl;
                cout << "norme du residu = " << residual << endl;
                cout << "norme du residu relatif =" << residual / b << endl << endl;
            }
        }
    }
    
    nb_points_patch.free();
    K.free();
    F.free();
    
    /// Construction des vecteurs E[ n ] pour chaque element n du maillage
    /// ------------------------------------------------------------------
    
    if ( debug_method )
        cout << "Construction des vecteurs E" << endl << endl;
    
    E.resize( m.elem_list.size() );
    
    for (unsigned n=0;n<m.elem_list.size();++n) {
        E[ n ].resize( nb_points_elem[ n ] * dim );
        E[ n ].set( 0. );
    }
    
    Calcul_Elem_Vector_E calcul_elem_vector_E;
    calcul_elem_vector_E.elem_list_vertex_node = &elem_list_vertex_node;
    calcul_elem_vector_E.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_elem_vector_E.patch_elem = &patch_elem;
    
    apply( m.elem_list, calcul_elem_vector_E, U, E );
    
    if ( debug_method ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "dimension du vecteur E associe a l'element " << n << " = " << nb_points_elem[ n ] * dim << endl;
            cout << "vecteur E associe a l'element " << n << " =" << endl;
            cout << E[ n ] << endl << endl;
        }
    }
    
    nb_points_elem.free();
    patch_elem.free();
    U.free();
    connect_node_to_vertex_node.free();
    elem_list_vertex_node.free();
    
    /// ------------------------------------------------------------------------------------------------------- ///
    /// Construction d'un champ de contrainte admissible par element et Calcul d'un estimateur d'erreur globale ///
    /// ------------------------------------------------------------------------------------------------------- ///
    
    cout << "Construction d'un champ de contrainte admissible par element et Calcul d'un estimateur d'erreur globale" << endl;
    cout << "-------------------------------------------------------------------------------------------------------" << endl << endl;
    
    theta = 0.;
    
    theta_elem.resize( m.elem_list.size() );
    theta_elem.set( 0. );

    Calcul_Elem_Error_Estimate_Init_SPET<T> calcul_elem_error_estimate_init_SPET;
    calcul_elem_error_estimate_init_SPET.E = &E;
    calcul_elem_error_estimate_init_SPET.theta_elem = &theta_elem;
    calcul_elem_error_estimate_init_SPET.theta_elem_init = &theta_elem_init;
    calcul_elem_error_estimate_init_SPET.theta_init = &theta_init;
    
    apply( m.elem_list, calcul_elem_error_estimate_init_SPET, m, f, theta );

    if ( debug_error_estimate or debug_method ) {
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
    cout << "theta / ||u_h|| = " << theta / m.norm_dep * 100. << " %" << endl;
    cout << "theta_init / ||u_h_init|| = " << theta_init / m.norm_dep_init * 100. << " %" << endl << endl;

    if ( want_global_discretization_error ) {
        m.eff_index_SPET = theta / m.discretization_error;
        cout << "indice d'efficacite global :" << endl;
        cout << "eta = theta / e" << endl;
        cout << "    = " << m.eff_index_SPET << endl << endl;
    }
    
    if ( want_local_discretization_error ) {
        Vec<T> eff_index_elem;
        eff_index_elem.resize( m.elem_list.size() );
        eff_index_elem.set( 0. );
        
        apply( m.elem_list, Calcul_Elem_Effectivity_Index(), method, eff_index_elem );
        
        if ( debug_local_effectivity_index or debug_method ) {
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
