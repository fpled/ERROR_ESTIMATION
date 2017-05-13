//
// C++ Interface: Construct_dep_hat
//
// Description: resolution des problemes locaux a partir de K_hat et F_hat avec choix du solveur
//              construction du vecteur dep_hat_patch sur chaque noeud sommet du maillage
//              construction du vecteur dep_hat sur chaque element du maillage
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_dep_hat_SPET_h
#define Construct_dep_hat_SPET_h

#include "SPET.h"

using namespace LMT;
using namespace std;

/// Resolution des problemes locaux K_hat[ j ] * dep_hat_patch[ j ] = F_hat[ j ] pour chaque noeud sommet j du maillage
/// Construction des vecteurs dep_hat_patch[ j ] pour chaque noeud sommet j du maillage
/// Construction des vecteurs dep_hat[ n ] pour chaque element n du maillage
/// -------------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class TVV, class TTVV, class TMatV>
void construct_dep_hat( const TM &m, const TF &f, const string &solver, const unsigned &nb_vertex_nodes, const Vec<unsigned> &connect_node_to_vertex_node, const Vec< Vec<unsigned> > &elem_list_vertex_node, const Vec< Vec< Vec<unsigned> > > &patch_elem, const Vec<unsigned> &nb_points_patch, const Vec<unsigned> &nb_points_elem, TMatV &K_hat, TVV &F_hat, TTVV &dep_hat, const bool verif_solver = false, const T tol_solver = 1e-6, const bool disp = false ) {
    
    static const unsigned dim = TM::dim;
    
    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    typedef Mat<T, Sym<>, SparseCholMod > TMatSparseCholMod;
    typedef Mat<T, Gen<>, SparseUMFPACK > TMatGenSparseUMFPACK;
    typedef Mat<T, Gen<>, SparseLU > TMatGenSparseLU;
    
    if ( disp ) {
        cout << "Resolution des pbs locaux auto-equilibres par patch" << endl;
        cout << "---------------------------------------------------" << endl << endl;
    }
    
    if ( disp ) {
        cout << "Resolution des pbs locaux K_hat * U_hat_patch = F_hat" << endl;
        cout << "Construction des vecteurs U_hat_patch" << endl << endl;
    }
    
    TTVV dep_hat_patch;
    dep_hat_patch.resize( nb_vertex_nodes );
    
    TicToc t;
    t.start();
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        dep_hat_patch[ j ].resize( nb_points_patch[ j ] * dim );
        dep_hat_patch[ j ].set( 0. );
    }
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        if ( solver == "CholMod" ) {
            #ifdef WITH_CHOLMOD
            TMatSymSparse K_hat_sym = K_hat[ j ];
            TMatSparseCholMod K_hat_CholMod = K_hat_sym;
            K_hat_CholMod.get_factorization();
            dep_hat_patch[ j ] = K_hat_CholMod.solve( F_hat[ j ] );
            K_hat_sym.clear();
            K_hat_CholMod.clear();
            #endif
        }
        else if ( solver == "LDL" ) {
            #ifdef WITH_LDL
            TMatSymSparse K_hat_LDL = K_hat[ j ];
            dep_hat_patch[ j ] = F_hat[ j ];
            LDL_solver ls;
            Vec< Vec<T> > Ker;
            Vec<int> Pivots;
            ls.get_factorization( K_hat_LDL, Ker, Pivots );
            ls.solve( dep_hat_patch[ j ] );
            Ker.free();
            Pivots.free();
            K_hat_LDL.clear();
            #endif
        }
        else if ( solver == "UMFPACK" ) {
            #ifdef WITH_UMFPACK
            TMatGenSparseUMFPACK K_hat_UMFPACK = K_hat[ j ];
            K_hat_UMFPACK.get_factorization();
            dep_hat_patch[ j ] = K_hat_UMFPACK.solve( F_hat[ j ] );
            K_hat_UMFPACK.clear();
            #endif
        }
        else if ( solver == "CholFactorize" ) {
            TMatSymSparse K_hat_Chol = K_hat[ j ];
            chol_factorize( K_hat_Chol );
            solve_using_chol_factorize( K_hat_Chol, F_hat[ j ], dep_hat_patch[ j ] );
            K_hat_Chol.clear();
        }
        else if ( solver == "LUFactorize" ) {
            TMatSymSparse K_hat_sym = K_hat[ j ];
            TMatGenSparseLU K_hat_LU = K_hat_sym;
//            Vec<int> vector_permutation;
            lu_factorize( K_hat_LU/*, vector_permutation*/ );
            solve_using_lu_factorize( K_hat_LU, /*vector_permutation,*/ F_hat[ j ], dep_hat_patch[ j ] );
            K_hat_sym.clear();
            K_hat_LU.clear();
        }
        else if ( solver == "Inv" ) {
            TMatSymSparse K_hat_Inv = K_hat[ j ];
            dep_hat_patch[ j ] = inv( K_hat_Inv ) * F_hat[ j ];
            K_hat_Inv.clear();
        }
        else {
            cerr << "Bing. Error : solveur " << solver << " pour la resolution des pbs locaux non implemente" << endl << endl;
        }
//        if ( not( disp ) and not( verif_solver ) ) {
//            K_hat[ j ].clear();
//            F_hat[ j ].free();
//        }
    }
    
    t.stop();
    cout << "temps de calcul de la resolution des pbs locaux auto-equilibres par patch = " << t.res << endl << endl;
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "dimension du vecteur U_hat_patch associee au noeud sommet " << j << " = " << nb_points_patch[ j ] * dim << endl;
            cout << "vecteur U_hat_patch associe au noeud sommet " << j << " =" << endl;
            cout << dep_hat_patch[ j ] << endl << endl;
        }
    }
    if ( verif_solver ) {
        if ( disp )
            cout << "Verification de la resolution des problemes locaux auto-equilibres par patch : tolerance = " << tol_solver << endl << endl;
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            T residual = norm_2( K_hat[ j ] * dep_hat_patch[ j ] - F_hat[ j ] );
            T b = norm_2( F_hat[ j ] );
            if ( residual / b > tol_solver ) {
                cout << "residu associe au noeud sommet " << j << " :" << endl;
//                cout << "K_hat * U_hat_patch - F_hat =" << endl;
//                cout << K_hat[ j ] * U_hat_patch[ j ] - F_hat[ j ] << endl;
                cout << "norme du residu = " << residual << endl;
                cout << "norme du residu relatif =" << residual / b << endl << endl;
            }
        }
    }
    
    /// Construction des vecteurs dep_hat[ n ] pour chaque element n du maillage
    /// ------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs U_hat" << endl << endl;
    
    dep_hat.resize( m.elem_list.size() );
    
    for (unsigned n=0;n<m.elem_list.size();++n) {
        dep_hat[ n ].resize( nb_points_elem[ n ] * dim );
        dep_hat[ n ].set( 0. );
    }
    
    Calcul_Elem_Vector_Dep_hat calcul_elem_vector_dep_hat;
    calcul_elem_vector_dep_hat.elem_list_vertex_node = &elem_list_vertex_node;
    calcul_elem_vector_dep_hat.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_elem_vector_dep_hat.patch_elem = &patch_elem;
    
    apply( m.elem_list, calcul_elem_vector_dep_hat, dep_hat_patch, dep_hat );
    
    dep_hat_patch.free();
    
    if ( disp ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "dimension du vecteur U_hat associe a l'element " << n << " = " << nb_points_elem[ n ] * dim << endl;
            cout << "vecteur U_hat associe a l'element " << n << " =" << endl;
            cout << dep_hat[ n ] << endl << endl;
        }
    }
}

#endif // Construct_dep_hat_SPET_h
