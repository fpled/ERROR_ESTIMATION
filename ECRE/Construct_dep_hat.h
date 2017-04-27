//
// C++ Interface: Construct_dep_hat
//
// Description: resolution des problemes locaux a partir de K_hat et F_hat avec choix du solveur
//              construction du vecteur dep_hat sur chaque element du maillage
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_dep_hat_h
#define Construct_dep_hat_h

#include "ECRE.h"

using namespace LMT;
using namespace std;

/// Construction des vecteurs dep_hat[ n ] pour chaque element n du maillage
/// ------------------------------------------------------------------------
template<class TM, class TF, class T, class TVV, class TTVV, class TMatV>
void construct_dep_hat( const TM &m, const TF &f, const string &solver, TMatV &K_hat, TVV &F_hat, TTVV &dep_hat, const bool verif_solver = false, const T tol_solver = 1e-6, const bool disp = false ) {
    
    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    typedef Mat<T, Sym<>, SparseCholMod > TMatSparseCholMod;
    typedef Mat<T, Gen<>, SparseUMFPACK > TMatGenSparseUMFPACK;
    typedef Mat<T, Gen<>, SparseLU > TMatGenSparseLU;
    
    if ( disp ) {
        cout << "Resolution des pbs locaux par element" << endl;
        cout << "-------------------------------------" << endl << endl;
    }
    
    if ( disp ) {
        cout << "Resolution des pbs locaux K_hat * U_hat = F_hat" << endl;
        cout << "Construction des vecteurs U_hat" << endl << endl;
    }
    
    TicToc t;
    t.start();
    
    dep_hat.resize( m.elem_list.size() );
    
    for (unsigned n=0;n<m.elem_list.size();++n) {
        if ( solver == "CholMod" ) {
            #ifdef WITH_CHOLMOD
            TMatSymSparse K_hat_sym = K_hat[ n ];
            TMatSparseCholMod K_hat_CholMod = K_hat_sym;
            K_hat_CholMod.get_factorization();
            dep_hat[ n ] = K_hat_CholMod.solve( F_hat[ n ] );
            #endif
        }
        else if ( solver == "LDL" ) {
            #ifdef WITH_LDL
            TMatSymSparse K_hat_LDL = K_hat[ n ];
            dep_hat[ n ] = F_hat[ n ];
            LDL_solver ls;
            Vec< Vec<T> > Ker;
            Vec<int> Pivots;
            ls.get_factorization( K_hat_LDL, Ker, Pivots );
            ls.solve( dep_hat[ n ] );
            #endif
        }
        else if ( solver == "UMFPACK" ) {
            #ifdef WITH_UMFPACK
            TMatGenSparseUMFPACK K_hat_UMFPACK = K_hat[ n ];
            K_hat_UMFPACK.get_factorization();
            dep_hat[ n ] = K_hat_UMFPACK.solve( F_hat[ n ] );
            #endif
        }
        else if ( solver == "CholFactorize" ) {
            TMatSymSparse K_hat_Chol = K_hat[ n ];
            chol_factorize( K_hat_Chol );
            solve_using_chol_factorize( K_hat_Chol, F_hat[ n ], dep_hat[ n ] );
        }
        else if ( solver == "LUFactorize" ) {
            TMatSymSparse K_hat_sym = K_hat[ n ];
            TMatGenSparseLU K_hat_LU = K_hat_sym;
//            Vec<int> vector_permutation;
            lu_factorize( K_hat_LU/*, vector_permutation*/ );
            solve_using_lu_factorize( K_hat_LU, /*vector_permutation,*/ F_hat[ n ], dep_hat[ n ] );
        }
        else if ( solver == "Inv" ) {
            TMatSymSparse K_hat_Inv = K_hat[ n ];
            dep_hat[ n ] = inv( K_hat_Inv ) * F_hat[ n ];
        }
        else {
            cerr << "Bing. Error : solveur " << solver << " pour la resolution des pbs locaux non implemente" << endl << endl;
        }
    }
    
    t.stop();
    cout << "temps de calcul de la resolution des pbs locaux par element = " << t.res << endl << endl;
    
    if ( disp ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "vecteur U_hat de l'element " << n << " =" << endl;
            cout << dep_hat[ n ] << endl << endl;
        }
    }
    if ( verif_solver ) {
        if ( disp )
            cout << "Verification de la resolution des pbs locaux par element : tolerance = " << tol_solver << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            T residual = norm_2( K_hat[ n ] * dep_hat[ n ] - F_hat[ n ] );
            T b = norm_2( F_hat[ n ] );
            if ( residual / b > tol_solver ) {
                cout << "residu associe a l'element " << n << " :" << endl;
//                cout << "K_hat * U_hat - F_hat = " << endl;
//                cout << K_hat[ n ] * U_hat[ n ] - F_hat[ n ] << endl;
                cout << "norme du residu = " << residual << endl;
                cout << "norme du residu relatif = " << residual / b << endl << endl;
            }
        }
    }
}

#endif // Construct_dep_hat_h
