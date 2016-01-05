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
template<class TM, class TF, class T>
void construct_dep_hat( const TM &m, const TF &f, const string &method, const string &solver, Vec< Mat<T, Sym<> > > &K_hat, Vec< Vec<T> > &F_hat, Vec< Vec<T> > &dep_hat, const bool verif_solver = false, const T tol_solver = 1e-6, const bool debug_method = false ) {

    cout << "Resolution des problemes locaux pour la technique " << method << endl;
    cout << "Construction des vecteurs dep_hat" << endl << endl;

    TicToc t_solve_local_EET_EESPT;
    t_solve_local_EET_EESPT.start();

    dep_hat.resize( m.elem_list.size() );

    for (unsigned n=0;n<m.elem_list.size();++n) {
        if ( solver == "CholMod" ) {
            #ifdef WITH_CHOLMOD
            Mat<T, Sym<>, SparseLine<> > K_hat_sym = K_hat[ n ];
            Mat<T, Sym<>, SparseCholMod > K_hat_CholMod = K_hat_sym;
            K_hat_CholMod.get_factorization();
            dep_hat[ n ] = K_hat_CholMod.solve( F_hat[ n ] );
            #endif
        }
        else if ( solver == "LDL" ) {
            #ifdef WITH_LDL
            Mat<T, Sym<>, SparseLine<> > K_hat_LDL = K_hat[ n ];
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
            Mat<T, Gen<>, SparseUMFPACK > K_hat_UMFPACK = K_hat[ n ];
            K_hat_UMFPACK.get_factorization();
            dep_hat[ n ] = K_hat_UMFPACK.solve( F_hat[ n ] );
            #endif
        }
        else if ( solver == "CholFactorize" ) {
            Mat<T, Sym<>, SparseLine<> > K_hat_Chol = K_hat[ n ];
            chol_factorize( K_hat_Chol );
            solve_using_chol_factorize( K_hat_Chol, F_hat[ n ], dep_hat[ n ] );
        }
        else if ( solver == "LUFactorize" ) {
            Mat<T, Sym<>, SparseLine<> > K_hat_sym = K_hat[ n ];
            Mat<T, Gen<>, SparseLU > K_hat_LU = K_hat_sym;
//             Vec<int> vector_permutation;
            lu_factorize( K_hat_LU/*, vector_permutation*/ );
            solve_using_lu_factorize( K_hat_LU, /*vector_permutation,*/ F_hat[ n ], dep_hat[ n ] );
        }
        else if ( solver == "Inv" ) {
            Mat<T, Sym<>, SparseLine<> > K_hat_Inv = K_hat[ n ];
            dep_hat[ n ] = inv( K_hat_Inv ) * F_hat[ n ];
        }
        else {
            cerr << "Bing. Error : solveur " << solver << " pour la resolution des problemes locaux non implemente" << endl << endl;
        }
    }

    t_solve_local_EET_EESPT.stop();
    cout << "Temps de calcul de la resolution des problemes locaux pour la technique " << method << " : " << t_solve_local_EET_EESPT.res << endl << endl;

    if ( debug_method ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "vecteur dep_hat de l'element " << n << " :" << endl;
            cout << dep_hat[ n ] << endl << endl;
        }
    }
    if ( verif_solver ) {
        cout << "Verification de la resolution des problemes locaux pour la technique " << method << " : tolerance = " << tol_solver << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            T residual = norm_2( K_hat[ n ] * dep_hat[ n ] - F_hat[ n ] );
            T b = norm_2( F_hat[ n ] );
            if ( residual / b > tol_solver ) {
                cout << "residu associe a l'element " << n << " :" << endl;
//                cout << "K_hat * dep_hat - F_hat = " << endl;
//                cout << K_hat[ n ] * dep_hat[ n ] - F_hat[ n ] << endl << endl;
                cout << "norme du residu = " << residual << endl << endl;
                cout << "norme du residu relatif = " << residual / b << endl << endl;
            }
        }
    }
}

#endif // Construct_dep_hat_h
