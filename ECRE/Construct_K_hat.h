//
// C++ Interface: Construct_K_hat
//
// Description: construction de la matrice K_hat sur chaque element du maillage pour les methodes basees sur la condition de prolongement (EESPT et EET)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_K_hat_h
#define Construct_K_hat_h

#include "ECRE.h"

using namespace LMT;
using namespace std;

/// Construction des matrices K_hat[ n ] pour chaque element n du maillage
/// ----------------------------------------------------------------------
template<class TM, class TF, class T>
void construct_K_hat( const TM &m, const TF &f, const string &method, Vec< Mat<T, Sym<> > > &K_hat, const bool debug_method = false ) {

    cout << "Construction des matrices K_hat" << endl << endl;

    TicToc t_construct_K_hat;
    t_construct_K_hat.start();

    K_hat.resize( m.elem_list.size() );

    apply( m.elem_list, Calcul_Elem_Matrix_K_hat(), m, f, K_hat );

    t_construct_K_hat.stop();
    cout << "Temps de calcul de remplissage des matrices K_hat associees aux problemes locaux pour la technique " << method << " : " << t_construct_K_hat.res << endl << endl;

    if ( debug_method ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "matrice K_hat de l'element " << n << " :" << endl;
            cout << K_hat[ n ] << endl << endl;
        }
    }
}

#endif // Construct_K_hat_h
