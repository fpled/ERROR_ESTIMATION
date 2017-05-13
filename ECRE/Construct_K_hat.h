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
template<class TM, class TF, class TMatV>
void construct_K_hat( const TM &m, const TF &f, TMatV &K_hat, const bool disp = false ) {
    
    if ( disp )
        cout << "Construction des matrices K_hat" << endl << endl;
    
    TicToc t;
    t.start();
    
    K_hat.resize( m.elem_list.size() );
    
    apply( m.elem_list, Calcul_Elem_Matrix_K_hat(), m, f, K_hat );
    
    t.stop();
    if ( disp )
        cout << "temps de calcul du remplissage des matrices associees aux pbs locaux par element = " << t.res << endl << endl;
    
    if ( disp ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "matrice K_hat de l'element " << n << " :" << endl;
            cout << K_hat[ n ] << endl;
        }
    }
}

#endif // Construct_K_hat_h
