//
// C++ Interface: Boundary_conditions_homog
//
// Description: creation des conditions aux limites
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2015
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Boundary_conditions_homog_h
#define Boundary_conditions_homog_h

#include "INTEREST_QUANTITY/Interest_quantity.h"
#include "CONNECTIVITY/Calcul_connectivity.h"

using namespace LMT;
using namespace std;

/// Creation des conditions aux limites
/// -----------------------------------
template<class TM>
void set_load_conditions_init( TM &m, const string &structure ) {

    static const unsigned dim = TM::dim;
    typedef typename TM::TNode::T T;

    m.update_skin();

    /// Carre 2D
    /// pre-deformation et pre-contrainte appliquees sur tous les elements
    /// ------------------------------------------------------------------
    if ( structure.find("square") != string::npos ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            m.elem_list[n]->set_field( "pre_epsilon_init", Vec<T,unsigned(dim*(dim+1)/2) >( 0., -1/sqrt(2.), 0. ) );
        }
    }
    /// Cube 3D
    /// pre-deformation et pre-contrainte appliquees sur tous les elements
    /// ------------------------------------------------------------------
    else if ( structure.find("cube") != string::npos ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            m.elem_list[n]->set_field( "pre_epsilon_init", Vec<T,unsigned(dim*(dim+1)/2) >( 0., -1/sqrt(2.), 0., 0., 0., 0. ) );
        }
    }
}

#endif // Boundary_conditions_homog_h
