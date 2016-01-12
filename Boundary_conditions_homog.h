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
#include "GEOMETRY/Calcul_connectivity.h"

using namespace LMT;
using namespace std;

/// Creation des conditions aux limites
/// -----------------------------------
template<class TF, class TM>
void set_boundary_conditions_init( TF &f, TM &m, const string &boundary_condition_D, const string &pb, const string &structure ) {

    static const unsigned dim = TM::dim;
    typedef typename TM::TNode::T T;

    if ( m.node_list.size() ) {
        T penalty_val;
        if ( boundary_condition_D == "lagrange" )
            penalty_val = 0;
        else if ( boundary_condition_D == "penalty" )
            penalty_val = 1e8;

        m.update_skin();
        if ( pb == "direct" ) {
            /// Carre 2D
            /// pre-deformation et pre-contrainte appliquees sur tous les elements
            /// -----------------------------------------------------------------
            if ( structure.find("square") != string::npos ) {
                Vec<T,unsigned(dim*(dim+1)/2) >  pre_eps_init;
                pre_eps_init.set( 0, 0. );
                pre_eps_init.set( 1, -1./2 );
                pre_eps_init.set( 2, 0. );
                for (unsigned n=0;n<m.elem_list.size();++n)
                    m.elem_list[n]->set_field( "pre_epsilon_init", pre_eps_init );
            }
        }
    }
}

#endif // Boundary_conditions_homog_h
