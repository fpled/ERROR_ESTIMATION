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

//#include "../CONNECTIVITY/Connectivity.h"

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
    /// pre-deformation appliquee sur tous les elements
    /// -----------------------------------------------
    if ( structure.find("square") != string::npos ) {
        T E_12 = 1;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            m.elem_list[n]->set_field( "pre_epsilon_init", Vec<T,unsigned(dim*(dim+1)/2) >( 0., -E_12/sqrt(2.), 0. ) );
        }
    }
    /// Hashin's coated shpere 3D
    /// pre-deformation appliquee sur tous les elements
    /// -----------------------------------------------
    else if ( structure.find("hashin") != string::npos ) {
        T E_v = 1;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            m.elem_list[n]->set_field( "pre_epsilon_init", Vec<T,unsigned(dim*(dim+1)/2) >( -E_v/3., 0., -E_v/3., 0., 0., -E_v/3. ) );
        }
    }
}

#endif // Boundary_conditions_homog_h
