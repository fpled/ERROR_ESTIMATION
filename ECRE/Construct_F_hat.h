//
// C++ Interface: Construct_F_hat
//
// Description: construction du vecteur F_hat sur chaque element du maillage pour les methodes basees sur la condition de prolongement (EESPT et EET)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_F_hat_h
#define Construct_F_hat_h

#include "ECRE.h"
#include "../GEOMETRY/Calcul_geometry.h"

using namespace LMT;
using namespace std;

/// Construction des vecteurs F_hat[ n ] pour chaque element n du maillage
/// ----------------------------------------------------------------------
template<class TM, class TF, class T>
void construct_F_hat( TM &m, const TF &f, const string &pb, const bool &balancing, const Vec<bool> &elem_flag_bal, const Vec<bool> &elem_flag_enh, const Vec< Vec< Vec<T> > > &vec_force_fluxes, Vec< Vec<T> > &F_hat, const bool want_local_enrichment = false, const bool debug_method = false, const bool debug_geometry = false ) {


    if ( debug_method ) {
        cout << "Construction des vecteurs F_hat";
        if ( balancing )
            cout << " avec procedure d'equilibrage";
        cout << endl << endl;
    }

    Vec<unsigned> node_cpt_face;
    Vec< Vec<unsigned> > node_list_face;
    construct_nodes_connected_to_face( m, node_cpt_face, node_list_face, debug_geometry );

    Vec<unsigned> elem_cpt_node;
    Vec< Vec<unsigned> > elem_list_node;
    construct_elems_connected_to_node( m, elem_cpt_node, elem_list_node, debug_geometry );

    elem_list_node.free();

    F_hat.resize( m.elem_list.size() );

    Calcul_Elem_Vector_F_hat<T> calcul_elem_vector_F_hat;
    calcul_elem_vector_F_hat.node_list_face = &node_list_face;
    calcul_elem_vector_F_hat.elem_cpt_node = &elem_cpt_node;
    calcul_elem_vector_F_hat.balancing = &balancing;
    calcul_elem_vector_F_hat.elem_flag_bal = &elem_flag_bal;
    calcul_elem_vector_F_hat.elem_flag_enh = &elem_flag_enh;
    calcul_elem_vector_F_hat.pb = &pb;
    calcul_elem_vector_F_hat.want_local_enrichment = &want_local_enrichment;
    calcul_elem_vector_F_hat.vec_force_fluxes = &vec_force_fluxes;

    apply( m.elem_list, calcul_elem_vector_F_hat, m, f, F_hat );

    if ( debug_method ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "vecteur F_hat de l'element " << n << " =" << endl;
            cout << F_hat[ n ] << endl << endl;
        }
    }
}

#endif // Construct_F_hat_h
