//
// C++ Interface: Verification
//
// Description: verification de l'equilibre
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Verification_h
#define Verification_h

#include "../GEOMETRY/Calcul_geometry.h"

using namespace LMT;
using namespace std;

/// Verification de l'equilibre du pb direct/adjoint
///-------------------------------------------------
template<class TF>
void check_equilibrium( const TF &f, const string &pb ) {
    typedef typename TF::ScalarType T;
    cout << "Verification de l'equilibre global du pb " << pb << " : erreur = " << norm_2( f.matrices(Number<0>()) * f.vectors[0] - f.sollicitation ) << ", erreur relative = " << norm_2( f.matrices(Number<0>()) * f.vectors[0] - f.sollicitation ) / norm_2( f.sollicitation ) << endl << endl;
}

/// Construction du vecteur de vecteurs residual_force_fluxes
///----------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class TTVVV, class TTWW, class S, class B, class TTVV>
void check_elem_eq_force_fluxes( const TE &elem, const TM &m, const TF &f, const TVV &node_list_face, const TV &elem_cpt_node, const TTVVV &vec_force_fluxes, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVV &residual_force_fluxes ) {}

template<class T>
struct Check_Elem_Eq_Force_Fluxes {
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec<unsigned>* elem_cpt_node;
    const Vec< Vec< Vec<T> > >* vec_force_fluxes;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec<T> > &residual_force_fluxes ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        check_elem_eq_force_fluxes( elem, m, f, *node_list_face, *elem_cpt_node, *vec_force_fluxes, f.vectors, ind, *pb, *want_local_enrichment, residual_force_fluxes );
    }
};

/// Verification de l'equilibre des densites d'effort pour les methodes EESPT et EET
///---------------------------------------------------------------------------------
template<class TM, class TF, class T>
void check_equilibrium_force_fluxes( TM &m, const TF &f, const string pb, const Vec< Vec< Vec<T> > > &vec_force_fluxes, const T tol_eq_force_fluxes = 1e-6, const bool want_local_enrichment = false, const bool debug_geometry = false ) {

    static const unsigned dim = TM::dim;

    Vec<unsigned> node_cpt_face;
    Vec< Vec<unsigned> > node_list_face;
    construct_nodes_connected_to_face( m, node_cpt_face, node_list_face, debug_geometry );

    Vec<unsigned> elem_cpt_node;
    Vec< Vec<unsigned> > elem_list_node;
    construct_elems_connected_to_node( m, elem_cpt_node, elem_list_node, debug_geometry );

    elem_list_node.free();

    Vec< Vec<T> > residual_force_fluxes;
    residual_force_fluxes.resize( m.elem_list.size() );

    for (unsigned n=0;n<m.elem_list.size();++n) {
        if ( dim == 1 ) {
            residual_force_fluxes[ n ].resize( dim );
        }
        else if ( dim == 2 ) {
            residual_force_fluxes[ n ].resize( dim + 1 );
        }
        else if ( dim == 3 ) {
            residual_force_fluxes[ n ].resize( 2 * dim );
        }
        residual_force_fluxes[ n ].set( 0. );
    }

    Check_Elem_Eq_Force_Fluxes<T> check_elem_eq_force_fluxes;
    check_elem_eq_force_fluxes.node_list_face = &node_list_face;
    check_elem_eq_force_fluxes.elem_cpt_node = &elem_cpt_node;
    check_elem_eq_force_fluxes.vec_force_fluxes = &vec_force_fluxes;
    check_elem_eq_force_fluxes.pb = &pb;
    check_elem_eq_force_fluxes.want_local_enrichment = &want_local_enrichment;

    apply( m.elem_list, check_elem_eq_force_fluxes, m, f, residual_force_fluxes );

    cout << "--------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Verification de l'equilibre des densites d'effort en resultante et en moment sur chaque element : tolerance = " << tol_eq_force_fluxes << endl;
    cout << "--------------------------------------------------------------------------------------------------------------------------" << endl << endl;

    for (unsigned n=0;n<m.elem_list.size();++n) {
        for (unsigned d=0;d<dim;++d) {
            if ( fabs( residual_force_fluxes[ n ][ d ] ) > tol_eq_force_fluxes ) {
                cout << "verification de l'equilibre des densites d'effort en resultante pour l'element " << n << " dans la direction " << d << " :" << endl;
                cout << residual_force_fluxes[ n ][ d ] << " != 0" << endl << endl;
            }
        }
        if ( dim == 1 or dim == 2 ) {
            if ( fabs( residual_force_fluxes[ n ][ dim ] ) > tol_eq_force_fluxes ) {
                cout << "verification de l'equilibre des densites d'effort en moment pour l'element " << n << " dans la direction " << dim << " :" << endl;
                cout << residual_force_fluxes[ n ][ dim ] << " != 0" << endl << endl;
            }
        }
        else if ( dim == 3 ) {
            for (unsigned d=0;d<dim;++d) {
                if ( fabs( residual_force_fluxes[ n ][ dim + d ] ) > tol_eq_force_fluxes ) {
                    cout << "verification de l'equilibre des densites d'effort en moment pour l'element " << n << " dans la direction " << d << " :" << endl;
                    cout << residual_force_fluxes[ n ][ dim + d ] << " != 0" << endl << endl;
                }
            }
        }
    }
}

#endif // Verification_h
