//
// C++ Interface: Construct_F_hat
//
// Description: construction du vecteur F_hat sur chaque noeud sommet du maillage pour la methode basee sur la partition de l'unite (SPET)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_F_hat_SPET_h
#define Construct_F_hat_SPET_h

#include "SPET.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"

using namespace LMT;
using namespace std;

/// Construction des vecteurs F_hat[ j ] pour chaque noeud sommet j du maillage
/// ---------------------------------------------------------------------------
template<class TM, class TF, class TVV>
void construct_F_hat( TM &m, const TF &f, const string &pb, const unsigned &nb_vertex_nodes, const Vec< Vec<unsigned> > &face_type, const Vec<unsigned> &connect_node_to_vertex_node, const Vec< Vec<unsigned> > &elem_list_vertex_node, const Vec< Vec< Vec<unsigned> > > &patch_elem, const Vec<unsigned> &nb_points_patch, const Vec< Vec< Vec<unsigned> > > &constrained_points_list_patch, TVV &F_hat, const bool want_local_enrichment = false, const bool disp = false ) {
    
    static const unsigned dim = TM::dim;
    
    if ( disp )
        cout << "Construction des vecteurs F_hat" << endl << endl;
    
    F_hat.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        F_hat[ j ].resize( nb_points_patch[ j ] * dim );
        F_hat[ j ].set( 0. );
    }
    
    Calcul_Vertex_Nodal_Vector_F_hat calcul_vertex_nodal_vector_F_hat;
    calcul_vertex_nodal_vector_F_hat.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_vertex_nodal_vector_F_hat.elem_list_vertex_node = &elem_list_vertex_node;
    calcul_vertex_nodal_vector_F_hat.face_type = &face_type;
    calcul_vertex_nodal_vector_F_hat.patch_elem = &patch_elem;
    calcul_vertex_nodal_vector_F_hat.pb = &pb;
    calcul_vertex_nodal_vector_F_hat.want_local_enrichment = &want_local_enrichment;
    
    apply( m.elem_list, calcul_vertex_nodal_vector_F_hat, m, f, F_hat );
    
    /// Prise en compte des conditions aux limites dans les vecteurs F_hat[ j ] pour chaque noeud sommet j du maillage
    /// --------------------------------------------------------------------------------------------------------------
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            for (unsigned p=0;p<constrained_points_list_patch[ j ][ d ].size();++p) {
                F_hat[ j ][ constrained_points_list_patch[ j ][ d ][ p ] * dim + d ] = 0;
            }
        }
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "dimension du vecteur F_hat associe au noeud sommet " << j << " = " << nb_points_patch[ j ] * dim << endl;
            cout << "vecteur F_hat associe au noeud sommet " << j << " =" << endl;
            cout << F_hat[ j ] << endl << endl;
        }
    }
}

#endif // Construct_F_hat_SPET_h
