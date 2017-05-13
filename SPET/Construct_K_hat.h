//
// C++ Interface: Construct_K_hat
//
// Description: construction de la matrice K_hat sur chaque noeud sommet du maillage pour la methode basee sur la partition de l'unite (SPET)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_K_hat_SPET_h
#define Construct_K_hat_SPET_h

#include "SPET.h"

using namespace LMT;
using namespace std;

/// Construction des matrices K_hat[ j ] pour chaque noeud sommet j du maillage
/// ---------------------------------------------------------------------------
template<class TM, class TF, class TMatV>
void construct_K_hat( const TM &m, const TF &f, const unsigned &nb_vertex_nodes, const Vec<unsigned> &connect_node_to_vertex_node, const Vec< Vec<unsigned> > &elem_list_vertex_node, const Vec< Vec< Vec<unsigned> > > &patch_elem, const Vec<unsigned> &nb_points_patch, const Vec< Vec< Vec<unsigned> > > &constrained_points_list_patch, TMatV &K_hat, const bool disp = false ) {
    
    static const unsigned dim = TM::dim;
    
    if ( disp )
        cout << "Construction des matrices K_hat" << endl << endl;
    
    TicToc t;
    t.start();
    
    K_hat.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        K_hat[ j ].resize( nb_points_patch[ j ] * dim );
        K_hat[ j ].set( 0. );
    }
    
    Calcul_Vertex_Nodal_Matrix_K_hat calcul_vertex_nodal_matrix_K_hat;
    calcul_vertex_nodal_matrix_K_hat.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_vertex_nodal_matrix_K_hat.elem_list_vertex_node = &elem_list_vertex_node;
    calcul_vertex_nodal_matrix_K_hat.patch_elem = &patch_elem;
    
    apply( m.elem_list, calcul_vertex_nodal_matrix_K_hat, m, f, K_hat );
    
    /// Prise en compte des conditions aux limites dans les matrices K_hat[ j ] pour chaque noeud sommet j du maillage
    /// --------------------------------------------------------------------------------------------------------------
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        Vec<unsigned> vec_unk_patch = range( nb_points_patch[ j ] * dim );
        Vec<unsigned> vec_null;
        vec_null.resize( nb_points_patch[ j ] * dim );
        vec_null.set( 0 );
        for (unsigned d=0;d<dim;++d) {
            for (unsigned p=0;p<constrained_points_list_patch[ j ][ d ].size();++p) {
                K_hat[ j ].row( constrained_points_list_patch[ j ][ d ][ p ] * dim + d ) = vec_null;
                K_hat[ j ].col( constrained_points_list_patch[ j ][ d ][ p ] * dim + d ) = vec_null;
                K_hat[ j ]( constrained_points_list_patch[ j ][ d ][ p ] * dim + d, constrained_points_list_patch[ j ][ d ][ p ] * dim + d ) = 1;
            }
        }
    }
    
    t.stop();
    if ( disp )
        cout << "temps de calcul du remplissage des matrices associees aux pbs locaux auto-equilibres par patch = " << t.res << endl << endl;
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "dimension de la matrice K_hat associee au noeud sommet " << j << " = ( " << nb_points_patch[ j ] * dim << ", " << nb_points_patch[ j ] * dim << " )" << endl;
            cout << "matrice K_hat associee au noeud sommet " << j << " =" << endl;
            cout << K_hat[ j ] << endl;
        }
    }
    
}

#endif // Construct_K_hat_SPET_h
