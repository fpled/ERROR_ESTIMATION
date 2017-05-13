//
// C++ Interface: Construct_kinematic_contraints_patch
//
// Description: construction des contraintes cinematiques sur chaque patch pour la methode basee sur la partition de l'unite (SPET)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_kinematic_contraints_patch_h
#define Construct_kinematic_contraints_patch_h

#include "SPET.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"

using namespace LMT;
using namespace std;

/// Construction des contraintes cinematiques sur chaque patch : nb_constraints_patch[ j ] pour chaque noeud sommet j du maillage
///                                                              constrained_points_list_patch[ j ][ d ][ i ] pour chaque noeud i du patch j dans la direction d
/// ------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM>
void construct_kinematic_constraints_patch( const TM &m, const unsigned &nb_vertex_nodes, const Vec< Vec<unsigned> > &face_type, const Vec< Vec<unsigned> > &face_list_patch, const Vec<unsigned> &nb_points_face, const Vec< Vec< Vec<unsigned> > > &patch_face, Vec< Vec< Vec<unsigned> > > &constrained_points_list_patch, Vec<unsigned> &nb_constraints_patch, const bool disp = false ) {
    
    static const unsigned dim = TM::dim;
    
    if ( disp )
        cout << "Construction des contraintes cinematiques" << endl << endl;
    
    constrained_points_list_patch.resize( nb_vertex_nodes );
    
    nb_constraints_patch.resize( nb_vertex_nodes );
    nb_constraints_patch.set( 0 );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        constrained_points_list_patch[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            for (unsigned k=0;k<face_list_patch[ j ].size();++k) {
                if ( face_type[ face_list_patch[ j ][ k ] ][ d ] == 1 ) { // si face_list_patch[ j ][ k ] (k eme face connectee au patch associe au noeud sommet j) est une face de Dirichlet (en deplacement) dans la direction d
                    for (unsigned p=0;p<nb_points_face[ face_list_patch[ j ][ k ] ];++p) {
                        constrained_points_list_patch[ j ][ d ].push_back( patch_face[ j ][ k ][ p ] );
                    }
                }
            }
            remove_doubles( constrained_points_list_patch[ j ][ d ] );
            nb_constraints_patch[ j ] += constrained_points_list_patch[ j ][ d ].size();
        }
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "liste des points bloques dans le patch associe au noeud sommet " << j << " dans la direction " << d << " = " << constrained_points_list_patch[ j ][ d ] << endl;
            }
            cout << "nb de contraintes sur le patch associe au noeud sommet " << j << " = " << nb_constraints_patch[ j ] << endl << endl;
        }
    }
}

#endif // Construct_kinematic_contraints_patch_h
