//
// C++ Interface: Construct_connectivity_patch
//
// Description: construction de la table de connectivite de chaque pacth pour la methode basee sur la partition de l'unite (SPET)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_connectivity_patch_h
#define Construct_connectivity_patch_h

#include "SPET.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"

using namespace LMT;
using namespace std;

/// Construction de la table de connectivite de chaque patch : nb_unk_patch[ j ] pour chaque noeud sommet j du maillage
///                                                            nb_unk_elem[ n ] pour chaque element n du maillage
///                                                            nb_unk_face[ n ] pour chaque face k du maillage 
///                                                            patch_elem[ j ][ n ][ i ] pour chaque noeud i de chaque element n du patch j
///                                                            patch_face[ j ][ n ][ i ] pour chaque noeud i de chaque face k du patch j
/// ---------------------------------------------------------------------------------------------------------------------------------------
template<class TM>
void construct_connectivity_patch( const TM &m, const unsigned &nb_vertex_nodes, Vec< Vec<unsigned> > &face_list_patch, const Vec<unsigned> &child_cpt, const Vec< Vec<unsigned> > &child_list, const Vec<unsigned> &elem_cpt_vertex_node, const Vec< Vec<unsigned> > &elem_list_vertex_node, Vec<unsigned> &nb_points_face, Vec<unsigned> &nb_points_elem, Vec<unsigned> &nb_points_patch, Vec< Vec< Vec<unsigned> > > &patch_face, Vec< Vec< Vec<unsigned> > > &patch_elem, const bool debug_method = false ) {
    
    static const unsigned dim = TM::dim;
    typedef typename TM::TNode::T T;
    
    if ( debug_method )
        cout << "Construction de la table de connectivite de chaque patch" << endl << endl;
    
    if ( debug_method )
        cout << "Construction de la liste des faces associees a chaque patch" << endl << endl;
    
    face_list_patch.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned n=0;n<elem_cpt_vertex_node[ j ];++n) {
            for (unsigned k=0;k<child_cpt[ elem_list_vertex_node[ j ][ n ] ];++k) {
                face_list_patch[ j ].push_back( child_list[ elem_list_vertex_node[ j ][ n ] ][ k ] );
            }
        }
        remove_doubles( face_list_patch[ j ] );
    }
    
    if ( debug_method )
        cout << "Construction du nombre d'inconnues associees a chaque element et a chaque face" << endl << endl;
    
    nb_points_elem.resize( m.elem_list.size(), 0 );
    
    nb_points_face.resize( m.sub_mesh(Number<1>()).elem_list.size(), 0 );
    
    apply( m.elem_list, Calc_Nb_Points_Elem(), nb_points_elem ); // nb_points_elem[ n ] contient le nb de points associees aux vecteurs e_i pour chaque element n du maillage
    
    apply( m.sub_mesh(Number<1>()).elem_list, Calc_Nb_Points_Face(), nb_points_face ); // nb_points_face[ k ] contient le nb de points associees aux vecteurs e_i pour chaque face k du maillage
    
    if ( debug_method ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "nb d'inconnues e_i et u_star associees a l'element " << n << " : " << nb_points_elem[ n ] << endl;
        }
        cout << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            cout << "nb d'inconnues e_i et u_star associees a la face " << k << " : " << nb_points_face[ k ] << endl;
        }
        cout << endl << endl;
    }
    
    if ( debug_method )
        cout << "Reperage de la position des points de chaque element et chaque face" << endl << endl;
    
    Vec< Vec< Vec<T> > > pos_elem;
    pos_elem.resize( m.elem_list.size() );
    
    for (unsigned n=0;n<m.elem_list.size();++n) {
        pos_elem[ n ].resize( nb_points_elem[ n ] );
        for (unsigned p=0;p<nb_points_elem[ n ];++p) {
            pos_elem[ n ][ p ].resize( dim, 0.0 );
        }
    }
    
    Vec< Vec< Vec<T> > > pos_face;
    pos_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        pos_face[ k ].resize( nb_points_face[ k ] );
        for (unsigned p=0;p<nb_points_face[ k ];++p) {
            pos_face[ k ][ p ].resize( dim, 0.0 );
        }
    }
    
    apply( m.elem_list, Calc_Pos_Elem(), pos_elem ); // pos_elem[ n ][ p ][ d ] contient la position du point p de l'element n du maillage dans la direction d
    
    apply( m.sub_mesh(Number<1>()).elem_list, Calc_Pos_Face(), pos_face ); // pos_face[ k ][ p ][ d ] contient la position du point p de la face k du maillage dans la direction d
    
    if ( debug_method ) {
        for (unsigned n=0;n<m.elem_list.size();++n) {
            for (unsigned p=0;p<nb_points_elem[ n ];++p) {
                cout << "position du point " << p << " de l'element " << n << " : " << pos_elem[ n ][ p ] << endl;
            }
            cout << endl << endl;
        }
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned p=0;p<nb_points_face[ k ];++p) {
                cout << "position du point " << p << " de la face " << k << " : " << pos_face[ k ][ p ] << endl;
            }
            cout << endl << endl;
        }
        cout << endl << endl;
    }
    
    if ( debug_method )
        cout << "Reperage de la position des points de chaque patch" << endl << endl;
    
    nb_points_patch.resize( nb_vertex_nodes, 0 );
    
    Vec< Vec< Vec<T> > > pos_patch;
    pos_patch.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned n=0;n<elem_cpt_vertex_node[ j ];++n) {
            for (unsigned p=0;p<nb_points_elem[ elem_list_vertex_node[ j ][ n ] ];++p) {
                pos_patch[ j ].push_back( pos_elem[ elem_list_vertex_node[ j ][ n ] ][ p ] );
            }
        }
        remove_doubles_approx( pos_patch[ j ], 1e-6 );
        nb_points_patch[ j ] = pos_patch[ j ].size();
    }
    
    if ( debug_method ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned p=0;p<nb_points_patch[ j ];++p) {
                cout << "position du point " << p << " du patch associee au noeud sommet " << j << " : " << pos_patch[ j ][ p ] << endl;
            }
            cout << endl << endl;
        }
    }
    
    if ( debug_method )
        cout << "Construction de la table de connectivite de chaque patch : patch_elem et patch_face" << endl << endl;
    
    patch_elem.resize( nb_vertex_nodes );
    
    patch_face.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        patch_elem[ j ].resize( elem_cpt_vertex_node[ j ] );
        patch_face[ j ].resize( face_list_patch[ j ].size() );
        for (unsigned n=0;n<elem_cpt_vertex_node[ j ];++n) {
            patch_elem[ j ][ n ].resize( nb_points_elem[ elem_list_vertex_node[ j ][ n ] ] );
            for (unsigned pe=0;pe<nb_points_elem[ elem_list_vertex_node[ j ][ n ] ];++pe) {
                for (unsigned pp=0;pp<nb_points_patch[ j ];++pp) {
                    if ( norm_2( pos_elem[ elem_list_vertex_node[ j ][ n ] ][ pe ] - pos_patch[ j ][ pp ] ) < 1e-6 ) {
                        patch_elem[ j ][ n ][ pe ] = pp;
                    }
                }
            }
        }
        for (unsigned k=0;k<face_list_patch[ j ].size();++k) {
            patch_face[ j ][ k ].resize( nb_points_face[ face_list_patch[ j ][ k ] ] );
            for (unsigned pf=0;pf<nb_points_face[ face_list_patch[ j ][ k ] ];++pf) {
                for (unsigned pp=0;pp<nb_points_patch[ j ];++pp) {
                    if ( norm_2( pos_face[ face_list_patch[ j ][ k ] ][ pf ] - pos_patch[ j ][ pp ] ) < 1e-6 ) {
                        patch_face[ j ][ k ][ pf ] = pp;
                    }
                }
            }
        }
    }
    
    if ( debug_method ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "Table de connectivite du patch " << j << endl << endl;
            for (unsigned n=0;n<elem_cpt_vertex_node[ j ];++n) {
                cout << "indice du " << n << " eme element ( element " << elem_list_vertex_node[ j ][ n ] << " du maillage ) du patch associee au noeud sommet " << j << " : " << patch_elem[ j ][ n ] << endl;
            }
            cout << endl;
            for (unsigned k=0;k<face_list_patch[ j ].size();++k) {
                cout << "indice de la " << k << " eme face ( face " << face_list_patch[ j ][ k ] << " du maillage ) du patch associee au noeud sommet " << j << " : " << patch_face[ j ][ k ] << endl;
            }
            cout << endl << endl;
        }
    }
}

#endif // Construct_connectivity_patch_h
