//
// C++ Interface: Calcul_geometry
//
// Description: calcul des informations relatives a la geometrie, au type de face et au type de noeud
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_geometry_h
#define Calcul_geometry_h

#include "Geometry.h"

using namespace LMT;
using namespace std;

/// Compteur du nb de faces connectees a l'element n du maillage : cpt_children[ n ]
/// Liste des faces connectees a l'element n du maillage  : list_children[ n ]
///---------------------------------------------------------------------------------
template<class TM>
void construct_children( TM &m, Vec<unsigned> &cpt_children, Vec< Vec<unsigned> > &list_children, const bool &debug_geometry ) {
    cpt_children.resize( m.elem_list.size() );
    cpt_children.set( 0 );
    list_children.resize( m.elem_list.size() );

    m.update_skin();
    apply( m.elem_list, Counter_Children(), m, cpt_children, list_children ); // cpt_children[ n ] contient le nb de faces connectees a l'element n du maillage // list_children[ n ] contient le numero des faces connectees a l'element n du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des faces connectees aux elements du maillage" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "nb de faces connectees a l'element " << n << " : " << cpt_children[ n ] << endl;
            cout << "liste des faces connectees a l'element " << n << " : " << list_children[ n ] << endl << endl;
        }
    }
}

/// Compteur du nb d'elements connectes au noeud i du maillage : cpt_elems_node[ i ]
/// Liste des elements connectes au noeud i du maillage  : list_elems_node[ i ]
///---------------------------------------------------------------------------------
template<class TM>
void construct_elems_connected_to_node( const TM &m, Vec<unsigned> &cpt_elems_node, Vec< Vec<unsigned> > &list_elems_node, const bool &debug_geometry ) {
    cpt_elems_node.resize( m.node_list.size() );
    cpt_elems_node.set( 0 );
    list_elems_node.resize( m.node_list.size() );

    apply( m.elem_list, Counter(), cpt_elems_node, list_elems_node ); // cpt_elems_node[ i ] contient le nb d'elements connectes au noeud i du maillage // list_elems_node[ i ] contient le numero des elements connectes au noeud i du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des elements connectes aux noeuds du maillage" << endl << endl;
        for (unsigned i=0;i<m.node_list.size();++i) {
            cout << "nb d'elements connectes au noeud " << i << " : " << cpt_elems_node[ i ] << endl;
            cout << "liste des elements connectes au noeud " << i << " : " << list_elems_node[ i ] << endl << endl;
        }
    }
}

/// Compteur du nb de faces connectees au noeud i du maillage : cpt_faces_node[ i ]
/// Liste des faces connectees au noeud i du maillage  : list_faces_node[ i ]
///--------------------------------------------------------------------------------
template<class TM>
void construct_faces_connected_to_node( TM &m, Vec<unsigned> &cpt_faces_node, Vec< Vec<unsigned> > &list_faces_node, const bool &debug_geometry ) {
    cpt_faces_node.resize( m.node_list.size() );
    cpt_faces_node.set( 0 );
    list_faces_node.resize( m.node_list.size() );

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, Counter(), cpt_faces_node, list_faces_node ); // cpt_faces_node[ i ] contient le nb de faces connectees au noeud i du maillage // list_faces_node[ i ] contient le numero des faces connectees au noeud i du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des faces connectees aux noeuds du maillage" << endl << endl;
        for (unsigned i=0;i<m.node_list.size();++i) {
            cout << "nb de faces connectees au noeud " << i << " : " << cpt_faces_node[ i ] << endl;
            cout << "liste des faces connectees au noeud " << i << " : " << list_faces_node[ i ] << endl << endl;
        }
    }
}

/// Compteur du nb de noeuds connectes a l'element n du maillage : cpt_nodes_elem[ n ]
/// Liste des noeuds connectes a l'element n du maillage  : list_nodes_elem[ n ]
///-----------------------------------------------------------------------------------
template<class TM>
void construct_nodes_connected_to_elem( const TM &m, Vec<unsigned> &cpt_nodes_elem, Vec< Vec<unsigned> > &list_nodes_elem, const bool &debug_geometry ) {
    cpt_nodes_elem.resize( m.elem_list.size() );
    cpt_nodes_elem.set( 0 );
    list_nodes_elem.resize( m.elem_list.size() );

    apply( m.elem_list, Counter_Nodes(), cpt_nodes_elem, list_nodes_elem ); // cpt_nodes_elem[ n ] contient le nb de noeuds connectes a l'element n du maillage // list_nodes_elem[ n ] contient le numero des noeuds connectes a l'element n du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des noeuds connectes aux elements du maillage" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "nb de noeuds connectes a l'element " << n << " : " << cpt_nodes_elem[ n ] << endl;
            cout << "liste des noeuds connectes a l'element " << n << " : " << list_nodes_elem[ n ] << endl << endl;
        }
    }
}

/// Compteur du nb de noeuds connectes a la face k du maillage : cpt_nodes_face[ k ]
/// Liste des noeuds connectes a la face k du maillage  : list_nodes_face[ k ]
///---------------------------------------------------------------------------------
template<class TM>
void construct_nodes_connected_to_face( TM &m, Vec<unsigned> &cpt_nodes_face, Vec< Vec<unsigned> > &list_nodes_face, const bool &debug_geometry ) {
    cpt_nodes_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    cpt_nodes_face.set( 0 );
    list_nodes_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, Counter_Nodes(), cpt_nodes_face, list_nodes_face ); // cpt_nodes_face[ k ] contient le nb de noeuds connectes a la face k du maillage // list_nodes_face[ k ] contient le numero des noeuds connectes a la face k du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des noeuds connectes aux faces du maillage" << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            cout << "nb de noeuds connectes a la face " << k << " : " << cpt_nodes_face[ k ] << endl;
            cout << "liste des noeuds connectes a la face " << k << " : " << list_nodes_face[ k ] << endl << endl;
        }
    }
}

/// Construction du vecteur ( analogue de m.node_list pour les noeuds ) contenant les numeros des noeuds sommets dans le maillage : connect_node_to_vertex_node[ i ] pour chaque noeud i du maillage
/// Compteur du nb de noeuds sommets dans le maillage : nb_vertex_nodes
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM>
unsigned match_node_to_vertex_node( const TM &m, Vec<bool> &correspondance_node_to_vertex_node, Vec<unsigned> &connect_node_to_vertex_node, const bool &debug_geometry ) {
    correspondance_node_to_vertex_node.resize( m.node_list.size() );
    correspondance_node_to_vertex_node.set( 0 );

    apply( m.elem_list, Construct_Correspondance_Node_To_Vertex_Node(), correspondance_node_to_vertex_node ); // si correspondance_node_to_vertex_node[ i ] == 1, i est un noeud sommet
//         PRINT( correspondance_node_to_vertex_node );

    connect_node_to_vertex_node.resize( m.node_list.size() );
    connect_node_to_vertex_node.set( 0 );
    unsigned nb_vertex_nodes = 0;
    for (unsigned i=0;i<m.node_list.size();++i) {
        if ( correspondance_node_to_vertex_node[ i ] )
            connect_node_to_vertex_node[ i ] = nb_vertex_nodes++; // nb_vertex_nodes contient le nb de noeuds sommets dans le maillage // connect_node_to_vertex_node[ i ] = j : le i eme noeud du maillage corrrepond au j eme noeud sommet si i est un noeud sommet, sinon connect_node_to_vertex_node[ i ] = 0
    }

    if ( debug_geometry ) {
        cout << "Correspondance entre noeuds et noeuds sommets du maillage" << endl << endl;
        cout << "nombre de noeuds sommets dans le maillage : " << nb_vertex_nodes << endl;
        cout << "correspondance entre noeuds et noeuds sommets du maillage : " << connect_node_to_vertex_node << endl << endl;
    }
    return nb_vertex_nodes;
}

/// Compteur du nb d'elements connectes au noeud sommet j du maillage : cpt_elems_vertex_node[ j ]
/// Liste des elements connectes au noeud i du maillage  : list_elems_vertex_node[ j ]
///-----------------------------------------------------------------------------------------------
template<class TM>
void construct_elems_connected_to_vertex_node( const TM &m, const unsigned &nb_vertex_nodes, const Vec<bool> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &cpt_elems_vertex_node, Vec< Vec<unsigned> > &list_elems_vertex_node, const bool &debug_geometry ) {
    cpt_elems_vertex_node.resize( nb_vertex_nodes );
    cpt_elems_vertex_node.set( 0 );
    list_elems_vertex_node.resize( nb_vertex_nodes );

    Counter_Vertex counter_vertex;
    counter_vertex.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;

    apply( m.elem_list, counter_vertex, connect_node_to_vertex_node, cpt_elems_vertex_node, list_elems_vertex_node ); // cpt_elems_vertex_node[ j ] contient le nb d'elements connectes au noeud sommet j du maillage // list_elems_vertex_node[ j ] contient le numero des elements connectes au noeud sommet j du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des elements connectes aux noeuds sommets du maillage" << endl << endl;
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "nb d'elements connectes au noeud sommet " << j << " : " << cpt_elems_vertex_node[ j ] << endl;
            cout << "liste des elements connectes au noeud sommet " << j << " : " << list_elems_vertex_node[ j ] << endl << endl;
        }
    }
}

/// Compteur du nb de faces connectees au noeud sommet j du maillage : cpt_faces_vertex_node[ j ]
/// Liste des faces connectees au noeud i du maillage  : list_faces_vertex_node[ j ]
///----------------------------------------------------------------------------------------------
template<class TM>
void construct_faces_connected_to_vertex_node( TM &m, const unsigned &nb_vertex_nodes, const Vec<bool> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &cpt_faces_vertex_node, Vec< Vec<unsigned> > &list_faces_vertex_node, const bool &debug_geometry ) {
    cpt_faces_vertex_node.resize( nb_vertex_nodes );
    cpt_faces_vertex_node.set( 0 );
    list_faces_vertex_node.resize( nb_vertex_nodes );

    Counter_Vertex counter_vertex;
    counter_vertex.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, counter_vertex, connect_node_to_vertex_node, cpt_faces_vertex_node, list_faces_vertex_node ); // cpt_faces_vertex_node[ j ] contient le nb de faces connectees au noeud sommet j du maillage // list_faces_vertex_node[ j ] contient le numero des faces connectees au noeud sommet j du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des faces connectees aux noeuds sommets du maillage" << endl << endl;
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "nb de faces connectees au noeud sommet " << j << " : " << cpt_faces_vertex_node[ j ] << endl;
            cout << "liste des faces connectees au noeud sommet " << j << " : " << list_faces_vertex_node[ j ] << endl << endl;
        }
    }
}

/// Compteur du nb de noeuds sommets connectes a l'element n du maillage : cpt_vertex_nodes_elem[ n ]
/// Liste des noeuds sommets connectes a l'element n du maillage  : list_vertex_nodes_elem[ n ]
///--------------------------------------------------------------------------------------------------
template<class TM>
void construct_vertex_nodes_connected_to_elem( const TM &m, const Vec<bool> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &cpt_vertex_nodes_elem, Vec< Vec<unsigned> > &list_vertex_nodes_elem, const bool &debug_geometry ) {
    cpt_vertex_nodes_elem.resize( m.elem_list.size() );
    cpt_vertex_nodes_elem.set( 0 );
    list_vertex_nodes_elem.resize( m.elem_list.size() );

    Counter_Vertex_Nodes counter_vertex_nodes;
    counter_vertex_nodes.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;

    apply( m.elem_list, counter_vertex_nodes, connect_node_to_vertex_node, cpt_vertex_nodes_elem, list_vertex_nodes_elem ); // cpt_vertex_nodes_elem[ n ] contient le nb de noeuds sommets connectes a l'element n du maillage // list_vertex_nodes_elem[ n ] contient le numero des noeuds sommets connectes a l'element n du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des noeuds sommets connectes aux elements du maillage" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "nb de noeuds sommets connectes a l'element " << n << " : " << cpt_vertex_nodes_elem[ n ] << endl;
            cout << "liste des noeuds sommets connectes a l'element " << n << " : " << list_vertex_nodes_elem[ n ] << endl << endl;
        }
    }
}

/// Compteur du nb de noeuds sommets connectes a la face k du maillage : cpt_vertex_nodes_face[ k ]
/// Liste des noeuds sommets connectes a la face k du maillage  : list_vertex_nodes_face[ k ]
///------------------------------------------------------------------------------------------------
template<class TM>
void construct_vertex_nodes_connected_to_face( TM &m, const Vec<bool> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &cpt_vertex_nodes_face, Vec< Vec<unsigned> > &list_vertex_nodes_face, const bool &debug_geometry ) {
    cpt_vertex_nodes_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    cpt_vertex_nodes_face.set( 0 );
    list_vertex_nodes_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );

    Counter_Vertex_Nodes counter_vertex_nodes;
    counter_vertex_nodes.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, counter_vertex_nodes, connect_node_to_vertex_node, cpt_vertex_nodes_face, list_vertex_nodes_face ); // cpt_vertex_nodes_face[ k ] contient le nb de noeuds sommets connectes a la face k du maillage // list_vertex_nodes_face[ k ] contient le numero des noeuds sommets connectes a la face k du maillage

    if ( debug_geometry ) {
        cout << "Compteur et liste des noeuds sommets connectes aux faces du maillage" << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            cout << "nb de noeuds sommets connectes a la face " << k << " : " << cpt_vertex_nodes_face[ k ] << endl;
            cout << "liste des noeuds sommets connectes a la face " << k << " : " << list_vertex_nodes_face[ k ] << endl << endl;
        }
    }
}

/// Construction du vecteur de vecteurs type_face : 
/// type_face[ k ][ d ] = 0 : face k interieur a Omega dans la direction d
/// type_face[ k ][ d ] = 1 : face k appartenant a un bord delta_1_Omega dans la direction d
/// type_face_dim[ k ][ d ] = 2 : face k appartenant a un bord delta_2_Omega dans la direction d
///---------------------------------------------------------------------------------------------
template<class TM, class TF>
void construct_type_face( TM &m, const TF &f, Vec< Vec<unsigned> > &type_face, const bool &debug_geometry ) {

    static const unsigned dim = TM::dim;

    type_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );

    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        type_face[ k ].resize( dim );
        type_face[ k ].set( 0 ) ;
    }

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, Construct_Type_Face(), m, f, type_face );

    if ( debug_geometry ) {
        cout << "Construction du vecteur de vecteurs type_face" << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "type de la face " << k << " dans la direction " << d << " : " << type_face[ k ][ d ] << endl;
            }
            cout << endl << endl;
        }
    }
}

/// Construction du vecteur de vecteurs type_node :
/// type_node[ i ][ d ] = 0 : noeud i interieur a Omega dans la direction d
/// type_node[ i ][ d ] = 1 : noeud i sur un bord delta_1_Omega dans la direction d
/// type_node[ i ][ d ] = 2 : noeud i sur un bord delta_2_Omega dans la direction d
/// type_node[ i ][ d ] = 12 : noeud i a la separation entre delta_1_Omega et delat_2_Omega dans la direction d
///------------------------------------------------------------------------------------------------------------
template<class TM, class TF>
void construct_type_node( TM &m, const TF &f, const Vec< Vec<unsigned> > &type_face, Vec< Vec<unsigned> > &type_node, const bool &debug_geometry ) {

    static const unsigned dim = TM::dim;

    type_node.resize( m.node_list.size() );
    for (unsigned i=0;i<m.node_list.size();++i) {
        type_node[ i ].resize( dim );
        type_node[ i ].set( 0 );
    }

    Vec<unsigned> cpt_faces_node;
    Vec< Vec<unsigned> > list_faces_node;
    construct_faces_connected_to_node( m, cpt_faces_node, list_faces_node, debug_geometry );

    Construct_Type_Node construct_type_node;
    construct_type_node.type_face = &type_face;
    construct_type_node.cpt_faces_node = &cpt_faces_node;
    construct_type_node.list_faces_node = &list_faces_node;

    m.update_skin();
    apply( m.skin.node_list, construct_type_node, m, f, type_node );

    if ( debug_geometry ) {
        cout << "Construction du vecteur de vecteurs type_node" << endl << endl;
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "type du noeud " << i << " dans la direction " << d << " : " << type_node[ i ][ d ] << endl;
            }
            cout << endl << endl;
        }
    }
}

/// Calcul des informations relatives a la geometrie
///-------------------------------------------------
template<class TM, class TF>
void calcul_display_geometry( TM &m, const TF &f, const bool &debug_geometry ) {

    if ( debug_geometry ) {
        cout << "--------------------------------------------------------------" << endl;
        cout << " Calcul et Affichage des informations relatives a la geometrie" << endl;
        cout << "--------------------------------------------------------------" << endl << endl;
        Vec<unsigned> cpt_children;
        Vec< Vec<unsigned> > list_children;
        construct_children( m, cpt_children, list_children, debug_geometry );

        Vec<unsigned> cpt_elems_node;
        Vec< Vec<unsigned> > list_elems_node;
        construct_elems_connected_to_node( m, cpt_elems_node, list_elems_node, debug_geometry );

        Vec<unsigned> cpt_faces_node;
        Vec< Vec<unsigned> > list_faces_node;
        construct_faces_connected_to_node( m, cpt_faces_node, list_faces_node, debug_geometry );

        Vec<unsigned> cpt_nodes_elem;
        Vec< Vec<unsigned> > list_nodes_elem;
        construct_nodes_connected_to_elem( m, cpt_nodes_elem, list_nodes_elem, debug_geometry );

        Vec<unsigned> cpt_nodes_face;
        Vec< Vec<unsigned> > list_nodes_face;
        construct_nodes_connected_to_face( m, cpt_nodes_face, list_nodes_face, debug_geometry );

        Vec<bool> correspondance_node_to_vertex_node;
        Vec<unsigned> connect_node_to_vertex_node;
        unsigned nb_vertex_nodes = match_node_to_vertex_node( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, debug_geometry );

        Vec<unsigned> cpt_elems_vertex_node;
        Vec< Vec<unsigned> > list_elems_vertex_node;
        construct_elems_connected_to_vertex_node( m, nb_vertex_nodes, correspondance_node_to_vertex_node, connect_node_to_vertex_node, cpt_elems_vertex_node, list_elems_vertex_node, debug_geometry );

        Vec<unsigned> cpt_faces_vertex_node;
        Vec< Vec<unsigned> > list_faces_vertex_node;
        construct_faces_connected_to_vertex_node( m, nb_vertex_nodes, correspondance_node_to_vertex_node, connect_node_to_vertex_node, cpt_faces_vertex_node, list_faces_vertex_node, debug_geometry );

        Vec<unsigned> cpt_vertex_nodes_elem;
        Vec< Vec<unsigned> > list_vertex_nodes_elem;
        construct_vertex_nodes_connected_to_elem( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, cpt_vertex_nodes_elem, list_vertex_nodes_elem, debug_geometry );

        Vec<unsigned> cpt_vertex_nodes_face;
        Vec< Vec<unsigned> > list_vertex_nodes_face;
        construct_vertex_nodes_connected_to_face( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, cpt_vertex_nodes_face, list_vertex_nodes_face, debug_geometry );

        Vec< Vec<unsigned> > type_face; // type_face[ k ][ d ] = 0 : face k interieur au domaine Omega dans la direction d
                                        // type_face[ k ][ d ] = 1 : face k appartenant a un bord de Dirichlet delta_1_Omega dans la direction d
                                        // type_face_dim[ k ][ d ] = 2 : face k appartenant a un bord de Neumann delta_2_Omega dans la direction d
        construct_type_face( m, f, type_face, debug_geometry );

        Vec< Vec<unsigned> > type_node; // type_node[ i ][ d ] = 0 : noeud i interieur au domaine Omega dans la direction d
                                        // type_node[ i ][ d ] = 1 : noeud i sur un bord de Dirichlet delta_1_Omega dans la direction d
                                        // type_node[ i ][ d ] = 2 : noeud i sur un bord de Neumann delta_2_Omega dans la direction d
                                        // type_node[ i ][ d ] = 12 : noeud i a la separation entre delta_1_Omega et delat_2_Omega dans la direction d
        construct_type_node( m, f, type_face, type_node, debug_geometry );

    }
}

#endif // Calcul_geometry_h
