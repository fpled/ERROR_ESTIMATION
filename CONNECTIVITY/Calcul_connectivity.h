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

#include "Connectivity.h"

using namespace LMT;
using namespace std;

/// Compteur du nb de faces connectees a l'element n du maillage : child_cpt[ n ]
/// Liste des faces connectees a l'element n du maillage : child_list[ n ]
/// --------------------------------------------------------------------------------
template<class TM>
void construct_child( TM &m, Vec<unsigned> &child_cpt, Vec< Vec<unsigned> > &child_list, const bool debug_mesh = false ) {
    child_cpt.resize( m.elem_list.size() );
    child_cpt.set( 0 );
    child_list.resize( m.elem_list.size() );

    m.update_skin();
    m.update_elem_children();
    apply( m.elem_list, Counter_Child(), m, child_cpt, child_list ); // child_cpt[ n ] contient le nb de faces connectees a l'element n du maillage // child_list[ n ] contient le numero des faces connectees a l'element n du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des faces connectees aux elements du maillage" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "nb de faces connectees a l'element " << n << " = " << child_cpt[ n ] << endl;
            cout << "liste des faces connectees a l'element " << n << " = " << child_list[ n ] << endl << endl;
        }
    }
}

/// Compteur du nb d'elements connectes au noeud i du maillage : elem_cpt_node[ i ]
/// Liste des elements connectes au noeud i du maillage : elem_list_node[ i ]
/// --------------------------------------------------------------------------------
template<class TM>
void construct_elems_connected_to_node( const TM &m, Vec<unsigned> &elem_cpt_node, Vec< Vec<unsigned> > &elem_list_node, const bool debug_mesh = false ) {
    elem_cpt_node.resize( m.node_list.size() );
    elem_cpt_node.set( 0 );
    elem_list_node.resize( m.node_list.size() );

    apply( m.elem_list, Counter(), elem_cpt_node, elem_list_node ); // elem_cpt_node[ i ] contient le nb d'elements connectes au noeud i du maillage // elem_list_node[ i ] contient le numero des elements connectes au noeud i du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des elements connectes aux noeuds du maillage" << endl << endl;
        for (unsigned i=0;i<m.node_list.size();++i) {
            cout << "nb d'elements connectes au noeud " << i << " = " << elem_cpt_node[ i ] << endl;
            cout << "liste des elements connectes au noeud " << i << " = " << elem_list_node[ i ] << endl << endl;
        }
    }
}

/// Compteur du nb de faces connectees au noeud i du maillage : face_cpt_node[ i ]
/// Liste des faces connectees au noeud i du maillage : face_list_node[ i ]
/// -------------------------------------------------------------------------------
template<class TM>
void construct_faces_connected_to_node( TM &m, Vec<unsigned> &face_cpt_node, Vec< Vec<unsigned> > &face_list_node, const bool debug_mesh = false ) {
    face_cpt_node.resize( m.node_list.size() );
    face_cpt_node.set( 0 );
    face_list_node.resize( m.node_list.size() );

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, Counter(), face_cpt_node, face_list_node ); // face_cpt_node[ i ] contient le nb de faces connectees au noeud i du maillage // face_list_node[ i ] contient le numero des faces connectees au noeud i du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des faces connectees aux noeuds du maillage" << endl << endl;
        for (unsigned i=0;i<m.node_list.size();++i) {
            cout << "nb de faces connectees au noeud " << i << " = " << face_cpt_node[ i ] << endl;
            cout << "liste des faces connectees au noeud " << i << " = " << face_list_node[ i ] << endl << endl;
        }
    }
}

/// Compteur du nb de noeuds connectes a l'element n du maillage : node_cpt_elem[ n ]
/// Liste des noeuds connectes a l'element n du maillage : node_list_elem[ n ]
/// ----------------------------------------------------------------------------------
template<class TM>
void construct_nodes_connected_to_elem( const TM &m, Vec<unsigned> &node_cpt_elem, Vec< Vec<unsigned> > &node_list_elem, const bool debug_mesh = false ) {
    node_cpt_elem.resize( m.elem_list.size() );
    node_cpt_elem.set( 0 );
    node_list_elem.resize( m.elem_list.size() );

    apply( m.elem_list, Counter_Node(), node_cpt_elem, node_list_elem ); // node_cpt_elem[ n ] contient le nb de noeuds connectes a l'element n du maillage // node_list_elem[ n ] contient le numero des noeuds connectes a l'element n du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des noeuds connectes aux elements du maillage" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "nb de noeuds connectes a l'element " << n << " = " << node_cpt_elem[ n ] << endl;
            cout << "liste des noeuds connectes a l'element " << n << " = " << node_list_elem[ n ] << endl << endl;
        }
    }
}

/// Compteur du nb de noeuds connectes a la face k du maillage : node_cpt_face[ k ]
/// Liste des noeuds connectes a la face k du maillage : node_list_face[ k ]
/// --------------------------------------------------------------------------------
template<class TM>
void construct_nodes_connected_to_face( TM &m, Vec<unsigned> &node_cpt_face, Vec< Vec<unsigned> > &node_list_face, const bool debug_mesh  = false ) {
    node_cpt_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    node_cpt_face.set( 0 );
    node_list_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, Counter_Node(), node_cpt_face, node_list_face ); // node_cpt_face[ k ] contient le nb de noeuds connectes a la face k du maillage // node_list_face[ k ] contient le numero des noeuds connectes a la face k du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des noeuds connectes aux faces du maillage" << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            cout << "nb de noeuds connectes a la face " << k << " = " << node_cpt_face[ k ] << endl;
            cout << "liste des noeuds connectes a la face " << k << " = " << node_list_face[ k ] << endl << endl;
        }
    }
}

/// Construction du vecteur ( analogue de m.node_list pour les noeuds ) contenant les numeros des noeuds sommets dans le maillage : connect_node_to_vertex_node[ i ] pour chaque noeud i du maillage
/// Compteur du nb de noeuds sommets dans le maillage : nb_vertex_nodes
/// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM>
unsigned match_node_to_vertex_node( const TM &m, Vec<bool> &correspondance_node_to_vertex_node, Vec<unsigned> &connect_node_to_vertex_node, const bool debug_mesh = false ) {
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

    if ( debug_mesh ) {
        cout << "Correspondance entre noeuds et noeuds sommets du maillage" << endl << endl;
        cout << "nombre de noeuds sommets dans le maillage = " << nb_vertex_nodes << endl;
        cout << "correspondance entre noeuds et noeuds sommets du maillage = " << connect_node_to_vertex_node << endl << endl;
    }
    return nb_vertex_nodes;
}

/// Compteur du nb d'elements connectes au noeud sommet j du maillage : elem_cpt_vertex_node[ j ]
/// Liste des elements connectes au noeud i du maillage : elem_list_vertex_node[ j ]
/// ----------------------------------------------------------------------------------------------
template<class TM>
void construct_elems_connected_to_vertex_node( const TM &m, const unsigned &nb_vertex_nodes, const Vec<bool> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &elem_cpt_vertex_node, Vec< Vec<unsigned> > &elem_list_vertex_node, const bool debug_mesh = false ) {
    elem_cpt_vertex_node.resize( nb_vertex_nodes );
    elem_cpt_vertex_node.set( 0 );
    elem_list_vertex_node.resize( nb_vertex_nodes );

    Counter_Vertex counter_vertex;
    counter_vertex.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;

    apply( m.elem_list, counter_vertex, connect_node_to_vertex_node, elem_cpt_vertex_node, elem_list_vertex_node ); // elem_cpt_vertex_node[ j ] contient le nb d'elements connectes au noeud sommet j du maillage // elem_list_vertex_node[ j ] contient le numero des elements connectes au noeud sommet j du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des elements connectes aux noeuds sommets du maillage" << endl << endl;
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "nb d'elements connectes au noeud sommet " << j << " = " << elem_cpt_vertex_node[ j ] << endl;
            cout << "liste des elements connectes au noeud sommet " << j << " = " << elem_list_vertex_node[ j ] << endl << endl;
        }
    }
}

/// Compteur du nb de faces connectees au noeud sommet j du maillage : face_cpt_vertex_node[ j ]
/// Liste des faces connectees au noeud i du maillage : face_list_vertex_node[ j ]
/// ---------------------------------------------------------------------------------------------
template<class TM>
void construct_faces_connected_to_vertex_node( TM &m, const unsigned &nb_vertex_nodes, const Vec<bool> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &face_cpt_vertex_node, Vec< Vec<unsigned> > &face_list_vertex_node, const bool debug_mesh = false ) {
    face_cpt_vertex_node.resize( nb_vertex_nodes );
    face_cpt_vertex_node.set( 0 );
    face_list_vertex_node.resize( nb_vertex_nodes );

    Counter_Vertex counter_vertex;
    counter_vertex.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, counter_vertex, connect_node_to_vertex_node, face_cpt_vertex_node, face_list_vertex_node ); // face_cpt_vertex_node[ j ] contient le nb de faces connectees au noeud sommet j du maillage // face_list_vertex_node[ j ] contient le numero des faces connectees au noeud sommet j du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des faces connectees aux noeuds sommets du maillage" << endl << endl;
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "nb de faces connectees au noeud sommet " << j << " = " << face_cpt_vertex_node[ j ] << endl;
            cout << "liste des faces connectees au noeud sommet " << j << " = " << face_list_vertex_node[ j ] << endl << endl;
        }
    }
}

/// Compteur du nb de noeuds sommets connectes a l'element n du maillage : vertex_node_cpt_elem[ n ]
/// Liste des noeuds sommets connectes a l'element n du maillage : vertex_node_list_elem[ n ]
/// -------------------------------------------------------------------------------------------------
template<class TM>
void construct_vertex_nodes_connected_to_elem( const TM &m, const Vec<bool> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &vertex_node_cpt_elem, Vec< Vec<unsigned> > &vertex_node_list_elem, const bool debug_mesh = false ) {
    vertex_node_cpt_elem.resize( m.elem_list.size() );
    vertex_node_cpt_elem.set( 0 );
    vertex_node_list_elem.resize( m.elem_list.size() );

    Counter_Vertex_Node counter_vertex_node;
    counter_vertex_node.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;

    apply( m.elem_list, counter_vertex_node, connect_node_to_vertex_node, vertex_node_cpt_elem, vertex_node_list_elem ); // vertex_node_cpt_elem[ n ] contient le nb de noeuds sommets connectes a l'element n du maillage // vertex_node_list_elem[ n ] contient le numero des noeuds sommets connectes a l'element n du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des noeuds sommets connectes aux elements du maillage" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "nb de noeuds sommets connectes a l'element " << n << " = " << vertex_node_cpt_elem[ n ] << endl;
            cout << "liste des noeuds sommets connectes a l'element " << n << " = " << vertex_node_list_elem[ n ] << endl << endl;
        }
    }
}

/// Compteur du nb de noeuds sommets connectes a la face k du maillage : vertex_node_cpt_face[ k ]
/// Liste des noeuds sommets connectes a la face k du maillage : vertex_node_list_face[ k ]
/// -----------------------------------------------------------------------------------------------
template<class TM>
void construct_vertex_nodes_connected_to_face( TM &m, const Vec<bool> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &vertex_node_cpt_face, Vec< Vec<unsigned> > &vertex_node_list_face, const bool debug_mesh = false ) {
    vertex_node_cpt_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    vertex_node_cpt_face.set( 0 );
    vertex_node_list_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );

    Counter_Vertex_Node counter_vertex_node;
    counter_vertex_node.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, counter_vertex_node, connect_node_to_vertex_node, vertex_node_cpt_face, vertex_node_list_face ); // vertex_node_cpt_face[ k ] contient le nb de noeuds sommets connectes a la face k du maillage // vertex_node_list_face[ k ] contient le numero des noeuds sommets connectes a la face k du maillage

    if ( debug_mesh ) {
        cout << "Compteur et liste des noeuds sommets connectes aux faces du maillage" << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            cout << "nb de noeuds sommets connectes a la face " << k << " = " << vertex_node_cpt_face[ k ] << endl;
            cout << "liste des noeuds sommets connectes a la face " << k << " = " << vertex_node_list_face[ k ] << endl << endl;
        }
    }
}

/// Construction du vecteur de vecteurs face_type :
/// face_type[ k ][ d ] = 0 : face k interieur a Omega dans la direction d
/// face_type[ k ][ d ] = 1 : face k appartenant a un bord delta_1_Omega dans la direction d
/// face_type_dim[ k ][ d ] = 2 : face k appartenant a un bord delta_2_Omega dans la direction d
/// --------------------------------------------------------------------------------------------
template<class TM, class TF>
void construct_face_type( TM &m, const TF &f, Vec< Vec<unsigned> > &face_type, const bool debug_mesh = false ) {

    static const unsigned dim = TM::dim;

    face_type.resize( m.sub_mesh(Number<1>()).elem_list.size() );

    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        face_type[ k ].resize( dim );
        face_type[ k ].set( 0 ) ;
    }

    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, Construct_Face_Type(), m, f, face_type );

    if ( debug_mesh ) {
        cout << "Construction du vecteur de vecteurs face_type" << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "type de la face " << k << " dans la direction " << d << " = " << face_type[ k ][ d ] << endl << endl;
            }
        }
    }
}

/// Construction du vecteur de vecteurs node_type :
/// node_type[ i ][ d ] = 0 : noeud i interieur a Omega dans la direction d
/// node_type[ i ][ d ] = 1 : noeud i sur un bord delta_1_Omega dans la direction d
/// node_type[ i ][ d ] = 2 : noeud i sur un bord delta_2_Omega dans la direction d
/// node_type[ i ][ d ] = 12 : noeud i a la separation entre delta_1_Omega et delat_2_Omega dans la direction d
/// -----------------------------------------------------------------------------------------------------------
template<class TM, class TF>
void construct_node_type( TM &m, const TF &f, const Vec< Vec<unsigned> > &face_type, Vec< Vec<unsigned> > &node_type, const bool debug_mesh = false ) {

    static const unsigned dim = TM::dim;

    node_type.resize( m.node_list.size() );
    for (unsigned i=0;i<m.node_list.size();++i) {
        node_type[ i ].resize( dim );
        node_type[ i ].set( 0 );
    }

    Vec<unsigned> face_cpt_node;
    Vec< Vec<unsigned> > face_list_node;
    construct_faces_connected_to_node( m, face_cpt_node, face_list_node, debug_mesh );

    Construct_Node_Type construct_node_type;
    construct_node_type.face_type = &face_type;
    construct_node_type.face_cpt_node = &face_cpt_node;
    construct_node_type.face_list_node = &face_list_node;

    m.update_skin();
    apply( m.skin.node_list, construct_node_type, m, f, node_type );

    if ( debug_mesh ) {
        cout << "Construction du vecteur de vecteurs node_type" << endl << endl;
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "type du noeud " << i << " dans la direction " << d << " = " << node_type[ i ][ d ] << endl << endl;
            }
        }
    }
}

/// Affichage des informations relatives a la connectivite du maillage
/// ------------------------------------------------------------------
template<class TM, class TF>
void display_mesh_connectivity( TM &m, const TF &f, const bool debug_mesh = false ) {

    cout << "----------------------------------------------------" << endl;
    cout << "Informations relatives a la connectivite du maillage" << endl;
    cout << "----------------------------------------------------" << endl << endl;
    Vec<unsigned> child_cpt;
    Vec< Vec<unsigned> > child_list;
    construct_child( m, child_cpt, child_list, debug_mesh );

    Vec<unsigned> elem_cpt_node;
    Vec< Vec<unsigned> > elem_list_node;
    construct_elems_connected_to_node( m, elem_cpt_node, elem_list_node, debug_mesh );

    Vec<unsigned> face_cpt_node;
    Vec< Vec<unsigned> > face_list_node;
    construct_faces_connected_to_node( m, face_cpt_node, face_list_node, debug_mesh );

    Vec<unsigned> node_cpt_elem;
    Vec< Vec<unsigned> > node_list_elem;
    construct_nodes_connected_to_elem( m, node_cpt_elem, node_list_elem, debug_mesh );

    Vec<unsigned> node_cpt_face;
    Vec< Vec<unsigned> > node_list_face;
    construct_nodes_connected_to_face( m, node_cpt_face, node_list_face, debug_mesh );

    Vec<bool> correspondance_node_to_vertex_node;
    Vec<unsigned> connect_node_to_vertex_node;
    unsigned nb_vertex_nodes = match_node_to_vertex_node( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, debug_mesh );

    Vec<unsigned> elem_cpt_vertex_node;
    Vec< Vec<unsigned> > elem_list_vertex_node;
    construct_elems_connected_to_vertex_node( m, nb_vertex_nodes, correspondance_node_to_vertex_node, connect_node_to_vertex_node, elem_cpt_vertex_node, elem_list_vertex_node, debug_mesh );

    Vec<unsigned> face_cpt_vertex_node;
    Vec< Vec<unsigned> > face_list_vertex_node;
    construct_faces_connected_to_vertex_node( m, nb_vertex_nodes, correspondance_node_to_vertex_node, connect_node_to_vertex_node, face_cpt_vertex_node, face_list_vertex_node, debug_mesh );

    Vec<unsigned> vertex_node_cpt_elem;
    Vec< Vec<unsigned> > vertex_node_list_elem;
    construct_vertex_nodes_connected_to_elem( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, vertex_node_cpt_elem, vertex_node_list_elem, debug_mesh );

    Vec<unsigned> vertex_node_cpt_face;
    Vec< Vec<unsigned> > vertex_node_list_face;
    construct_vertex_nodes_connected_to_face( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, vertex_node_cpt_face, vertex_node_list_face, debug_mesh );

    Vec< Vec<unsigned> > face_type; // face_type[ k ][ d ] = 0 : face k interieur au domaine Omega dans la direction d
    // face_type[ k ][ d ] = 1 : face k appartenant a un bord de Dirichlet delta_1_Omega dans la direction d
    // face_type_dim[ k ][ d ] = 2 : face k appartenant a un bord de Neumann delta_2_Omega dans la direction d
    construct_face_type( m, f, face_type, debug_mesh );

    Vec< Vec<unsigned> > node_type; // node_type[ i ][ d ] = 0 : noeud i interieur au domaine Omega dans la direction d
    // node_type[ i ][ d ] = 1 : noeud i sur un bord de Dirichlet delta_1_Omega dans la direction d
    // node_type[ i ][ d ] = 2 : noeud i sur un bord de Neumann delta_2_Omega dans la direction d
    // node_type[ i ][ d ] = 12 : noeud i au bord entre delta_1_Omega et delat_2_Omega dans la direction d
    construct_node_type( m, f, face_type, node_type, debug_mesh );

}

#endif // Calcul_geometry_h
