//
// C++ Interface: Geometry
//
// Description: informations relatives a la geometrie
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef Geometry_h
#define Geometry_h

#include "../LMT/include/mesh/triangle.h"
#include "../LMT/include/mesh/triangle_6.h"
#include "../LMT/include/mesh/quad.h"
#include "../LMT/include/mesh/tetra.h"
#include "../LMT/include/mesh/tetra_10.h"

using namespace LMT;
using namespace std;

/// Compteur du nb de faces connectees a l'element n du maillage : cpt
/// Liste des faces connectees a l'element n du maillage : list
///-------------------------------------------------------------------
struct Counter_Children {
    template<class TE, class TM> void operator()( const TE &elem, const TM &m, Vec<unsigned> &cpt, Vec< Vec<unsigned> > &list ) const {
        cpt[ elem.number ] = NbChildrenElement<typename TE::NE,1>::res;
        for (unsigned k=0;k<NbChildrenElement<typename TE::NE,1>::res;++k) {
            list[ elem.number ].push_back( m.get_children_of( elem,Number<1>() )[ k ]->number );
        }
    }
};

/// Compteur du nb d'elements ( ou de faces ) connectes au noeud i du maillage : cpt
/// Liste des elements (ou des face ) connectes au noeud i du maillage : list
///---------------------------------------------------------------------------------
struct Counter {
    template<class TE> void operator()( const TE &elem, Vec<unsigned> &cpt, Vec< Vec<unsigned> > &list ) const {
        for (unsigned i=0;i<elem.nb_nodes;++i) {
            cpt[ elem.node( i )->number ]++;
            list[ elem.node( i )->number ].push_back( elem.number );
        }
    }
};

/// Construction de la correspondance entre noeuds et noeuds sommets du maillage
///-----------------------------------------------------------------------------
template<class TE, class BV> 
void construct_correspondance_node_to_vertex_node( const TE &elem, BV &correspondance_node_to_vertex_node ) {}

struct Construct_Correspondance_Node_To_Vertex_Node {
    template<class TE> void operator()( const TE &elem, Vec<bool> &correspondance_node_to_vertex_node ) const {
        construct_correspondance_node_to_vertex_node( elem, correspondance_node_to_vertex_node );
    }
};

/// Compteur du nb d'elements ( ou de faces ) connectes au noeud sommet j du maillage : cpt
/// Liste des elements (ou des face ) connectes au noeud sommet j du maillage : list
///----------------------------------------------------------------------------------------
struct Counter_Vertex {
    const Vec<bool>* correspondance_node_to_vertex_node;
    template<class TE> void operator()( const TE &elem, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &cpt, Vec< Vec<unsigned> > &list ) const {
        for (unsigned i=0;i<elem.nb_nodes;++i) {
            if ( (*correspondance_node_to_vertex_node)[ elem.node( i )->number ] == true ) {
                cpt[ connect_node_to_vertex_node[ elem.node( i )->number ] ]++;
                list[ connect_node_to_vertex_node[ elem.node( i )->number ] ].push_back( elem.number );
            }
        }
    }
};

/// Compteur du nb de noeuds sommets connectes a l'element e ( ou a la face face ) du maillage : cpt
/// Liste des noeuds sommets connectes a l'element e ( ou a la face face ) du maillage : list
///-------------------------------------------------------------------------------------------------
struct Counter_Vertex_Nodes {
    const Vec<bool>* correspondance_node_to_vertex_node;
    template<class TE> void operator()( const TE &elem, const Vec<unsigned> &connect_node_to_vertex_node, Vec<unsigned> &cpt, Vec< Vec<unsigned> > &list ) const {
        unsigned elem_nb_nodes = TE::nb_nodes; // nb de noeuds de l'element elem
        for (unsigned i=0;i<elem_nb_nodes;++i) {
            if ( (*correspondance_node_to_vertex_node)[ elem.node( i )->number ] == true ) {
                cpt[ elem.number ]++;
                list[ elem.number ].push_back( connect_node_to_vertex_node[ elem.node( i )->number ] );
            }
        }
    }
};

/// Compteur du nb de noeuds connectes a l'element e ( ou a la face face ) du maillage : cpt
/// Liste des noeuds connectes a l'element e ( ou a la face face ) du maillage : list
///-----------------------------------------------------------------------------------------
struct Counter_Nodes {
    template<class TE> void operator()( const TE &elem, Vec<unsigned> &cpt, Vec< Vec<unsigned> > &list ) const {
        unsigned elem_nb_nodes = TE::nb_nodes; // nb de noeuds de l'element elem
        for (unsigned i=0;i<elem_nb_nodes;++i) {
            cpt[ elem.number ]++;
            list[ elem.number ].push_back( elem.node( i )->number );
//            PRINT( elem.node( i )->number );
//            PRINT( elem.node( i )->pos );
        }
    }
};

/// Construction de la correspondance entre noeuds sommets et noeuds sommets appartenant a delta_Omega : correspondance_vertex_node_to_skin_vertex_node
///----------------------------------------------------------------------------------------------------------------------------------------------------
struct Constr_Correspondance_Vertex_Node_To_Skin_Vertex_Node {
    Vec<bool>* correspondance_node_to_vertex_node;
    template<class TE> void operator()( const TE &skin_node, const Vec<unsigned> &correspondance_node_to_vertex_node, const Vec<unsigned> &connect_node_to_vertex_node, Vec<bool> &correspondance_vertex_node_to_skin_vertex_node ) const {
        if ( correspondance_node_to_vertex_node[ skin_node.number ] == true ) { // si le noeud skin_node est un noeud sommet du maillage
            correspondance_vertex_node_to_skin_vertex_node[ connect_node_to_vertex_node[ skin_node.number ] ] = 0;
        }
    }
};

/// Construction du vecteur type_face
///----------------------------------
struct Construct_Type_Face {
    template<class TE, class TM, class TF> void operator()( const TE &child_elem, const TM &m, const TF &f, Vec< Vec<unsigned> > &type_face ) const {
        if ( m.sub_mesh(Number<1>()).get_parents_of( child_elem ).size() != 2 ) {
            for (unsigned d=0;d< TE::dim;++d) {
                unsigned cpt = 0;
                for (unsigned i=0;i<child_elem.nb_nodes;++i) {
                    if ( f.constrained_nodes_in_dim( d )[ child_elem.node( i )->number ] )
                        cpt++;
                }
                if ( cpt == child_elem.nb_nodes )
                    type_face[ child_elem.number ][ d ] = 1;
                else
                    type_face[ child_elem.number ][ d ] = 2;
            }
        }
    }
};

/// Construction du vecteur type_node
///----------------------------------
struct Construct_Type_Node {
    const Vec< Vec<unsigned> >* type_face;
    const Vec<unsigned>* cpt_faces_node;
    const Vec< Vec<unsigned> >* list_faces_node;
    template<class TE, class TM, class TF> void operator()( const TE &skin_node, const TM &m, const TF &f, Vec< Vec<unsigned> > &type_node ) const {
//         PRINT( skin_node.number ); // numero du noeud de bord skin_node dans le maillage global 
        for (unsigned d =0;d<TE::dim;++d) {
            unsigned cpt_type_face_1 = 0;
            unsigned cpt_type_face_2 = 0;
            for (unsigned k=0;k<(*cpt_faces_node)[ skin_node.number ];++k) {
                if ( (*type_face)[ (*list_faces_node)[ skin_node.number ][ k ] ][ d ] == 1 ) {
                    cpt_type_face_1++;
                }
                if ( (*type_face)[ (*list_faces_node)[ skin_node.number ][ k ] ][ d ] == 2 ) {
                    cpt_type_face_2++;
                }
            }
            if ( cpt_type_face_1 == 0 and cpt_type_face_2 != 0 ) {
                type_node[ skin_node.number ][ d ] = 2;
            }
            if ( cpt_type_face_1 != 0 and cpt_type_face_2 == 0 ) {
                type_node[ skin_node.number ][ d ] = 1;
            }
            if ( cpt_type_face_1 != 0 and cpt_type_face_2 != 0 ) {
                type_node[ skin_node.number ][ d ] = 12;
            }
        }
    }
};

/// Construction du maillage de ref pour un maillage comportant des Triangle : subdivision de chaque Triangle du maillage en 4 Triangle dans le maillage de ref
///------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE,class TM>
void divide_element_Triangle(TE &e, const TM &m, TM &m_ref) {}

template<class TN,class TNG,class TD,unsigned NET,class TM>
void divide_element_Triangle(Element<Triangle,TN,TNG,TD,NET> &e, const TM &m, TM &m_ref) {
    typedef typename TM::TNode TNode;
    BestialNodeAdder<TM> ban; ban.m = &m_ref; ban.prec = 1e-6;
    Vec<TNode*,3> node_Triangle;
    for (unsigned i=0;i<3;++i) {
        node_Triangle[i] = ban.get_node( e.node(i)->pos );
    }
    Vec<TNode*,3> node_center_Bar;
    for (unsigned i=0;i<3;++i) {
        node_center_Bar[i] = ban.get_node( center(*m.get_children_of( e, Number<1>() )[i]) );
    }
    DM::copy( e, *m_ref.add_element( Triangle(), TN(), node_Triangle[0], node_center_Bar[0], node_center_Bar[2] ) );
    DM::copy( e, *m_ref.add_element( Triangle(), TN(), node_center_Bar[0], node_Triangle[1], node_center_Bar[1] ) );
    DM::copy( e, *m_ref.add_element( Triangle(), TN(), node_center_Bar[2], node_center_Bar[1], node_Triangle[2] ) );
    DM::copy( e, *m_ref.add_element( Triangle(), TN(), node_center_Bar[0], node_center_Bar[1], node_center_Bar[2] ) );
}

struct Construction_Ref_Mesh_Triangle {
    template<class TE, class TM> void operator()( TE &elem, TM &m, TM &m_ref ) const {
        m.update_elem_children();
        m.update_elem_children( Number<1>() );
        divide_element_Triangle( elem, m, m_ref );
    }
};

/// Construction du maillage de ref pour un maillage comportant des Quad : subdivision de chaque Quad du maillage en 4 Quad dans le maillage de ref
///------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE,class TM>
void divide_element_Quad(TE &e, const TM &m, TM &m_ref) {}

template<class TN,class TNG,class TD,unsigned NET,class TM>
void divide_element_Quad(Element<Quad,TN,TNG,TD,NET> &e, const TM &m, TM &m_ref) {
    typedef typename TM::TNode TNode;
    BestialNodeAdder<TM> ban; ban.m = &m_ref; ban.prec = 1e-6;
    TNode *node_center_Quad = ban.get_node( center(e) );
    Vec<TNode*,4> node_Quad;
    for (unsigned i=0;i<4;++i) {
        node_Quad[i] = ban.get_node( e.node(i)->pos );
    }
    Vec<TNode*,4> node_center_Bar;
    for (unsigned i=0;i<4;++i) {
        node_center_Bar[i] = ban.get_node( center(*m.get_children_of( e, Number<1>() )[i]) );
    }
    DM::copy( e, *m_ref.add_element( Quad(), TN(), node_Quad[0], node_center_Bar[0], node_center_Quad, node_center_Bar[3] ) );
    DM::copy( e, *m_ref.add_element( Quad(), TN(), node_center_Bar[0], node_Quad[1], node_center_Bar[1], node_center_Quad ) );
    DM::copy( e, *m_ref.add_element( Quad(), TN(), node_center_Quad, node_center_Bar[1], node_Quad[2], node_center_Bar[2] ) );
    DM::copy( e, *m_ref.add_element( Quad(), TN(), node_center_Bar[3], node_center_Quad, node_center_Bar[2], node_Quad[3] ) );
}

struct Construction_Ref_Mesh_Quad {
    template<class TE, class TM> void operator()( TE &elem, TM &m, TM &m_ref ) const {
        m.update_elem_children();
        m.update_elem_children( Number<1>() );
        divide_element_Quad( elem, m, m_ref );
    }
};

/// Construction du maillage de ref pour un maillage comportant des Tetra : subdivision de chaque Tetra du maillage en 8 Tetra dans le maillage de ref
///---------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE,class TM>
void divide_element_Tetra(TE &e, const TM &m, TM &m_ref) {}

template<class TN,class TNG,class TD,unsigned NET,class TM>
void divide_element_Tetra(Element<Tetra,TN,TNG,TD,NET> &e, const TM &m, TM &m_ref) {
    typedef typename TM::TNode TNode;
    BestialNodeAdder<TM> ban; ban.m = &m_ref; ban.prec = 1e-6;
    Vec<TNode*,4> node_Tetra;
    for (unsigned i=0;i<4;++i) {
        node_Tetra[i] = ban.get_node( e.node(i)->pos );
    }
    Vec<TNode*,6> node_center_Bar;
    for (unsigned i=0;i<6;++i) {
        node_center_Bar[i] = ban.get_node( center(*m.get_children_of( e, Number<2>() )[i]) );
    }
    DM::copy( e, *m_ref.add_element( Tetra(), TN(), node_center_Bar[2], node_center_Bar[4], node_center_Bar[5], node_Tetra[3] ) );
    DM::copy( e, *m_ref.add_element( Tetra(), TN(), node_Tetra[0], node_center_Bar[0], node_center_Bar[1], node_center_Bar[2] ) );
    DM::copy( e, *m_ref.add_element( Tetra(), TN(), node_center_Bar[1], node_center_Bar[3], node_Tetra[2], node_center_Bar[5] ) );
    DM::copy( e, *m_ref.add_element( Tetra(), TN(), node_center_Bar[0], node_Tetra[1], node_center_Bar[3], node_center_Bar[4] ) );
    DM::copy( e, *m_ref.add_element( Tetra(), TN(), node_center_Bar[2], node_center_Bar[0], node_center_Bar[1], node_center_Bar[4] ) );
    DM::copy( e, *m_ref.add_element( Tetra(), TN(), node_center_Bar[1], node_center_Bar[3], node_center_Bar[5], node_center_Bar[4] ) );
    DM::copy( e, *m_ref.add_element( Tetra(), TN(), node_center_Bar[2], node_center_Bar[1], node_center_Bar[5], node_center_Bar[4] ) );
    DM::copy( e, *m_ref.add_element( Tetra(), TN(), node_center_Bar[0], node_center_Bar[3], node_center_Bar[1], node_center_Bar[4] ) );
}

struct Construction_Ref_Mesh_Tetra {
    template<class TE, class TM> void operator()( TE &elem, TM &m, TM &m_ref ) const {
        m.update_elem_children();
        m.update_elem_children( Number<2>() );
        divide_element_Tetra( elem, m, m_ref );
    }
};

/// Construction du maillage de ref pour un maillage comportant des Hexa : subdivision de chaque Hexa du maillage en 8 Hexa dans le maillage de ref
///------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE,class TM>
void divide_element_Hexa(TE &e, const TM &m, TM &m_ref) {}

template<class TN,class TNG,class TD,unsigned NET,class TM>
void divide_element_Hexa(Element<Hexa,TN,TNG,TD,NET> &e, const TM &m, TM &m_ref) {
    typedef typename TM::TNode TNode;
    BestialNodeAdder<TM> ban; ban.m = &m_ref; ban.prec = 1e-6;
    TNode *node_center_Hexa = ban.get_node( center(e) );
    Vec<TNode*,8> node_Hexa;
    for (unsigned i=0;i<8;++i) {
        node_Hexa[i] = ban.get_node( e.node(i)->pos );
    }
    Vec<TNode*,6> node_center_Quad;
    for (unsigned i=0;i<6;++i) {
        node_center_Quad[i] = ban.get_node( center(*m.get_children_of( e, Number<1>() )[i]) );
    }
    Vec<TNode*,12> node_center_Bar;
    for (unsigned i=0;i<12;++i) {
        node_center_Bar[i] = ban.get_node( center(*m.get_children_of( e, Number<2>() )[i]) );
    }
    DM::copy( e, *m_ref.add_element( Hexa(), TN(), node_Hexa[0], node_center_Bar[0], node_center_Quad[0], node_center_Bar[3], node_center_Bar[8], node_center_Quad[2], node_center_Hexa, node_center_Quad[4] ) );
    DM::copy( e, *m_ref.add_element( Hexa(), TN(), node_center_Bar[0], node_Hexa[1], node_center_Bar[1], node_center_Quad[0], node_center_Quad[2], node_center_Bar[9], node_center_Quad[5], node_center_Hexa ) );
    DM::copy( e, *m_ref.add_element( Hexa(), TN(), node_center_Bar[3], node_center_Quad[0], node_center_Bar[2], node_Hexa[3], node_center_Quad[4], node_center_Hexa, node_center_Quad[3], node_center_Bar[11] ) );
    DM::copy( e, *m_ref.add_element( Hexa(), TN(), node_center_Quad[0], node_center_Bar[1], node_Hexa[2], node_center_Bar[2], node_center_Hexa, node_center_Quad[5], node_center_Bar[10], node_center_Quad[3] ) );
    DM::copy( e, *m_ref.add_element( Hexa(), TN(), node_center_Bar[8], node_center_Quad[2], node_center_Hexa, node_center_Quad[4], node_Hexa[4], node_center_Bar[4], node_center_Quad[1], node_center_Bar[7] ) );
    DM::copy( e, *m_ref.add_element( Hexa(), TN(), node_center_Quad[2], node_center_Bar[9], node_center_Quad[5], node_center_Hexa, node_center_Bar[4], node_Hexa[5], node_center_Bar[5], node_center_Quad[1] ) );
    DM::copy( e, *m_ref.add_element( Hexa(), TN(), node_center_Quad[4], node_center_Hexa, node_center_Quad[3], node_center_Bar[11], node_center_Bar[7], node_center_Quad[1], node_center_Bar[6], node_Hexa[7] ) );
    DM::copy( e, *m_ref.add_element( Hexa(), TN(), node_center_Hexa, node_center_Quad[5], node_center_Bar[10], node_center_Quad[3], node_center_Quad[1], node_center_Bar[5], node_Hexa[6], node_center_Bar[6] ) );
}

struct Construction_Ref_Mesh_Hexa {
    template<class TE, class TM> void operator()( TE &elem, TM &m, TM &m_ref ) const {
        m.update_elem_children();
        m.update_elem_children( Number<2>() );
        divide_element_Hexa( elem, m, m_ref );
    }
};

/// Construction du maillage de ref local par subdivision des elements dans list_elems en 4 pour Triangle et Quad, et en 8 pour Tetra et Hexa
///------------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class T>
struct Divide_element_locally {
    Divide_element_locally( TM &m, const Vec<T> &list_elems ) : ptr_m( &m ), ptr_list_elems( &list_elems ) {} 
    
    template< class NameBehavior, class TNode_, class TData, unsigned num_in_elem_list_ >
    bool operator() ( Element< Triangle, NameBehavior, TNode_, TData, num_in_elem_list_ > &e ) {
       if ( find( ptr_list_elems, _1 == e.number ) ) {
            typedef typename Element< Triangle, NameBehavior, TNode_, TData, num_in_elem_list_ >::TNode TNode;
            Vec< TNode*, 3 > vn;
            ptr_m->update_elem_children();
            BestialNodeAdder<TM> ban; ban.m = ptr_m; ban.prec = 1e-6;
            Vec<TNode*,3> node_center_Bar;
            for (unsigned i=0;i<3;++i) {
                node_center_Bar[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<1>())[i]) ) );
            }
            /// Triangle 0
            vn[ 0 ] = e.node( 0 );
            vn[ 1 ] = node_center_Bar[ 0 ];
            vn[ 2 ] = node_center_Bar[ 2 ];
            permutation_if_jac_neg ( Triangle(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Triangle(), NameBehavior(), vn.ptr() ) );
            /// Triangle 1
            vn[ 0 ] = node_center_Bar[ 0 ];
            vn[ 1 ] = e.node( 1 );
            vn[ 2 ] = node_center_Bar[ 1 ];
            permutation_if_jac_neg ( Triangle(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Triangle(), NameBehavior(), vn.ptr() ) );
            /// Triangle 2
            vn[ 0 ] = node_center_Bar[ 2 ];
            vn[ 1 ] = node_center_Bar[ 1 ];
            vn[ 2 ] = e.node( 2 );
            permutation_if_jac_neg ( Triangle(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Triangle(), NameBehavior(), vn.ptr() ) );
            /// Triangle 3
            vn[ 0 ] = node_center_Bar[ 0 ];
            vn[ 1 ] = node_center_Bar[ 1 ];
            vn[ 2 ] = node_center_Bar[ 2 ];
            permutation_if_jac_neg ( Triangle(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Triangle(), NameBehavior(), vn.ptr() ) );
            return true;
       }
       else
           return false;
    }
    
    template< class NameBehavior, class TNode_, class TData, unsigned num_in_elem_list_ >
    bool operator() ( Element< Quad, NameBehavior, TNode_, TData, num_in_elem_list_ > &e ) {
        if ( find( ptr_list_elems, _1 == e.number ) ) {
            typedef typename Element< Quad, NameBehavior, TNode_, TData, num_in_elem_list_ >::TNode TNode;
            Vec< TNode*, 4 > vn;
            ptr_m->update_elem_children();
            BestialNodeAdder<TM> ban; ban.m = ptr_m; ban.prec = 1e-6;
            TNode* node_center_Quad = ban.get_node( center(e) );
            Vec<TNode*,4> node_center_Bar;
            for (unsigned i=0;i<4;++i) {
                node_center_Bar[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<1>())[i]) ) );
            }
            /// Quad 0
            vn[ 0 ] = e.node( 0 );
            vn[ 1 ] = node_center_Bar[ 0 ];
            vn[ 2 ] = node_center_Quad;
            vn[ 3 ] = node_center_Bar[ 3 ];
            permutation_if_jac_neg ( Quad(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Quad(), NameBehavior(), vn.ptr() ) );
            /// Quad 1
            vn[ 0 ] = node_center_Bar[ 0 ];
            vn[ 1 ] = e.node( 1 );
            vn[ 2 ] = node_center_Bar[ 1 ];
            vn[ 3 ] = node_center_Quad;
            permutation_if_jac_neg ( Quad(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Quad(), NameBehavior(), vn.ptr() ) );
            /// Quad 2
            vn[ 0 ] = node_center_Quad;
            vn[ 1 ] = node_center_Bar[ 1 ];
            vn[ 2 ] = e.node( 2 );
            vn[ 3 ] = node_center_Bar[ 2 ];
            permutation_if_jac_neg ( Quad(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Quad(), NameBehavior(), vn.ptr() ) );
            /// Quad 3
            vn[ 0 ] = node_center_Bar[ 3 ];
            vn[ 1 ] = node_center_Quad;
            vn[ 2 ] = node_center_Bar[ 2 ];
            vn[ 3 ] = e.node( 3 );
            permutation_if_jac_neg ( Quad(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Quad(), NameBehavior(), vn.ptr() ) );
            return true;
       }
       else
           return false;
    }
    
    template< class NameBehavior, class TNode_, class TData, unsigned num_in_elem_list_ >
    bool operator() ( Element< Tetra, NameBehavior, TNode_, TData, num_in_elem_list_ > &e ) {
        if ( find( ptr_list_elems, _1 == e.number ) ) {
            typedef typename Element< Tetra, NameBehavior, TNode_, TData, num_in_elem_list_ >::TNode TNode;
            Vec< TNode*, 4 > vn;
            ptr_m->update_elem_children();
            ptr_m->update_elem_children( Number<2>() );
            BestialNodeAdder<TM> ban; ban.m = ptr_m; ban.prec = 1e-6;
            Vec<TNode*,6> node_center_Bar;
            for (unsigned i=0;i<6;++i) {
                node_center_Bar[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<2>())[i]) ) );
            }
            /// Tetra 0
            vn[ 0 ] = node_center_Bar[ 2 ];
            vn[ 1 ] = node_center_Bar[ 4 ];
            vn[ 2 ] = node_center_Bar[ 5 ];
            vn[ 3 ] = e.node( 3 );
            permutation_if_jac_neg ( Tetra(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
            /// Tetra 1
            vn[ 0 ] = e.node( 0 );
            vn[ 1 ] = node_center_Bar[ 0 ];
            vn[ 2 ] = node_center_Bar[ 1 ];
            vn[ 3 ] = node_center_Bar[ 2 ];
            permutation_if_jac_neg ( Tetra(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
            /// Tetra 2
            vn[ 0 ] = node_center_Bar[ 1 ];
            vn[ 1 ] = node_center_Bar[ 3 ];
            vn[ 2 ] = e.node( 2 );
            vn[ 3 ] = node_center_Bar[ 5 ];
            permutation_if_jac_neg ( Tetra(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
            /// Tetra 3
            vn[ 0 ] = node_center_Bar[ 0 ];
            vn[ 1 ] = e.node( 1 );
            vn[ 2 ] = node_center_Bar[ 3 ];
            vn[ 3 ] = node_center_Bar[ 4 ];
            permutation_if_jac_neg ( Tetra(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
            /// Tetra 4
            vn[ 0 ] = node_center_Bar[ 2 ];
            vn[ 1 ] = node_center_Bar[ 0 ];
            vn[ 2 ] = node_center_Bar[ 1 ];
            vn[ 3 ] = node_center_Bar[ 4 ];
            permutation_if_jac_neg ( Tetra(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
            /// Tetra 5
            vn[ 0 ] = node_center_Bar[ 1 ];
            vn[ 1 ] = node_center_Bar[ 3 ];
            vn[ 2 ] = node_center_Bar[ 5 ];
            vn[ 3 ] = node_center_Bar[ 4 ];
            permutation_if_jac_neg ( Tetra(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
            /// Tetra 6
            vn[ 0 ] = node_center_Bar[ 2 ];
            vn[ 1 ] = node_center_Bar[ 1 ];
            vn[ 2 ] = node_center_Bar[ 5 ];
            vn[ 3 ] = node_center_Bar[ 4 ];
            permutation_if_jac_neg ( Tetra(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
            /// Tetra 7
            vn[ 0 ] = node_center_Bar[ 0 ];
            vn[ 1 ] = node_center_Bar[ 3 ];
            vn[ 2 ] = node_center_Bar[ 1 ];
            vn[ 3 ] = node_center_Bar[ 4 ];
            permutation_if_jac_neg ( Tetra(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
            return true;
       }
       else
           return false;
    }

    template< class NameBehavior, class TNode_, class TData, unsigned num_in_elem_list_ >
    bool operator() ( Element< Hexa, NameBehavior, TNode_, TData, num_in_elem_list_ > &e ) {
        if ( find( ptr_list_elems, _1 == e.number ) ) {
            typedef typename Element< Hexa, NameBehavior, TNode_, TData, num_in_elem_list_ >::TNode TNode;
            Vec< TNode*, 8 > vn;
            ptr_m->update_elem_children();
            ptr_m->update_elem_children( Number<2>() );
            BestialNodeAdder<TM> ban; ban.m = ptr_m; ban.prec = 1e-6;
            TNode* node_center_Hexa = ban.get_node( center(e) );
            Vec<TNode*,6> node_center_Quad;
            for (unsigned i=0;i<6;++i) {
                node_center_Quad[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<1>())[i]) ) );
            }
            Vec<TNode*,12> node_center_Bar;
            for (unsigned i=0;i<12;++i) {
                node_center_Bar[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<2>())[i]) ) );
            }
            /// Hexa 0
            vn[ 0 ] = e.node( 0 );
            vn[ 1 ] = node_center_Bar[ 0 ];
            vn[ 2 ] = node_center_Quad[ 0 ];
            vn[ 3 ] = node_center_Bar[ 3 ];
            vn[ 4 ] = node_center_Bar[ 8 ];
            vn[ 5 ] = node_center_Quad[ 2 ];
            vn[ 6 ] = node_center_Hexa;
            vn[ 7 ] = node_center_Quad[ 4 ];
            permutation_if_jac_neg ( Hexa(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
            /// Hexa 1
            vn[ 0 ] = node_center_Bar[ 0 ];
            vn[ 1 ] = e.node( 1 );
            vn[ 2 ] = node_center_Bar[ 1 ];
            vn[ 3 ] = node_center_Quad[ 0 ];
            vn[ 4 ] = node_center_Quad[ 2 ];
            vn[ 5 ] = node_center_Bar[ 9 ];
            vn[ 6 ] = node_center_Quad[ 5 ];
            vn[ 7 ] = node_center_Hexa;
            permutation_if_jac_neg ( Hexa(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
            /// Hexa 2
            vn[ 0 ] = node_center_Bar[ 3 ];
            vn[ 1 ] = node_center_Quad[ 0 ];
            vn[ 2 ] = node_center_Bar[ 2 ];
            vn[ 3 ] = e.node( 3 );
            vn[ 4 ] = node_center_Quad[ 4 ];
            vn[ 5 ] = node_center_Hexa;
            vn[ 6 ] = node_center_Quad[ 3 ];
            vn[ 7 ] = node_center_Bar[ 11 ];
            permutation_if_jac_neg ( Hexa(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
            /// Hexa 3
            vn[ 0 ] = node_center_Quad[ 0 ];
            vn[ 1 ] = node_center_Bar[ 1 ];
            vn[ 2 ] = e.node( 2 );
            vn[ 3 ] = node_center_Bar[ 2 ];
            vn[ 4 ] = node_center_Hexa;
            vn[ 5 ] = node_center_Quad[ 5 ];
            vn[ 6 ] = node_center_Bar[ 10 ];
            vn[ 7 ] = node_center_Quad[ 3 ];
            permutation_if_jac_neg ( Hexa(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
            /// Hexa 4
            vn[ 0 ] = node_center_Bar[ 8 ];
            vn[ 1 ] = node_center_Quad[ 2 ];
            vn[ 2 ] = node_center_Hexa;
            vn[ 3 ] = node_center_Quad[ 4 ];
            vn[ 4 ] = e.node( 4 );
            vn[ 5 ] = node_center_Bar[ 4 ];
            vn[ 6 ] = node_center_Quad[ 1 ];
            vn[ 7 ] = node_center_Bar[ 7 ];
            permutation_if_jac_neg ( Hexa(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
            /// Hexa 5
            vn[ 0 ] = node_center_Quad[ 2 ];
            vn[ 1 ] = node_center_Bar[ 9 ];
            vn[ 2 ] = node_center_Quad[ 5 ];
            vn[ 3 ] = node_center_Hexa;
            vn[ 4 ] = node_center_Bar[ 4 ];
            vn[ 5 ] = e.node( 5 );
            vn[ 6 ] = node_center_Bar[ 5 ];
            vn[ 7 ] = node_center_Quad[ 1 ];
            permutation_if_jac_neg ( Hexa(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
            /// Hexa 6
            vn[ 0 ] = node_center_Quad[ 4 ];
            vn[ 1 ] = node_center_Hexa;
            vn[ 2 ] = node_center_Quad[ 3 ];
            vn[ 3 ] = node_center_Bar[ 11 ];
            vn[ 4 ] = node_center_Bar[ 7 ];
            vn[ 5 ] = node_center_Quad[ 1 ];
            vn[ 6 ] = node_center_Bar[ 6 ];
            vn[ 7 ] = e.node( 7 );
            permutation_if_jac_neg ( Hexa(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
            /// Hexa 7
            vn[ 0 ] = node_center_Hexa;
            vn[ 1 ] = node_center_Quad[ 5 ];
            vn[ 2 ] = node_center_Bar[ 10 ];
            vn[ 3 ] = node_center_Quad[ 3 ];
            vn[ 4 ] = node_center_Quad[ 1 ];
            vn[ 5 ] = node_center_Bar[ 5 ];
            vn[ 6 ] = e.node( 6 );
            vn[ 7 ] = node_center_Bar[ 6 ];
            permutation_if_jac_neg ( Hexa(), vn.ptr() );
            DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
            return true;
       }
       else
           return false;
    }
    
    /// pour tous les autres éléments autres que \a Triangle, \a Quad, \a Tetra, \a Hexa, on renvoie faux
    template<class TE> 
    bool operator() ( TE &e ) { 
        return false;
    }
    TM* ptr_m;
    const Vec<T>* ptr_list_elems;
};

template <class TM, class T>
void divide_element( TM &m, const Vec<T> &list_elems ) {
    Divide_element_locally<TM,T> de( m, list_elems );
    m.remove_elements_if( de );
}

template<class TM>
struct Divide_element {
    Divide_element( TM &m ) : ptr_m( &m ) {}
    
    template< class NameBehavior, class TNode_, class TData, unsigned num_in_elem_list_ >
    bool operator() ( Element< Triangle, NameBehavior, TNode_, TData, num_in_elem_list_ > &e ) {
        typedef typename Element< Triangle, NameBehavior, TNode_, TData, num_in_elem_list_ >::TNode TNode;
        Vec< TNode*, 3 > vn;
        ptr_m->update_elem_children();
        BestialNodeAdder<TM> ban; ban.m = ptr_m; ban.prec = 1e-6;
        Vec<TNode*,3> node_center_Bar;
        for (unsigned i=0;i<3;++i) {
            node_center_Bar[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<1>())[i]) ) );
        }
        /// Triangle 0
        vn[ 0 ] = e.node( 0 );
        vn[ 1 ] = node_center_Bar[ 0 ];
        vn[ 2 ] = node_center_Bar[ 2 ];
        permutation_if_jac_neg ( Triangle(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Triangle(), NameBehavior(), vn.ptr() ) );
        /// Triangle 1
        vn[ 0 ] = node_center_Bar[ 0 ];
        vn[ 1 ] = e.node( 1 );
        vn[ 2 ] = node_center_Bar[ 1 ];
        permutation_if_jac_neg ( Triangle(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Triangle(), NameBehavior(), vn.ptr() ) );
        /// Triangle 2
        vn[ 0 ] = node_center_Bar[ 2 ];
        vn[ 1 ] = node_center_Bar[ 1 ];
        vn[ 2 ] = e.node( 2 );
        permutation_if_jac_neg ( Triangle(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Triangle(), NameBehavior(), vn.ptr() ) );
        /// Triangle 3
        vn[ 0 ] = node_center_Bar[ 0 ];
        vn[ 1 ] = node_center_Bar[ 1 ];
        vn[ 2 ] = node_center_Bar[ 2 ];
        permutation_if_jac_neg ( Triangle(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Triangle(), NameBehavior(), vn.ptr() ) );
        return true;
    }
    
    template< class NameBehavior, class TNode_, class TData, unsigned num_in_elem_list_ >
    bool operator() ( Element< Quad, NameBehavior, TNode_, TData, num_in_elem_list_ > &e ) {
        typedef typename Element< Quad, NameBehavior, TNode_, TData, num_in_elem_list_ >::TNode TNode;
        Vec< TNode*, 4 > vn;
        ptr_m->update_elem_children();
        BestialNodeAdder<TM> ban; ban.m = ptr_m; ban.prec = 1e-6;
        TNode* node_center_Quad = ban.get_node( center(e) );
        Vec<TNode*,4> node_center_Bar;
        for (unsigned i=0;i<4;++i) {
            node_center_Bar[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<1>())[i]) ) );
        }
        /// Quad 0
        vn[ 0 ] = e.node( 0 );
        vn[ 1 ] = node_center_Bar[ 0 ];
        vn[ 2 ] = node_center_Quad;
        vn[ 3 ] = node_center_Bar[ 3 ];
        permutation_if_jac_neg ( Quad(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Quad(), NameBehavior(), vn.ptr() ) );
        /// Quad 1
        vn[ 0 ] = node_center_Bar[ 0 ];
        vn[ 1 ] = e.node( 1 );
        vn[ 2 ] = node_center_Bar[ 1 ];
        vn[ 3 ] = node_center_Quad;
        permutation_if_jac_neg ( Quad(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Quad(), NameBehavior(), vn.ptr() ) );
        /// Quad 2
        vn[ 0 ] = node_center_Quad;
        vn[ 1 ] = node_center_Bar[ 1 ];
        vn[ 2 ] = e.node( 2 );
        vn[ 3 ] = node_center_Bar[ 2 ];
        permutation_if_jac_neg ( Quad(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Quad(), NameBehavior(), vn.ptr() ) );
        /// Quad 3
        vn[ 0 ] = node_center_Bar[ 3 ];
        vn[ 1 ] = node_center_Quad;
        vn[ 2 ] = node_center_Bar[ 2 ];
        vn[ 3 ] = e.node( 3 );
        permutation_if_jac_neg ( Quad(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Quad(), NameBehavior(), vn.ptr() ) );
        return true;
    }
    
    template< class NameBehavior, class TNode_, class TData, unsigned num_in_elem_list_ >
    bool operator() ( Element< Tetra, NameBehavior, TNode_, TData, num_in_elem_list_ > &e ) {
        typedef typename Element< Tetra, NameBehavior, TNode_, TData, num_in_elem_list_ >::TNode TNode;
        Vec< TNode*, 4 > vn;
        ptr_m->update_elem_children();
        ptr_m->update_elem_children( Number<2>() );
        BestialNodeAdder<TM> ban; ban.m = ptr_m; ban.prec = 1e-6;
        Vec<TNode*,6> node_center_Bar;
        for (unsigned i=0;i<6;++i) {
            node_center_Bar[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<2>())[i]) ) );
        }
        /// Tetra 0
        vn[ 0 ] = node_center_Bar[ 2 ];
        vn[ 1 ] = node_center_Bar[ 4 ];
        vn[ 2 ] = node_center_Bar[ 5 ];
        vn[ 3 ] = e.node( 3 );
        permutation_if_jac_neg ( Tetra(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
        /// Tetra 1
        vn[ 0 ] = e.node( 0 );
        vn[ 1 ] = node_center_Bar[ 0 ];
        vn[ 2 ] = node_center_Bar[ 1 ];
        vn[ 3 ] = node_center_Bar[ 2 ];
        permutation_if_jac_neg ( Tetra(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
        /// Tetra 2
        vn[ 0 ] = node_center_Bar[ 1 ];
        vn[ 1 ] = node_center_Bar[ 3 ];
        vn[ 2 ] = e.node( 2 );
        vn[ 3 ] = node_center_Bar[ 5 ];
        permutation_if_jac_neg ( Tetra(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
        /// Tetra 3
        vn[ 0 ] = node_center_Bar[ 0 ];
        vn[ 1 ] = e.node( 1 );
        vn[ 2 ] = node_center_Bar[ 3 ];
        vn[ 3 ] = node_center_Bar[ 4 ];
        permutation_if_jac_neg ( Tetra(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
        /// Tetra 4
        vn[ 0 ] = node_center_Bar[ 2 ];
        vn[ 1 ] = node_center_Bar[ 0 ];
        vn[ 2 ] = node_center_Bar[ 1 ];
        vn[ 3 ] = node_center_Bar[ 4 ];
        permutation_if_jac_neg ( Tetra(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
        /// Tetra 5
        vn[ 0 ] = node_center_Bar[ 1 ];
        vn[ 1 ] = node_center_Bar[ 3 ];
        vn[ 2 ] = node_center_Bar[ 5 ];
        vn[ 3 ] = node_center_Bar[ 4 ];
        permutation_if_jac_neg ( Tetra(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
        /// Tetra 6
        vn[ 0 ] = node_center_Bar[ 2 ];
        vn[ 1 ] = node_center_Bar[ 1 ];
        vn[ 2 ] = node_center_Bar[ 5 ];
        vn[ 3 ] = node_center_Bar[ 4 ];
        permutation_if_jac_neg ( Tetra(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
        /// Tetra 7
        vn[ 0 ] = node_center_Bar[ 0 ];
        vn[ 1 ] = node_center_Bar[ 3 ];
        vn[ 2 ] = node_center_Bar[ 1 ];
        vn[ 3 ] = node_center_Bar[ 4 ];
        permutation_if_jac_neg ( Tetra(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Tetra(), NameBehavior(), vn.ptr() ) );
        return true;
    }

    template< class NameBehavior, class TNode_, class TData, unsigned num_in_elem_list_ >
    bool operator() ( Element< Hexa, NameBehavior, TNode_, TData, num_in_elem_list_ > &e ) {
        typedef typename Element< Hexa, NameBehavior, TNode_, TData, num_in_elem_list_ >::TNode TNode;
        Vec< TNode*, 8 > vn;
        ptr_m->update_elem_children();
        ptr_m->update_elem_children( Number<2>() );
        BestialNodeAdder<TM> ban; ban.m = ptr_m; ban.prec = 1e-6;
        TNode* node_center_Hexa = ban.get_node( center(e) );
        Vec<TNode*,6> node_center_Quad;
        for (unsigned i=0;i<6;++i) {
            node_center_Quad[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<1>())[i]) ) );
        }
        Vec<TNode*,12> node_center_Bar;
        for (unsigned i=0;i<12;++i) {
            node_center_Bar[i] = ban.get_node( center( *(ptr_m->get_children_of(e,Number<2>())[i]) ) );
        }
        /// Hexa 0
        vn[ 0 ] = e.node( 0 );
        vn[ 1 ] = node_center_Bar[ 0 ];
        vn[ 2 ] = node_center_Quad[ 0 ];
        vn[ 3 ] = node_center_Bar[ 3 ];
        vn[ 4 ] = node_center_Bar[ 8 ];
        vn[ 5 ] = node_center_Quad[ 2 ];
        vn[ 6 ] = node_center_Hexa;
        vn[ 7 ] = node_center_Quad[ 4 ];
        permutation_if_jac_neg ( Hexa(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
        /// Hexa 1
        vn[ 0 ] = node_center_Bar[ 0 ];
        vn[ 1 ] = e.node( 1 );
        vn[ 2 ] = node_center_Bar[ 1 ];
        vn[ 3 ] = node_center_Quad[ 0 ];
        vn[ 4 ] = node_center_Quad[ 2 ];
        vn[ 5 ] = node_center_Bar[ 9 ];
        vn[ 6 ] = node_center_Quad[ 5 ];
        vn[ 7 ] = node_center_Hexa;
        permutation_if_jac_neg ( Hexa(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
        /// Hexa 2
        vn[ 0 ] = node_center_Bar[ 3 ];
        vn[ 1 ] = node_center_Quad[ 0 ];
        vn[ 2 ] = node_center_Bar[ 2 ];
        vn[ 3 ] = e.node( 3 );
        vn[ 4 ] = node_center_Quad[ 4 ];
        vn[ 5 ] = node_center_Hexa;
        vn[ 6 ] = node_center_Quad[ 3 ];
        vn[ 7 ] = node_center_Bar[ 11 ];
        permutation_if_jac_neg ( Hexa(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
        /// Hexa 3
        vn[ 0 ] = node_center_Quad[ 0 ];
        vn[ 1 ] = node_center_Bar[ 1 ];
        vn[ 2 ] = e.node( 2 );
        vn[ 3 ] = node_center_Bar[ 2 ];
        vn[ 4 ] = node_center_Hexa;
        vn[ 5 ] = node_center_Quad[ 5 ];
        vn[ 6 ] = node_center_Bar[ 10 ];
        vn[ 7 ] = node_center_Quad[ 3 ];
        permutation_if_jac_neg ( Hexa(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
        /// Hexa 4
        vn[ 0 ] = node_center_Bar[ 8 ];
        vn[ 1 ] = node_center_Quad[ 2 ];
        vn[ 2 ] = node_center_Hexa;
        vn[ 3 ] = node_center_Quad[ 4 ];
        vn[ 4 ] = e.node( 4 );
        vn[ 5 ] = node_center_Bar[ 4 ];
        vn[ 6 ] = node_center_Quad[ 1 ];
        vn[ 7 ] = node_center_Bar[ 7 ];
        permutation_if_jac_neg ( Hexa(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
        /// Hexa 5
        vn[ 0 ] = node_center_Quad[ 2 ];
        vn[ 1 ] = node_center_Bar[ 9 ];
        vn[ 2 ] = node_center_Quad[ 5 ];
        vn[ 3 ] = node_center_Hexa;
        vn[ 4 ] = node_center_Bar[ 4 ];
        vn[ 5 ] = e.node( 5 );
        vn[ 6 ] = node_center_Bar[ 5 ];
        vn[ 7 ] = node_center_Quad[ 1 ];
        permutation_if_jac_neg ( Hexa(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
        /// Hexa 6
        vn[ 0 ] = node_center_Quad[ 4 ];
        vn[ 1 ] = node_center_Hexa;
        vn[ 2 ] = node_center_Quad[ 3 ];
        vn[ 3 ] = node_center_Bar[ 11 ];
        vn[ 4 ] = node_center_Bar[ 7 ];
        vn[ 5 ] = node_center_Quad[ 1 ];
        vn[ 6 ] = node_center_Bar[ 6 ];
        vn[ 7 ] = e.node( 7 );
        permutation_if_jac_neg ( Hexa(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
        /// Hexa 7
        vn[ 0 ] = node_center_Hexa;
        vn[ 1 ] = node_center_Quad[ 5 ];
        vn[ 2 ] = node_center_Bar[ 10 ];
        vn[ 3 ] = node_center_Quad[ 3 ];
        vn[ 4 ] = node_center_Quad[ 1 ];
        vn[ 5 ] = node_center_Bar[ 5 ];
        vn[ 6 ] = e.node( 6 );
        vn[ 7 ] = node_center_Bar[ 6 ];
        permutation_if_jac_neg ( Hexa(), vn.ptr() );
        DM::copy( e, *ptr_m->add_element( Hexa(), NameBehavior(), vn.ptr() ) );
        return true;
    }

    /// pour tous les autres éléments autres que \a Triangle, \a Quad, \a Tetra, \a Hexa, on renvoie faux
    template<class TE> 
    bool operator() ( TE &e ) { 
        return false;
    }
    TM* ptr_m;
};

template <class TM>
void divide_element( TM &m ) {
    Divide_element<TM> de( m );
    m.remove_elements_if( de );
}

/// Construction de la liste des elements list_elems_local_ref du maillage de reference contenus dans la liste des elements list_elems du maillage
///-----------------------------------------------------------------------------------------------------------------------------------------------
struct Construct_List_Elems_Local_Ref {
    template<class TE, class TM> void operator()( const TE &elem, TM &m, const Vec<unsigned> &list_elems, Vec<unsigned> &list_elems_local_ref ) const {
        for (unsigned n=0;n<list_elems.size();++n) {
            Vec<Vec<typename TE::T,TE::dim>, TE::nb_nodes > pos_nodes;
            for(unsigned i=0;i<(m.elem_list[ list_elems[ n ] ]->nb_nodes_virtual());++i)
                pos_nodes[i] = m.elem_list[ list_elems[ n ] ]->node_virtual(i)->pos;
            if ( is_inside_linear( typename TE::NE(), pos_nodes, center( elem ) ) ) {
                list_elems_local_ref.push_back( elem.number );
            }
        }
    }
};

/// Construction du noeud node_local_ref du maillage de reference correspondant au noeud node du maillage
///------------------------------------------------------------------------------------------------------
struct Construct_Node_Local_Ref {
    template<class TN, class TM> void operator()( const TN &node, TM &m, const unsigned &node_number, unsigned &node_number_local_ref ) const {
        if ( node.pos == m.node_list[ node_number ].pos ) {
            node_number_local_ref = node.number ;
        }
    }
};

/// Criteres de rafinement du maillage adjoint
///-------------------------------------------
/*!
    Objectif :
        Ce foncteur est conçu pour la fonction \a refinement () . Il permet de raffiner localement un maillage.
        
    Attributs :
        * <strong> c </strong> le centre de la zone que l'on veut raffiner. c n'est pas forcément un point dans le maillage.
        * <strong> l_min </strong> la longueur minimale des côtés des éléments du maillage.
        * <strong> k </strong> le coefficient d'augmentation de la longueur maximale des côtés des éléments en fonction de la distance au point c.
        * <strong> id </strong> le nom de l'attribut nodal qui compte le nombre de découpe. Remarque : il faut que le MeshCarac du maillage contienne une classe tag_refinement_DM.
        
    Description :
        on décide de couper le côté d'un élément ( i.e. une \a Bar ) si sa longueur est supérieure à d * k + l_min où d est la distance entre le milieu du côté et le centre c.
*/
template<class T, class Pvec>
struct Local_refinement_point_w_id {
    Local_refinement_point_w_id( T length_min, T _k, Pvec _c ) : l_min( length_min ), k( _k ), c( _c ), id( 1 ) {}

    template<class TE> 
    bool operator()( TE &e ) const {
        T l = length( e.node( 1 )->pos - e.node( 0 )->pos );
        T v = length( center( e ) - c ) * k + l_min;
        if ( l > v ) {
            e.node( 0 )->tag_refinement = id;
            e.node( 1 )->tag_refinement = id;
            return true;
        } else
            return false;
    }

    T l_min, k;
    Pvec c; /// centre
    unsigned id;
};

/*!
    Objectif :
        Ce foncteur est conçu pour la fonction \a refinement () . Il permet de raffiner localement un maillage.
        
    Attributs :
        * <strong> c </strong> le centre de la zone (cercle) autour duquel on veut raffiner. c n'est pas forcément un point dans le maillage.
        * <strong> R </strong> le rayon de la zone (cercle) autour duquel on veut raffiner.
        * <strong> l_min </strong> la longueur minimale des côtés des éléments du maillage.
        * <strong> k </strong> le coefficient d'augmentation de la longueur maximale des côtés des éléments en fonction de la distance au cercle.
        * <strong> id </strong> le nom de l'attribut nodal qui compte le nombre de découpe. Remarque : il faut que le MeshCarac du maillage contienne une classe tag_refinement_DM.
        
    Description :
        on décide de couper le côté d'un élément ( i.e. une \a Bar ) si sa longueur est supérieure à d * k + l_min où d est la distance entre le milieu du côté et le cercle de centre c et de rayon R.
*/
template<class T, class Pvec>
struct Local_refinement_circle_w_id {
    Local_refinement_circle_w_id( T length_min, T _k, Pvec _c, T _R ) : l_min( length_min ), k( _k ), c( _c ), R( _R ), id( 1 ) {}

    template<class TE> 
    bool operator()( TE &e ) const {
        T l = length( e.node( 1 )->pos - e.node( 0 )->pos );
        T v = fabs( R - length( center( e ) - c ) ) * k + l_min;
        if ( l > v ) {
            e.node( 0 )->tag_refinement = id;
            e.node( 1 )->tag_refinement = id;
            return true;
        } else
            return false;
    }

    T l_min, k;
    T R; /// rayon du cercle
    Pvec c; /// centre du cercle
    unsigned id;
};

/// Construction de la correspondance entre la liste des elements du maillage extrait (m_extracted) et la liste des elements du maillage original (m)
///--------------------------------------------------------------------------------------------------------------------------------------------------
struct Construct_Correspondance_Elem_Mesh_Extracted_To_Elem_Mesh {
    Vec<unsigned>* correspondance_elem_m_extracted_to_elem_m;
    template<class TE_extracted, class TE> void operator()( const TE_extracted &elem_extracted, unsigned i, const TE &elem, unsigned j ) const {
        Vec<Vec<typename TE::T,TE::dim>, TE::nb_nodes > pos_nodes;
        for(unsigned n=0;n<elem.nb_nodes;++n)
            pos_nodes[n] = elem.pos(n);
        if ( is_inside_linear( typename TE::NE(), pos_nodes, center( elem_extracted ) ) ) {
            (*correspondance_elem_m_extracted_to_elem_m)[ elem_extracted.number ] = elem.number;
        }
    }
};

#endif // Geometry_h
