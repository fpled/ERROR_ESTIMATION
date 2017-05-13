//
// C++ Interface: Set_patch
//
// Description: construction de la table de connectivite de chaque pacth pour la methode basee sur la partition de l'unite (SPET)
//              construction des contraintes cinematiques sur chaque patch pour la methode basee sur la partition de l'unite (SPET)
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Set_patch_h
#define Set_patch_h

#include "SPET.h"
#include "Construct_connectivity_patch.h"
#include "Construct_kinematic_constraints_patch.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"

using namespace LMT;
using namespace std;

/// Construction de la table de connectivite et des contraintes cinematiques sur chaque patch
/// -----------------------------------------------------------------------------------------
template<class TM, class TF>
void set_patch( TM &m, const TF &f, unsigned &nb_vertex_nodes, Vec< Vec<unsigned> > &face_type, Vec<unsigned> &connect_node_to_vertex_node, Vec< Vec<unsigned> > &elem_list_vertex_node, Vec< Vec< Vec<unsigned> > > &patch_elem, Vec<unsigned> &nb_points_elem, Vec<unsigned> &nb_points_patch, Vec< Vec< Vec<unsigned> > > &constrained_points_list_patch,  Vec<unsigned> &nb_constraints_patch, const bool disp = false ) {
    
    Vec<unsigned> child_cpt;
    Vec< Vec<unsigned> > child_list;
    construct_child( m, child_cpt, child_list );
    
    Vec<bool> correspondance_node_to_vertex_node;
    nb_vertex_nodes = match_node_to_vertex_node( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node );
    
    Vec<unsigned> elem_cpt_vertex_node;
    construct_elems_connected_to_vertex_node( m, nb_vertex_nodes, correspondance_node_to_vertex_node, connect_node_to_vertex_node, elem_cpt_vertex_node, elem_list_vertex_node );
    
    correspondance_node_to_vertex_node.free();
    
    construct_face_type( m, f, face_type );
    
    /// Construction de la table de connectivite de chaque patch
    /// --------------------------------------------------------
    
    Vec< Vec<unsigned> > face_list_patch;
    Vec<unsigned> nb_points_face;
    Vec< Vec< Vec<unsigned> > > patch_face;
    
    construct_connectivity_patch( m, nb_vertex_nodes, face_type, face_list_patch, child_cpt, child_list, elem_cpt_vertex_node, elem_list_vertex_node, nb_points_face, nb_points_elem, nb_points_patch, patch_face, patch_elem );
    
    child_cpt.free();
    child_list.free();
    elem_cpt_vertex_node.free();
    
    /// Construction des contraintes cinematiques pour chaque noeud sommet j du maillage et chaque direction d
    /// ------------------------------------------------------------------------------------------------------
    
    construct_kinematic_constraints_patch( m, nb_vertex_nodes, face_type, face_list_patch, nb_points_face, patch_face, constrained_points_list_patch, nb_constraints_patch, disp );
    
    face_list_patch.free();
    nb_points_face.free();
    patch_face.free();
    nb_constraints_patch.free();
}

#endif // Set_patch_h
