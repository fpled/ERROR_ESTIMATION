//
// C++ Interface: SPET
//
// Description: construction de champs admissibles, methode SPET
//
//
// Author: Pled Florent,,,, These 2009 <pled@lmt.ens-cachan.fr>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SPET_h
#define SPET_h

using namespace LMT;
using namespace std;

/// Calcul du nb d'inconnues associes aux vecteurs u[ j ] et v[ j ] : nb_points_elem[ n ] pour chaque element n du maillage
///                                                                   nb_points_face[ n ] pour chaque face k du maillage
/// -----------------------------------------------------------------------------------------------------------------------
template<class TE, class TV>
void calc_nb_points_elem( const TE &elem, TV &nb_points_elem ) {}

struct Calc_Nb_Points_Elem {
    template<class TE> void operator()( const TE &elem, Vec<unsigned> &nb_points_elem ) const {
        calc_nb_points_elem( elem, nb_points_elem );
    }
};

template<class TE, class TV>
void calc_nb_points_face( const TE &child_elem, TV &nb_points_face ) {}

struct Calc_Nb_Points_Face {
    template<class TE> void operator()( const TE &child_elem, Vec<unsigned> &nb_points_face ) const {
        calc_nb_points_face( child_elem, nb_points_face );
    }
};

/// Construction de la position des inconnues associes aux vecteurs u[ j ] et v[ j ] : pos_elem[ n ][ p ][ n ] position du point p de l'element n du maillage dans la direction d
///                                                                                    pos_face[ k ][ p ][ n ] position du point p de la face k du maillage dans la direction d
/// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TTVVV>
void calc_pos_elem( const TE &elem, TTVVV &pos_elem ) {}

struct Calc_Pos_Elem {
    template<class TE, class TVVV> void operator()( const TE &elem, TVVV &pos_elem ) const {
        calc_pos_elem( elem, pos_elem );
    }
};

template<class TE, class TTVVV>
void calc_pos_face( const TE &elem, TTVVV &pos_face ) {}

struct Calc_Pos_Face {
    template<class TE, class TVVV> void operator()( const TE &elem, TVVV &pos_face ) const {
        calc_pos_face( elem, pos_face );
    }
};

/// Construction de la table de connectivite de chaque patch : nb_points_patch[ j ] pour chaque noeud sommet j du maillage
///                                                            patch_elem[ j ][ n ][ i ] pour chaque noeud i de chaque element n du patch j
///                                                            patch_face[ j ][ n ][ i ] pour chaque noeud i de chaque face k du patch j
/// ---------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class TVVV>
void construct_patch( const TE &elem, const TM &m, const TF &f, const TVV &elem_list_vertex_node, const TVV &face_list_patch, const TV &connect_node_to_vertex_node, TV &nb_points_patch, TVVV &patch_elem, TVVV &patch_face ) {}

struct Construct_Patch {
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec< Vec<unsigned> >* face_list_patch;
    const Vec<unsigned>* connect_node_to_vertex_node;
    Vec<unsigned>* nb_points_patch;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<unsigned> > > &patch_face, Vec< Vec< Vec<unsigned> > > &patch_elem ) const {
        construct_patch( elem, m, f, *elem_list_vertex_node, *face_list_patch, *connect_node_to_vertex_node, *nb_points_patch, patch_face, patch_elem );
    }
};

/// Construction des matrices K_hat[ j ] pour chaque noeud sommet j du maillage
/// ---------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class TVVV, class TMatV>
void calc_vertex_nodal_matrix_K_hat( const TE &elem, const TM &m, const TF &f, const TVV &elem_list_vertex_node, const TV &connect_node_to_vertex_node, const TVVV &patch_elem, const TMatV &K_hat ) {}

struct Calcul_Vertex_Nodal_Matrix_K_hat {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* patch_elem;
    template<class TE, class TM, class TF, class TMatV> void operator()( const TE &elem, const TM &m, const TF &f, TMatV &K_hat ) const {
        calc_vertex_nodal_matrix_K_hat( elem, m, f, *connect_node_to_vertex_node, *elem_list_vertex_node, *patch_elem, K_hat );
    }
};

/// Construction des vecteurs F_hat[ j ] pour chaque noeud sommet j du maillage
/// ---------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class TVVV, class TTWW, class S, class B, class TTVV>
void calc_vertex_nodal_vector_F_hat( const TE &elem, const TM &m, const TF &f, const TVV &elem_list_vertex_node, const TVV &face_type, const TV &connect_node_to_vertex_node, const TVVV &patch_elem, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVV &F_hat ) {}

struct Calcul_Vertex_Nodal_Vector_F_hat {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* patch_elem;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class TVV> void operator()( const TE &elem, const TM &m, const TF &f, TVV &F_hat ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_vertex_nodal_vector_F_hat( elem, m, f, *connect_node_to_vertex_node, *elem_list_vertex_node, *face_type, *patch_elem, f.vectors, ind, *pb, *want_local_enrichment, F_hat );
    }
};

template<class TVV>
struct Calcul_Vertex_Nodal_Vector_F_hat_PGD {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* patch_elem;
    const string* pb;
    const bool* want_local_enrichment;
    const TVV* dep;
    const Vec< Vec<unsigned> >* elem_group;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, const TF &f, TVV &F_hat ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        Vec<Vec<typename TE::T>,TF::nb_vectors> vectors;
        for (unsigned g=0;g<(*elem_group).size();++g) {
            if ( find( (*elem_group)[g], _1 == elem.number ) )
                vectors[0] = (*dep)[g];
        }
        calc_vertex_nodal_vector_F_hat( elem, m, f, *connect_node_to_vertex_node, *elem_list_vertex_node, *face_type, *patch_elem, vectors, ind, *pb, *want_local_enrichment, F_hat );
    }
};

/// Construction des vecteurs dep_hat[ n ] pour chaque element n du maillage
/// ------------------------------------------------------------------------
template<class TE, class TVV, class TV, class TVVV, class TTVV>
void calc_elem_vector_dep_hat( const TE &elem, const TV &connect_node_to_vertex_node, const TVV &elem_list_vertex_node, const TVVV &patch_elem, const TTVV &dep_hat_patch, TTVV &dep_hat ) {}

struct Calcul_Elem_Vector_Dep_hat {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* patch_elem;
    template<class TE, class TVV> void operator()( const TE &elem, const TVV &dep_hat_patch, TVV &dep_hat ) const {
        calc_elem_vector_dep_hat( elem, *connect_node_to_vertex_node, *elem_list_vertex_node, *patch_elem, dep_hat_patch, dep_hat );
    }
};

/// Construction de la matrice sigma_hat & Calcul d'un estimateur d'erreur globale au carre theta
/// ---------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TTVV, class TTWW, class TTV, class TT>
void calc_elem_error_estimate_SPET( TE &elem, const TM &m, const TF &f, const TTVV &dep_hat, const TTWW &vectors, const Vec<unsigned> &indices, TTV &theta_elem, TT &theta ) {}

template<class TV, class TVV>
struct Calcul_Elem_Error_Estimate_SPET {
    const TVV* dep_hat;
    TV* theta_elem;
    template<class TE, class TM, class TF, class T> void operator()( TE &elem, const TM &m, const TF &f, T &theta ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_error_estimate_SPET( elem, m, f, *dep_hat, f.vectors, ind, *theta_elem, theta );
    }
};

template<class TE, class TM, class TF, class TTVV, class TTWW, class TTV, class TT>
void calc_elem_error_estimate_init_SPET( TE &elem, const TM &m, const TF &f, const TTVV &dep_hat, const TTWW &vectors, const Vec<unsigned> &indices, TTV &theta_elem, TTV &theta_elem_init, TT &theta, TT &theta_init ) {}

template<class T, class TV, class TVV>
struct Calcul_Elem_Error_Estimate_Init_SPET {
    const TVV* dep_hat;
    TV* theta_elem;
    TV* theta_elem_init;
    T* theta_init;
    template<class TE, class TM, class TF> void operator()( TE &elem, const TM &m, const TF &f, T &theta ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_error_estimate_init_SPET( elem, m, f, *dep_hat, f.vectors, ind, *theta_elem, *theta_elem_init, theta, *theta_init );
    }
};

#endif // SPET_h
