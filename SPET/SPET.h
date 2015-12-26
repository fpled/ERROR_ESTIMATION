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
/// ------------------------------------------------------------------------------------------------------------------------
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
    template<class TE, class T> void operator()( const TE &elem, Vec< Vec< Vec<T> > > &pos_elem ) const {
        calc_pos_elem( elem, pos_elem );
    }
};

template<class TE, class TTVVV>
void calc_pos_face( const TE &elem, TTVVV &pos_face ) {}

struct Calc_Pos_Face {
    template<class TE, class T> void operator()( const TE &elem, Vec< Vec< Vec<T> > > &pos_face ) const {
        calc_pos_face( elem, pos_face );
    }
};

/// Construction de la table de connectivite de chaque patch : nb_unk_patch[ j ] pour chaque noeud sommet j du maillage 
///                                                            patch_elem[ j ][ n ][ i ] pour chaque noeud i de chaque element n du patch j
///                                                            patch_face[ j ][ n ][ i ] pour chaque noeud i de chaque face k du patch j
/// ---------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class TVVV>
void construct_patch( const TE &elem, const TM &m, const TF &f, const TVV &elem_list_vertex_node, const TVV &face_list_patch, const TV &connect_node_to_vertex_node, TV &nb_unk_patch, TVVV &patch_elem, TVVV &patch_face ) {}

struct Construct_Patch {
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec< Vec<unsigned> >* face_list_patch;
    const Vec<unsigned>* connect_node_to_vertex_node;
    Vec<unsigned>* nb_unk_patch;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<unsigned> > > &patch_face, Vec< Vec< Vec<unsigned> > > &patch_elem ) const {
        construct_patch( elem, m, f, *elem_list_vertex_node, *face_list_patch, *connect_node_to_vertex_node, *nb_unk_patch, patch_face, patch_elem );
    }
};

/// Construction des matrices K[ j ] pour chaque noeud sommet j du maillage
/// -----------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class TVVV, class TTMV> 
void calc_vertex_nodal_matrix_K( const TE &elem, const TM &m, const TF &f, const TVV &elem_list_vertex_node, const TV &connect_node_to_vertex_node, const TVVV &patch_elem, const TTMV &K ) {}

struct Calcul_Vertex_Nodal_Matrix_K {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* patch_elem;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Mat<T, Sym<> > > &K ) const {
        calc_vertex_nodal_matrix_K( elem, m, f, *elem_list_vertex_node, *connect_node_to_vertex_node, *patch_elem, K );
    }
};

/// Construction des vecteurs F[ j ] pour chaque noeud sommet j du maillage
/// -----------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class TVVV, class TTWW, class S, class B, class TTVV> 
void calc_vertex_nodal_vector_F( const TE &elem, const TM &m, const TF &f, const TVV &elem_list_vertex_node, const TVV &face_type, const TV &connect_node_to_vertex_node, const TVVV &patch_elem, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVV &F ) {}

struct Calcul_Vertex_Nodal_Vector_F {
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec< Vec<unsigned> >* face_type;
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* patch_elem;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec<T> > &F ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_vertex_nodal_vector_F( elem, m, f, *elem_list_vertex_node, *face_type, *connect_node_to_vertex_node, *patch_elem, f.vectors, ind, *pb, *want_local_enrichment, F );
    }
};

/// Construction des vecteurs E[ n ] pour chaque element n du maillage
/// ------------------------------------------------------------------
template<class TE, class TVV, class TV, class TVVV, class TTVV> 
void calc_elem_vector_E( const TE &elem, const TVV &elem_list_vertex_node, const TV &connect_node_to_vertex_node, const TVVV &patch_elem, const TTVV &U, TTVV &E ) {}

struct Calcul_Elem_Vector_E {
    const Vec< Vec<unsigned> >* elem_list_vertex_node;
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* patch_elem;
    template<class TE, class T> void operator()( const TE &elem, const Vec< Vec<T> > &U, Vec< Vec<T> > &E ) const {
        calc_elem_vector_E( elem, *elem_list_vertex_node, *connect_node_to_vertex_node, *patch_elem, U, E );
    }
};

/// Construction de la matrice sigma_hat et Calcul d'un estimateur d'erreur globale au carre theta
/// ----------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TTVV, class TTWW, class TTV, class TT> 
void calc_elem_error_estimate_SPET( TE &elem, const TM &m, const TF &f, const TTVV &E, const TTWW &vectors, const Vec<unsigned> &indices, TTV &theta_elem, TT &theta ) {}

template<class T>
struct Calcul_Elem_Error_Estimate_SPET {
    const Vec< Vec<T> >* E;
    Vec<T>* theta_elem;
    template<class TE, class TM, class TF> void operator()( TE &elem, const TM &m, const TF &f, T &theta ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_error_estimate_SPET( elem, m, f, *E, f.vectors, ind, *theta_elem, theta );
    }
};

template<class TE, class TM, class TF, class TTVV, class TTWW, class TTV, class TT>
void calc_elem_error_estimate_init_SPET( TE &elem, const TM &m, const TF &f, const TTVV &E, const TTWW &vectors, const Vec<unsigned> &indices, TTV &theta_elem, TTV &theta_elem_init, TT &theta, TT &theta_init ) {}

template<class T>
struct Calcul_Elem_Error_Estimate_Init_SPET {
    const Vec< Vec<T> >* E;
    Vec<T>* theta_elem;
    Vec<T>* theta_elem_init;
    T* theta_init;
    template<class TE, class TM, class TF> void operator()( TE &elem, const TM &m, const TF &f, T &theta ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_error_estimate_init_SPET( elem, m, f, *E, f.vectors, ind, *theta_elem, *theta_elem_init, theta, *theta_init );
    }
};

#endif // SPET_h
