//
// C++ Interface: EESPT
//
// Description: construction de champs admissibles, methode EESPT : construction standard des densites d'effort
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef EESPT_h
#define EESPT_h

using namespace LMT;
using namespace std;

/// Reperage pour chaque face k et chaque direction d de l'indice de debut de ligne dans les matrices A[ j ][ d ] : face_ind[ k ][ d ]
/// Calcul du nb de lignes de la matrice A[ j ][ d ] : nb_unk[ j ][ d ]
/// ----------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TVV, class TVVV, class TV> 
void calc_face_ind_EESPT( const TE &child_elem, TVVV &face_ind, TVV &nb_unk, const TV &connect_node_to_vertex_node ) {}

struct Calcul_Face_Ind_EESPT {
    template<class TE> void operator()( const TE &child_elem, Vec< Vec< Vec<unsigned> > > &face_ind, Vec< Vec<unsigned> > &nb_unk, const Vec<unsigned> &connect_node_to_vertex_node ) const {
        calc_face_ind_EESPT( child_elem, face_ind, nb_unk, connect_node_to_vertex_node );
    }
};

/// Reperage pour chaque element n et chaque direction d de l'indice de debut de colonne dans les matrices A[ j ][ d ] et de debut de ligne dans les vecteurs R[ j ][ d ] : vertex_nodal_ind[ n ][ d ]
/// Calcul du nb de colonnes de la matrice A[ j ][ d ] : nb_eq[ j ][ d ]
/// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TVV, class TVVV, class TV> 
void calc_vertex_nodal_ind_EESPT( const TE &elem, TVVV &vertex_nodal_ind, TVV &nb_eq, const TV &connect_node_to_vertex_node ) {}

struct Calcul_Vertex_Nodal_Ind_EESPT {
    template<class TE> void operator()( const TE &elem, Vec< Vec< Vec<unsigned> > > &vertex_nodal_ind, Vec< Vec<unsigned> > &nb_eq, const Vec<unsigned> &connect_node_to_vertex_node ) const {
        calc_vertex_nodal_ind_EESPT( elem, vertex_nodal_ind, nb_eq, connect_node_to_vertex_node );
    }
};

/// Construction des matrices A[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
/// --------------------------------------------------------------------------------------------------
template<class TE, class TM, class TV, class TVVV, class TVV, class TTMVV> 
void calc_vertex_nodal_matrix_A( const TE &elem, const TM &m, const TV &connect_node_to_vertex_node, const TVVV &face_ind, const TVVV &vertex_nodal_ind, const TVV &vertex_node_list_elem, const TVV &node_list_face, TTMVV &A ) {}

struct Calcul_Vertex_Nodal_Matrix_A {
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec< Vec<unsigned> > >* vertex_nodal_ind;
    const Vec< Vec<unsigned> >* vertex_node_list_elem;
    const Vec< Vec<unsigned> >* node_list_face;
    template<class TE, class TM, class T> void operator()( const TE &elem, const TM &m, const Vec<unsigned> &connect_node_to_vertex_node, Vec< Vec< Mat<T, Gen<>, SparseLine<> > > > &A ) const {
        calc_vertex_nodal_matrix_A( elem, m, connect_node_to_vertex_node, *face_ind, *vertex_nodal_ind, *vertex_node_list_elem, *node_list_face, A );
    }
};

/// Construction des vecteurs R[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
/// --------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TV, class TVVV, class TVV, class TTWW, class S, class B, class TTVVV> 
void calc_vertex_nodal_vector_R( const TE &elem, const TM &m, const TF &f, const TV &connect_node_to_vertex_node, const TVVV &vertex_nodal_ind, const TVV &node_list_face, const TVV &vertex_node_list_elem, const TV &elem_cpt_node, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVVV &R ) {}

struct Calcul_Vertex_Nodal_Vector_R {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* vertex_nodal_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec< Vec<unsigned> >* vertex_node_list_elem;
    const Vec<unsigned>* elem_cpt_node;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<T> > > &R ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_vertex_nodal_vector_R( elem, m, f, *connect_node_to_vertex_node, *vertex_nodal_ind, *node_list_face, *vertex_node_list_elem, *elem_cpt_node, f.vectors, ind, *pb, *want_local_enrichment, R );
    }
};

/// Suppression du noyau des matrices A[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
/// Stockage des inconnues non bloques : eq_indep[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
/// --------------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TV, class TVVV, class TVV, class BVV> 
void remove_kernel( const TE &elem , const TM &m, const TV &connect_node_to_vertex_node, const TVVV &vertex_nodal_ind, const TVV &node_list_vertex_node, const TVV &edge_node_list_vertex_node, BVV &node_flag, TVV &elem_flag, TVVV &eq_indep ) {}

struct Remove_Kernel {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* vertex_nodal_ind;
    const Vec< Vec<unsigned> >* node_list_vertex_node;
    const Vec< Vec<unsigned> >* edge_node_list_vertex_node;
    Vec< Vec<bool> >* node_flag;
    Vec< Vec<unsigned> >* elem_flag;
    template<class TE, class TM> void operator()( const TE &elem, const TM &m, Vec< Vec< Vec<unsigned> > > &eq_indep ) const {
        remove_kernel( elem, m, *connect_node_to_vertex_node, *vertex_nodal_ind, *node_list_vertex_node, *edge_node_list_vertex_node, *node_flag, *elem_flag, eq_indep );
    }
};

/// Calcul du degre p de l'analyse elements finis
/// ---------------------------------------------
template<class TE, class TV> 
void get_elem_degree( const TE &elem , TV &degree ) {}

struct Get_Elem_Degree {
    template<class TE> void operator()( const TE &elem, Vec<unsigned> &degree ) const {
        get_elem_degree( elem, degree );
    }
};

/// Construction des vecteurs lambda_F[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
/// ------------------------------------------------------------------------------------------------------------
    /// Cas p = 1
    /// ---------
        /// Construction des matrices B[ k ][ d ] pour chaque face k du maillage et chaque direction d
        /// ------------------------------------------------------------------------------------------
template<class TE, class TM, class TTMVV> 
void calc_skin_elem_matrix_B_p_1( const TE &child_elem, const TM &m, TTMVV &B ) {}

struct Calcul_Skin_Elem_Matrix_B_p_1 {
    template<class TE, class TM, class T> void operator()( const TE &child_elem, const TM &m, Vec< Vec< Mat< T, Gen<>, SparseUMFPACK > > > &B ) const {
        calc_skin_elem_matrix_B_p_1( child_elem, m, B );
    }
};

        /// Construction des vecteurs Q[ k ][ d ] pour chaque face k du maillage et chaque direction d
        /// ------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TV, class TVV, class TVVV, class TTWW, class S, class B, class TTVVV> 
void calc_skin_elem_vector_Q_p_1( const TE &elem , const TM &m, const TF &f, const TV &connect_node_to_vertex_node, const TVV &face_type, const TVVV &face_ind, const TVV &node_list_face, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVVV &Q ) {}

struct Calcul_Skin_Elem_Vector_Q_p_1 {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<T> > > &Q ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_skin_elem_vector_Q_p_1( elem, m, f, *connect_node_to_vertex_node, *face_type, *face_ind, *node_list_face, f.vectors, ind, *pb, *want_local_enrichment, Q );
    }
};

template<class TE, class TV, class TVVV, class TTVVV>
void calc_vertex_nodal_vector_lambda_F_p_1( const TE &child_elem, const TV &connect_node_to_vertex_node, const TVVV &face_ind, const TTVVV &lambda_F_face, TTVVV &lambda_F ) {}

struct Calcul_Vertex_Nodal_Vector_lambda_F_p_1 {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    template<class TE, class T> void operator()( const TE &child_elem, const Vec< Vec< Vec<T> > > &lambda_F_face, Vec< Vec< Vec<T> > > &lambda_F ) const {
        calc_vertex_nodal_vector_lambda_F_p_1( child_elem, *connect_node_to_vertex_node, *face_ind, lambda_F_face, lambda_F );
    }
};
    /// Cas p >= 2
    /// ----------
template<class TE, class TM, class TF, class TV, class TVV, class TVVV, class TTWW, class S, class B, class TTVVV>
void calc_vertex_nodal_vector_lambda_F_p_2( const TE &elem , const TM &m, const TF &f, const TV &connect_node_to_vertex_node, const TVV &face_type, const TVVV &face_ind, const TVV &node_list_face, const TVV &vertex_node_list_face, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVVV &lambda_F ) {}

struct Calcul_Vertex_Nodal_Vector_lambda_F_p_2 {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec< Vec<unsigned> >* vertex_node_list_face;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<T> > > &lambda_F ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_vertex_nodal_vector_lambda_F_p_2( elem, m, f, *connect_node_to_vertex_node, *face_type, *face_ind, *node_list_face, *vertex_node_list_face, f.vectors, ind, *pb, *want_local_enrichment, lambda_F );
    }
};

/// Construction des matrices de minimisation M[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
/// ------------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TV, class TVV, class TVVV, class T, class TT, class BVV, class TTMVV>
void calc_vertex_nodal_matrix_M( const TE &child_elem , const TM &m, const TF &f, const TV &connect_node_to_vertex_node, const TVV face_type, const TVVV &face_ind, const T &cost_function, const TT &penalty_val_N, const BVV &minimisation, TTMVV &M ) {}

template<class T>
struct Calcul_Vertex_Nodal_Matrix_M {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const unsigned* cost_function;
    const T* penalty_val_N;
    const Vec< Vec<bool> >* minimisation;
    template<class TE, class TM, class TF> void operator()( const TE &child_elem, const TM &m, const TF &f, Vec< Vec< Mat< T, Diag<> > > > &M ) const {
        calc_vertex_nodal_matrix_M( child_elem, m, f, *connect_node_to_vertex_node, *face_type, *face_ind, *cost_function, *penalty_val_N, *minimisation, M );
    }
};

/// Modification des vecteurs lambda_F_hat[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d (si amelioration)
/// ----------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TV, class B, class BV, class TVVV, class TTVVV> 
void reset_vertex_nodal_vector_lambda_F_hat( const TE &child_elem , const TM &m, const TV &connect_node_to_vertex_node, const B &enhancement, const BV &flag_face_enh, const TVVV &face_ind, TTVVV &lambda_F_hat ) {}

struct Reset_Vertex_Nodal_Vector_lambda_F_hat {
	const Vec<unsigned>* connect_node_to_vertex_node;
    const bool* enhancement;
    const Vec<bool>* flag_face_enh;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    template<class TE, class TM, class T> void operator()( const TE &child_elem, const TM &m, Vec< Vec< Vec<T> > > &lambda_F_hat ) const {
        reset_vertex_nodal_vector_lambda_F_hat( child_elem, m, *connect_node_to_vertex_node, *enhancement, *flag_face_enh, *face_ind, lambda_F_hat );
    }
};

/// Construction des vecteurs force_fluxes[ k ] pour chaque face k du maillage
/// --------------------------------------------------------------------------
template<class TE, class TV, class TVVV, class TTVVV>
void calc_skin_elem_force_fluxes( const TE &child_elem, const TV &connect_node_to_vertex_node, const TVVV &face_ind, const TTVVV &lambda_F_hat, TTVVV &force_fluxes ) {}

struct Calcul_Skin_Elem_Force_Fluxes {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    template<class TE, class T> void operator()( const TE &child_elem, const Vec< Vec< Vec<T> > > &lambda_F_hat, Vec< Vec< Vec<T> > > &force_fluxes ) const {
        calc_skin_elem_force_fluxes( child_elem, *connect_node_to_vertex_node, *face_ind, lambda_F_hat, force_fluxes );
    }
};

#endif // EESPT_h
