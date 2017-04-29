//
// C++ Interface: EET
//
// Description: construction de champs admissibles, methode EET : construction standard des densites d'effort
//
//
// Author: Pled Florent These 2009 <pled@lmt.ens-cachan.fr>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef EET_h
#define EET_h

using namespace LMT;
using namespace std;

/// Reperage pour chaque face k et chaque direction d de l'indice de debut de ligne dans les vecteurs b[ i ][ d ] et de debut de colonne dans les matrices B[ i ][ d ] : face_ind[ k ][ d ]
/// Calcul du nb de lignes du vecteur b[ i ][ d ] et de colonnes de la matrice B[ i ][ d ] : nb_unk[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TVVV, class TVV>
void calc_face_ind_EET( const TE &child_elem, TVVV &face_ind, TVV &nb_unk_b ) {}

struct Calcul_Face_Ind_EET {
    template<class TE> void operator()( const TE &child_elem, Vec< Vec< Vec<unsigned> > > &face_ind, Vec< Vec<unsigned> > &nb_unk ) const {
        calc_face_ind_EET( child_elem, face_ind, nb_unk );
    }
};

/// Reperage pour chaque element n et chaque direction d de l'indice de debut de ligne dans les vecteurs r[ i ][ d ] et dans les matrices B[ i ][ d ] : elem_ind[ n ][ d ]
/// Calcul du nb de lignes du vecteur r[ i ][ d ] et de la matrice B[ i ][ d ] : nb_eq[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TVVV, class TVV>
void calc_elem_ind_EET( const TE &elem, TVVV &elem_ind, TVV &nb_eq ) {}

struct Calcul_Elem_Ind_EET {
    template<class TE> void operator()( const TE &elem, Vec< Vec< Vec<unsigned> > > &elem_ind, Vec< Vec<unsigned> > &nb_eq ) const {
        calc_elem_ind_EET( elem, elem_ind, nb_eq );
    }
};

/// Construction des vecteurs r[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// -------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVVV, class TVV, class TV, class TTWW, class S, class B, class TTVVV>
void calc_nodal_vector_r( const TE &elem, const TM &m, const TF &f, const TVVV &elem_ind, const TVV &node_list_face, const TV &elem_cpt_node, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVVV &r ) {}

struct Calcul_Nodal_Vector_r {
    const Vec< Vec< Vec<unsigned> > >* elem_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec<unsigned>* elem_cpt_node;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class TVVV> void operator()( const TE &elem, const TM &m, const TF &f, TVVV &r ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_nodal_vector_r( elem, m, f, *elem_ind, *node_list_face, *elem_cpt_node, f.vectors, ind, *pb, *want_local_enrichment, r );
    }
};

template<class TVV>
struct Calcul_Nodal_Vector_r_PGD {
    const Vec< Vec< Vec<unsigned> > >* elem_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec<unsigned>* elem_cpt_node;
    const string* pb;
    const bool* want_local_enrichment;
    const TVV* dep;
    const Vec< Vec<unsigned> >* elem_group;
    template<class TE, class TM, class TF, class TVVV> void operator()( const TE &elem, const TM &m, TF &f, TVVV &r ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        for (unsigned g=0;g<(*elem_group).size();++g) {
            if ( find( (*elem_group)[g], _1 == elem.number ) )
                f.vectors[0] = (*dep)[g];
        }
        calc_nodal_vector_r( elem, m, f, *elem_ind, *node_list_face, *elem_cpt_node, f.vectors, ind, *pb, *want_local_enrichment, r );
    }
};

/// Construction des matrices B[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// -------------------------------------------------------------------------------------------
template<class TE, class TM, class TVVV, class TVV, class TMatVV>
void calc_nodal_matrix_B( const TE &elem, const TM &m, const TVVV &elem_ind, const TVVV &face_ind, const TVV &node_list_face, TMatVV &B ) {}

struct Calcul_Nodal_Matrix_B {
    const Vec< Vec< Vec<unsigned> > >* elem_ind;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    template<class TE, class TM, class TMatVV> void operator()( const TE &elem, const TM &m, TMatVV &B ) const {
        calc_nodal_matrix_B( elem, m, *elem_ind, *face_ind, *node_list_face, B );
    }
};

/// Reperage pour chaque face k connectee au noeud i et chaque direction d de l'indice de debut de ligne dans les matrices C[ i ][ d ] et dans les vecteurs q[ i ][ d ] : nodal_ind[ i ][ d ][ k ]
/// Calcul du nb de lignes de la matrice C[ i ][ d ] et du vecteur q[ i ][ d ] : nb_eq_imp[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Calcul_Nodal_Ind {
    const Vec<unsigned>* face_cpt_node;
    const Vec< Vec<unsigned> >* face_list_node;
    const Vec< Vec<unsigned> >* node_type;
    const Vec< Vec<unsigned> >* face_type;
    template<class TN> void operator()( const TN &node, Vec< Vec<unsigned> > &nb_eq_imp, Vec< Vec< Vec<unsigned> > > &nodal_ind ) const {
        for (unsigned d=0;d<TN::dim;++d) {
            if ( (*node_type)[ node.number ][ d ] == 2 or (*node_type)[ node.number ][ d ] == 12 ) {
                for (unsigned k=0;k<(*face_cpt_node)[ node.number ];++k) {
                    if ( (*face_type)[ (*face_list_node)[ node.number ][ k ] ][ d ] == 2 ) {
                        nodal_ind[ node.number ][ d ][ k ] += nb_eq_imp[ node.number ][ d ];
                        nb_eq_imp[ node.number ][ d ]++;
                    }
                }
            }
        }
    }
};

/// Construction des matrices C[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// -------------------------------------------------------------------------------------------
template<class TE, class TM, class TVV, class TVVV, class TMatVV>
void calc_nodal_matrix_C( const TE &child_elem, const TM &m, const TVV &face_type, const TVVV &nodal_ind, const TVVV &face_ind, const TVV &face_list_node, TMatVV &C ) {}

struct Calcul_Nodal_Matrix_C {
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* nodal_ind;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* face_list_node;
    template<class TE, class TM, class TMatVV> void operator()( const TE &child_elem, const TM &m, TMatVV &C ) const {
        calc_nodal_matrix_C( child_elem, m, *face_type, *nodal_ind, *face_ind, *face_list_node, C );
    }
};

/// Construction des vecteurs q[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// -------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TVVV, class TTWW, class TTVVV>
void calc_nodal_vector_q( const TE &elem, const TM &m, const TF &f, const TVV &face_type, const TVVV &nodal_ind, const TVV &node_list_face, const TVV &face_list_node, const TTWW &vectors, const Vec<unsigned> &indices, TTVVV &q ) {}

struct Calcul_Nodal_Vector_q {
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* nodal_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec< Vec<unsigned> >* face_list_node;
    template<class TE, class TM, class TF, class TVVV> void operator()( const TE &elem, const TM &m, const TF &f, TVVV &q ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_nodal_vector_q( elem, m, f, *face_type, *nodal_ind, *node_list_face, *face_list_node, f.vectors, ind, q );
    }
};

template<class TVV>
struct Calcul_Nodal_Vector_q_PGD {
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* nodal_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec< Vec<unsigned> >* face_list_node;
    const TVV* dep;
    const Vec< Vec<unsigned> >* elem_group;
    template<class TE, class TM, class TF, class TVVV> void operator()( const TE &elem, const TM &m, TF &f, TVVV &q ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        for (unsigned g=0;g<(*elem_group).size();++g) {
            if ( find( (*elem_group)[g], _1 == elem.number ) )
                f.vectors[0] = (*dep)[g];
        }
        calc_nodal_vector_q( elem, m, f, *face_type, *nodal_ind, *node_list_face, *face_list_node, f.vectors, ind, q );
    }
};

/// Construction des matrices de minimisation M[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// -----------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVVV, class T, class BVV, class TTDMVV>
void calc_nodal_matrix_M( const TE &child_elem, const TM &m, const TF &f, const TVVV &face_ind, const T &cost_function, const BVV &minimisation, TTDMVV &M ) {}

struct Calcul_Nodal_Matrix_M {
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const unsigned* cost_function;
    const Vec< Vec<bool> >* minimisation;
    template<class TE, class TM, class TF, class TMatVV> void operator()( const TE &child_elem, const TM &m, const TF &f, TMatVV &M ) const {
        calc_nodal_matrix_M( child_elem, m, f, *face_ind, *cost_function, *minimisation, M );
    }
};

/// Construction des vecteurs de minimisation b[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
/// -----------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class BVV, class TVV, class TVVV, class TTWW, class S, class B, class TTVVV>
void calc_nodal_vector_b( const TE &elem, const TM &m, const TF &f, const BVV &minimisation, const TVV &face_type, const TVVV &face_ind, const TVV &node_list_face, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVVV &b ) {}

struct Calcul_Nodal_Vector_b {
    const Vec< Vec<bool> >* minimisation;
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class TVVV> void operator()( const TE &elem, const TM &m, const TF &f, TVVV &b ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_nodal_vector_b( elem, m, f, *minimisation, *face_type, *face_ind, *node_list_face, f.vectors, ind, *pb, *want_local_enrichment, b );
    }
};

template<class TVV>
struct Calcul_Nodal_Vector_b_PGD {
    const Vec< Vec<bool> >* minimisation;
    const Vec< Vec<unsigned> >* face_type;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* node_list_face;
    const string* pb;
    const bool* want_local_enrichment;
    const TVV* dep;
    const Vec< Vec<unsigned> >* elem_group;
    template<class TE, class TM, class TF, class TVVV> void operator()( const TE &elem, const TM &m, TF &f, TVVV &b ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        for (unsigned g=0;g<(*elem_group).size();++g) {
            if ( find( (*elem_group)[g], _1 == elem.number ) )
                f.vectors[0] = (*dep)[g];
        }
        calc_nodal_vector_b( elem, m, f, *minimisation, *face_type, *face_ind, *node_list_face, f.vectors, ind, *pb, *want_local_enrichment, b );
    }
};

/// Modification des vecteurs b_hat[ i ][ d ] pour chaque noeud i du maillage et chaque direction d (si amelioration)
/// -----------------------------------------------------------------------------------------------------------------

template<class TE, class TM, class B, class BV, class TVVV, class TTVVV>
void reset_nodal_vector_b_hat( const TE &child_elem, const TM &m, const B &enhancement, const BV &flag_face_enh, const TVVV &face_ind, TTVVV &b_hat ) {}

struct Reset_Nodal_Vector_b_hat {
    const bool* enhancement;
    const Vec<bool>* flag_face_enh;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    template<class TE, class TM, class TVVV> void operator()( const TE &child_elem, const TM &m, TVVV &b_hat ) const {
        reset_nodal_vector_b_hat( child_elem, m, *enhancement, *flag_face_enh, *face_ind, b_hat );
    }
};

/// Construction des vecteurs de projection b_face[ k ] pour chaque face k du maillage
/// ----------------------------------------------------------------------------------
template<class TE, class TVVV, class TTVVV>
void calc_skin_elem_vector_b_face( const TE &child_elem, const TVVV &face_ind, const TTVVV &b_hat, TTVVV &b_face ) {}

struct Calcul_Skin_Elem_Vector_b_face {
    template<class TE, class TVVV> void operator()( const TE &child_elem, const Vec< Vec< Vec<unsigned> > > &face_ind, const TVVV &b_hat, TVVV &b_face ) const {
        calc_skin_elem_vector_b_face( child_elem, face_ind, b_hat, b_face );
    }
};

/// Construction des matrices K_face[ k ] pour chaque face k du maillage
/// --------------------------------------------------------------------
template<class TE, class TMatVV>
void calc_skin_elem_matrix_K_face( const TE &child_elem, TMatVV &K_face ) {}

struct Calcul_Skin_Elem_Matrix_K_face {
    template<class TE, class TMatVV> void operator()( const TE &child_elem, TMatVV &K_face ) const {
        calc_skin_elem_matrix_K_face( child_elem, K_face );
    }
};

#endif // EET_h
