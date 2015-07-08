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
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TVVV, class TVV>
void calc_face_ind_EET( const TE &child_elem, TVVV &face_ind, TVV &nb_unk_b ) {}

struct Calcul_Face_Ind_EET {
    template<class TE> void operator()( const TE &child_elem, Vec< Vec< Vec<unsigned> > > &face_ind, Vec< Vec<unsigned> > &nb_unk ) const {
        calc_face_ind_EET( child_elem, face_ind, nb_unk );
    }
};

/// Reperage pour chaque element n et chaque direction d de l'indice de debut de ligne dans les vecteurs r[ i ][ d ] et dans les matrices B[ i ][ d ] : elem_ind[ n ][ d ]
/// Calcul du nb de lignes du vecteur r[ i ][ d ] et de la matrice B[ i ][ d ] : nb_eq[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TVVV, class TVV>
void calc_elem_ind_EET( const TE &elem, TVVV &elem_ind, TVV &nb_eq ) {}

struct Calcul_Elem_Ind_EET {
    template<class TE> void operator()( const TE &elem, Vec< Vec< Vec<unsigned> > > &elem_ind, Vec< Vec<unsigned> > &nb_eq ) const {
        calc_elem_ind_EET( elem, elem_ind, nb_eq );
    }
};

/// Construction des vecteurs r[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
///--------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVVV, class TVV, class TV, class TTWW, class S, class B, class TTVVV> 
void calc_nodal_vector_r( const TE &elem, const TM &m, const TF &f, const TVVV &elem_ind, const TVV &list_nodes_face, const TV &cpt_elems_node, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVVV &r ) {}

struct Calcul_Nodal_Vector_r {
    const Vec< Vec< Vec<unsigned> > >* elem_ind;
    const Vec< Vec<unsigned> >* list_nodes_face;
    const Vec<unsigned>* cpt_elems_node;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<T> > > &r ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_nodal_vector_r( elem, m, f, *elem_ind, *list_nodes_face, *cpt_elems_node, f.vectors, ind, *pb, *want_local_enrichment, r );
    }
};

template<class T, class TMAT>
struct Calcul_Nodal_Vector_r_PGD {
    const Vec< Vec< Vec<unsigned> > >* Q_i_ind;
    const Vec< Vec<unsigned> >* list_nodes_face;
    const Vec<unsigned>* cpt_elems_node;
    const string* pb;
    const bool* want_local_enrichment;
    const Vec< Vec<T> >* dep_psi;
    const Vec< Vec<T> >* dep_lambda;
    const Vec<T>* dep_part;
    const Vec<T>* kappa;
    const TMAT* K_k_p;
    const TMAT* K_unk_p;
    const Vec<unsigned>* list_elems_PGD_unknown_parameter;
    const unsigned* mode;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, TF &f, Vec< Vec< Vec<T> > > &Q_i ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        f.vectors[0] = - dot( (*dep_lambda)[ *mode ], (*K_k_p) * (*kappa) ) * (*dep_part);
        for (unsigned i=0;i<(*mode)+1;++i) {
            if ( find( *list_elems_PGD_unknown_parameter, _1 == elem.number ) )
                f.vectors[0] += dot( (*dep_lambda)[ *mode ], (*K_unk_p) * (*dep_lambda)[ i ] ) * (*dep_psi)[ i ];
            else
                f.vectors[0] += dot( (*dep_lambda)[ *mode ], (*K_k_p) * (*dep_lambda)[ i ] ) * (*dep_psi)[ i ];
        }
        calc_nodal_vector_r( elem, m, f, *Q_i_ind, *list_nodes_face, *cpt_elems_node, f.vectors, ind, *pb, *want_local_enrichment, Q_i );
    }
};

/// Construction des matrices B[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
///--------------------------------------------------------------------------------------------
template<class TE, class TM, class TVVV, class TVV, class TTMVV> 
void calc_nodal_matrix_B( const TE &elem, const TM &m, const TVVV &elem_ind, const TVVV &face_ind, const TVV &list_nodes_face, TTMVV &B ) {}

struct Calcul_Nodal_Matrix_B {
    const Vec< Vec< Vec<unsigned> > >* elem_ind;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* list_nodes_face;
    template<class TE, class TM, class T> void operator()( const TE &elem, const TM &m, Vec< Vec< Mat<T, Gen<>, SparseLine<> > > > &B ) const {
        calc_nodal_matrix_B( elem, m, *elem_ind, *face_ind, *list_nodes_face, B );
    }
};

/// Reperage pour chaque face k connectee au noeud i et chaque direction d de l'indice de debut de ligne dans les matrices C[ i ][ d ] et dans les vecteurs q[ i ][ d ] : nodal_ind[ i ][ d ][ k ]
/// Calcul du nb de lignes de la matrice C[ i ][ d ] et du vecteur q[ i ][ d ] : nb_eq_imp[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
///-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Calcul_Nodal_Ind {
    const Vec<unsigned>* cpt_faces_node;
    const Vec< Vec<unsigned> >* list_faces_node;
    const Vec< Vec<unsigned> >* type_node;
    const Vec< Vec<unsigned> >* type_face;
    template<class TN> void operator()( const TN &node, Vec< Vec<unsigned> > &nb_eq_imp, Vec< Vec< Vec<unsigned> > > &nodal_ind ) const {
        for (unsigned d=0;d<TN::dim;++d) {
            if ( (*type_node)[ node.number ][ d ] == 2 or (*type_node)[ node.number ][ d ] == 12 ) {
                for (unsigned k=0;k<(*cpt_faces_node)[ node.number ];++k) {
                    if ( (*type_face)[ (*list_faces_node)[ node.number ][ k ] ][ d ] == 2 ) {
                        nodal_ind[ node.number ][ d ][ k ] += nb_eq_imp[ node.number ][ d ];
                        nb_eq_imp[ node.number ][ d ]++;
                    }
                }
            }
        }
    }
};

/// Construction des matrices C[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
///--------------------------------------------------------------------------------------------
template<class TE, class TM, class TVV, class TVVV, class TTMVV> 
void calc_nodal_matrix_C( const TE &child_elem, const TM &m, const TVV &type_face, const TVVV &nodal_ind, const TVVV &face_ind, const TVV &list_faces_node, TTMVV &C ) {}

struct Calcul_Nodal_Matrix_C {
    const Vec< Vec<unsigned> >* type_face;
    const Vec< Vec< Vec<unsigned> > >* nodal_ind;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* list_faces_node;
    template<class TE, class TM, class T> void operator()( const TE &child_elem, const TM &m, Vec< Vec< Mat<T, Gen<>, SparseLine<> > > > &C ) const {
        calc_nodal_matrix_C( child_elem, m, *type_face, *nodal_ind, *face_ind, *list_faces_node, C );
    }
};

/// Construction des vecteurs q[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
///--------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TVVV, class TTWW, class TTVVV> 
void calc_nodal_vector_q( const TE &elem, const TM &m, const TF &f, const TVV &type_face, const TVVV &nodal_ind, const TVV &list_nodes_face, const TVV &list_faces_node, const TTWW &vectors, const Vec<unsigned> &indices, TTVVV &q ) {}

struct Calcul_Nodal_Vector_q {
    const Vec< Vec<unsigned> >* type_face;
    const Vec< Vec< Vec<unsigned> > >* nodal_ind;
    const Vec< Vec<unsigned> >* list_nodes_face;
    const Vec< Vec<unsigned> >* list_faces_node;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<T> > > &q ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_nodal_vector_q( elem, m, f, *type_face, *nodal_ind, *list_nodes_face, *list_faces_node, f.vectors, ind, q );
    }
};

template<class T, class TMAT>
struct Calcul_Nodal_Vector_q_PGD {
    const Vec< Vec<unsigned> >* type_face;
    const Vec< Vec< Vec<unsigned> > >* C_i_ind;
    const Vec< Vec<unsigned> >* list_nodes_face;
    const Vec< Vec<unsigned> >* list_faces_node;
    const Vec< Vec<T> >* dep_psi;
    const Vec< Vec<T> >* dep_lambda;
    const Vec<T>* dep_part;
    const Vec<T>* kappa;
    const TMAT* K_k_p;
    const TMAT* K_unk_p;
    const Vec<unsigned>* list_elems_PGD_unknown_parameter;
    const unsigned* mode;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, TF &f, Vec< Vec< Vec<T> > > &q_i ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        f.vectors[0] = - dot( (*dep_lambda)[ *mode ], (*K_k_p) * (*kappa) ) * (*dep_part);
        for (unsigned i=0;i<(*mode)+1;++i) {
            if ( find( *list_elems_PGD_unknown_parameter, _1 == elem.number ) )
                f.vectors[0] += dot( (*dep_lambda)[ *mode ], (*K_unk_p) * (*dep_lambda)[ i ] ) * (*dep_psi)[ i ];
            else
                f.vectors[0] += dot( (*dep_lambda)[ *mode ], (*K_k_p) * (*dep_lambda)[ i ] ) * (*dep_psi)[ i ];
        }
        calc_nodal_vector_q( elem, m, f, *type_face, *C_i_ind, *list_nodes_face, *list_faces_node, f.vectors, ind, q_i );
    }
};

/// Construction des matrices de minimisation M[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
///------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVVV, class T, class BVV, class TTDMVV> 
void calc_nodal_matrix_M( const TE &child_elem, const TM &m, const TF &f, const TVVV &face_ind, const T &cost_function, const BVV &minimisation, TTDMVV &M ) {}

struct Calcul_Nodal_Matrix_M {
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const unsigned *cost_function;
    const Vec< Vec<bool> >* minimisation;
    template<class TE, class TM, class TF, class T> void operator()( const TE &child_elem, const TM &m, const TF &f, Vec< Vec< Mat< T, Diag<> > > > &M ) const {
        calc_nodal_matrix_M( child_elem, m, f, *face_ind, *cost_function, *minimisation, M );
    }
};

/// Construction des vecteurs de minimisation b[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
///------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class BVV, class TVV, class TVVV, class TTWW, class S, class B, class TTVVV> 
void calc_nodal_vector_b( const TE &elem, const TM &m, const TF &f, const BVV &minimisation, const TVV &type_face, const TVVV &face_ind, const TVV &list_nodes_face, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, TTVVV &b ) {}

struct Calcul_Nodal_Vector_b {
    const Vec< Vec<bool> >* minimisation;
    const Vec< Vec<unsigned> >* type_face;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    const Vec< Vec<unsigned> >* list_nodes_face;
    const string* pb;
    const bool* want_local_enrichment;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<T> > > &b ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_nodal_vector_b( elem, m, f, *minimisation, *type_face, *face_ind, *list_nodes_face, f.vectors, ind, *pb, *want_local_enrichment, b );
    }
};

template<class T, class TMAT>
struct Calcul_Nodal_Vector_b_PGD {
    const Vec< Vec<bool> >* minimisation;
    const Vec< Vec<unsigned> >* type_face;
    const Vec< Vec< Vec<unsigned> > >* b_i_ind;
    const Vec< Vec<unsigned> >* list_nodes_face;
    const string* pb;
    const bool* want_local_enrichment;
    const Vec< Vec<T> >* dep_psi;
    const Vec< Vec<T> >* dep_lambda;
    const Vec<T>* dep_part;
    const Vec<T>* kappa;
    const TMAT* K_k_p;
    const TMAT* K_unk_p;
    const Vec<unsigned>* list_elems_PGD_unknown_parameter;
    const unsigned* mode;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, TF &f, Vec< Vec< Vec<T> > > &b_i ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        f.vectors[0] = - dot( (*dep_lambda)[ *mode ], (*K_k_p) * (*kappa) ) * (*dep_part);
        for (unsigned i=0;i<(*mode)+1;++i) {
            if ( find( *list_elems_PGD_unknown_parameter, _1 == elem.number ) )
                f.vectors[0] += dot( (*dep_lambda)[ *mode ], (*K_unk_p) * (*dep_lambda)[ i ] ) * (*dep_psi)[ i ];
            else
                f.vectors[0] += dot( (*dep_lambda)[ *mode ], (*K_k_p) * (*dep_lambda)[ i ] ) * (*dep_psi)[ i ];
        }
        calc_nodal_vector_b( elem, m, f, *minimisation, *type_face, *b_i_ind, *list_nodes_face, f.vectors, ind, *pb, *want_local_enrichment, b_i );
    }
};

/// Modification des vecteurs b_hat[ i ][ d ] pour chaque noeud i du maillage et chaque direction d (si amelioration)
///------------------------------------------------------------------------------------------------------------------

template<class TE, class TM, class B, class BV, class TVVV, class TTVVV> 
void reset_nodal_vector_b_hat( const TE &child_elem, const TM &m, const B &enhancement, const BV &flag_face_enh, const TVVV &face_ind, TTVVV &b_hat ) {}

struct Reset_Nodal_Vector_b_hat {
    const bool* enhancement;
    const Vec<bool>* flag_face_enh;
    const Vec< Vec< Vec<unsigned> > >* face_ind;
    template<class TE, class TM, class T> void operator()( const TE &child_elem, const TM &m, Vec< Vec< Vec<T> > > &b_hat ) const {
        reset_nodal_vector_b_hat( child_elem, m, *enhancement, *flag_face_enh, *face_ind, b_hat );
    }
};

/// Construction des vecteurs de projection b_face[ k ] pour chaque face k du maillage
///-----------------------------------------------------------------------------------
template<class TE, class TVVV, class TTVVV> 
void calc_skin_elem_vector_b_face( const TE &child_elem, const TVVV &face_ind, const TTVVV &b_hat, TTVVV &b_face ) {}

struct Calcul_Skin_Elem_Vector_b_face {
    template<class TE, class T> void operator()( const TE &child_elem, const Vec< Vec< Vec<unsigned> > > &face_ind, const Vec< Vec< Vec<T> > > &b_hat, Vec< Vec< Vec<T> > > &b_face ) const {
        calc_skin_elem_vector_b_face( child_elem, face_ind, b_hat, b_face );
    }
};

/// Construction des matrices K_face[ k ] pour chaque face k du maillage
///---------------------------------------------------------------------
template<class TE, class TTMVV> 
void calc_skin_elem_matrix_K_face( const TE &child_elem, TTMVV &K_face ) {}

struct Calcul_Skin_Elem_Matrix_K_face {
    template<class TE, class T> void operator()( const TE &child_elem, Vec< Vec< Mat<T, Gen<>, SparseUMFPACK > > > &K_face ) const {
        calc_skin_elem_matrix_K_face( child_elem, K_face );
    }
};

#endif // EET_h
