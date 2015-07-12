//
// C++ Interface: Enhancement_EET_EESPT
//
// Description: construction de champs admissibles, methodes EET et EESPT : amelioration des densites d'effort
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef Enhancement_EET_EESPT_h
#define Enhancement_EET_EESPT_h

using namespace LMT;
using namespace std;

/// Construction des vecteurs F_hat_enh[ n ] pour chaque element ameliore n du maillage
///------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TTWW, class BV, class TV, class TTVVV>
void calc_elem_vector_F_hat_enh( const TE &elem , const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, const BV &elem_flag_enh, const TV &elem_list_enh, TV &nb_unk_local_enh, TTVVV &F_hat_enh ) {}

struct Calcul_Elem_Vector_F_hat_enh {
    const Vec<bool>* elem_flag_enh;
    const Vec<unsigned>* elem_list_enh;
    Vec<unsigned>* nb_unk_local_enh;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec< Vec<T> > > &F_hat_enh ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_vector_F_hat_enh( elem, m, f, f.vectors, ind, *elem_flag_enh, *elem_list_enh, *nb_unk_local_enh, F_hat_enh );
    }
};

/// Construction des matrices A_local_enh[ n ] pour chaque element ameliore n du maillage
///--------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class BV, class TV, class TTVVV, class TTMV>
void calc_elem_matrix_A_enh( const TE &elem , const TM &m, const TF &f, const BV &elem_flag_enh, const TV &elem_list_enh, const TTVVV &dep_hat_enh, TTMV &A_local_enh ) {}

template<class T>
struct Calcul_Elem_Matrix_A_enh {
    const Vec<bool>* elem_flag_enh;
    const Vec<unsigned>* elem_list_enh;
    const Vec< Vec< Vec<T> > >* dep_hat_enh;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Mat<T, Gen<>, SparseLine<> > > &A_local_enh ) const {
        calc_elem_matrix_A_enh( elem, m, f, *elem_flag_enh, *elem_list_enh, *dep_hat_enh, A_local_enh );
    }
};

/// Construction des vecteurs d_local_enh[ n ] pour chaque element ameliore n du maillage
///--------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class BV, class TV, class TTVV, class TTVVV, class TTWW>
void calc_elem_vector_d_enh( const TE &elem , const TM &m, const TF &f, const BV &elem_flag_enh, const TV &elem_list_enh, const TTVV &dep_hat, const TTVVV &dep_hat_enh, const TTWW &vectors, const Vec<unsigned> &indices, TTVV &d_local_enh ) {}

template<class T>
struct Calcul_Elem_Vector_d_enh {
    const Vec<bool>* elem_flag_enh;
    const Vec<unsigned>* elem_list_enh;
    const Vec< Vec<T> >* dep_hat;
    const Vec< Vec< Vec<T> > >* dep_hat_enh;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec<T> > &d_local_enh ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_vector_d_enh( elem, m, f, *elem_flag_enh, *elem_list_enh, *dep_hat, *dep_hat_enh, f.vectors, ind, d_local_enh );
    }
};

/// Construction des matrices L_local_enh[ n ] pour chaque element ameliore n du maillage
///--------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class BV, class TV, class TTMV>
void calc_elem_matrix_L_enh( const TE &elem , const TM &m, const TF &f, const BV &face_flag_enh, const BV &elem_flag_bal, const TV &elem_list_bal, TV &nb_unk_local_bal, TV &nb_eq_f_vol_local_enh, TTMV &L_local_enh ) {}

struct Calcul_Elem_Matrix_L_enh {
    const Vec<bool>* face_flag_enh;
    const Vec<bool>* elem_flag_bal;
    const Vec<unsigned>* elem_list_bal;
    Vec<unsigned>* nb_unk_local_bal;
    Vec<unsigned>* nb_eq_f_vol_local_enh;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Mat<T, Gen<>, SparseLine<> > > &L_local_enh ) const {
        calc_elem_matrix_L_enh( elem, m, f, *face_flag_enh, *elem_flag_bal, *elem_list_bal, *nb_unk_local_bal, *nb_eq_f_vol_local_enh, L_local_enh );
    }
};

/// Construction des vecteurs b_local_enh[ n ] pour chaque element ameliore n du maillage
///--------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TTVVV, class BV, class TV, class TTVV>
void calc_elem_vector_b_enh( const TE &elem , const TM &m, const TF &f, const TVV &node_list_face, const TV &elem_cpt_node, const TTVVV &vec_force_fluxes, const BV &elem_flag_bal, const TV &elem_list_bal, TTVV &b_local_enh ) {}

template<class T>
struct Calcul_Elem_Vector_b_enh {
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec<unsigned>* elem_cpt_node;
    const Vec< Vec< Vec<T> > >* vec_force_fluxes;
    const Vec<bool>* elem_flag_bal;
    const Vec<unsigned>* elem_list_bal;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m, const TF &f, Vec< Vec<T> > &b_local_enh ) const {
        calc_elem_vector_b_enh( elem, m, f, *node_list_face, *elem_cpt_node, *vec_force_fluxes, *elem_flag_bal, *elem_list_bal, b_local_enh );
    }
};

/// Construction de la matrice globale A_enh
///-----------------------------------------
template<class TE, class TM, class TVV, class BV, class TV, class TTMV, class TTM>
void calc_glob_matrix_A_enh( const TE &elem , const TM &m, const TVV &node_list_face, const BV &elem_flag_enh, const TV &elem_list_enh, const TV &face_list_enh, const TTMV &A_local_enh, TTM &A_enh ) {}

struct Calcul_Global_Matrix_A_enh {
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec<bool>* elem_flag_enh;
    const Vec<unsigned>* elem_list_enh;
    const Vec<unsigned>* face_list_enh;
    template<class TE, class TM, class T> void operator()( const TE &elem, const TM &m, const Vec< Mat<T, Gen<>, SparseLine<> > > &A_local_enh, Mat<T, Gen<>, SparseLine<> > &A_enh ) const {
        calc_glob_matrix_A_enh( elem, m, *node_list_face, *elem_flag_enh, *elem_list_enh, *face_list_enh, A_local_enh, A_enh );
    }
};

/// Construction du vecteur global d_enh
///-------------------------------------
template<class TE, class TM, class TVV, class BV, class TV, class TTVV, class TTV>
void calc_glob_vector_d_enh( const TE &elem , const TM &m, const TVV &node_list_face, const BV &elem_flag_enh, const TV &elem_list_enh, const TV &face_list_enh, const TTVV &d_local_enh, TTV &d_enh ) {}

struct Calcul_Global_Vector_d_enh {
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec<bool>* elem_flag_enh;
    const Vec<unsigned>* elem_list_enh;
    const Vec<unsigned>* face_list_enh;
    template<class TE, class TM, class T> void operator()( const TE &elem, const TM &m, const Vec< Vec<T> > &d_local_enh, Vec<T> &d_enh ) const {
        calc_glob_vector_d_enh( elem, m, *node_list_face, *elem_flag_enh, *elem_list_enh, *face_list_enh, d_local_enh, d_enh );
    }
};

/// Construction de la matrice globale L_enh
///-----------------------------------------
template<class TE, class TM, class TVV, class BV, class TV, class TTMV, class TTM>
void calc_glob_matrix_L_enh( const TE &elem , const TM &m, const TVV &node_list_face, const BV &elem_flag_bal, const TV &elem_list_bal, const BV &face_flag_enh, const TV &face_list_enh, const TV &nb_eq_f_vol_local_enh, const TTMV &L_local_enh, TTM &L_enh ) {}

struct Calcul_Global_Matrix_L_enh {
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec<bool>* elem_flag_bal;
    const Vec<unsigned>* elem_list_bal;
    const Vec<bool>* face_flag_enh;
    const Vec<unsigned>* face_list_enh;
    const Vec<unsigned>* nb_eq_f_vol_local_enh;
    template<class TE, class TM, class T> void operator()( const TE &elem, const TM &m, const Vec< Mat<T, Gen<>, SparseLine<> > > &L_local_enh, Mat<T, Gen<>, SparseLine<> > &L_enh ) const {
        calc_glob_matrix_L_enh( elem, m, *node_list_face, *elem_flag_bal, *elem_list_bal, *face_flag_enh, *face_list_enh, *nb_eq_f_vol_local_enh, L_local_enh, L_enh );
    }
};

/// Construction du vecteur global b_enh
///-------------------------------------
template<class TE, class TM, class BV, class TV, class TTVV, class TTV>
void calc_glob_vector_b_enh( const TE &elem , const TM &m, const BV &elem_flag_bal, const TV &elem_list_bal, const TV &nb_eq_f_vol_local_enh, const TTVV &b_local_enh, TTV &b_enh ) {}

struct Calcul_Global_Vector_b_enh {
    const Vec<bool>* elem_flag_bal;
    const Vec<unsigned>* elem_list_bal;
    const Vec<unsigned>* nb_eq_f_vol_local_enh;
    template<class TE, class TM, class T> void operator()( const TE &elem, const TM &m, const Vec< Vec<T> > &b_local_enh, Vec<T> &b_enh ) const {
        calc_glob_vector_b_enh( elem, m, *elem_flag_bal, *elem_list_bal, *nb_eq_f_vol_local_enh, b_local_enh, b_enh );
    }
};

/// Construction de la matrice globale C_enh
///-----------------------------------------
template<class TE, class TVV, class BV, class TV, class T, class TTM>
void calc_glob_matrix_C_enh( const TE &child_elem , const TVV &face_type, const BV &face_flag_enh, const TV &face_list_enh, T &cpt_eq_f_surf_enh, TTM &C_enh ) {}

struct Calcul_Global_Matrix_C_enh {
    const Vec< Vec<unsigned> >* face_type;
    const Vec<bool>* face_flag_enh;
    const Vec<unsigned>* face_list_enh;
    template<class TE, class T> void operator()( const TE &child_elem, unsigned &cpt_eq_f_surf_enh, Mat<T, Gen<>, SparseLine<> > &C_enh ) const {
        calc_glob_matrix_C_enh( child_elem, *face_type, *face_flag_enh, *face_list_enh, cpt_eq_f_surf_enh, C_enh );
    }
};

/// Construction du vecteur global q_enh
///-------------------------------------
template<class TE, class TM, class TF, class TVV, class TTWW, class BV, class T, class TTV>
void calc_glob_vector_q_enh( const TE &elem , const TM &m, const TF &f, const TVV &face_type, const TTWW &vectors, const Vec<unsigned> &indices, const BV &elem_flag_enh, T cpt_eq_f_surf_enh, TTV &q_enh ) {}

struct Calcul_Global_Vector_q_enh {
    const Vec< Vec<unsigned> >* face_type;
    const Vec<bool>* elem_flag_enh;
    unsigned* cpt_eq_f_surf_enh;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, Vec<T> &q_enh ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_glob_vector_q_enh( elem, m, f, *face_type, f.vectors, ind, *elem_flag_enh, *cpt_eq_f_surf_enh, q_enh );
    }
};

/// Construction de la matrice globale P_enh
///-----------------------------------------
template<class TE, class BV, class TV, class T, class TTM>
void calc_glob_matrix_P_enh( const TE &child_elem, const BV &face_flag_enh, const TV &face_list_enh, T &cpt_eq_proj_f_surf_enh, TTM &P_enh ) {}

struct Calcul_Global_Matrix_P_enh {
    const Vec<bool>* face_flag_enh;
    const Vec<unsigned>* face_list_enh;
    template<class TE, class T> void operator()( const TE &child_elem, unsigned &cpt_eq_proj_f_surf_enh, Mat<T, Gen<>, SparseLine<> > &P_enh ) const {
        calc_glob_matrix_P_enh( child_elem, *face_flag_enh, *face_list_enh, cpt_eq_proj_f_surf_enh, P_enh );
    }
};

#endif // Enhancement_EET_EESPT_h
