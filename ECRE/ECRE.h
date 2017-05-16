//
// C++ Interface: ECRE
//
// Description: construction de champs admissibles, methodes EESPT et EET : calcul d'un champ admissible et d'un estimateur d'erreur globale
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef ECRE_h
#define ECRE_h

using namespace LMT;
using namespace std;

/// Construction des matrices K_hat[ n ] pour chaque element n du maillage
/// ----------------------------------------------------------------------
template<class TE, class TM, class TF, class TTWW, class TMatV>
void calc_elem_matrix_K_hat( TE &elem, const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, TMatV &K_hat ) {}

struct Calcul_Elem_Matrix_K_hat {
    template<class TE, class TM, class TF, class TMatV> void operator()( TE &elem, const TM &m, const TF &f, TMatV &K_hat ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_matrix_K_hat( elem, m, f, f.vectors, ind, K_hat );
    }
};

/// Construction des vecteurs F_hat[ n ] pour chaque element n du maillage
/// ----------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class B, class BV, class TTWW, class S, class TTVVV, class TTVV>
void calc_elem_vector_F_hat( TE &elem, const TM &m, const TF &f, const TVV &node_list_face, const TV &elem_cpt_node, const B &balancing, const BV &elem_flag_bal, const BV &elem_flag_enh, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, const TTVVV &force_fluxes, TTVV &F_hat ) {}

template<class TVVV>
struct Calcul_Elem_Vector_F_hat {
    const Vec< Vec<unsigned> >* node_list_face;
    const Vec<unsigned>* elem_cpt_node;
    const bool* balancing;
    const Vec<bool>* elem_flag_bal;
    const Vec<bool>* elem_flag_enh;
    const string* pb;
    const bool* want_local_enrichment;
    const TVVV* force_fluxes;
    template<class TE, class TM, class TF, class TVV> void operator()( TE &elem, const TM &m, const TF &f, TVV &F_hat ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_vector_F_hat( elem, m, f, *node_list_face, *elem_cpt_node, *balancing, *elem_flag_bal, *elem_flag_enh, f.vectors, ind, *pb, *want_local_enrichment, *force_fluxes, F_hat );
    }
};

/// Construction d'un champ admissible sigma_hat & Calcul d'un estimateur d'erreur globale au carre theta
/// -----------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TTWW, class TTVV, class TTV, class TT>
void calc_elem_error_estimate( TE &elem, const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, const TTVV &dep_hat, TTV &theta_elem, TT &theta ) {}

template<class TV, class TVV>
struct Calc_Elem_Error_Estimate {
    const TVV* dep_hat;
    TV* theta_elem;
    template<class TE, class TM, class TF, class T> void operator()( TE &elem, const TM &m, const TF &f, T &theta ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        elem.ecre_elem = 0.0;
        elem.error_estimate_elem = 0.0;
        calc_elem_error_estimate( elem, m, f, f.vectors, ind, *dep_hat, *theta_elem, theta );
    }
};

template<class TE, class TM, class TF, class TTWW, class TTV, class TT>
void calc_elem_error_estimate_PGD( TE &elem, const TM &m, const TF &f, const TTWW &vectors, const TTWW &FE_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &FE_indices, TTV &theta_elem_PGD, TT &theta_PGD ) {}

template<class TV, class TVV>
struct Calc_Elem_Error_Estimate_PGD {
    const TVV* dep;
    const Vec< Vec<unsigned> >* elem_group;
    TV* theta_elem_PGD;
    template<class TE, class TM, class TF, class T> void operator()( TE &elem, const TM &m, const TF &f, T &theta_PGD ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        Vec<Vec<typename TE::T>,TF::nb_vectors> vectors;
        for (unsigned g=0;g<(*elem_group).size();++g) {
            if ( find( (*elem_group)[g], _1 == elem.number ) )
                vectors[0] = (*dep)[g];
        }
        elem.ecre_elem_PGD = 0.0;
        elem.error_estimate_elem_PGD = 0.0;
        calc_elem_error_estimate_PGD( elem, m, f, f.vectors, vectors, ind, ind, *theta_elem_PGD, theta_PGD );
    }
};

struct Calc_Elem_Error_Estimate_Dis {
    template<class TE, class TV> void operator()( TE &elem, TV &theta_elem_dis ) const {
        elem.ecre_elem_dis = elem.ecre_elem - elem.ecre_elem_PGD;
        elem.error_estimate_elem_dis = elem.error_estimate_elem - elem.error_estimate_elem_PGD;
        theta_elem_dis[ elem.number ] = elem.ecre_elem_dis;
    }
};

template<class TE, class TM, class TF, class TTWW, class TTVV, class TTV, class TT>
void calc_elem_error_estimate_init( TE &elem, const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, const TTVV &dep_hat, TTV &theta_elem, TTV &theta_elem_init, TT &theta, TT &theta_init ) {}

template<class T, class TV, class TVV>
struct Calc_Elem_Error_Estimate_Init {
    const TVV* dep_hat;
    TV* theta_elem;
    TV* theta_elem_init;
    T* theta_init;
    template<class TE, class TM, class TF> void operator()( TE &elem, const TM &m, const TF &f, T &theta ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        elem.ecre_elem = 0.0;
        elem.ecre_elem_init = 0.0;
        elem.error_estimate_elem = 0.0;
        elem.error_estimate_elem_init = 0.0;
        calc_elem_error_estimate_init( elem, m, f, f.vectors, ind, *dep_hat, *theta_elem, *theta_elem_init, theta, *theta_init );
    }
};

#endif // ECRE_h
