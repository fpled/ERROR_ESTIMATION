//
// C++ Interface: ECRE
//
// Description: construction de champs admissibles, methodes EESPT et EET : calcul d'un champ admissible et d'un estimateur theta d'erreur globale
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
///-----------------------------------------------------------------------
template<class TE, class TM, class TF, class TTWW, class TTMV>
void calc_elem_matrix_K_hat( TE &elem, const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, TTMV &K_hat ) {}

struct Calcul_Elem_Matrix_K_hat {
    template<class TE, class TM, class TF, class T> void operator()( TE &elem, const TM &m, const TF &f, Vec< Mat<T, Sym<> > > &K_hat ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_matrix_K_hat( elem, m, f, f.vectors, ind, K_hat );
    }
};

/// Construction des vecteurs F_hat[ n ] pour chaque element n du maillage
///-----------------------------------------------------------------------
template<class TE, class TM, class TF, class TVV, class TV, class B, class BV, class TTWW, class S, class TTVVV, class TTVV>
void calc_elem_vector_F_hat( TE &elem, const TM &m, const TF &f, const TVV &list_nodes_face, const TV &cpt_elems_node, const B &balancing, const BV &flag_elem_bal, const BV &flag_elem_enh, const TTWW &vectors, const Vec<unsigned> &indices, const S &pb, const B &want_local_enrichment, const TTVVV &vec_force_fluxes, TTVV &F_hat ) {}

template<class T>
struct Calcul_Elem_Vector_F_hat {
    const Vec< Vec<unsigned> >* list_nodes_face;
    const Vec<unsigned>* cpt_elems_node;
    const bool* balancing;
    const Vec<bool>* flag_elem_bal;
    const Vec<bool>* flag_elem_enh;
    const string* pb;
    const bool* want_local_enrichment;
    const Vec< Vec< Vec<T> > >* vec_force_fluxes;
    template<class TE, class TM, class TF> void operator()( TE &elem, const TM &m, const TF &f, Vec< Vec<T> > &F_hat ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        calc_elem_vector_F_hat( elem, m, f, *list_nodes_face, *cpt_elems_node, *balancing, *flag_elem_bal, *flag_elem_enh, f.vectors, ind, *pb, *want_local_enrichment, *vec_force_fluxes, F_hat );
    }
};

/// Construction de la matrice sigma_hat et Calcul d'un estimateur d'erreur globale au carre theta
///-----------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TTWW, class TTVV, class S, class TTV, class TT>
void calc_elem_error_estimate_EET_EESPT( TE &elem, const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, const TTVV &dep_hat, const S &method, TTV &theta_elem, TT &theta ) {}

template<class T, class TMAT>
struct Calc_Elem_Error_Estimate_EET_EESPT {
    const Vec< Vec<T> >* dep_hat;
    const string* method;
    Vec<T>* theta_elem;
    const Vec< Vec<T> >* dep_psi;
    const Vec< Vec<T> >* dep_lambda;
    const Vec<T>* dep_part;
    const Vec<T>* kappa;
    const TMAT* K_k_p;
    const TMAT* K_unk_p;
    const Vec<unsigned>* list_elems_PGD_unknown_parameter;
    const unsigned* nb_modes;
    template<class TE, class TM, class TF> void operator()( TE &elem, const TM &m, TF &f, T &theta ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        f.vectors[0] = (*dep_part);
        for (unsigned n=0;n<(*nb_modes);++n) {
            f.vectors[0] = - dot( (*dep_lambda)[ n ], (*K_k_p) * (*kappa) ) * (*dep_part);
            for (unsigned i=0;i<n+1;++i) {
                if ( find( *list_elems_PGD_unknown_parameter, _1 == elem.number ) )
                    f.vectors[0] += dot( (*dep_lambda)[ n ], (*K_unk_p) * (*dep_lambda)[ i ] ) * (*dep_psi)[ i ];
                else
                    f.vectors[0] += dot( (*dep_lambda)[ n ], (*K_k_p) * (*dep_lambda)[ i ] ) * (*dep_psi)[ i ];
            }
        }
        if ( *method == "EESPT" ) {
            elem.ecre_elem_EESPT = 0.0;
            elem.theta_elem_EESPT = 0.0;
            elem.theta_elem_rel_EESPT = 0.0;
        }
        if ( *method == "EET" ) {
            elem.ecre_elem_EET = 0.0;
            elem.theta_elem_EET = 0.0;
            elem.theta_elem_rel_EET = 0.0;
        }
        calc_elem_error_estimate_EET_EESPT( elem, m, f, f.vectors, ind, *dep_hat, *method, *theta_elem, theta );
    }
};

#endif // ECRE_h
