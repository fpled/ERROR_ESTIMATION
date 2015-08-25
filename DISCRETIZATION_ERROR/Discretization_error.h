//
// C++ Interface: Discretization_error
//
// Description: calcul de la mesure de l'erreur de discretisation globale et locale
//
//
// Author: Pled Florent These 2009 <pled@lmt.ens-cachan.fr>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef Discretization_error_h
#define Discretization_error_h

using namespace LMT;
using namespace std;

/// Calcul de la norme au carre du champ de deplacement
/// ---------------------------------------------------
template<class TE, class TM, class TF, class TTWW, class TT>
void add_elem_norm_dep( const TE &elem, const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, TT &norm_dep ) {}

struct Add_Elem_Norm_Dep {
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, T &norm_dep ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
        add_elem_norm_dep( elem, m, f, f.vectors, ind, norm_dep );
    }
};

/// Calcul de l'erreur de discretisation au carre
/// ---------------------------------------------
template<class TE, class TE_REF, class TM, class TF, class TTVV, class TTV>
void calc_elem_discretization_error( const TE &elem, const TE_REF &elem_ref, const TM &m, const TM &m_ref, const TF &f, const TF &f_ref, const TTVV &vectors, const TTVV &ref_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &ref_indices, TTV &discretization_error_elem ) {}

template<class TM, class TF, class T>
struct Calcul_Elem_Discretization_Error {
    Vec<T>* discretization_error_elem;
    const TM* m;
    const TM* m_ref;
    const TF* f;
    const TF* f_ref;
    template<class TE, class TE_REF> void operator()( const TE &elem, unsigned i, const TE_REF &elem_ref, unsigned j ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = (*f).indices_for_element( elem );
        Vec<unsigned,TE_REF::nb_nodes+1+TF::nb_global_unknowns> ind_ref = (*f_ref).indices_for_element( elem_ref );
        calc_elem_discretization_error( elem, elem_ref, *m, *m_ref, *f, *f_ref, (*f).vectors, (*f_ref).vectors, ind, ind_ref, *discretization_error_elem );
    }
};

struct Set_Elem_Discretization_Error {
    template<class TE, class T> void operator()( TE &elem, const Vec<T> &discretization_error_elem ) const {
        elem.discretization_error_elem = discretization_error_elem[ elem.number ];
        elem.discretization_error_elem_rel = discretization_error_elem[ elem.number ] * pow( elem.measure_virtual(), -1 );
    }
};

/// Calcul de l'indice d'efficacite local
/// -------------------------------------
struct Calcul_Elem_Effectivity_Index {
    template<class TE, class T> void operator()( TE &elem, const string &method, Vec<T> &eff_index_elem ) const {
        if ( method == "EET" ) {
            elem.eff_index_elem_EET = sqrt( elem.theta_elem_EET / elem.discretization_error_elem );
            eff_index_elem[ elem.number ] = elem.eff_index_elem_EET;
        }
        else if ( method == "SPET" ) {
            elem.eff_index_elem_SPET = sqrt( elem.theta_elem_SPET / elem.discretization_error_elem );
            eff_index_elem[ elem.number ] = elem.eff_index_elem_SPET;
        }
        else if ( method == "EESPT" ) {
            elem.eff_index_elem_EESPT = sqrt( elem.theta_elem_EESPT / elem.discretization_error_elem );
            eff_index_elem[ elem.number ] = elem.eff_index_elem_EESPT;
        } 
    }
};

#endif // Discretization_error_h
