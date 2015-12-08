//
// C++ Interface: Material_properties
//
// Description: creation des proprietes materiaux
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Material_properties_h
#define Material_properties_h

using namespace LMT;
using namespace std;

template<class TM, unsigned n, class T,class DM> void set_field_alternativeontype( TM &m, const Number<n> &, const T &val, const DM &dm ) {}

/// Number<0> : champ global
/// Number<1> : champ elementaire
/// -----------------------------
template<class TM, class T, class DM> void set_field_alternativeontype( TM &m, Number<0>, const T &val, const DM &dm ) {
    ExtractDM<DM> ed;
    ed( m ) = val;
}

template<class TM, class T, class DM> void set_field_alternativeontype( TM &m, Number<1>, const T &val, const DM &dm ) {
    for (unsigned n=0;n<m.elem_list.size();++n) {
        m.elem_list[n]->set_field( DM::name(), val );
    }
}

/// Number<0> : champ global
/// Number<1> : champ elementaire
/// -----------------------------
template<class TM, class TF, unsigned n1, unsigned n2> void calc_material_coefficients_alternativeontype( TM &m, TF &f, const Number<n1> &, const Number<n2> & );

template<class TM, class TF> void calc_material_coefficients_alternativeontype( TM &m, TF &f, Number<0>, Number<0> ) {
    static const unsigned dim = TM::dim;
    calc_material_coefficients( m, f, Number<dim>() );
}

struct Calcul_Elem_Material_Coefficients {
    template<class TE, class TM, class TF> void operator()( TE &elem, const TM &m, const TF &f ) const {
        calc_elem_material_coefficients( elem, m, f );
    }
};

template<class TM, class TF> void calc_material_coefficients_alternativeontype( TM &m, TF &f, Number<1>, Number<0> ) {
    apply( m.elem_list, Calcul_Elem_Material_Coefficients(), m, f );
}

template<class TM, class TF> void calc_material_coefficients_alternativeontype( TM &m, TF &f, Number<0>, Number<1> ) {
    apply( m.elem_list, Calcul_Elem_Material_Coefficients(), m, f );
}

template<class TM, class TF> void calc_material_coefficients_alternativeontype( TM &m, TF &f, Number<1>, Number<1> ) {
    apply( m.elem_list, Calcul_Elem_Material_Coefficients(), m, f );
}

/// Creation des proprietes materiaux
/// ---------------------------------
template<class TF, class TM>
void set_material_properties( TF &f, TM &m, const string &structure ) {

    typedef typename TM::TNode::T T;

    if ( m.node_list.size() ) {
        T young = 1;
        T poisson = 0.3;
        T density = 1;
        /// Plaque fissuree 2D
        /// ------------------
        if ( structure == "plate_crack" ) {
            young = 2.1e11;
            poisson = 0.3;
            density = 7820;
        }
        /// SAP 3D
        /// ------
        else if ( structure == "SAP" ) {
            young = 2.1e11;
            poisson = 0.29;
            density = 7820;
        }
        /// Carre 2D
        /// --------
        else if ( structure.find("square") != string::npos ) {
            T mu = 0.9;
            poisson = 0.3;
            young = 2*(1+poisson)*mu;
        }

        set_field_alternativeontype( m, Number< AreSameType< typename ExtractDM<young_DM>::ReturnType<TM>::T, void >::res >(), young, young_DM() );

        set_field_alternativeontype( m, Number< AreSameType< typename ExtractDM<poisson_DM>::ReturnType<TM>::T, void >::res >(), poisson, poisson_DM() );

        set_field_alternativeontype( m, Number< AreSameType< typename ExtractDM<density_DM>::ReturnType<TM>::T, void >::res >(), density, density_DM() );

        calc_material_coefficients_alternativeontype( m, f, Number< AreSameType< typename ExtractDM<young_DM>::ReturnType<TM>::T, void >::res >(), Number< AreSameType< typename ExtractDM<poisson_DM>::ReturnType<TM>::T, void >::res >() );
    }
}

#endif // Material_properties_h
