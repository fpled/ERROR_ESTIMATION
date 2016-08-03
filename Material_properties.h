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

#include "LMT/include/util/Hdf.h"

using namespace LMT;
using namespace std;

template<class TM, unsigned n, class T,class DM>
void set_field_alternativeontype( TM &m, const Number<n> &, const T &val, const DM &dm ) {}

/// Number<0> : champ global
/// Number<1> : champ elementaire
/// -----------------------------
template<class TM, class T, class DM>
void set_field_alternativeontype( TM &m, Number<0>, const T &val, const DM &dm ) {
    ExtractDM<DM> ed;
    ed( m ) = val;
}

template<class TM, class T, class DM> void set_field_alternativeontype( TM &m, Number<1>, const T &val, const DM &dm ) {
    for (unsigned n=0;n<m.elem_list.size();++n)
        m.elem_list[n]->set_field( dm.name(), val ); // m.elem_list[n]->set_field( DM::name(), val );
}

//struct Set_Elem_Field {
//    template<class TE, class T, class DM>
//    void operator()( TE& elem, const T &val, const DM &dm ) const {
//        ExtractDM<DM> ed;
//        ed( elem ) = val;
//    }
//};

//template<class TM, class T, class DM>
//void set_field_alternativeontype( TM &m, Number<1>, const T &val, const DM &dm ) {
//    apply( m.elem_list, Set_Elem_Field(), val, dm );
//}

/// Number<0> : champ global
/// Number<1> : champ elementaire
/// -----------------------------
template<class TM, class TF, unsigned n1, unsigned n2>
void calc_material_coefficients_alternativeontype( TM &m, TF &f, const Number<n1> &, const Number<n2> & );

template<class TM, class TF>
void calc_material_coefficients_alternativeontype( TM &m, TF &f, Number<0>, Number<0> ) {
    static const unsigned dim = TM::dim;
    calc_material_coefficients( m, f, Number<dim>() );
}

struct Calcul_Elem_Material_Coefficients {
    template<class TE, class TM, class TF>
    void operator()( TE &elem, const TM &m, const TF &f ) const {
        calc_elem_material_coefficients( elem, m, f );
    }
};

template<class TM, class TF>
void calc_material_coefficients_alternativeontype( TM &m, TF &f, Number<1>, Number<0> ) {
    apply( m.elem_list, Calcul_Elem_Material_Coefficients(), m, f );
}

template<class TM, class TF>
void calc_material_coefficients_alternativeontype( TM &m, TF &f, Number<0>, Number<1> ) {
    apply( m.elem_list, Calcul_Elem_Material_Coefficients(), m, f );
}

template<class TM, class TF>
void calc_material_coefficients_alternativeontype( TM &m, TF &f, Number<1>, Number<1> ) {
    apply( m.elem_list, Calcul_Elem_Material_Coefficients(), m, f );
}

/// Creation des proprietes materiaux
/// ---------------------------------
template<class TF, class TM>
void set_material_properties( TF &f, TM &m, const string &structure ) {

    static const unsigned dim = TM::dim;
    typedef typename TM::TNode::T T;

    if ( m.node_list.size() ) {
        T young = 1;
        T poisson = 0.3;
        T density = 1;
        /// Dimension 2
        /// -----------
        if ( dim == 2 ) {
            /// Plaque fissuree 2D
            /// ------------------
            if ( structure == "plate_crack" ) {
                young = 2.1e11;
                poisson = 0.3;
                density = 7820;
            }
            /// Carre 2D
            /// --------
            else if ( structure.find("square") != string::npos ) {
                T mu_m = 1;
                if ( structure.find("init") == string::npos ) {
                    T mu = 0.9*mu_m;
                    poisson = 0.3;
                    young = 2*(1+poisson)*mu;
                }
                else {
                    for (unsigned n=0;n<m.elem_list.size();++n) {
                        T mu = mu_m;
                        poisson = 0.3;
                        if ( center( *m.elem_list[n] )[0] < 0.5 and center( *m.elem_list[n] )[1] < 0.5 ) { // x < 0.5 and y < 0.5
                            mu *= 100;
                            poisson = 0.2;
                        }
                        young = 2*(1+poisson)*mu;
                        m.elem_list[n]->set_field( "young", young );
                        m.elem_list[n]->set_field( "poisson", poisson );
                    }
                }
            }
        }
        /// Dimension 3
        /// -----------
        else if ( dim == 3 ) {
            /// Barre rectangulaire trouee 3D
            /// -----------------------------
            if ( structure == "beam_hole" ) {
                young = 100;
                poisson = 0.3;
                density = 80;
            }
            /// Quart de conduite 3D
            /// blocage des noeuds situes en x = 0 dans la direction x
            /// blocage des noeuds situes en y = 0.528 dans la direction y
            /// blocage des noeuds situes en z = 0.229 dans la direction z
            ///-----------------------------------------------------------
            else if ( structure == "pipe" ) {
                young = 209e3;
                poisson = 0.3;
                density = 80;
            }
            /// SAP 3D
            /// ------
            else if ( structure == "SAP" ) {
                young = 2.1e11;
                poisson = 0.29;
                density = 7820;
            }
            /// Eprouvette 3D
            /// -------------
            else if ( structure.find("test_specimen") != string::npos ) {
                young = 1.58e5;
                poisson = 0.28;
                T sigma_y = 289e6;
                T K = 1300e6;
                T n = 0.44;
            }
            /// Hashin's coated shpere 3D
            /// -------------------------
            else if ( structure.find("hashin") != string::npos ) {
                size_t offset = structure.rfind( "_" )+1;
                const string str = structure.substr( offset );
                istringstream buffer(str);
                unsigned N; buffer >> N;
                Hdf hdf("DATA/hashin-" + str + "x" + str + "x" + str + ".hdf5");

                hdf.read_tag( "/", "nu", poisson );
                if ( structure.find("init") == string::npos ) {
                    T k0;
                    hdf.read_tag( "/", "k0", k0 );
                    T mu = 3/5.*k0;
                    young = 2*(1+poisson)*mu;
                }
                else {
                    T k1, k2, k3;
                    hdf.read_tag( "/", "k1", k1 );
                    hdf.read_tag( "/", "k2", k2 );
                    hdf.read_tag( "/", "k3", k3 );
                    Tens3<double> f1, f2, f3;
                    hdf.read( "/f1", f1 );
                    hdf.read( "/f2", f2 );
                    f3.resize(N);
                    f3.set(1.);
                    f3 -= f1 + f2;
                    for (unsigned n=0;n<m.elem_list.size();++n) {
                        unsigned i = unsigned(center( *m.elem_list[n] )[0]*N-1/2);
                        unsigned j = unsigned(center( *m.elem_list[n] )[1]*N-1/2);
                        unsigned k = unsigned(center( *m.elem_list[n] )[2]*N-1/2);
                        T kappa = f1( k, j, i ) * k1 + f2( k, j, i ) * k2 + f3( k, j, i ) * k3;
                        T mu = 3/5.*kappa;
                        young = 2*(1+poisson)*mu;
                        m.elem_list[n]->set_field( "young", young );
                        m.elem_list[n]->set_field( "poisson", poisson );
                    }
                }
            }
        }

        if ( ( structure.find("square") == string::npos and structure.find("hashin") == string::npos ) or ( ( structure.find("square") != string::npos or structure.find("hashin") != string::npos ) and structure.find("init") == string::npos ) ) {
            set_field_alternativeontype( m, Number< AreSameType< typename ExtractDM<young_DM>::ReturnType<TM>::T, void >::res >(), young, young_DM() );
            set_field_alternativeontype( m, Number< AreSameType< typename ExtractDM<poisson_DM>::ReturnType<TM>::T, void >::res >(), poisson, poisson_DM() );
        }
        calc_material_coefficients_alternativeontype( m, f, Number< AreSameType< typename ExtractDM<young_DM>::ReturnType<TM>::T, void >::res >(), Number< AreSameType< typename ExtractDM<poisson_DM>::ReturnType<TM>::T, void >::res >() );

        set_field_alternativeontype( m, Number< AreSameType< typename ExtractDM<density_DM>::ReturnType<TM>::T, void >::res >(), density, density_DM() );
    }
}

#endif // Material_properties_h
