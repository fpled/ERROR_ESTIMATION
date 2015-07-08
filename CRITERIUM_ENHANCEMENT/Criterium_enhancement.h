//
// C++ Interface: Criterium_enhancement
//
// Description: construction du critere d'amelioration geometrique et/ou sur l'estimateur
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef Criterium_enhancement_h
#define Criterium_enhancement_h

#include "../LMT/include/mesh/triangle.h"
#include "../LMT/include/mesh/triangle_6.h"
#include "../LMT/include/mesh/quad.h"
#include "../LMT/include/mesh/quad_8.h"
#include "../LMT/include/mesh/quad_9.h"
#include "../LMT/include/mesh/tetra.h"
#include "../LMT/include/mesh/tetra_10.h"
#include "../LMT/include/mesh/hexa.h"
#include "../LMT/include/mesh/hexa_20.h"

using namespace LMT;
using namespace std;

/// Construction du vecteur de vecteurs circum_center et du vecteur circum_radius
///------------------------------------------------------------------------------
template<class TE, class TTVV, class TTV>
void update_circum_center( const TE &e, TTVV &circum_center, TTV &circum_radius ) {}

struct Calcul_Circum_Center_Radius {
    template<class TE, class T> void operator()( const TE &elem, Vec< Vec<T> > &circum_center, Vec<T> &circum_radius ) const {
        update_circum_center( elem, circum_center[ elem.number ], circum_radius[ elem.number ] );
    }
};

struct Construct_Circum_Center_Radius_List_Elems {
    const Vec<unsigned>* list_elems;
    template<class TE, class T> void operator()( const TE &elem, Vec< Vec<T> > &circum_center, Vec<T> &circum_radius ) const {
        if ( find( *list_elems, _1 == elem.number ) ) {
            Vec<unsigned> ind_in_list_elems = find_with_index( *list_elems == elem.number );
            update_circum_center( elem, circum_center[ ind_in_list_elems[ 0 ] ], circum_radius[ ind_in_list_elems[ 0 ] ] );
        }
    }
};

/// Construction du vecteur de vecteurs in_center et du vecteur in_radius
///----------------------------------------------------------------------
template<class TE, class TTVV, class TTV>
void update_in_center( const TE &e, TTVV &in_center, TTV &in_radius ) {}

struct Calcul_In_Center_Radius {
    template<class TE, class T> void operator()( const TE &elem, Vec< Vec<T> > &in_center, Vec<T> &in_radius ) const {
        update_in_center( elem, in_center[ elem.number ], in_radius[ elem.number ] );
    }
};

/// Construction du vecteur radius_ratio
///-------------------------------------
template<class TE, class TTV>
void update_radius_ratio( const TE &e, TTV &radius_ratio ) {}

struct Calcul_Radius_Ratio {
    template<class TE, class T> void operator()( const TE &elem, Vec<T> &radius_ratio ) const {
        update_radius_ratio( elem, radius_ratio[ elem.number ] );
    }
};

/// Construction du vecteur edge_ratio
///-----------------------------------
template<class TE, class TM, class TTV>
void update_edge_ratio( const TE &e, TM &m, TTV &edge_ratio ) {}

struct Calcul_Edge_Ratio {
    template<class TE, class TM, class T> void operator()( const TE &elem, TM &m, Vec<T> &edge_ratio ) const {
        update_edge_ratio( elem, m, edge_ratio[ elem.number ] );
    }
};

/// Construction du vecteur geometric_ratio
///----------------------------------------
struct Calcul_Geometric_Ratio {
    template<class TE, class TM, class T> void operator()( const TE &elem, const TM &m, const string &geometric_criterium, Vec<T> &geometric_ratio ) const {
        if ( geometric_criterium == "radius_ratio" and ( elem.name() == "Triangle" or elem.name() == "Triangle_6" or elem.name() == "Tetra" or elem.name() == "Tetra_10" ) ) {
            update_radius_ratio( elem, geometric_ratio[ elem.number ] );
        }
        else if ( geometric_criterium == "edge_ratio" ) {
            update_edge_ratio( elem, m, geometric_ratio[ elem.number ] );
        }
        else {
            std::cout << "type de critere geometrique pour l'amelioration non implemente..." << std::endl;
        }
    }
};

/// Construction du vecteur estimator_ratio
///----------------------------------------
struct Calcul_Estimator_Ratio {
    template<class TE, class TM, class T> void operator()( const TE &elem, const TM &m, const Vec<T> &theta_2_elem, Vec<T> &estimator_ratio ) const {
		T theta_2_elem_min, theta_2_elem_max;
		get_min_max( theta_2_elem, theta_2_elem_min, theta_2_elem_max);
        estimator_ratio[ elem.number ] = theta_2_elem[ elem.number ] / theta_2_elem_max;
    }
};

/// Construction des vecteurs flag_elem_enh, flag_face_enh, flag_node_enh
/// Construction des vecteurs list_elems_enh, list_faces_enh, list_nodes_enh
///-------------------------------------------------------------------------
template<class TE, class TM, class S, class B, class TTV, class T, class BV, class TV>
void apply_criterium_enhancement( TE &e, const TM &m, const S &method, const B &enhancement_with_estimator_criterium, const B &enhancement_with_geometric_criterium, const TTV &estimator_ratio, const TTV &geometric_ratio, const T &val_estimator_criterium, const T &val_geometric_criterium, BV &flag_elem_enh, BV &flag_face_enh, BV &flag_node_enh, TV &list_elems_enh, TV &list_faces_enh, TV &list_nodes_enh ) {}

template<class T>
struct Apply_Criterium_Enhancement {
    const string* method;
    const bool* enhancement_with_estimator_criterium;
    const bool* enhancement_with_geometric_criterium;
    const Vec<T>* estimator_ratio;
    const Vec<T>* geometric_ratio;
    const T* val_estimator_criterium;
    const T* val_geometric_criterium;
    Vec<bool>* flag_elem_enh;
    Vec<bool>* flag_face_enh;
    Vec<bool>* flag_node_enh;
    Vec<unsigned>* list_nodes_enh;
    template<class TE, class TM> void operator()( TE &elem, const TM &m, Vec<unsigned> &list_elems_enh, Vec<unsigned> &list_faces_enh ) const {
        if ( ( *enhancement_with_estimator_criterium and *enhancement_with_geometric_criterium == 0 and (*estimator_ratio)[ elem.number ] >= *val_estimator_criterium ) or ( *enhancement_with_estimator_criterium == 0 and *enhancement_with_geometric_criterium and (*geometric_ratio)[ elem.number ] <= *val_geometric_criterium ) or ( *enhancement_with_estimator_criterium and *enhancement_with_geometric_criterium and (*geometric_ratio)[ elem.number ] <= *val_geometric_criterium and (*estimator_ratio)[ elem.number ] >= *val_estimator_criterium ) ) {
            if ( *method == "EET" ) {
                elem.enhancement_EET = 1;
            }
            if ( *method == "EESPT" ) {
                elem.enhancement_EESPT = 1;
            }
            (*flag_elem_enh)[ elem.number ] = 1;
            list_elems_enh.push_back( elem.number );
            for (unsigned k=0;k<NbChildrenElement<typename TE::NE,1>::res;++k) {
                (*flag_face_enh)[ m.get_children_of( elem, Number<1>() )[k]->number ] = 1;
                list_faces_enh.push_back( m.get_children_of( elem, Number<1>() )[k]->number );
            }
            for (unsigned i=0;i<TE::nb_nodes;++i) {
                (*flag_node_enh)[ elem.node( i )->number ] = 1;
                (*list_nodes_enh).push_back( elem.node( i )->number );
            }
        }
    }
};

/// Construction du vecteur flag_elem_bal
/// Construction du vecteur list_elems_bal
///---------------------------------------
struct Construct_Balancing {
    const Vec<bool>* flag_face_enh;
    template<class TE, class TM> void operator()( const TE &elem, const TM &m, Vec<bool> &flag_elem_bal, Vec<unsigned> &list_elems_bal ) const {
        for (unsigned k=0;k<NbChildrenElement<typename TE::NE,1>::res;++k) {
            if ( (*flag_face_enh)[ m.get_children_of( elem, Number<1>() )[k]->number ] ) {
                flag_elem_bal[ elem.number ] = 1;
            }
        }
        if ( flag_elem_bal[ elem.number ] ) {
            list_elems_bal.push_back( elem.number );
        }
    }
};

#endif // Criterium_enhancement_h
