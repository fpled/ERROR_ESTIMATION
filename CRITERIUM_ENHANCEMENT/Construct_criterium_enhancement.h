//
// C++ Interface: Construction_criterium_enhnacement
//
// Description: contruction d'un critere d'amelioration : 
//              - geometrique : si geometric_ratio[ elem ] <= val_geometric_criterium, alors amelioration de la construction des densites d'effort sur elem
//              - sur l'estimateur d'erreur theta : si estimator_ratio[ elem ] >= val_estimator_criterium, alors amelioration de la construction des densites d'effort sur elem
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construction_criterium_enhnacement_h
#define Construction_criterium_enhnacement_h

#include "../GEOMETRY/Geometry.h"
#include "Criterium_enhancement.h"

using namespace LMT;
using namespace std;

/// Construction d'un critere geometrique pour le choix des elements dont les densites d'effort doivent etre ameliorees
///--------------------------------------------------------------------------------------------------------------------
template<class TM, class T>
void construct_geometric_criterium( TM &m, const string &geometric_criterium, Vec<T> &geometric_ratio, const bool &debug_criterium_enhancement ) {

    static const unsigned dim = TM::dim;

    cout << "--------------------------------------------------" << endl;
    cout << "Construction du critere d'amelioration geometrique" << endl;
    cout << "--------------------------------------------------" << endl << endl;

    /// Construction du vecteur de vecteurs circum_center et du vecteur circum_radius
    /// circum_center[ n ] : position du centre du cercle/sphere circonscrit a l'element n
    /// circum_radius[ n ] : rayon du cercle/sphere circonscrit a l'element n
    ///-----------------------------------------------------------------------------------

//     Vec< Vec<T> > circum_center;
//     circum_center.resize( m.elem_list.size() );
// 
//     for (unsigned n=0;n<m.elem_list.size();++n) {
//         circum_center[ n ].resize( dim );
//         circum_center[ n ].set( 0. );
//     }
// 
//     Vec<T> circum_radius;
//     circum_radius.resize( m.elem_list.size() );
//     circum_radius.set( 0. );
// 
//     apply( m.elem_list, Calcul_Circum_Center_Radius(), circum_center, circum_radius );
// 
//     if ( debug_criterium_enhancement ) {
//         cout << "Construction du vecteur de vecteurs circum_center et du vecteur circum_radius" << endl << endl;
//         for (unsigned n=0;n<m.elem_list.size();++n) {
//             if ( dim == 2 ) {
//                 cout << "position du centre du cercle circonscrit a l'element " << n << " : " << circum_center[ n ] << endl;
//                 cout << "rayon du cercle circonscrit a l'element " << n << " : " << circum_radius[ n ] << endl;
//             }
//             else if ( dim == 3 ) {
//                 cout << "position du centre de la sphere circonscrite a l'element " << n << " : " << circum_center[ n ] << endl;
//                 cout << "rayon de la sphere circonscrite a l'element " << n << " : " << circum_radius[ n ] << endl;
//             }
//             cout << endl << endl;
//         }
//     }

    /// Construction du vecteur de vecteurs in_center et du vecteur in_radius
    /// in_center[ n ] : position du centre du cercle/sphere inscrit dans l'element n
    /// in_radius[ n ] : rayon du cercle/sphere inscrit dans l'element n
    ///------------------------------------------------------------------------------

//     Vec< Vec<T> > in_center;
//     in_center.resize( m.elem_list.size() );
// 
//     for (unsigned n=0;n<m.elem_list.size();++n) {
//         in_center[ n ].resize( dim );
//         in_center[ n ].set( 0. );
//     }
// 
//     Vec<T> in_radius;
//     in_radius.resize( m.elem_list.size() );
//     in_radius.set( 0. );
// 
//     apply( m.elem_list, Calcul_In_Center_Radius(), in_center, in_radius );
// 
//     if ( debug_criterium_enhancement ) {
//         cout << "Construction du vecteur de vecteurs in_center et du vecteur in_radius" << endl << endl;
//         for (unsigned n=0;n<m.elem_list.size();++n) {
//             if ( dim == 2 ) {
//                 cout << "position du centre du cercle inscrit dans l'element " << n << " : " << in_center[ n ] << endl;
//                 cout << "rayon du cercle inscrit dans l'element " << n << " : " << in_radius[ n ] << endl;
//             }
//             else if ( dim == 3 ) {
//                 cout << "position du centre de la sphere inscrite dans l'element " << n << " : " << in_center[ n ] << endl;
//                 cout << "rayon de la sphere inscrite dans l'element " << n << " : " << in_radius[ n ] << endl;
//             }
//             cout << endl << endl;
//         }
//     }

    /// Construction du vecteur radius_ratio
    /// radius_ratio[ n ] : rapport du rayon du cercle/sphere circonscrit sur le rayon du cercle/sphere inscrit pour l'element n
    ///-------------------------------------------------------------------------------------------------------------------------

//     Vec<T> radius_ratio;
//     radius_ratio.resize( m.elem_list.size() );
//     radius_ratio.set( 0. );
// 
//     apply( m.elem_list, Calcul_Radius_Ratio(), radius_ratio );
// 
//     if ( debug_criterium_enhancement ) {
//         cout << "Construction du vecteur radius_ratio" << endl << endl;
//         for (unsigned n=0;n<m.elem_list.size();++n) {
//             if ( dim == 2 ) {
//                 cout << "rapport du rayon du cercle inscrit sur le rayon du cercle circonscrit pour l'element " << n << " : " << radius_ratio[ n ] << endl;
//             }
//             else if ( dim == 3 ) {
//                 cout << "rapport du rayon de la sphere inscrite sur le rayon de la sphere circonscrite pour l'element " << n << " : " << radius_ratio[ n ] << endl;
//             }
//             cout << endl << endl;
//         }
//     }

    /// Construction du vecteur edge_ratio
    /// edge_ratio[ n ] : rapport de la face la plus petite sur la face la plus grande pour l'element n
    ///------------------------------------------------------------------------------------------------

//     Vec<T> edge_ratio;
//     edge_ratio.resize( m.elem_list.size() );
//     edge_ratio.set( 0. );
// 
//     apply( m.elem_list, Calcul_Edge_Ratio(), m, edge_ratio );
// 
//     if ( debug_criterium_enhancement ) {
//         cout << "Construction du vecteur edge_ratio" << endl << endl;
//         for (unsigned n=0;n<m.elem_list.size();++n) {
//             cout << "rapport de la face la plus petite sur la face la plus grande pour l'element " << n << " : " << edge_ratio[ n ] << endl << endl;
//         }
//     }

    /// Construction du vecteur geometric_ratio
    /// geometric_ratio[ n ] : radius_ratio[ n ] ou edge_ratio[ n ] selon geometric_criterium pour l'element n
    ///-------------------------------------------------------------------------------------------------------

    cout << "Construction du vecteur geometric_ratio" << endl << endl;

    geometric_ratio.resize( m.elem_list.size() );
    geometric_ratio.set( 0. );

    apply( m.elem_list, Calcul_Geometric_Ratio(), m, geometric_criterium, geometric_ratio );
	
    if ( debug_criterium_enhancement ) {
        cout << "Construction du vecteur geometric_ratio" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "rapport geometrique pour l'element " << n << " : " << geometric_ratio[ n ] << endl << endl;
        }
    }
//     T geometric_ratio_max, geometric_ratio_min;
//     get_min_max( geometric_ratio, geometric_ratio_min, geometric_ratio_max );
//     cout << "plus petit rapport geometrique : " << geometric_ratio_min << endl << endl;
//     cout << "plus grand rapport geometrique : " << geometric_ratio_max << endl << endl;
	
}

/// Construction d'un critere sur l'estimateur d'erreur theta pour le choix des elements dont les densites d'effort doivent etre ameliorees
///----------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class T>
void construct_estimator_criterium( TM &m, Vec<T> &estimator_ratio, const Vec<T> &theta_2_elem, const bool &debug_criterium_enhancement ) {

    cout << "----------------------------------------------------------------" << endl;
    cout << "Construction du critere d'amelioration sur l'estimateur d'erreur" << endl;
    cout << "----------------------------------------------------------------" << endl << endl;

    /// Construction du vecteur estimator_ratio
    ///----------------------------------------

    cout << "Construction du vecteur estimator_ratio" << endl << endl;

    estimator_ratio.resize( m.elem_list.size() );
    estimator_ratio.set( 0. );

    apply( m.elem_list, Calcul_Estimator_Ratio(), m, theta_2_elem, estimator_ratio );
	
    if ( debug_criterium_enhancement ) {
        cout << "Construction du vecteur estimator_ratio" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "rapport sur l'estimateur d'erreur theta pour l'element " << n << " : " << estimator_ratio[ n ] << endl << endl;
        }
    }
//     T estimator_ratio_max, estimator_ratio_min;
//     get_min_max( estimator_ratio, estimator_ratio_min, estimator_ratio_max);
//     cout << "plus petit rapport sur l'estimateur d'erreur theta : " << estimator_ratio_min << endl << endl;
//     cout << "plus grand rapport sur l'estimateur d'erreur theta : " << estimator_ratio_max << endl << endl;
	
}

/// Application du critere d'amelioration geometrique et/ou sur l'estimateur d'erreur theta
///----------------------------------------------------------------------------------------
template<class TM, class T>
void apply_criterium_enhancement( TM &m, const string &method, const bool &enhancement_with_estimator_criterium, const bool &enhancement_with_geometric_criterium, const Vec<T> &estimator_ratio, const Vec<T> &geometric_ratio, const T &val_estimator_criterium, const T &val_geometric_criterium, Vec<bool> &flag_elem_enh, Vec<bool> &flag_face_enh, Vec<bool> &flag_elem_bal, Vec<unsigned> &list_elems_enh, Vec<unsigned> &list_faces_enh, Vec<unsigned> &list_elems_bal, const bool &debug_criterium_enhancement ) {
	
    Vec<bool> flag_node_enh;
    Vec<unsigned> list_nodes_enh;
	
    flag_elem_enh.resize( m.elem_list.size() );
    flag_elem_enh.set( 0 );
    flag_face_enh.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    flag_face_enh.set( 0 );
    flag_node_enh.resize( m.node_list.size() );
    flag_node_enh.set( 0 );
    flag_elem_bal.resize( m.elem_list.size() );
    flag_elem_bal.set( 0 );
	
    /// Construction des vecteurs flag_elem_enh, flag_face_enh, flag_node_enh
    /// Construction des vecteurs list_elems_enh, list_faces_enh, list_nodes_enh
    ///-------------------------------------------------------------------------

    cout << "Construction des vecteurs flag_elem_enh, flag_face_enh et flag_node_enh" << endl;
    cout << "Construction des vecteurs list_elems_enh, list_faces_enh et list_nodes_enh" << endl << endl;

    Apply_Criterium_Enhancement<T> apply_criterium_enhancement;
	apply_criterium_enhancement.method = &method;
	apply_criterium_enhancement.enhancement_with_estimator_criterium = &enhancement_with_estimator_criterium;
    apply_criterium_enhancement.enhancement_with_geometric_criterium = &enhancement_with_geometric_criterium;
    apply_criterium_enhancement.estimator_ratio = &estimator_ratio;
    apply_criterium_enhancement.geometric_ratio = &geometric_ratio;
	apply_criterium_enhancement.val_estimator_criterium = &val_estimator_criterium;
    apply_criterium_enhancement.val_geometric_criterium = &val_geometric_criterium;
    apply_criterium_enhancement.flag_elem_enh = &flag_elem_enh;
    apply_criterium_enhancement.flag_face_enh = &flag_face_enh;
    apply_criterium_enhancement.flag_node_enh = &flag_node_enh;
    apply_criterium_enhancement.list_nodes_enh = &list_nodes_enh;

    apply( m.elem_list, apply_criterium_enhancement, m, list_elems_enh, list_faces_enh );

    remove_doubles( list_faces_enh );
    remove_doubles( list_nodes_enh );

    /// Construction du vecteur flag_elem_bal
    /// Construction du vecteur list_elems_bal
    ///---------------------------------------

    cout << "Construction des vecteurs flag_elem_bal et list_elems_bal" << endl << endl;

    Construct_Balancing construct_balancing;
    construct_balancing.flag_face_enh = &flag_face_enh;

    apply( m.elem_list, construct_balancing, m, flag_elem_bal, list_elems_bal );

    if ( debug_criterium_enhancement ) {
        cout << "Construction du vecteur flag_elem_enh" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "participation de l'element " << n << " a l'amelioration : " << flag_elem_enh[ n ] << endl;
        }
        cout << endl << endl;
        cout << "Construction du vecteur flag_face_enh" << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            cout << "participation de la face " << k << " a l'amelioration : " << flag_face_enh[ k ] << endl;
        }
        cout << endl << endl;
        cout << "Construction du vecteur flag_node_enh" << endl << endl;
        for (unsigned i=0;i<m.node_list.size();++i) {
            cout << "participation du noeud " << i << " a l'amelioration : " << flag_node_enh[ i ] << endl;
        }
        cout << endl << endl;
        cout << "Construction du vecteur flag_elem_bal" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "participation de l'element " << n << " a l'equilibrage des densites d'effort ameliorees : " << flag_elem_bal[ n ] << endl;
        }
        cout << endl << endl;
        cout << "liste des elements participant a l'amelioration de la construction des densites d'effort : " << list_elems_enh << endl << endl;
        cout << "liste des faces participant a l'amelioration de la construction des densites d'effort : " << list_faces_enh << endl << endl;
        cout << "liste des noeuds participant a l'amelioration de la construction des densites d'effort : " << list_nodes_enh << endl << endl;
        cout << "liste des elements participant a l'equilibrage des densites d'effort ameliorees : " << list_elems_bal << endl << endl;
    }
    
    if ( enhancement_with_estimator_criterium )
        cout << "valeur du critere d'amelioration sur l'estimateur d'erreur : " << val_estimator_criterium << endl << endl;
    if ( enhancement_with_geometric_criterium )
	    cout << "valeur du critere d'amelioration geometrique : " << val_geometric_criterium << endl << endl;
    cout << "nb d'elements participant a l'amelioration de la construction des densites d'effort : " << list_elems_enh.size() << endl << endl;
    cout << "nb de faces participant a l'amelioration de la construction des densites d'effort : " << list_faces_enh.size() << endl << endl;
    cout << "nb de noeuds participant a l'amelioration de la construction des densites d'effort : " << list_nodes_enh.size() << endl << endl;
    cout << "nb d'elements participant a l'equilibrage des densites d'effort ameliorees : " << list_elems_bal.size() << endl << endl;
	
    if (list_elems_enh.size() == 0) {
        cerr << "Arret brutal, car il n'y a aucune densite d'effort a ameliorer pour le critere choisi..." << endl << endl;
        throw "Baleinou sous caillou...";
    }

}

#endif // Construction_criterium_enhancement_h
