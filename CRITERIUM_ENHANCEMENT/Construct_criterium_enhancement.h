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
void construct_geometric_criterium( TM &m, const string &geometric_criterium, Vec<T> &geometric_ratio, const bool debug_criterium_enhancement = false ) {

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
void construct_estimator_criterium( TM &m, Vec<T> &estimator_ratio, const Vec<T> &theta_2_elem, const bool debug_criterium_enhancement = false ) {

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
void apply_criterium_enhancement( TM &m, const string &method, const bool &enhancement_with_estimator_criterium, const bool &enhancement_with_geometric_criterium, const Vec<T> &estimator_ratio, const Vec<T> &geometric_ratio, const T &val_estimator_criterium, const T &val_geometric_criterium, Vec<bool> &elem_flag_enh, Vec<bool> &face_flag_enh, Vec<bool> &elem_flag_bal, Vec<unsigned> &elem_list_enh, Vec<unsigned> &face_list_enh, Vec<unsigned> &elem_list_bal, const bool debug_criterium_enhancement = false ) {
	
    Vec<bool> node_flag_enh;
    Vec<unsigned> node_list_enh;
	
    elem_flag_enh.resize( m.elem_list.size() );
    elem_flag_enh.set( 0 );
    face_flag_enh.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    face_flag_enh.set( 0 );
    node_flag_enh.resize( m.node_list.size() );
    node_flag_enh.set( 0 );
    elem_flag_bal.resize( m.elem_list.size() );
    elem_flag_bal.set( 0 );
	
    /// Construction des vecteurs elem_flag_enh, face_flag_enh, node_flag_enh
    /// Construction des vecteurs elem_list_enh, face_list_enh, node_list_enh
    ///-------------------------------------------------------------------------

    cout << "Construction des vecteurs elem_flag_enh, face_flag_enh et node_flag_enh" << endl;
    cout << "Construction des vecteurs elem_list_enh, face_list_enh et node_list_enh" << endl << endl;

    Apply_Criterium_Enhancement<T> apply_criterium_enhancement;
	apply_criterium_enhancement.method = &method;
	apply_criterium_enhancement.enhancement_with_estimator_criterium = &enhancement_with_estimator_criterium;
    apply_criterium_enhancement.enhancement_with_geometric_criterium = &enhancement_with_geometric_criterium;
    apply_criterium_enhancement.estimator_ratio = &estimator_ratio;
    apply_criterium_enhancement.geometric_ratio = &geometric_ratio;
	apply_criterium_enhancement.val_estimator_criterium = &val_estimator_criterium;
    apply_criterium_enhancement.val_geometric_criterium = &val_geometric_criterium;
    apply_criterium_enhancement.elem_flag_enh = &elem_flag_enh;
    apply_criterium_enhancement.face_flag_enh = &face_flag_enh;
    apply_criterium_enhancement.node_flag_enh = &node_flag_enh;
    apply_criterium_enhancement.node_list_enh = &node_list_enh;

    apply( m.elem_list, apply_criterium_enhancement, m, elem_list_enh, face_list_enh );

    remove_doubles( face_list_enh );
    remove_doubles( node_list_enh );

    /// Construction du vecteur elem_flag_bal
    /// Construction du vecteur elem_list_bal
    ///---------------------------------------

    cout << "Construction des vecteurs elem_flag_bal et elem_list_bal" << endl << endl;

    Construct_Balancing construct_balancing;
    construct_balancing.face_flag_enh = &face_flag_enh;

    apply( m.elem_list, construct_balancing, m, elem_flag_bal, elem_list_bal );

    if ( debug_criterium_enhancement ) {
        cout << "Construction du vecteur elem_flag_enh" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "participation de l'element " << n << " a l'amelioration : " << elem_flag_enh[ n ] << endl;
        }
        cout << endl << endl;
        cout << "Construction du vecteur face_flag_enh" << endl << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            cout << "participation de la face " << k << " a l'amelioration : " << face_flag_enh[ k ] << endl;
        }
        cout << endl << endl;
        cout << "Construction du vecteur node_flag_enh" << endl << endl;
        for (unsigned i=0;i<m.node_list.size();++i) {
            cout << "participation du noeud " << i << " a l'amelioration : " << node_flag_enh[ i ] << endl;
        }
        cout << endl << endl;
        cout << "Construction du vecteur elem_flag_bal" << endl << endl;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            cout << "participation de l'element " << n << " a l'equilibrage des densites d'effort ameliorees : " << elem_flag_bal[ n ] << endl;
        }
        cout << endl << endl;
        cout << "liste des elements participant a l'amelioration de la construction des densites d'effort : " << elem_list_enh << endl << endl;
        cout << "liste des faces participant a l'amelioration de la construction des densites d'effort : " << face_list_enh << endl << endl;
        cout << "liste des noeuds participant a l'amelioration de la construction des densites d'effort : " << node_list_enh << endl << endl;
        cout << "liste des elements participant a l'equilibrage des densites d'effort ameliorees : " << elem_list_bal << endl << endl;
    }
    
    if ( enhancement_with_estimator_criterium )
        cout << "valeur du critere d'amelioration sur l'estimateur d'erreur : " << val_estimator_criterium << endl << endl;
    if ( enhancement_with_geometric_criterium )
	    cout << "valeur du critere d'amelioration geometrique : " << val_geometric_criterium << endl << endl;
    cout << "nb d'elements participant a l'amelioration de la construction des densites d'effort : " << elem_list_enh.size() << endl << endl;
    cout << "nb de faces participant a l'amelioration de la construction des densites d'effort : " << face_list_enh.size() << endl << endl;
    cout << "nb de noeuds participant a l'amelioration de la construction des densites d'effort : " << node_list_enh.size() << endl << endl;
    cout << "nb d'elements participant a l'equilibrage des densites d'effort ameliorees : " << elem_list_bal.size() << endl << endl;
	
    if (elem_list_enh.size() == 0) {
        cerr << "Arret brutal, car il n'y a aucune densite d'effort a ameliorer pour le critere choisi..." << endl << endl;
        throw "Baleinou sous caillou...";
    }

}

#endif // Construction_criterium_enhancement_h
