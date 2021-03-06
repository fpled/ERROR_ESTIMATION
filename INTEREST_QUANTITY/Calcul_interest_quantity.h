//
// C++ Interface: Calcul_interest_quantity
//
// Description: calcul de la quantite d'interet locale
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Calcul_interest_quantity_h
#define Calcul_interest_quantity_h

#include "Interest_quantity.h"
#include "../LMT/include/containers/apply_ij.h"
#include "../Mesh.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"
#include "../Display.h"
#include "../LMT/include/mesh/gauss_point_for_Bar.h"
#include "../LMT/include/mesh/gauss_point_for_simplex_Triangle.h"
#include "../LMT/include/mesh/gauss_point_for_Quad.h"
#include "../LMT/include/mesh/gauss_point_for_simplex_Tetra.h"
#include "../LMT/include/mesh/gauss_point_for_Hexa.h"

using namespace LMT;
using namespace std;

/// Definition de l'extracteur
/// --------------------------
template<class TM, class TF, class T, class Pvec>
void define_extractor( TM &m, TM &m_crown, const TF &f, TF &f_crown, const string &interest_quantity, const string &direction_extractor, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list_interest_quantity, const unsigned &node_interest_quantity, const Pvec &pos_interest_quantity, const Pvec &pos_crack_tip, const T &angle_crack, const T &radius_Ri, const T &radius_Re, const bool want_local_enrichment = false, const bool display_vtu_crown = false, const string &filename = "paraview" ) {
    
    static const unsigned dim = TM::dim;
    
    if ( interest_quantity == "mean_sigma" or interest_quantity == "mean_epsilon" ) {
        T mes = 0.;
        for (unsigned n=0;n<m.elem_list.size();++n) {
            if ( find( elem_list_interest_quantity, _1 == n ) ) {
                mes += m.elem_list[n]->measure_virtual();
            }
        }
        for (unsigned n=0;n<m.elem_list.size();++n) {
            if ( find( elem_list_interest_quantity, _1 == n ) ) {
                Vec<T,unsigned(dim*(dim+1)/2)> extractor;
                extractor.set( 0. );
                if ( dim == 1 ) {
                    if ( direction_extractor == "xx" )
                        extractor[0] = 1. / mes;
                    else {
                        cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                        throw "Anguille sous coquille...";
                    }
                }
                else if ( dim == 2 ) {
                    if ( direction_extractor == "xx" )
                        extractor[0] = 1. / mes;
                    else if ( direction_extractor == "xy" or direction_extractor == "yx" )
                        extractor[1] = 1. / ( 2 * mes );
                    else if ( direction_extractor == "yy" )
                        extractor[2] = 1. / mes;
                    else {
                        cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                        throw "Anguille sous coquille...";
                    }
                }
                else if ( dim == 3 ) {
                    if ( direction_extractor == "xx" )
                        extractor[0] = 1. / mes;
                    else if ( direction_extractor == "xy" or direction_extractor == "yx" )
                        extractor[1] = 1. / ( 2 * mes );
                    else if ( direction_extractor == "yy" )
                        extractor[2] = 1. / mes;
                    else if ( direction_extractor == "xz" or direction_extractor == "zx" )
                        extractor[3] = 1. / ( 2 * mes );
                    else if ( direction_extractor == "yz" or direction_extractor == "zy" )
                        extractor[4] = 1. / ( 2 * mes );
                    else if ( direction_extractor == "zz" )
                        extractor[5] = 1. / mes;
                    else {
                        cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                        throw "Anguille sous coquille...";
                    }
                }
                else {
                    cerr << "Arret brutal, car la dimension " << dim << " pour la quantite d'interet " << interest_quantity << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                if ( interest_quantity == "mean_epsilon" )
                    m.elem_list[n]->set_field( "pre_sigma", extractor );
                else if ( interest_quantity == "mean_sigma" )
                    m.elem_list[n]->set_field( "pre_epsilon", extractor );
//                Mat<T,Sym<dim> > pre_eps = m.elem_list[n]->get_field( "pre_epsilon", StructForType<Mat<T,Sym<dim> > >() );
//                cout << pre_eps << endl;
//                Mat<T,Sym<dim> > pre_sig = m.elem_list[n]->get_field( "pre_sigma", StructForType<Mat<T,Sym<dim> > >() );
//                cout << pre_sig << endl;
            }
        }
    }
    else if ( interest_quantity == "pointwise_dep" ) {
        if ( dim == 1 ) {
            if ( direction_extractor == "x" )
                m.node_list[node_interest_quantity].pre_f_nodal[0] = 1.;
            else {
                cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                throw "Anguille sous coquille...";
            }
        }
        else if ( dim == 2 ) {
            if ( direction_extractor == "x" )
                m.node_list[node_interest_quantity].pre_f_nodal[0] = 1.;
            else if ( direction_extractor == "y" )
                m.node_list[node_interest_quantity].pre_f_nodal[1] = 1.;
            else {
                cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                throw "Anguille sous coquille...";
            }
        }
        else if ( dim == 3 ) {
            if ( direction_extractor == "x" )
                m.node_list[node_interest_quantity].pre_f_nodal[0] = 1.;
            else if ( direction_extractor == "y" )
                m.node_list[node_interest_quantity].pre_f_nodal[1] = 1.;
            else if ( direction_extractor == "z" )
                m.node_list[node_interest_quantity].pre_f_nodal[2] = 1.;
            else {
                cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                throw "Anguille sous coquille...";
            }
        }
    }
    else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
        if ( dim == 2 ) {
            f_crown.allocate_matrices();
            f_crown.shift();
            f_crown.assemble();
            f_crown.update_variables();
            f_crown.call_after_solve();
            
            Define_Extractor_SIF<T, Pvec> define_extractor_SIF;
            define_extractor_SIF.interest_quantity = &interest_quantity;
            define_extractor_SIF.direction_extractor = &direction_extractor;
            define_extractor_SIF.pos_crack_tip = &pos_crack_tip;
            define_extractor_SIF.angle_crack = &angle_crack;
            define_extractor_SIF.radius_Ri = &radius_Ri;
            define_extractor_SIF.radius_Re = &radius_Re;
            
            apply( m_crown.elem_list, define_extractor_SIF, m_crown, f_crown );
            
            const string filename_crown = filename + "_crown";
            if ( display_vtu_crown )
                display( m_crown, filename_crown );
            else
                save( m_crown, filename_crown );
        }
        else {
            cerr << "Arret brutal, car la dimension " << dim << " pour la quantite d'interet " << interest_quantity << " n'est pas implementee..." << endl << endl;
            throw "Anguille sous coquille...";
        }
    }
    else {
        cerr << "Arret brutal, car la quantite d'interet " << interest_quantity << " n'est pas implementee..." << endl << endl;
        throw "Anguille sous coquille...";
    }
    
    
    if ( want_local_enrichment and interest_quantity.find("pointwise") != string::npos ) {
        if ( pointwise_interest_quantity == "node" )
            m.pos_handbook = m.node_list[ node_interest_quantity ].pos;
        else if ( pointwise_interest_quantity == "pos" )
            m.pos_handbook = pos_interest_quantity;
        else {
            cerr << "Arret brutal, car la defintion " << pointwise_interest_quantity << " de la quantite d'interet ponctuelle " << interest_quantity << " n'est pas implementee..." << endl << endl;
            throw "Anguille sous coquille...";
        }
    }
}

/// Calcul de la quantite d'interet locale I
/// ----------------------------------------
template<class TM, class TF, class T, class Pvec>
void calcul_interest_quantity( const TM &m, const TM &m_crown, const TF &f, const TF &f_crown, const string &pb, const string &interest_quantity, const string &direction_extractor, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list_interest_quantity, const unsigned &node_interest_quantity, const Pvec &pos_interest_quantity, const Pvec &pos_crack_tip, const T &angle_crack, const T &radius_Ri, const T &radius_Re, T &I_h, const bool disp = false ) {
    
    I_h = 0.;
    
    static const unsigned dim = TM::dim;
    
    if ( disp ) {
        if ( pb == "direct" )
            cout << "Calcul de la quantite d'interet locale approchee I_h" << endl << endl;
        else if ( pb == "reference" )
            cout << "Calcul de la quantite d'interet locale (quasi-)exacte I_ex" << endl << endl;
    }
    
    TicToc t;
    t.start();
    
    if ( interest_quantity == "mean_sigma" or interest_quantity == "mean_epsilon" ) {
        Calcul_Interest_Quantity_Mean_Sigma_Epsilon calcul_interest_quantity_mean_sigma_epsilon;
        calcul_interest_quantity_mean_sigma_epsilon.elem_list_interest_quantity = &elem_list_interest_quantity;
        
        apply( m.elem_list, calcul_interest_quantity_mean_sigma_epsilon, m, f, I_h );
    }
    else if ( interest_quantity.find("pointwise") != string::npos ) {
        if ( pointwise_interest_quantity == "pos" ) {
            Calcul_Interest_Quantity_Pointwise_Dep_Sigma_Epsilon<Pvec> calcul_interest_quantity_pointwise_dep_sigma_epsilon;
            calcul_interest_quantity_pointwise_dep_sigma_epsilon.interest_quantity = &interest_quantity;
            calcul_interest_quantity_pointwise_dep_sigma_epsilon.direction_extractor = &direction_extractor;
            calcul_interest_quantity_pointwise_dep_sigma_epsilon.pos_interest_quantity = &pos_interest_quantity;
            
            apply( m.elem_list, calcul_interest_quantity_pointwise_dep_sigma_epsilon, m, f, I_h );
        }
        else if ( pointwise_interest_quantity == "node" and interest_quantity == "pointwise_dep" ) {
            Calcul_Interest_Quantity_Pointwise_Dep_Node calcul_interest_quantity_pointwise_dep_node;
            calcul_interest_quantity_pointwise_dep_node.interest_quantity = &interest_quantity;
            
            apply( m.node_list, calcul_interest_quantity_pointwise_dep_node, direction_extractor, node_interest_quantity, I_h );
        }
        else {
            cerr << "Arret brutal, car la defintion " << pointwise_interest_quantity << " de la quantite d'interet ponctuelle " << interest_quantity << " n'est pas implementee..." << endl << endl;
            throw "Anguille sous coquille...";
        }
    }
    else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
        Calcul_Interest_Quantity_SIF<T, Pvec> calcul_interest_quantity_SIF;
        calcul_interest_quantity_SIF.direction_extractor = &direction_extractor;
        calcul_interest_quantity_SIF.pos_crack_tip = &pos_crack_tip;
        calcul_interest_quantity_SIF.angle_crack = &angle_crack;
        calcul_interest_quantity_SIF.radius_Ri = &radius_Ri;
        calcul_interest_quantity_SIF.radius_Re = &radius_Re;
        
        apply( m_crown.elem_list, calcul_interest_quantity_SIF, m_crown, f_crown, I_h );
    }
    else {
        cerr << "Arret brutal, car la quantite d'interet " << interest_quantity << " n'est pas implementee..." << endl << endl;
        throw "Anguille sous coquille...";
    }
    
    if ( pb == "direct" ) {
        cout << "quantite d'interet locale approchee I_h :" << endl;
        cout << "I_h = " << I_h << endl << endl;
    }
    else if ( pb == "reference" ) {
        cout << "quantite d'interet locale (quasi-)exacte I_ex :" << endl;
        cout << "I_ex = " << I_h << endl << endl;
    }
    
    t.stop();
    if ( disp ) {
        if ( pb == "direct" )
            cout << "temps de calcul de la quantite d'interet approchee I_h = " << t.res << " s" << endl << endl;
        else if ( pb == "reference" )
            cout << "temps de calcul de la quantite d'interet (quasi-)exacte I_ex = " << t.res << " s" << endl << endl;
    }
}

/// Calcul de la correction I_hh sur la quantite d'interet locale I (avec ou sans introduction de sigma_hat_m)
/// ----------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class TV, class TVV>
void calcul_correction_interest_quantity( TM &m, TM &m_adjoint, const TF &f, const TF &f_adjoint, const string &interest_quantity, const string &method, const string &method_adjoint, const T &theta, const T &theta_adjoint, const TV &theta_elem_adjoint, const Vec<unsigned> &correspondance_elem_m_adjoint_to_elem_m, const TVV &dep_hat, const TVV &dep_adjoint_hat, const T &I_h, T &I_hh, const bool want_local_enrichment = false, const bool want_introduction_sigma_hat_m = true, const bool disp = false ) {
    
    I_hh = 0.;
    
    if ( disp )
        cout << "Calcul de la correction I_hh sur la quantite d'interet locale I" << endl << endl;
    
    TicToc t;
    t.start();
    
    if ( not want_introduction_sigma_hat_m ) {
        Calcul_Correction_Interest_Quantity_w_sigma_hat_m<TM, TF, T, TVV> calcul_correction_interest_quantity_w_sigma_hat_m;
        calcul_correction_interest_quantity_w_sigma_hat_m.m = &m;
        calcul_correction_interest_quantity_w_sigma_hat_m.m_adjoint = &m_adjoint;
        calcul_correction_interest_quantity_w_sigma_hat_m.f = &f;
        calcul_correction_interest_quantity_w_sigma_hat_m.f_adjoint = &f_adjoint;
        calcul_correction_interest_quantity_w_sigma_hat_m.interest_quantity = &interest_quantity;
        calcul_correction_interest_quantity_w_sigma_hat_m.method = &method;
        calcul_correction_interest_quantity_w_sigma_hat_m.method_adjoint = &method_adjoint;
        calcul_correction_interest_quantity_w_sigma_hat_m.want_local_enrichment = &want_local_enrichment;
        calcul_correction_interest_quantity_w_sigma_hat_m.correspondance_elem_m_adjoint_to_elem_m = &correspondance_elem_m_adjoint_to_elem_m;
        calcul_correction_interest_quantity_w_sigma_hat_m.dep_hat = &dep_hat;
        calcul_correction_interest_quantity_w_sigma_hat_m.dep_adjoint_hat = &dep_adjoint_hat;
        calcul_correction_interest_quantity_w_sigma_hat_m.I_hh = &I_hh;
        
        apply( m_adjoint.elem_list, calcul_correction_interest_quantity_w_sigma_hat_m );
    }
    else {
        Calcul_Correction_Interest_Quantity_wo_sigma_hat_m<TM, TF, T, TVV> calcul_correction_interest_quantity_wo_sigma_hat_m;
        calcul_correction_interest_quantity_wo_sigma_hat_m.m = &m;
        calcul_correction_interest_quantity_wo_sigma_hat_m.m_adjoint = &m_adjoint;
        calcul_correction_interest_quantity_wo_sigma_hat_m.f = &f;
        calcul_correction_interest_quantity_wo_sigma_hat_m.f_adjoint = &f_adjoint;
        calcul_correction_interest_quantity_wo_sigma_hat_m.interest_quantity = &interest_quantity;
        calcul_correction_interest_quantity_wo_sigma_hat_m.method = &method;
        calcul_correction_interest_quantity_wo_sigma_hat_m.want_local_enrichment = &want_local_enrichment;
        calcul_correction_interest_quantity_wo_sigma_hat_m.correspondance_elem_m_adjoint_to_elem_m = &correspondance_elem_m_adjoint_to_elem_m;
        calcul_correction_interest_quantity_wo_sigma_hat_m.dep_hat = &dep_hat;
        calcul_correction_interest_quantity_wo_sigma_hat_m.I_hh = &I_hh;
        
        apply( m_adjoint.elem_list, calcul_correction_interest_quantity_wo_sigma_hat_m );
    }
    
    cout << "correction I_hh sur la quantite d'interet :" << endl;
    cout << "I_hh = " << I_hh << endl << endl;
    
    cout << "nouvelle approximation I_h + I_hh de la quantite d'interet :" << endl;
    cout << "I_h + I_hh = " << I_h + I_hh << endl << endl;
    
    t.stop();
    if ( disp )
        cout << "temps de calcul de la correction I_hh sur la quantite d'interet = " << t.res << " s" << endl << endl;
}

/// Construction d'un champ de contrainte admissible & Calcul d'un estimateur d'erreur globale sur la structure extraite, Affichage de l'estimateur
/// -----------------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class TVV>
void calcul_error_estimate_lambda( const TM &m, TM &m_lambda, const TF &f, TF &f_lambda, const string &pb, const string &method, const string &shape, const T &k, T &theta_lambda, const TVV &dep_hat, TVV &dep_hat_lambda, const bool disp = false ) {
    
    typedef Vec<T> TV;
    
    /// -------------------------------------------------------------------------------------------------------------------------------- ///
    /// Construction d'un champ de contrainte admissible par element & Calcul d'un estimateur d'erreur globale sur la structure extraite ///
    /// -------------------------------------------------------------------------------------------------------------------------------- ///
    
    /// Construction de la correspondance entre maillages extraits et maillages initiaux direct/adjoint
    /// -----------------------------------------------------------------------------------------------
    Vec<unsigned> correspondance_elem_m_lambda_to_elem_m;
    correspondance_elem_m_lambda_to_elem_m.resize( m_lambda.elem_list.size() );
    
    Construct_Correspondance_Elem_Mesh_Extracted_To_Elem_Mesh construct_correspondance_elem_m_extracted_to_elem_m;
    construct_correspondance_elem_m_extracted_to_elem_m.correspondance_elem_m_extracted_to_elem_m = &correspondance_elem_m_lambda_to_elem_m;
    apply_ij( m_lambda.elem_list, m.elem_list, construct_correspondance_elem_m_extracted_to_elem_m );
    
    TicToc t;
    t.start();
    
    /// Calcul d'un estimateur d'erreur globale sur les maillages extraits des maillages initiaux direct/adjoint
    /// --------------------------------------------------------------------------------------------------------
    f_lambda.allocate_matrices();
    f_lambda.shift();
    f_lambda.assemble();
    f_lambda.update_variables();
    f_lambda.call_after_solve();
    
    theta_lambda = 0.;
    TV theta_lambda_elem;
    theta_lambda_elem.resize( m.elem_list.size() );
    theta_lambda_elem.set( 0. );
    
    Calcul_Error_Estimate_Lambda<TM, TF, T, TV, TVV> calcul_error_estimate_lambda;
    calcul_error_estimate_lambda.m = &m;
    calcul_error_estimate_lambda.m_lambda = &m_lambda;
    calcul_error_estimate_lambda.f = &f;
    calcul_error_estimate_lambda.f_lambda = &f_lambda;
    calcul_error_estimate_lambda.method = &method;
    calcul_error_estimate_lambda.correspondance_elem_m_lambda_to_elem_m = &correspondance_elem_m_lambda_to_elem_m;
    calcul_error_estimate_lambda.dep_hat = &dep_hat;
    calcul_error_estimate_lambda.dep_hat_lambda = &dep_hat_lambda;
    calcul_error_estimate_lambda.theta_elem = &theta_lambda_elem;
    calcul_error_estimate_lambda.theta = &theta_lambda;
    
    apply( m_lambda.elem_list, calcul_error_estimate_lambda );
    
//    for (unsigned n=0;n<m_lambda.elem_list.size();++n) {
//        cout << "contribution a l'estimateur d'erreur globale au carre sur la structure extraite associee au pb " << pb << " de type " << shape << " et de taille " << k << " de l'element " << n << " :" << endl;
//        cout << "theta_elem^2 = " << theta_lambda_elem[ n ] << endl;
//    }
    
    theta_lambda = sqrt( theta_lambda );
    m_lambda.ecre = theta_lambda / sqrt(2.);
    m_lambda.error_estimate = theta_lambda;
    
    if ( disp ) {
        cout << "estimateur d'erreur globale sur la structure extraite associee au pb " << pb << " de type " << shape << " et de taille " << k << " :" << endl;
        cout << "theta = " << theta_lambda << endl << endl;
    }
    
    t.stop();
    if ( disp )
        cout << "temps de calcul de l'estimateur d'erreur globale sur la structure extraite associee au pb " << pb << " de type " << shape << " et de taille " << k << " = " << t.res << " s" << endl << endl;
}

/// Calcul d'un estimateur pondere d'erreur globale weighted_theta
/// --------------------------------------------------------------
template<class TM, class TF, class T, class TVV, class Pvec>
void calcul_weighted_error_estimate_lambda( const TM &m, TM &m_lambda, const TF &f, TF &f_lambda, const string &pb, const string &method, const string &shape, const T &h, const Pvec &domain_center, const Vec<T> &domain_length, const T &k_min, T &weighted_theta_lambda, const TVV &dep_hat, const bool disp = false ) {
    
    /// ---------------------------------------------------------------------------------------------------------------------------------------- ///
    /// Construction d'un champ de contrainte admissible par element & Calcul d'un estimateur pondere d'erreur globale sur la structure extraite ///
    /// ---------------------------------------------------------------------------------------------------------------------------------------- ///
    
    /// Construction de la correspondance entre maillages extraits et maillages initiaux direct/adjoint
    /// -----------------------------------------------------------------------------------------------
    Vec<unsigned> correspondance_elem_m_lambda_to_elem_m;
    correspondance_elem_m_lambda_to_elem_m.resize( m_lambda.elem_list.size() );
    
    Construct_Correspondance_Elem_Mesh_Extracted_To_Elem_Mesh construct_correspondance_elem_m_extracted_to_elem_m;
    construct_correspondance_elem_m_extracted_to_elem_m.correspondance_elem_m_extracted_to_elem_m = &correspondance_elem_m_lambda_to_elem_m;
    apply_ij( m_lambda.elem_list, m.elem_list, construct_correspondance_elem_m_extracted_to_elem_m );
    
    TicToc t;
    t.start();
    
    /// Calcul d'un estimateur pondere weighted_theta de l'erreur globale sur les maillages extraits des maillages initiaux direct/adjoint
    /// ----------------------------------------------------------------------------------------------------------------------------------
    f_lambda.allocate_matrices();
    f_lambda.shift();
    f_lambda.assemble();
    f_lambda.update_variables();
    f_lambda.call_after_solve();
    
    weighted_theta_lambda = 0.;
    
    Calcul_Weighted_Error_Estimate_Lambda<TM, TF, T, TVV, Pvec> calcul_weighted_error_estimate_lambda;
    calcul_weighted_error_estimate_lambda.m = &m;
    calcul_weighted_error_estimate_lambda.m_lambda = &m_lambda;
    calcul_weighted_error_estimate_lambda.f = &f;
    calcul_weighted_error_estimate_lambda.f_lambda = &f_lambda;
    calcul_weighted_error_estimate_lambda.method = &method;
    calcul_weighted_error_estimate_lambda.correspondance_elem_m_lambda_to_elem_m = &correspondance_elem_m_lambda_to_elem_m;
    calcul_weighted_error_estimate_lambda.dep_hat = &dep_hat;
    calcul_weighted_error_estimate_lambda.h = &h;
    calcul_weighted_error_estimate_lambda.domain_length = &domain_length;
    calcul_weighted_error_estimate_lambda.domain_center = &domain_center;
    calcul_weighted_error_estimate_lambda.k_min = &k_min;
    calcul_weighted_error_estimate_lambda.weighted_theta = &weighted_theta_lambda;
    
    apply( m_lambda.elem_list, calcul_weighted_error_estimate_lambda );
    
//    for (unsigned n=0;n<m_lambda.elem_list.size();++n) {
//        cout << "contribution a l'estimateur pondere d'erreur globale au carre sur la structure extraite associee au pb " << pb << " de type " << shape << " et de taille " << k_min << " de l'element " <<  n << " :" << endl;
//        cout << "weighted_theta_elem^2 = " << weighted_theta_lambda_elem[ n ] << endl;
//    }
    
    
    if ( disp ) {
        cout << "estimateur pondere d'erreur globale au carre sur la structure extraite associee au pb " << pb << " de type " << shape << " :" << endl;
        cout << "weighted_theta^2 = " << weighted_theta_lambda << endl << endl;
    }
    
    t.stop();
    if ( disp )
        cout << "temps de calcul de l'estimateur pondere d'erreur globale au carre (technique " << method << ") sur la structure extraite associee au pb " << pb << " de type " << shape << " = " << t.res << " s" << endl << endl;
}

/// Construction d'un champ de contrainte admissible & Calcul d'un estimateur d'erreur globale sur le bord de la structure extraite, Affichage de l'estimateur
/// ----------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class TVV, class Pvec>
void calcul_error_estimate_lambda_boundary( const TM &m, TM &m_lambda, const TF &f, TF &f_lambda, const string &pb, const string &method, const string &shape, const Pvec &domain_center, const T &k, T &theta_lambda_boundary, const TVV &dep_hat, const bool disp = false ) {
    
    typedef Vec<T> TV;
    
    /// ------------------------------------------------------------------------------------------------------------------------------------------- ///
    /// Construction d'un champ de contrainte admissible par element & Calcul d'un estimateur d'erreur globale sur le bord de la structure extraite ///
    /// ------------------------------------------------------------------------------------------------------------------------------------------- ///
    
    /// Construction de la correspondance entre maillages extraits et maillages initiaux direct/adjoint
    /// -----------------------------------------------------------------------------------------------
    Vec<unsigned> correspondance_elem_m_lambda_to_elem_m;
    correspondance_elem_m_lambda_to_elem_m.resize( m_lambda.elem_list.size() );
    
    Construct_Correspondance_Elem_Mesh_Extracted_To_Elem_Mesh construct_correspondance_elem_m_extracted_to_elem_m;
    construct_correspondance_elem_m_extracted_to_elem_m.correspondance_elem_m_extracted_to_elem_m = &correspondance_elem_m_lambda_to_elem_m;
    apply_ij( m_lambda.elem_list, m.elem_list, construct_correspondance_elem_m_extracted_to_elem_m );
    
    TicToc t;
    t.start();
    
    /// Calcul d'un estimateur d'erreur globale sur le bord des maillages extraits des maillages initiaux direct/adjoint
    /// ----------------------------------------------------------------------------------------------------------------
    f_lambda.allocate_matrices();
    f_lambda.shift();
    f_lambda.assemble();
    f_lambda.update_variables();
    f_lambda.call_after_solve();
    
    theta_lambda_boundary = 0.;
    TV theta_face_lambda_boundary;
    
    Calcul_Error_Estimate_Lambda_Boundary<TM, TF, T, TV, TVV, Pvec> calcul_error_estimate_lambda_boundary;
    calcul_error_estimate_lambda_boundary.m = &m;
    calcul_error_estimate_lambda_boundary.m_lambda = &m_lambda;
    calcul_error_estimate_lambda_boundary.f = &f;
    calcul_error_estimate_lambda_boundary.f_lambda = &f_lambda;
    calcul_error_estimate_lambda_boundary.method = &method;
    calcul_error_estimate_lambda_boundary.correspondance_elem_m_lambda_to_elem_m = &correspondance_elem_m_lambda_to_elem_m;
    calcul_error_estimate_lambda_boundary.dep_hat = &dep_hat;
    calcul_error_estimate_lambda_boundary.domain_center = &domain_center;
    calcul_error_estimate_lambda_boundary.theta_boundary_face = &theta_face_lambda_boundary;
    calcul_error_estimate_lambda_boundary.theta_boundary = &theta_lambda_boundary;
    
    m_lambda.update_skin();
    apply( m_lambda.elem_list, calcul_error_estimate_lambda_boundary );
    
//    unsigned cpt_face = 0;
//    for (unsigned n=0;n<m_lambda.sub_mesh(Number<1>()).elem_list.size();++n) {
//        if ( m_lambda.sub_mesh(Number<1>()).get_parents_of_EA( m_lambda.sub_mesh(Number<1>()).elem_list[ n ] ).size() == 1 ) {
//            cout << "contribution a l'estimateur d'erreur globale au carre sur le bord de la structure extraite associee au pb " << pb << " de type " << shape << " et de taille " << k << " de la face " << m_lambda.sub_mesh(Number<1>()).elem_list[ n ]->number << " :" << endl;
//            cout << "theta_face^2 = " << theta_face_lambda_boundary[ cpt_face ] << endl;
//            cpt_face++;
//        }
//    }
    
    theta_lambda_boundary = sqrt( theta_lambda_boundary );
    
    if ( disp ) {
        cout << "estimateur d'erreur globale sur le bord de la structure extraite associee au pb " << pb << " de type " << shape << " et de taille " << k << " :" << endl;
        cout << "theta = " << theta_lambda_boundary << endl << endl;
    }
    
    t.stop();
    if ( disp )
        cout << "temps de calcul de l'estimateur d'erreur globale sur la structure extraite associee au pb " << pb << " de type " << shape << " et de taille " << k << " = " << t.res << " s" << endl << endl;
}

/// Calcul du terme gamma pour le calcul des bornes locales ameliorees
/// ------------------------------------------------------------------
template<class TM, class TF, class T, class TVV, class Pvec>
void calcul_gamma( TM &m, TM m_adjoint, TM &m_adjoint_lambda_, const TF &f, const TF &f_adjoint, TF &f_adjoint_lambda_, const string &method, const string &local_improvement, const string &shape, const T &k_min, const T &k_max, const T &k_opt, const T &theta_lambda_min, const T &theta_lambda_max, const T &h, const Pvec &domain_center, const Vec<T> &domain_length, const bool &spread_cut, const TVV &dep_hat, const TVV &dep_adjoint_hat, const string &integration_k, const unsigned &integration_nb_points, T &gamma, const bool disp = false ) {
    
    typedef Vec<T> TV;
    
    gamma = 0.;
    
    static const unsigned dim = TM::dim;
    
    if ( disp )
        cout << "Calcul de la quantite gamma pour la technique " << local_improvement << endl << endl;
    
    TicToc t;
    t.start();
    
//    if ( disp ) {
//        cout << "Construction d'un champ de contrainte admissible par element" << endl;
//        if ( local_improvement == "steklov" )
//            cout << "Calcul d'un estimateur d'erreur globale sur un ensemble de structures extraites associees au pb direct de type " << shape << " pour un parametre lambda compris entre " << k_min << " et " << k_max << endl << endl;
//        else if ( local_improvement == "rayleigh" )
//            cout << "Calcul d'un estimateur d'erreur globale sur le bord d'un ensemble de structures extraites associees au pb adjoint de type " << shape << " pour un parametre lambda compris entre 0 et " << k_opt << endl << endl;
//    }
    
    if ( local_improvement == "steklov" ) {
        TV theta_lambda;
        Vec<T> func_lambda;
        Vec<T> lambda;
        if ( integration_k == "trapeze" ) {
            for (unsigned i=0; i<integration_nb_points; ++i) {
                T k = ( k_max - k_min ) * i / ( integration_nb_points - 1 ) + k_min;
                TM m_lambda;
                set_mesh_domain( m_lambda, m, shape, k, domain_length, domain_center, spread_cut );
                TF f_lambda( m_lambda );
                T theta_k = 0.;
                Vec< Vec<T> > dep_hat_lambda;
                dep_hat_lambda.resize( m_lambda.elem_list.size() );
                calcul_error_estimate_lambda( m, m_lambda, f, f_lambda, "direct", method, shape, k, theta_k, dep_hat, dep_hat_lambda, disp );
                theta_lambda.push_back( pow( theta_k, 2 ) );
                func_lambda.push_back( pow( k / k_min, -1./h ) * pow( theta_k, 2 ) / ( h * k ) );
                lambda.push_back( k );
            }
            for (unsigned i=0; i<integration_nb_points-1; ++i)
                gamma += ( lambda[ i + 1 ] - lambda[ i ] ) * ( func_lambda[ i + 1 ] + func_lambda[ i ] ) / 2;
        }
        else if ( integration_k == "gauss" ) {
            Vec<T> poids;
            Vec<Vec<T, 1> > valeurs;
            gauss_points( Bar(), 13, poids, valeurs );
            for (unsigned n=0; n<poids.size(); ++n) {
                T k = ( k_max - k_min ) * valeurs[ n ][ 0 ] + k_min;
                TM m_lambda;
                set_mesh_domain( m_lambda, m, shape, k, domain_length, domain_center, spread_cut );
                TF f_lambda( m_lambda );
                T theta_k = 0.;
                Vec< Vec<T> > dep_hat_lambda;
                dep_hat_lambda.resize( m_lambda.elem_list.size() );
                calcul_error_estimate_lambda( m, m_lambda, f, f_lambda, "direct", method, shape, k, theta_k, dep_hat, dep_hat_lambda, disp );
                theta_lambda.push_back( pow( theta_k, 2 ) );
                func_lambda.push_back( pow( k / k_min, -1./h ) * pow( theta_k, 2 ) / ( h * k ) );
                lambda.push_back( k );
                gamma += poids[ n ] * pow( k / k_min, -1./h ) * pow( theta_k, 2 ) / ( h * k ) * ( k_max - k_min );
            }
        }
        else if ( integration_k == "IPP" ) {
            TM m_crown_lambda;
            T lambda_min = k_min * domain_length[ 0 ];
            T lambda_max = k_max * domain_length[ 0 ];
            set_mesh_crown( m_crown_lambda, m, domain_center, lambda_min, lambda_max, spread_cut );
            TF f_crown_lambda( m_crown_lambda );
            T weighted_theta_crown_lambda_2 = 0.;
            calcul_weighted_error_estimate_lambda( m, m_crown_lambda, f, f_crown_lambda, "direct", method, shape, h, domain_center, domain_length, k_min, weighted_theta_crown_lambda_2, dep_hat, disp );
            gamma = pow( theta_lambda_min, 2 ) - pow( k_min / k_max, 1./h ) * pow( theta_lambda_max, 2 ) + weighted_theta_crown_lambda_2;
        }
    }
    else if ( local_improvement == "rayleigh" ) {
        T theta_adjoint_lambda_boundary = 0.;
        calcul_error_estimate_lambda_boundary( m_adjoint, m_adjoint_lambda_, f_adjoint, f_adjoint_lambda_, "adjoint", method, shape, domain_center, k_opt, theta_adjoint_lambda_boundary, dep_adjoint_hat, disp );
        unsigned n = dim - 1;
        gamma = 2 * sqrt( h ) / ( h + n + 1 ) * theta_adjoint_lambda_boundary;
    }
    
    cout << "gamma = " << gamma << endl << endl;
    
    t.stop();
    if ( disp )
        cout << "temps de calcul de la quantite gamma = " << t.res << " s" << endl << endl;
}

/// Calcul de la correction I_hhh sur la quantite d'interet locale I pour le calcul des bornes locales ameliorees
/// -------------------------------------------------------------------------------------------------------------
template<class TM, class TF, class T, class TVV>
void calcul_correction_interest_quantity_lambda( const TM &m_lambda_min, const TM &m_adjoint_lambda_min, const TF &f_lambda_min, const TF &f_adjoint_lambda_min, const string &interest_quantity, const string &method, const string &method_adjoint, const TVV &dep_hat_lambda_min, const TVV &dep_adjoint_hat_lambda_min, T &I_hhh, const bool want_introduction_sigma_hat_m = true, const bool disp = false ) {
    
    I_hhh = 0.;
    
    if ( disp )
        cout << "Calcul de la correction I_hhh sur la quantite d'interet locale I" << endl << endl;
    
    TicToc t;
    t.start();
    
    Calcul_Correction_Interest_Quantity_Lambda<TM, TF, T, TVV> calcul_correction_interest_quantity_lambda;
    calcul_correction_interest_quantity_lambda.m = &m_lambda_min;
    calcul_correction_interest_quantity_lambda.m_adjoint = &m_adjoint_lambda_min;
    calcul_correction_interest_quantity_lambda.f = &f_lambda_min;
    calcul_correction_interest_quantity_lambda.f_adjoint = &f_adjoint_lambda_min;
    calcul_correction_interest_quantity_lambda.interest_quantity = &interest_quantity;
    calcul_correction_interest_quantity_lambda.method = &method;
    calcul_correction_interest_quantity_lambda.method_adjoint = &method_adjoint;
    calcul_correction_interest_quantity_lambda.dep_hat = &dep_hat_lambda_min;
    calcul_correction_interest_quantity_lambda.dep_adjoint_hat = &dep_adjoint_hat_lambda_min;
    calcul_correction_interest_quantity_lambda.I_hhh = &I_hhh;
    
    apply_ij( m_lambda_min.elem_list, m_adjoint_lambda_min.elem_list, calcul_correction_interest_quantity_lambda );
    
    if ( not want_introduction_sigma_hat_m )
        I_hhh *= 2;
    
    cout << "correction I_hhh sur la quantite d'interet :" << endl;
    cout << "I_hhh = " << I_hhh << endl;
    
    t.stop();
    if ( disp )
        cout << "temps de calcul de la correction I_hhh sur la quantite d'interet = " << t.res << " s" << endl << endl;
}

/// Calcul de la constante h pour le calcul des bornes locales ameliorees
/// ---------------------------------------------------------------------
template<class TM, unsigned n1, unsigned n2, class T, class S> void calcul_h_alternativeontype( const TM &m, const Number<n1> &, const Number<n2> &, T &h, const S &shape );

template<class TM, class T, class S> void calcul_h_alternativeontype( const TM &m, Number<0>, Number<0>, T &h, const S &shape ) {
    static const unsigned dim = TM::dim;
    if ( dim == 2 ) {
        if ( shape == "circle" )
            h = ( 2 * m.mu + m.la ) / ( 2 * ( m.mu + m.la ) );
        else if ( shape == "circle_crack" )
            h = ( 2 * m.mu + m.la ) / ( 2 * ( m.mu + m.la ) ) + m.mu / ( 6 * M_PI * ( m.mu + m.la ) );
        else if ( shape == "square" or shape == "rectangle" )
            h = ( 7 * m.mu + 3 * m.la ) / ( 6 * ( m.mu + m.la ) );
    }
    else if ( dim == 3 ) {
        if ( shape == "sphere" )
            h = ( 2 * m.mu + m.la ) / ( 2 * m.mu + 3 * m.la );
        else if ( shape == "cube" or shape == "cuboid" )
            h = ( 8 * m.mu + 3 * m.la ) / ( 6 * m.mu + 9 * m.la );
    }
    else
        cerr << "forme " << shape << " en dim " << dim << " non implementee pour amelioration locale basee sur steklov..." << endl << endl;
};

template<class TM, class T, class S> void calcul_h_alternativeontype( const TM &m, Number<1>, Number<0>, T &h, const S &shape ) {
    static const unsigned dim = TM::dim;
    for (unsigned n=0;n<m.elem_list.size();++n) {
        T la = m.elem_list[n]->get_field("la", StructForType<T>());
        if ( dim == 2 ) {
            if ( shape == "cercle" )
                h = max( h, ( 2 * m.mu + la ) / ( 2 * ( m.mu + la ) ) );
            else if ( shape == "circle_crack" )
                h = max( h, ( 2 * m.mu + la ) / ( 2 * ( m.mu + la ) ) + m.mu / ( 6 * M_PI * ( m.mu + la ) ) );
            else if ( shape == "square" or shape == "rectangle" )
                h = max( h, ( 7 * m.mu + 3 * la ) / ( 6 * ( m.mu + la ) ) );
        }
        else if ( dim == 3 ) {
            if ( shape == "sphere" )
                h = max( h, ( 2 * m.mu + la ) / ( 2 * m.mu + 3 * la ) );
            else if ( shape == "cube" or shape == "cuboid" )
                h = max( h, ( 8 * m.mu + 3 * la ) / ( 6 * m.mu + 9 * la ));
        }
        else
            cerr << "forme " << shape << " en dim " << dim << " non implementee pour amelioration locale basee sur steklov..." << endl << endl;
    }
};

template<class TM, class T, class S> void calcul_h_alternativeontype( const TM &m, Number<0>, Number<1>, T &h, const S &shape ) {
    static const unsigned dim = TM::dim;
    for (unsigned n=0;n<m.elem_list.size();++n) {
        T mu = m.elem_list[n]->get_field("mu", StructForType<T>());
        if ( dim == 2 ) {
            if ( shape == "circle" )
                h = max( h, ( 2 * mu + m.la ) / ( 2 * ( mu + m.la ) ) );
            else if ( shape == "circle_crack" )
                h = max( h, ( 2 * mu + m.la ) / ( 2 * ( mu + m.la ) ) + mu / ( 6 * M_PI * ( mu + m.la ) ) );
            else if ( shape == "square" or shape == "rectangle" )
                h = max( h, ( 7 * mu + 3 * m.la ) / ( 6 * ( mu + m.la ) ) );
        }
        else if ( dim == 3 ) {
            if ( shape == "sphere" )
                h = max( h, ( 2 * mu + m.la ) / ( 2 * mu + 3 * m.la ) );
            else if ( shape == "cube" or shape == "cuboid" )
                h = max( h, ( 8 * mu + 3 * m.la ) / ( 6 * mu + 9 * m.la ));
        }
        else
            cerr << "forme " << shape << " en dim " << dim << " non implementee pour amelioration locale basee sur steklov..." << endl << endl;
    }
};

template<class TM, class T, class S> void calcul_h_alternativeontype( const TM &m, Number<1>, Number<1>, T &h, const S &shape ) {
    static const unsigned dim = TM::dim;
    for (unsigned n=0;n<m.elem_list.size();++n) {
        T la = m.elem_list[n]->get_field("la", StructForType<T>());
        T mu = m.elem_list[n]->get_field("mu", StructForType<T>());
        if ( dim == 2 ) {
            if ( shape == "circle" )
                h = max( h, ( 2 * mu + la ) / ( 2 * ( mu + la ) ) );
            else if ( shape == "circle_crack" )
                h = max( h, ( 2 * mu + la ) / ( 2 * ( mu + la ) ) + mu / ( 6 * M_PI * ( mu + la ) ) );
            else if ( shape == "square" or shape == "rectangle" )
                h = max( h, ( 7 * mu + 3 * la ) / ( 6 * ( mu + la ) ) );
        }
        else if ( dim == 3 ) {
            if ( shape == "sphere" )
                h = max( h, ( 2 * mu + la ) / ( 2 * mu + 3 * la ) );
            else if ( shape == "cube" or shape == "cuboid" )
                h = max( h, ( 8 * mu + 3 * la ) / ( 6 * mu + 9 * la ));
        }
        else
            cerr << "forme " << shape << " en dim " << dim << " non implementee pour amelioration locale basee sur steklov..." << endl << endl;
    }
}

/// Calcul du deplacement total
/// ---------------------------
template<class TM>
void calcul_dep_tot_after_solve( TM &m ) {
    for (unsigned i=0;i<m.node_list.size();++i) {
        m.node_list[i].dep_tot = m.node_list[i].dep + m.node_list[i].dep_handbook * m.node_list[i].phi_nodal_handbook;
    }
}

#endif // Calcul_interest_quantity_h
