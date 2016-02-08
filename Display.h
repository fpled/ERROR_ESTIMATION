//
// C++ Interface: Display
//
// Description: Affichage
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Display_h
#define Display_h

using namespace LMT;
using namespace std;

/// Display dimension, structure and mesh parameters
/// --------------------------------------------------------
void display_pb( const unsigned &dim, const string &structure, const unsigned &deg_p  ) {

    cout << "----------------------------------------------------" << endl;
    cout << "dimension = " << dim << endl;
    cout << "structure = " << structure << endl;
    cout << "degre p   = " << deg_p << endl;
    cout << "----------------------------------------------------" << endl << endl;
    
}

/// Display method for constructing admissible fields (EET,SPET,EESPT) with/without enhancement, cost-function and solver
/// ---------------------------------------------------------------------------------------------------------------------
void display_method( const string &pb, const string &method, const unsigned &cost_function, const bool &enhancement_with_geometric_criterium, const bool &enhancement_with_estimator_criterium, const string &solver, const string &solver_minimisation ) {
    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << "Estimation d'erreur globale du pb " << pb << endl;
    cout << "Methode de construction de champs admissibles = " << method;
    if ( method.find("EET") != string::npos or method.find("EESPT") != string::npos ) {
        if ( enhancement_with_geometric_criterium == 0 and enhancement_with_estimator_criterium == 0 )
            cout << " standard";
        else if ( enhancement_with_geometric_criterium and enhancement_with_estimator_criterium == 0 )
            cout << " avec amelioration geometrique";
        else if ( enhancement_with_geometric_criterium == 0 and enhancement_with_estimator_criterium )
            cout << " avec amelioration sur l'estimateur";
        else
            cout << "avec amelioration geometrique et sur l'estimateur";
    }
    cout << endl;
    cout << "-----------------------------------------------------------------------------------------" << endl << endl;

    if ( method == "EET" or method == "EESPT" ) {
        cout << "fonction-cout = " << cost_function << endl;
        cout << "solveur pour la resolution des pbs de minimisation    = " << solver_minimisation << endl;
        cout << "solveur pour la resolution des pbs locaux par element = " << solver << endl;
    }
    else if ( method == "SPET" )
        cout << "solveur pour la resolution des pbs locaux auto-equilibres par patch = " << solver << endl;
    cout << endl;
}

/// Display quantity of interest and zone of interest
/// -------------------------------------------------
template<class T, class Pvec>
void display_interest_quantity( const string &interest_quantity, const string &direction_extractor, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list, const unsigned &node, const Pvec &pos, const Pvec &pos_crack_tip, const T &angle_crack, const T &radius_Ri, const T &radius_Re ) {

    cout << "-------------------------------------------------------" << endl;
    cout << "Estimation d'erreur locale" << endl;
    cout << "Quantite d'interet : " << interest_quantity << " dans la direction " << direction_extractor << endl;
    cout << "-------------------------------------------------------" << endl << endl;

    cout << "zone d'interet :" << endl;
    if ( interest_quantity.find("mean") != string::npos ) {
        cout << "liste des elements = " << elem_list << endl;
    }
    else if ( interest_quantity.find("pointwise") != string::npos ) {
        if (pointwise_interest_quantity == "node") {
            cout << "noeud = " << node << endl;
        }
        else if (pointwise_interest_quantity == "pos") {
            cout << "position du point = " << pos << endl;
        }
    }
    else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
        cout << "position de la pointe de fissure = " << pos_crack_tip << endl;
        cout << "angle de la fissure = " << angle_crack << " rad" << endl;
        cout << "                    = " << angle_crack * 180 / M_PI << " deg" << endl;
        cout << "rayon interieur de la couronne entourant la pointe de fissure : Ri = " << radius_Ri << endl;
        cout << "rayon exterieur de la couronne entourant la pointe de fissure : Re = " << radius_Re << endl;
    }
    cout << endl;
}

/// Display adjoint mesh parameters
/// -------------------------------
template<class T>
void display_params_adjoint( const bool &want_local_refinement, const T &l_min, const T &k, const bool &spread_cut, const bool &want_local_enrichment, const unsigned &nb_layers_nodes_enrichment, const Vec<unsigned> &elem_list_adjoint_enrichment_zone_1, const Vec<unsigned> &elem_list_adjoint_enrichment_zone_2, const Vec<unsigned> &face_list_adjoint_enrichment_zone_12, const Vec<unsigned> &node_list_adjoint_enrichment, const bool &want_local_improvement, const string &local_improvement, const string &shape, const T &k_min, const T &k_max, const T &k_opt ) {

    if ( want_local_refinement ) {
        cout << "Parametres associes au raffinement local du maillage adjoint" << endl;
        cout << "------------------------------------------------------------" << endl << endl;
        cout << "longueur minimale des cotes des elements du maillage adjoint : l_min = " << l_min << endl;
        cout << "coefficient d'augmentation de la longueur maximale des cotes des elements : k = " << k << endl;
        cout << "propagation du raffinement au reste du maillage = " << spread_cut << endl << endl;
    }
    if ( want_local_enrichment ) {
        cout << "Parametres associes a l'enrichissement local avec fonctions handbook" << endl;
        cout << "--------------------------------------------------------------------" << endl << endl;
        cout << "nb de couches/rangées de noeuds enrichis = " << nb_layers_nodes_enrichment << endl;
        cout << "nb de noeuds enrichis = " << node_list_adjoint_enrichment.size() << endl;
        cout << "nb d'elements enrichis dans Omega_1 = " << elem_list_adjoint_enrichment_zone_1.size() << endl;
        cout << "nb d'elements enrichis dans Omega_2 = " << elem_list_adjoint_enrichment_zone_2.size() << endl;
        cout << "nb de faces enrichies a l'interface entre Omega_1 et Omega_2 = " << face_list_adjoint_enrichment_zone_12.size() << endl;
        cout << "liste des noeuds enrichis = " << node_list_adjoint_enrichment << endl;
        cout << "liste des elements enrichis dans Omega_1 = " << elem_list_adjoint_enrichment_zone_1 << endl;
        cout << "liste des elements enrichis dans Omega_2 = " << elem_list_adjoint_enrichment_zone_2 << endl;
        cout << "liste des faces enrichies a l'interface entre Omega_1 et Omega_2 = " << face_list_adjoint_enrichment_zone_12 << endl << endl;
    }
    if ( want_local_improvement ) {
        cout << "Technique d'amelioration de l'erreur locale : " << local_improvement << endl;
        cout << "------------------------------------------------------" << endl << endl;
        cout << "forme geometrique des domaines homothetiques = " << shape << endl;
        if ( local_improvement == "steklov" ) {
            cout << "parametres des domaines homothetiques : k_min = " << k_min << endl << " et k_max = " << k_max << endl << endl;
        }
        else if ( local_improvement == "rayleigh" )
            cout << "parametre du domaine homothetique : k_opt = " << k_opt << endl << endl;
    }
}

/// Display/Save vtu and pvd files
/// ------------------------------
template<class TM>
string define_prefix( TM &m, const string &pb, const string &structure, const string &loading, const string &mesh_size ) {

    static const unsigned dim = TM::dim;

    string prefix = structure;
    if ( structure == "plate_crack" or structure == "spot_weld" or (structure.find("test_specimen") != string::npos and dim == 3) )
        prefix += "_" + loading;
    if ( structure == "plate_hole" or structure == "plate_crack" or structure == "structure_crack" or (structure == "test_specimen" and dim == 2) or structure == "weight_sensor" or structure == "circle" or structure == "beam_hole" or structure == "plate_hole_full" or structure == "spot_weld" or structure == "reactor_head" or structure == "door_seal" or structure == "sphere" or structure == "sphere_center" or structure == "SAP" )
        prefix += "_" + mesh_size;
    for (unsigned n=0;n<m.type_elements().size();++n)
        prefix += "_" + m.type_elements()[n];
    prefix += "_" + pb;
    return prefix;
}

template<class TM, class T, class Pvec>
void display_vtu_pvd( TM &m, TM &m_ref, TM &m_lambda_min, TM &m_lambda_max, TM &m_lambda_opt, TM &m_crown, const string &pb, const string &method, const string &structure, const string &loading, const string &mesh_size, const bool &enhancement_with_geometric_criterium, const bool &enhancement_with_estimator_criterium, const T &val_geometric_criterium, const T &val_estimator_criterium, const string &geometric_criterium, const unsigned &refinement_level_ref, const bool &want_global_discretization_error, const bool &want_local_discretization_error, const bool &want_global_estimation, const bool &want_local_estimation, const bool &want_local_improvement, const string &interest_quantity, const string &direction_extractor, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list_interest_quantity, const unsigned &node_interest_quantity, const Pvec &pos_interest_quantity, const Pvec &pos_crack_tip, const T &radius_Ri, const T &radius_Re, const string &local_improvement, const string &shape, const T &k_min, const T &k_max, const T &k_opt, const bool &want_local_enrichment, const unsigned &nb_layers_nodes_enrichment, const bool save_vtu = true, const bool display_vtu = false, const bool save_pvd = false, const bool display_pvd = false, const bool save_vtu_ref = false, const bool display_vtu_ref = false, const bool save_vtu_lambda = true, const bool display_vtu_lambda = false, const bool save_vtu_crown = true, const bool display_vtu_crown = false ) {
    
    static const unsigned dim = TM::dim;
    
    string prefix = define_prefix( m, pb, structure, loading, mesh_size );
    string prefix_ref = define_prefix( m_ref, pb, structure, loading, mesh_size );
    
    if ( want_global_estimation or want_local_estimation ) {
        prefix += "_" + method;
        if ( ( enhancement_with_geometric_criterium or enhancement_with_estimator_criterium ) and ( method.find("EET") != string::npos or method.find("EESPT") != string::npos ) ) {
            prefix += "_enhancement";
            if ( enhancement_with_estimator_criterium )
                prefix += "_estimate_ratio_" + to_string( val_estimator_criterium );
            if ( enhancement_with_geometric_criterium )
                prefix += "_" + geometric_criterium + "_ratio_" + to_string( val_geometric_criterium );
        }
    }
    if ( pb == "direct" and ( want_global_discretization_error or want_local_discretization_error ) ) {
        if ( want_global_discretization_error )
            prefix += "_global";
        if ( want_local_discretization_error )
            prefix += "_local";
        prefix += "_discretization_error";
    }
    if ( want_local_estimation ) {
        prefix += "_local_estimation_" + interest_quantity + "_" + direction_extractor;
        prefix_ref += "_local_estimation_" + interest_quantity + "_" + direction_extractor;
        if ( interest_quantity.find("mean") != string::npos ) {
            prefix += "_elem";
            prefix_ref += "_elem";
            for (unsigned n=0;n<elem_list_interest_quantity.size();++n) {
                prefix += "_" + to_string( elem_list_interest_quantity[ n ] );
                prefix_ref += "_" + to_string( elem_list_interest_quantity[ n ] );
            }
        }
        else if ( interest_quantity.find("pointwise") != string::npos ) {
            if ( pointwise_interest_quantity == "node" ) {
                prefix += "_node_" + to_string( node_interest_quantity );
                prefix_ref += "_node_" + to_string( node_interest_quantity );
            }
            else if ( pointwise_interest_quantity == "pos" ) {
                prefix += "_pos";
                prefix_ref += "_pos";
                for (unsigned d=0;d<dim;++d) {
                    prefix += "_" + to_string( pos_interest_quantity[ d ] );
                    prefix_ref += "_" + to_string( pos_interest_quantity[ d ] );
                }
            }
        }
        else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
            prefix += "_pos_crack_tip";
            prefix_ref += "_pos_crack_tip";
            for (unsigned d=0;d<dim;++d) {
                prefix += "_" + to_string( pos_crack_tip[ d ] );
                prefix_ref += "_" + to_string( pos_crack_tip[ d ] );
            }
            prefix += "_Ri_" + to_string( radius_Ri ) + "_Re_" + to_string( radius_Re );
            prefix_ref += "_Ri_" + to_string( radius_Ri ) + "_Re_" + to_string( radius_Re );
        }
        if ( pb == "adjoint" )
            prefix += "_" + to_string( m.node_list.size() ) + "_nodes_" + to_string( m.elem_list.size() ) + "_elems";
    }
    prefix_ref += "_ref_" + to_string( refinement_level_ref );
    
    string prefix_crown = prefix;
    
    if ( want_local_estimation and ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) )
        prefix_crown += "_crown";
    
    string prefix_lambda_min = prefix;
    string prefix_lambda_max = prefix;
    string prefix_lambda_opt = prefix;
    
    if ( want_local_estimation and want_local_improvement ) {
        prefix += "_local_improvement_" + local_improvement + "_shape_" + shape;
        if ( local_improvement == "steklov" ) {
            prefix += "_lambda_min_" + to_string( k_min ) + "_max_" + to_string( k_max );
            prefix_lambda_min += "_lambda_min_" + to_string( k_min );
            prefix_lambda_max += "_lambda_max_" + to_string( k_max );
        }
        else if ( local_improvement == "rayleigh" ) {
            prefix += "_lambda_opt_" + to_string( k_opt );
            prefix_lambda_opt += "_lambda_opt_" + to_string( k_opt );
        }
    }
    
    if ( want_local_estimation and want_local_enrichment ) {
        prefix += "_local_enrichement_nb_layers_" + to_string( nb_layers_nodes_enrichment );
        if ( want_local_improvement ) {
            if ( local_improvement == "steklov" ) {
                prefix_lambda_min += "_local_enrichement_nb_layers_" + to_string( nb_layers_nodes_enrichment );
                prefix_lambda_max += "_local_enrichement_nb_layers_" + to_string( nb_layers_nodes_enrichment );
            }
            else if ( local_improvement == "rayleigh" ) {
                prefix_lambda_opt += "_local_enrichement_nb_layers_" + to_string( nb_layers_nodes_enrichment );
            }
        }
    }
    
    const string prefix_ = prefix;
    const string prefix_crown_ = prefix_crown;
    const string prefix_lambda_min_ = prefix_lambda_min;
    const string prefix_lambda_max_ = prefix_lambda_max;
    const string prefix_lambda_opt_ = prefix_lambda_opt;
    const string prefix_ref_ = prefix_ref;
    
    if ( display_vtu )
        display( m, prefix_ );
    else if ( save_vtu )
        save( m, prefix_ );
    if ( display_pvd or save_pvd ) {
        DisplayParaview dp;
        dp.add_mesh( m, prefix_ );
        if ( display_pvd )
            dp.exec( prefix_ );
        else
            dp.make_pvd_file( prefix_ );
    }
    if ( display_vtu_ref )
        display( m_ref, prefix_ref_ );
    else if ( save_vtu_ref )
        save( m_ref, prefix_ref_ );
    if ( want_local_estimation and ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) ) {
        if ( display_vtu_crown )
            display( m_crown, prefix_crown_ );
        else if ( save_vtu_crown )
            save( m_crown, prefix_crown_ );
    }
    if ( want_local_estimation and want_local_improvement ) {
        if ( local_improvement == "steklov" ) {
            if ( display_vtu_lambda ) {
                display( m_lambda_min, prefix_lambda_min_ );
                display( m_lambda_max, prefix_lambda_max_ );
            }
            else if ( save_vtu_lambda ) {
                save( m_lambda_min, prefix_lambda_min_ );
                save( m_lambda_max, prefix_lambda_max_ );
            }
        }
        else if ( local_improvement == "rayleigh" ) {
            if ( display_vtu_lambda )
                display( m_lambda_opt, prefix_lambda_opt_ );
            else if ( save_vtu_lambda )
                save( m_lambda_opt, prefix_lambda_opt_ );
        }
    }
}

#endif // Display_h
