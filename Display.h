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

/// Display dimension, degree, structure and mesh parameters
/// --------------------------------------------------------
template<class TM>
void display_structure( TM &m, TM &m_ref, const string &pb, const string &structure, const unsigned &deg_p, const bool want_ref = false ) {
    
    static const unsigned dim = TM::dim;
    typedef typename TM::Pvec Pvec;
    
    cout << "-------------" << endl;
    cout << "Dimension : " << dim << endl;
    cout << "-------------" << endl << endl;
    
    cout << "-----------" << endl;
    cout << "Degre p : " << deg_p << endl;
    cout << "-----------" << endl << endl;
    
    cout << "----------------------------------------------------" << endl;
    cout << "Structure : " << structure << endl;
    cout << "----------------------------------------------------" << endl << endl;
    
    cout << "----------------------------------------------------" << endl;
    cout << "Construction de la solution du pb " << pb << endl;
    cout << "----------------------------------------------------" << endl << endl;
    
    if ( m.node_list.size() ) {
        if ( remove_lonely_nodes( m ) )
            cerr << "Des noeuds seuls ont ete retires du maillage " << pb << "..." << endl << endl;
        cout << "nb de ddl du pb " << pb << " : " << m.node_list.size() * dim << endl << endl;
        cout << "nb de noeuds du pb " << pb << " : " << m.node_list.size() << endl << endl;
        cout << "nb d'elements du pb " << pb << " : " << m.elem_list.size() << endl << endl;
    }
    
    if ( want_ref ) {
        if ( m_ref.node_list.size() ) {
            if ( remove_lonely_nodes( m_ref ) )
                cerr << "Des noeuds seuls ont ete retires du maillage de reference associe au pb " << pb << "..." << endl << endl;
            cout << "nb de ddl du pb de reference associe au pb " << pb << " : " << m_ref.node_list.size() * dim << endl << endl;
            cout << "nb de noeuds du pb de reference associe au pb " << pb << " : " << m_ref.node_list.size() << endl << endl;
            cout << "nb d'elements du pb de reference associe au pb " << pb << " : " << m_ref.elem_list.size() << endl << endl;
        }
    }
}

/// Display quantity of interest and zone of interest
/// -------------------------------------------------
template<class T, class Pvec>
void display_interest_quantity( const string &interest_quantity, const string &direction_extractor, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list, const unsigned &node, const Pvec &pos, const Pvec &pos_crack_tip, const T &angle_crack, const T &radius_Ri, const T &radius_Re ) {

    cout << "-------------------------------------------------------" << endl;
    cout << "Quantite d'interet : " << interest_quantity << " dans la direction " << direction_extractor << endl;
    cout << "-------------------------------------------------------" << endl << endl;

    cout << "--------------" << endl;
    cout << "Zone d'interet" << endl;
    cout << "--------------" << endl << endl;
    if ( interest_quantity.find("mean") != string::npos ) {
        cout << "liste des elements du maillage direct :" << endl << elem_list << endl << endl;
    }
    else if ( interest_quantity.find("pointwise") != string::npos ) {
        if (pointwise_interest_quantity == "node") {
            cout << "numero du noeud dans le maillage direct :" << endl << node << endl << endl;
        }
        else if (pointwise_interest_quantity == "pos") {
            cout << "position dans le maillage direct :" << endl << pos << endl << endl;
        }
    }
    else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
        cout << "position de la pointe de fissure :" << endl << pos_crack_tip << endl;
        cout << "angle de la fissure :" << endl << angle_crack << " rad = " << angle_crack * 180 / M_PI << " deg" << endl;
        cout << "rayon du cercle interieur a la couronne omega (qui entoure la pointe de fissure) :" << endl << "Ri = " << radius_Ri << endl;
        cout << "rayon du cercle interieur a la couronne omega (qui entoure la pointe de fissure) :" << endl << "Re = " << radius_Re << endl << endl;
    }
}

/// Display adjoint mesh parameters
/// -------------------------------
template<class T>
void display_params_adjoint( const bool &want_local_refinement, const T &l_min, const T &k, const bool &spread_cut, const bool &want_local_enrichment, const unsigned &nb_layers_nodes_enrichment, const Vec<unsigned> &elem_list_adjoint_enrichment_zone_1, const Vec<unsigned> &elem_list_adjoint_enrichment_zone_2, const Vec<unsigned> &face_list_adjoint_enrichment_zone_12, const Vec<unsigned> &node_list_adjoint_enrichment, const bool &want_local_improvement, const string &local_improvement, const string &shape, const T &k_min, const T &k_max, const T &k_opt ) {

    if ( want_local_refinement ) {
        cout << "------------------------------------------------------------" << endl;
        cout << "Parametres associes au raffinement local du maillage adjoint" << endl;
        cout << "------------------------------------------------------------" << endl << endl;
        cout << "longueur minimale des cotes des elements du maillage adjoint = " << l_min << endl << endl;
        cout << "coefficient d'augmentation de la longueur maximale des cotes des elements = " << k << endl << endl;
        cout << "propagation du raffinement au reste du maillage : " << spread_cut << endl << endl;
    }
    if ( want_local_enrichment ) {
        cout << "--------------------------------------------------------------------" << endl;
        cout << "Parametres associes a l'enrichissement local avec fonctions handbook" << endl;
        cout << "--------------------------------------------------------------------" << endl << endl;
        cout << "nb de couches/rangées de noeuds enrichis par la PUM sur le pb direct = " << nb_layers_nodes_enrichment << endl << endl;
        cout << "nb de noeuds enrichis par la PUM = " << node_list_adjoint_enrichment.size() << endl << endl;
        cout << "liste des noeuds enrichis par la PUM :" << endl << node_list_adjoint_enrichment << endl << endl;
        cout << "nb d'elements enrichis dans la zone Omega_1 = " << elem_list_adjoint_enrichment_zone_1.size() << endl << endl;
        cout << "liste des elements enrichis dans la zone Omega_1 :" << endl << elem_list_adjoint_enrichment_zone_1 << endl << endl;
        cout << "nb d'elements enrichis dans la zone Omega_2 = " << elem_list_adjoint_enrichment_zone_2.size() << endl << endl;
        cout << "liste des elements enrichis dans la zone Omega_2 :" << endl << elem_list_adjoint_enrichment_zone_2 << endl << endl;
        cout << "nb de faces enrichies a l'interface entre les zones Omega_1 et Omega_2 = " << face_list_adjoint_enrichment_zone_12.size() << endl << endl;
        cout << "liste des faces enrichies a l'interface entre les zones Omega_1 et Omega_2 :" << endl << face_list_adjoint_enrichment_zone_12 << endl << endl;
    }
    if ( want_local_improvement ) {
        cout << "------------------------------------------------------" << endl;
        cout << "Technique d'amelioration de l'erreur locale : " << local_improvement << endl;
        cout << "------------------------------------------------------" << endl << endl;
        cout << "forme geometrique des domaines homothetiques : " << shape << endl << endl;
        if ( local_improvement == "steklov" ) {
            cout << "parametres des domaines homothetiques : " << k_min << " et " << k_max << endl << endl;
        }
        else if ( local_improvement == "rayleigh" ) {
            cout << "parametre du domaine homothetique : " << k_opt << endl << endl;
        }
    }
}

/// Display method for constructing admissible fields (EET,SPET,EESPT) with/without enhancement, cost-function and solver
/// ---------------------------------------------------------------------------------------------------------------------
void display_method( const string &pb, const string &method, const unsigned &cost_function, const bool &enhancement_with_geometric_criterium, const bool &enhancement_with_estimator_criterium, const string &solver, const string &solver_minimisation ) {
    cout << "-----------------------------------------------------------------------------------------" << endl;
    cout << "Methode de construction de champs admissibles associee au pb " << pb << " : " << method;
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
        cout << "Fonction-cout : " << cost_function << endl << endl;
        cout << "Solveur pour la resolution des problemes de minimisation : " << solver_minimisation << endl << endl;
    }
    cout << "Solveur pour la resolution des problemes locaux : " << solver << endl << endl;
}

/// Display/Save vtu and pvd files
/// ------------------------------
template<class TM>
string define_prefix( TM &m, const string &pb, const string &structure, const string &loading, const string &mesh_size ) {
    string prefix = structure;
    if ( structure == "plate_crack" or structure == "spot_weld" )
        prefix += "_" + loading;
    if ( structure == "plate_hole" or structure == "plate_crack" or structure == "structure_crack" or structure == "test_specimen" or structure == "weight_sensor" or structure == "circle" or structure == "beam_hole" or structure == "plate_hole_full" or structure == "spot_weld" or structure == "reactor_head" or structure == "door_seal" or structure == "sphere" or structure == "sphere_center" or structure == "SAP" )
        prefix += "_" + mesh_size;
    prefix += "_" + m.type_elements()[0] + "_" + pb;
    return prefix;
}

template<class TM, class T, class Pvec>
void display_vtu_pvd( TM &m, TM &m_ref, TM &m_lambda_min, TM &m_lambda_max, TM &m_lambda_opt, TM &m_crown, const string &pb, const string &method, const string &structure, const string &loading, const string &mesh_size, const unsigned &cost_function, const bool &enhancement_with_geometric_criterium, const bool &enhancement_with_estimator_criterium, const T &val_geometric_criterium, const T &val_estimator_criterium, const string &geometric_criterium, const unsigned &deg_k, const unsigned &refinement_deg_ref, const bool &want_global_discretization_error, const bool &want_local_discretization_error, const bool &want_global_estimation, const bool &want_local_estimation, const bool &want_local_improvement, const string &interest_quantity, const string &direction_extractor, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list_interest_quantity, const unsigned &node_interest_quantity, const Pvec &pos_interest_quantity, const Pvec &pos_crack_tip, const T &radius_Ri, const T &radius_Re, const string &local_improvement, const string &shape, const T &k_min, const T &k_max, const T &k_opt, const bool &want_local_enrichment, const unsigned &nb_layers_nodes_enrichment, const bool save_vtu = true, const bool display_vtu = false, const bool save_pvd = false, const bool display_pvd = false, const bool save_vtu_ref = false, const bool display_vtu_ref = false, const bool save_vtu_lambda = true, const bool display_vtu_lambda = false, const bool save_vtu_crown = true, const bool display_vtu_crown =false ) {
    
    static const unsigned dim = TM::dim;
    
    string prefix = define_prefix( m, pb, structure, loading, mesh_size );
    string prefix_ref = define_prefix( m_ref, pb, structure, loading, mesh_size );
    
    if ( ( want_global_estimation or want_local_estimation ) and ( method.find("EET") != string::npos or method.find("SPET") != string::npos or method.find("EESPT") != string::npos ) ) {
        prefix += "_" + method + "_k_" + to_string( deg_k ) + "_J" + to_string( cost_function );
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
    prefix_ref += "_ref_" + to_string( refinement_deg_ref );
    
    string prefix_crown = prefix;
    
    if ( want_local_estimation and ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) ) {
        prefix_crown += "_crown";
    }
    
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
    
    if ( display_vtu ) {
        display( m, prefix_ );
    }
    else if ( save_vtu ) {
        DisplayParaview dp;
        dp.set_mesh( m, prefix_ );
    }
    if ( save_pvd or display_pvd ) {
        DisplayParaview dp;
        dp.add_mesh( m, prefix_ );
        if ( display_pvd )
            dp.exec( prefix_ );
        else
            dp.make_pvd_file( prefix_ );
    }
    if ( display_vtu_ref ) {
        display( m_ref, prefix_ref_ );
    }
    else if ( save_vtu_ref ) {
        DisplayParaview dp;
        dp.set_mesh( m_ref, prefix_ref_ );
    }
    if ( want_local_estimation and ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) ) {
        if ( display_vtu_crown ) {
            display( m_crown, prefix_crown_ );
        }
        else if ( ( save_vtu_crown ) ) {
            DisplayParaview dp;
            dp.set_mesh( m_crown, prefix_crown_ );
        }
    }
    if ( want_local_estimation and want_local_improvement ) {
        if ( local_improvement == "steklov" ) {
            if ( display_vtu_lambda ) {
                display( m_lambda_min, prefix_lambda_min_ );
                display( m_lambda_max, prefix_lambda_max_ );
            }
            else if ( save_vtu_lambda ) {
                DisplayParaview dp_min;
                dp_min.set_mesh( m_lambda_min, prefix_lambda_min_ );
                DisplayParaview dp_max;
                dp_max.set_mesh( m_lambda_max, prefix_lambda_max_ );
            }
        }
        else if ( local_improvement == "rayleigh" ) {
            if ( display_vtu_lambda )
                display( m_lambda_opt, prefix_lambda_opt_ );
            else if ( save_vtu_lambda ) {
                DisplayParaview dp_opt;
                dp_opt.set_mesh( m_lambda_opt, prefix_lambda_opt_ );
            }
        }
    }
}

#endif // Display_h
