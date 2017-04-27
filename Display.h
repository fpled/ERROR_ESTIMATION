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
    
    cout << "Estimation d'erreur globale du pb " << pb << endl;
    cout << "-----------------------------------------" << endl << endl;
    
    cout << "methode de construction de champs admissibles = " << method;
    if ( method.find("EET") != string::npos or method.find("EESPT") != string::npos ) {
        if ( not( enhancement_with_geometric_criterium ) and not( enhancement_with_estimator_criterium ) )
            cout << " standard";
        else if ( enhancement_with_geometric_criterium and not( enhancement_with_estimator_criterium ) )
            cout << " avec amelioration geometrique";
        else if ( not( enhancement_with_geometric_criterium ) and enhancement_with_estimator_criterium )
            cout << " avec amelioration sur l'estimateur";
        else
            cout << "avec amelioration geometrique et sur l'estimateur";
    }
    cout << endl;
    if ( method == "EET" or method == "EESPT" ) {
        cout << "fonction-cout = " << cost_function << endl;
        cout << "solveur pour la resolution des pbs de minimisation = " << solver_minimisation << endl;
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
    
    cout << "Quantite d'interet locale" << endl;
    cout << "-------------------------" << endl << endl;
    
    cout << "quantite d'interet = " << interest_quantity << " dans la direction " << direction_extractor << endl;
    cout << "zone d'interet : " << endl;
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

/// Display/Save vtu and pvd files
/// ------------------------------
template<class TM, class T, class Pvec>
string define_prefix( TM &m, const string &pb, const string &structure, const string loading = "pull", const string mesh_size = "fine", const string method = "EET", const bool enhancement_with_geometric_criterium = false, const bool enhancement_with_estimator_criterium = false, const T val_geometric_criterium = 0.34, const T val_estimator_criterium = 0.8, const string geometric_criterium = "radius_ratio", const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool want_global_estimation = false, const bool want_local_estimation = false, const string interest_quantity = "mean_sigma", const string direction_extractor = "xx", const string pointwise_interest_quantity = "node", const Vec<unsigned> elem_list_interest_quantity = Vec<unsigned>( 4886 ), const unsigned node_interest_quantity = 661, const Pvec pos_interest_quantity = Pvec( 49.5, 135.5 ), const Pvec pos_crack_tip = Pvec( 109., 105. ), const T radius_Ri = 6, const T radius_Re = 8, const bool want_local_improvement = false, const string local_improvement = false, const string shape = "circle", const T k_min = 2.5, const T k_max = 7., const T k_opt = 4.4, const bool want_local_enrichment = false, const unsigned nb_layers_nodes_enrichment = 2 ) {
    
    static const unsigned dim = TM::dim;
    
    string prefix = structure;
    if ( structure == "plate_crack" or structure == "spot_weld" or (structure.find("test_specimen") != string::npos and dim == 3) )
        prefix += "_" + loading;
    if ( structure == "plate_hole" or structure == "plate_crack" or structure == "structure_crack" or (structure == "test_specimen" and dim == 2) or structure == "weight_sensor" or structure == "circle" or structure == "beam_hole" or structure == "plate_hole_full" or structure == "spot_weld" or structure == "reactor_head" or structure == "door_seal" or structure == "sphere" or structure == "sphere_center" or structure == "SAP" )
        prefix += "_" + mesh_size;
    for (unsigned n=0;n<m.type_elements().size();++n)
        prefix += "_" + m.type_elements()[n];
    prefix += "_" + pb;
    
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
        if ( interest_quantity.find("mean") != string::npos ) {
            prefix += "_elem";
            for (unsigned n=0;n<elem_list_interest_quantity.size();++n)
                prefix += "_" + to_string( elem_list_interest_quantity[ n ] );
        }
        else if ( interest_quantity.find("pointwise") != string::npos ) {
            if ( pointwise_interest_quantity == "node" )
                prefix += "_node_" + to_string( node_interest_quantity );
            else if ( pointwise_interest_quantity == "pos" ) {
                prefix += "_pos";
                for (unsigned d=0;d<dim;++d)
                    prefix += "_" + to_string( pos_interest_quantity[ d ] );
            }
        }
        else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
            prefix += "_pos_crack_tip";
            for (unsigned d=0;d<dim;++d) {
                prefix += "_" + to_string( pos_crack_tip[ d ] );
            }
            prefix += "_Ri_" + to_string( radius_Ri ) + "_Re_" + to_string( radius_Re );
        }
        if ( pb == "adjoint" )
            prefix += "_" + to_string( m.node_list.size() ) + "_nodes_" + to_string( m.elem_list.size() ) + "_elems";
        if ( want_local_improvement ) {
            prefix += "_local_improvement_" + local_improvement + "_shape_" + shape;
            if ( local_improvement == "steklov" )
                prefix += "_lambda_min_" + to_string( k_min ) + "_max_" + to_string( k_max );
            else if ( local_improvement == "rayleigh" )
                prefix += "_lambda_opt_" + to_string( k_opt );
        }
        if ( want_local_enrichment )
            prefix += "_local_enrichement_" + to_string( nb_layers_nodes_enrichment ) + "_layers";
    }
    return prefix;
}

#endif // Display_h
