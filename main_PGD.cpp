//
// C++ Implementation: main_cpp
//
// Description: Global/Goal-oriented error estimation methods
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2009
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "build/problem_space/all_in_one.h" // sert a forcer le logiciel scons a generer le repertoire build et ses codes sources .h et .cpp correspondant a la formulation
#include "build/problem_parameter/all_in_one.h" // sert a forcer le logiciel scons a generer le repertoire build et ses codes sources .h et .cpp correspondant a la formulation
#include "Structure.h"
#include "Material_properties.h"
#include "Boundary_conditions.h"
#include "GEOMETRY/Calcul_geometry.h"
#include "GEOMETRY/Geometry.h"
#include "DISCRETIZATION_ERROR/Calcul_discretization_error.h"
#include "Display.h"
#include "Calcul_global_error_estimation.h"
#include "Calcul_goal_oriented_error_estimation.h"
#include "PGD/PGD.h"

using namespace LMT;
using namespace std;

int main( int argc, char **argv ) {
    TicToc t_total;
    t_total.start();
    static const unsigned dim = 2;
    static const bool wont_add_nz = true; // pour completer plus vite la matrice de rigidite si ajout de termes non nuls. 
                                          // si wont_add_nz == true et " -DWITH_CHOLMOD " dans les parametres de compilation, pas de solveur iteratif
    typedef Mesh<Mesh_carac_space<double,dim> > TM; // definit un synonyme du type Mesh<Mesh_carac_space<double,dim> >
    typedef Formulation<TM,FormulationElasticity,DefaultBehavior,double,wont_add_nz> TF; // definit un synonyme du type Formulation<TM,FormulationElasticity>
    typedef TM::Pvec Pvec; // definit un synonyme du type TM::Pvec servant a representer la position d'un point
    typedef TM::TNode::T T;
    typedef Mesh<Mesh_carac_parameter<double,1> > TM_param;
    typedef Formulation<TM_param,FormulationUnknownParam,DefaultBehavior,double,wont_add_nz> TF_unknown_param;
    typedef Formulation<TM_param,FormulationKnownParam,DefaultBehavior,double,wont_add_nz> TF_known_param;
    typedef TM_param::Pvec Pvec_param;
    static const string structure = "weight_sensor"; // structure 2D : plate_traction, plate_flexion, plate_hole, plate_crack, structure_crack, eprouvette, weight_sensor, circle
                                                     // structure 3D : beam_traction, beam_flexion, beam_hole, plate_hole, plate_hole_full, hub_rotor_helico, reactor_head, door_seal, spot_weld, blade, pipe, SAP, spherical_inclusions, sphere, sphere_center, sphere_hollow
    static const string mesh_size = "fine"; // mesh_size pour les structures plate_hole (2D ou  3D), plate_crack, structure_crack, test_specimen, weigth_sensor, spot_weld (3D), reactor_head (3D) : coarse, fine
    static const string loading = "pull"; // loading pour la structure spot_weld (3D) : pull, shear, peeling et pour la structure plate_crack (2D) : pull, shear
                                          // loading pour square_... (2D) : pre_epsilon, pre_sigma
    static const unsigned deg_p = 1; // degre de l'analyse elements finis : 1, 2, ...
    static const unsigned deg_k = 3; // degre supplementaire : 1, 2 , 3, ...
    static const string boundary_condition_D = "penalisation"; // methodes de prise en compte des conditions aux limites de Dirichlet (en deplacement) pour le pb direct : lagrange, penalisation
    static const bool display_constraints = 0; // affichage des contraintes cinematiques
    
    /// Global discretization error
    ///----------------------------
    static const bool want_global_discretization_error = 0; // calcul de l'erreur de discretisation globale associee au pb direct
    static const bool want_local_discretization_error = 0; // calcul de l'erreur de discretisation locale associee au pb direct
    static const bool want_solve_ref = 0; // calcul d'une solution de reference sur un maillage de reference (tres fin)
    static const unsigned refinement_degree_ref = 2; // degre du h-refinement pour la construction du maillage de reference associe au pb de reference :
                                                     // 1 -> sous-decoupage en 4/8 elements en 2D/3D
                                                     // 2 -> sous-decoupage en 16/64 elements en 2D/3D
                                                     // 3 -> sous-decoupage en 64/512 elements en 2D/3D
                                                     // 4 -> sous-decoupage en 256/4096 elements en 2D/3D
                                                     // 5 -> sous-decoupage en 1024/32768 elements en 2D/3D
                                                     // 6 -> sous-decoupage en 4096/32768 elements en 2D/3D
                                                     // 7 -> sous-decoupage en 16384/262144 elements en 2D/3D
                                                     // 8 -> sous-decoupage en 65536/2097152 elements en 2D/3D

    /// Global error estimation method
    ///-------------------------------
    static const bool want_global_estimation = 0; // calcul de l'erreur globale au sens de la norme energetique usuelle
    static const string method = "EET"; // methodes de construction de champs admissibles pour le pb direct : EET, EESPT, SPET, EET_SPET_EESPT, EET_SPET, EET_EESPT, SPET_EESPT
    static const string method_adjoint = "EET"; // methodes de construction de champs admissibles pour le pb adjoint : EET, EESPT, SPET, EET_SPET_EESPT, EET_SPET, EET_EESPT, SPET_EESPT

    static const unsigned cost_function = 0; // fonction-cout pour les methodes EET, EESPT :
                                             // 0 : norme matricielle M sans coefficient de ponderation : M = Identite
                                             // 1 : norme matricielle M avec coeff de ponderation en 1/mes(face)^2
                                             // 2 : norme energetique
    static const T pen_N = 1e5; // coefficient de penalisation pour la prise en compte des conditions aux limites de Neumann (en effort) pour la methode EESPT
    static const string solver = "LDL"; // types de solveur pour la resolution des problemes locaux avec blocage auto du noyau : CholMod (sym, def, pos), LDL (sym) // types de solveur sans blocage auto du noyau (-> ne marche pas!) : CholFactorize (sym, def, pos), LUFactorize, Inv, UMFPACK
    static const string solver_minimisation = "UMFPACK"; // types de solveur pour la resolution des problemes de minimisation : LDL (sym), UMFPACK, LUFactorize, Inv
    
    /// Enhanced technique
    ///-------------------
    static const bool enhancement_with_geometric_criterium = 0; // amelioration de la construction des densites d'effort (pour les methodes EET, EESPT) basee sur un critere geometrique
    static const string geometric_criterium = "radius_ratio"; // choix du type de critere d'amelioration geometrique : radius_ratio, edge_ratio
    static const T val_geometric_criterium = 0.34; // valeur du critere d'amelioration geometrique pour le choix des elements dont les densites d'effort doivent etre ameliorees
    
    static const bool enhancement_with_estimator_criterium = 0; // amelioration de la construction des densites d'effort (pour les methodes EET, EESPT) basee sur un critere sur l'estimateur d'erreur theta
    static const T val_estimator_criterium = 0.8; // valeur du critere d'amelioration sur l'estimateur d'erreur theta (% a la contribution max au carre a l'estimateur d'erreur theta) pour le choix des elements dont les densites d'effort doivent etre ameliorees
    
    /// Goal-oriented error estimation method
    ///--------------------------------------
    static const bool want_local_estimation = 0; // calcul de l'erreur sur une quantite d'interet locale I
    static const bool want_interest_quantity_only = 0; // calcul de la quantite d'interet locale I uniquement
    static const bool want_introduction_sigma_hat_m = 1; // introduction de sigma_hat_m pour le calcul de l'erreur sur une quantite d'interet locale I
    static const bool want_local_refinement_adjoint = 1; // raffinement local du mailage adjoint
    static const bool want_local_enrichment = 0; // enrichissement local avec fonctions handbook
    static const bool want_handbook_only = 0; // calcul de la solution handbook uniquement
    static const bool want_local_improvement = 0; // amelioration des bornes pour le calcul de l'erreur sur une quantite d'interet locale I
    static const bool want_compute_eig_local_improvement = 0; // resolution du pb aux valeurs propres generalisees pour determiner la constante impliquee dans l'amelioration des bornes d'erreur locale
    static const bool use_mask_eig_local_improvement = 0; // utilisation d'un masque (image) pour definir le maillage associe au pb aux valeurs propres generalisees
    static const bool want_compute_local_ref = 0; // calcul de la quantite d'interet sur le maillage de reference
    static const string interest_quantity = "pointwise_dep"; // quantite d'interet : mean_sigma, mean_epsilon, pointwise_dep, pointwise_sigma, pointwise_epsilon, SIF (stress intensity factor)
    static const string direction_extractor = "x"; // direction de l'operateur d'extraction pour quantite d'interet mean_sigma, mean_epsilon, pointwise_sigma, pointwise_epsilon : xx, yy, xy, zz, xz, yz
                                                   // direction de l'operateur d'extraction pour quantite d'interet pointwise_dep : x, y, z
                                                   // direction de l'operateur d'extraction pour quantite d'interet SIF (stress intensity factor) : I, II, III
    
    /// Zone of interest
    ///-----------------
    static const Vec<unsigned> list_elems_interest_quantity( 4886 ); // liste des elements definissant la zone d'interet pour quantite d'interet mean_sigma, mean_epsilon
    static const string pointwise_interest_quantity = "node"; // definition de la quantite d'interet ponctuelle : node, pos
    static const unsigned node_interest_quantity( 661 ); // noeud definissant la zone d'interet pour quantite d'interet pointwise_dep, pointwise_sigma, pointwise_epsilon
    static const Pvec pos_interest_quantity( 49.5, 135.5 ); // position definissant la zone d'interet pour quantite d'interet pointwise_dep, pointwise_sigma, pointwise_epsilon
    static const Pvec pos_crack_tip( 109, 105. ); // position de la pointe de fissure : ( 3.5, 0. ) pour plate_crack et ( 109., 105. ) pour structure_crack
    static const T angle_crack = atan2( -17, -3 ); // angle de la fissure (en rad) : 0. pour plate_crack et atan2( -17, -3 ) pour structure_crack
    static const T radius_R1 = 6; // rayon du cercle C_1 interieur a la couronne omega (qui entoure la pointe de fissure) : 1.6 pour plate_crack et 6 pour structure_crack
    static const T radius_R2 = 8; // rayon du cercle C_2 exterieur a la couronne omega (qui entoure la pointe de fissure) : 3.4 pour plate_crack et 8 pour structure_crack
    
    /// Local refinement parameters for adjoint problem
    ///------------------------------------------------
    // on decide de couper le cote d'un element ( i.e. une Bar ) si sa longueur est superieure à d * k + l_min ou d est la distance entre le milieu du cote et le centre
    static const T local_refinement_adjoint_l_min = 1.0; // la longueur minimale des cotes des elements du maillage adjoint
    static const T local_refinement_adjoint_k = 1.0; // le coefficient d'augmentation de la longueur maximale des cotes des elements en fonction de la distance au point, au cercle, ... autour duquel on souhaite raffiner le maillage
    static const bool spread_cut = true; // on propage le raffinement au reste du maillage, i.e. on étend la coupe si l'arête coupée n'est pas la plus longue de l'élément
    
    /// Local enrichment with handbook functions
    ///--------------------------------------------------------
    static const unsigned nb_layers_nodes_enrichment = 2; // nombre de couches/rangées de noeuds enrichis par la PUM sur le pb direct
    
    /// Improved goal-oriented error estimation methods based on Steklov/Rayleigh constants
    ///------------------------------------------------------------------------------------
    static const string local_improvement = "rayleigh"; // choix du type d'amelioration pour le calcul de l'erreur sur une quantite d'interet locale I, basee sur la constante de Steklov ou sur le quotient de Rayleigh : steklov, rayleigh
    static const string shape = "circle"; // forme geometrique des domaines homothetiques associee au calcul ameliore des bornes d'erreur locales
    static const T k_min = 2.5; 
    static const T k_max = 7.; // parametres associes aux domaines homothetiques servant au calcul ameliore des bornes d'erreur locales avec amelioration steklov
                                // pour un cercle de rayon lambda = k * Rc : facteur multiplicatif devant le rayon Rc du cercle (sphere) circonscrit(e) au triangle (tetra)
    static const T k_opt = 4.4; // parametre associe domaine homothetique servant au calcul ameliore des bornes d'erreur locales avec amelioration rayleigh
    static const string integration_k = "trapeze"; // pour amelioration steklov : type d'integration sur le parametre k : gauss, trapeze, IPP
    static const unsigned integration_nb_steps = 100; // pour amelioration steklov : nb d'intervalles pour l'integration type trapeze sur le parametre k
    
    /// Parameters PGD
    ///---------------
    static const bool want_PGD = 1; // methode de reduction de modele : PGD
    static const bool want_normalization_space = 1; // normalisation des fonctions psi en espace
    static Vec<unsigned> list_elems_PGD_unknown_parameter; // liste des elements definissant la zone avec parametre materiau inconnu
    static const unsigned nb_modes_max = 4; // nb max de modes dans la decomposition
    static unsigned nb_modes; // nb de modes dans la decomposition (critere d'arret global sur le residu au sens faible associe a la solution a variables separees avec precision tol_global_convergence_criterium)
    static const unsigned nb_iterations_max = 10; // nb d'iterations max de l'algorithme de type point fixe
    static Vec<unsigned> nb_iterations; // nb d'iterations de l'algorithme de type point fixe pour chaque mode (critere d'arret local sur la stationnarite du produit des fonctions espace-parametre avec precision tol_local_convergence_criterium)
    static Vec<T> residual_mode; // residu au sens faible associe a la solution a variables separees pour chaque mode
    static T residual = 0.; // residu au sens faible associe a la solution a variables separees
    static const T tol_global_convergence_criterium = 1e-4; // precision pour critere d'arret global (boucle sur les modes)
    static const T tol_local_convergence_criterium = 1e-8; // precision pour critere d'arret local (processus iteratif)
    static const T initial_val_param = 0.0; // valeur initiale sur l'intervalle des parametres
    static const T final_val_param = 1.0; // valeur finale sur l'intervalle des parametres
    static const unsigned nb_steps_param = 100; // nb d'elements du maillage parametrique
    static const bool want_verif_kinematic_PGD = 1; // verification de la decomposition PGD cinematique
    static const unsigned nb_vals_param_verif = 5; // nb de valeurs des parametres pris aleatoirement servant a la verification de la decomposition PGD
    
    /// Parameters for iterative solver
    ///--------------------------------
    static const bool want_iterative_solver = 0; // resolution du pb direct avec solveur iteratif
    static const bool want_iterative_solver_ref = 0; // resolution du pb de reference associe au pb direct avec solveur iteratif
    static const bool want_iterative_solver_adjoint = 0; // resolution du pb adjoint avec solveur iteratif
    static const bool want_iterative_solver_adjoint_ref = 0; // resolution du pb de reference associe au pb adjoint avec solveur iteratif
    static const T criterium_iterative_solver = 1e-3; // critere de resolution pour le solveur iteratif du pb direct : residu en norme inf
    static const T criterium_iterative_solver_ref = 1e-3; // critere de resolution pour le solveur iteratif du pb de reference associe au pb direct : residu en norme inf
    static const T criterium_iterative_solver_adjoint = 1e-3; // critere de resolution pour le solveur iteratif du pb adjoint : residu en norme inf
    static const T criterium_iterative_solver_adjoint_ref = 1e-3; // critere de resolution pour le solveur iteratif du pb de reference associe au pb adjoint : residu en norme inf
    
    /// Verification equilibre / solveur
    ///---------------------------------
    static const bool verif_eq = 1; // verification de l'equilibre global elements finis
    static const bool verif_compatibility_conditions = 1; // verification des conditions de compatibilite dans la methode EET correspondant a l'equilibre elements finis
    static const T tol_compatibility_conditions = 1e-6; // tolerance pour la verification des conditions de compatibilite dans la methode EET correspondant a l'equilibre elements finis
    static const bool verif_eq_force_fluxes = 1; // verification de l'equilibre des densites d'effort pour les methodes EET, EESPT
    static const T tol_eq_force_fluxes = 1e-6; // tolerance pour la verification de l'equilibre des densites d'effort pour les methodes EET, EESPT

    static const bool verif_solver = 1; // verification du solveur pour la resolution des problemes locaux des methodes EET, SPET, EESPT
    static const T tol_solver = 1e-6; // tolerance pour la verification du solveur pour la resolution des problemes locaux des methodes EET, SPET, EESPT
    static const bool verif_solver_enhancement = 0; // verification du solveur pour la resolution des problemes locaux pour l'amelioration des methodes EET EESPT
    static const T tol_solver_enhancement = 1e-6; // tolerance pour la verification du solveur pour la resolution des problemes locaux pour l'amelioration des methodes EET EESPT
    static const bool verif_solver_minimisation = 1; // verification du solver pour la resolution des problemes de minimisation des methodes EET, EESPT
    static const T tol_solver_minimisation = 1e-6; // tolerance pour la verification du solver pour la resolution des problemes de minimisation des methodes EET, EESPT
    static const bool verif_solver_minimisation_enhancement = 0; // verification du solver pour la resolution des problemes de minimisation pour l'amelioration des methodes EET, EESPT
    static const T tol_solver_minimisation_enhancement = 1e-6; // tolerance pour la verification du solver pour la resolution des problemes de minimisation pour l'amelioration des methodes EET, EESPT
    
    /// Debug
    ///------
    static const bool debug_method = 0; // debug total des methodes EET, SPET, EESPT pour le pb direct
    static const bool debug_method_adjoint = 0; // debug total des methodes EET, SPET, EESPT pour le pb adjoint
    static const bool debug_method_enhancement = 0; // debug total de l'amelioration des methodes EET, EESPT du pb direct
    static const bool debug_method_enhancement_adjoint = 0; // debug total de l'amelioration des methodes EET, EESPT du pb adjoint
    static const bool debug_criterium_enhancement = 0; // debug du critere d'amelioration du pb direct
    static const bool debug_criterium_enhancement_adjoint = 0; // debug du critere d'amelioration du pb adjoint
    static const bool debug_geometry = 0; // debug de la geometrie du pb direct
    static const bool debug_geometry_adjoint = 0; // debug de la geometrie du pb adjoint
    static const bool debug_discretization_error = 0; // debug de la contribution a l'erreur de discretisation du pb direct
    static const bool debug_discretization_error_adjoint = 0; // debug de la contribution a l'erreur de discretisation du pb adjoint
    static const bool debug_force_fluxes = 0; // debug des densites d'effort pour les methodes EET et EESPT du pb direct
    static const bool debug_force_fluxes_adjoint = 0; // debug des densites d'effort pour les methodes EET et EESPT du pb adjoint
    static const bool debug_force_fluxes_enhancement = 0; // debug des densites d'effort ameliorees pour les methodes EET et EESPT du pb direct
    static const bool debug_force_fluxes_enhancement_adjoint = 0; // debug des densites d'effort ameliorees pour les methodes EET et EESPT du pb adjoint
    static const bool debug_error_estimate = 0; // debug de l'estimateur d'erreur globale du pb direct
    static const bool debug_error_estimate_adjoint = 0; // debug de l'estimateur d'erreur globale du pb adjoint
    static const bool debug_local_effectivity_index = 0; // debug de l'indice d'efficacite local du pb direct
    static const bool debug_local_effectivity_index_adjoint = 0; // debug de l'indice d'efficacite local du pb adjoint
    
    /// Sauvegarde / Affichage
    ///-----------------------
    static const bool save_vtu = 1;
    static const bool display_vtu = 0;
    static const bool save_pvd = 0;
    static const bool display_pvd = 0;
    static const bool save_vtu_ref = 0;
    static const bool display_vtu_ref = 0;
    
    static const bool save_vtu_adjoint = 1;
    static const bool display_vtu_adjoint = 0;
    static const bool save_pvd_adjoint = 0;
    static const bool display_adjoint_pvd = 0;
    static const bool save_vtu_local_ref = 0;
    static const bool display_vtu_local_ref = 0;
    
    static const bool save_vtu_lambda = 1;
    static const bool display_vtu_lambda = 0;
    
    static const bool save_vtu_adjoint_lambda = 1;
    static const bool display_vtu_adjoint_lambda = 0;
    
    static const bool save_vtu_crown = 1;
    static const bool display_vtu_crown = 0;
    
    static const bool save_pvd_PGD_space = 1;
    static const bool display_pvd_PGD_space = 0;
    static const bool save_pvd_PGD_parameter = 1;
    static const bool display_pvd_PGD_parameter = 0;
    static const bool save_plot_PGD_parameter = 1;
    static const bool save_pvd_PGD_space_verif = 1;
    static const bool display_pvd_PGD_space_verif = 0;
    
    ///--------------------------------------------///
    /// Construction de la solution elements finis ///
    ///--------------------------------------------///
    
    /// Maillage du pb direct
    ///----------------------
    TM m; // declaration d'un maillage de type TM
    TM m_ref;
    create_structure( m, m_ref, "direct", structure, mesh_size, loading, deg_p, want_global_discretization_error, want_local_discretization_error, refinement_degree_ref, want_solve_ref );
    display_structure( m, m_ref, "direct", structure, deg_p, want_solve_ref );
    TM_param m_param;
    create_structure_param( m_param, initial_val_param, final_val_param, nb_steps_param );
    
    /// Formulation du pb direct
    ///-------------------------
    TF f( m ); // creation d'une formulation du type TF avec le maillage m
    TF f_ref( m_ref );
    TF_unknown_param f_unknown_param( m_param );
    TF_known_param f_known_param( m_param );
    
    /// Proprietes materiaux et Conditions aux limites du pb direct
    ///------------------------------------------------------------
    create_material_properties( f, m, structure, loading );
    create_boundary_conditions( f, m, boundary_condition_D, "direct", structure, loading, mesh_size );
    define_unknown_parameter_zone( f, m, structure, list_elems_PGD_unknown_parameter );
    if ( want_solve_ref ) {
        create_material_properties( f_ref, m_ref, structure, loading );
        create_boundary_conditions( f_ref, m_ref, boundary_condition_D, "direct", structure, loading, mesh_size );
    }
    
    /// Defintion des fonctions a variables separees
    ///---------------------------------------------
    Vec< Vec<T> > dep_psi;
    Vec< Vec<T> > dep_lambda;
    Vec< Vec< Vec<T>, nb_iterations_max+1 >, nb_modes_max > psi;
    Vec< Vec< Vec<T>, nb_iterations_max+1 >, nb_modes_max > lambda;
    nb_iterations.resize( nb_modes_max );
    nb_iterations.set( 1 );
    residual_mode.resize( nb_modes_max );
    residual_mode.set( 0. );
    Vec<T> vals_param;
    for (unsigned i=0;i<m_param.node_list.size();++i)
        vals_param.push_back( m_param.node_list[i].pos[0] );
    Vec< DisplayParaview, nb_modes_max > dp_space, dp_parameter;
    Vec< DisplayParaview, nb_vals_param_verif > dp_space_verif;
    string prefix_space = define_prefix( m, "direct", structure, loading, mesh_size ) + "_space_mesh";
    string prefix_parameter = define_prefix( m, "direct", structure, loading, mesh_size ) + "_parameter_mesh";
    Vec<string> lp_space;
    lp_space.push_back( "dep" );
    lp_space.push_back( "young_eff" );
    Vec<string> lp_parameter;
    lp_parameter.push_back( "dep" );
    
    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    TMatSymSparse K_s, K_unk_s, K_k_s, K_unk_p, K_k_p;
    Vec<T> F_s, F_p;
    
    /// Verification des conditions cinematiques
    ///-----------------------------------------
    verification_kinematic_conditions( f, verif_kinematic_conditions );
    
    /// Resolution du pb direct
    ///------------------------
    TicToc t;
    t.start();
    if ( want_PGD == 0 ) {
        if ( want_iterative_solver == 0 )
            f.solve();
        else
            f.solve( criterium_iterative_solver );
    }
    else {
        f.allocate_matrices();
        f.shift();
        f.assemble();
        K_s = f.matrices(Number<0>());
        F_s = f.sollicitation;
        for (unsigned i=0;i<m.elem_list.size();++i) {
            if ( find( list_elems_PGD_unknown_parameter, _1 == i ) )
                m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", 1. );
            else
                m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", 0. );
        }
        f.assemble( true, false );
        K_unk_s = f.matrices(Number<0>());
        for (unsigned i=0;i<m.elem_list.size();++i) {
            if ( find( list_elems_PGD_unknown_parameter, _1 == i ) )
                m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", 0. );
            else
                m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", 1. );
        }
        f.assemble( true, false );
        K_k_s = f.matrices(Number<0>());
        f_unknown_param.allocate_matrices();
        f_unknown_param.shift();
        f_unknown_param.assemble();
        K_unk_p = f_unknown_param.matrices(Number<0>());
        F_p = f_unknown_param.sollicitation;
        f_known_param.allocate_matrices( true, false );
        f_known_param.shift();
        f_known_param.assemble( true, false );
        K_k_p = f_known_param.matrices(Number<0>());
        unsigned n = 0;
        while ( true ) {
//             cout << "Mode n = " << n << endl << endl;
            string prefix_mode = "mode_" + to_string( n );
            
            /// Initialisation
            ///---------------
            unsigned k = 0;
            lambda[ n ][ 0 ].resize( f_unknown_param.vectors[0].size() );
            lambda[ n ][ 0 ].set( 1. );
            f_unknown_param.vectors[0] = lambda[ n ][ 0 ];
            f_unknown_param.update_variables();
            f_unknown_param.call_after_solve();
            psi[ n ][ 0 ].resize( f.vectors[0].size() );
            solve_PGD_space( m, f, n, k, nb_iterations, F_s, F_p, K_unk_p, K_k_p, want_iterative_solver_EF, criterium_iterative_solver_EF, list_elems_PGD_unknown_parameter, lambda, psi );
            if ( want_normalization_space ) {
                psi[ n ][ k ] /= sqrt( dot( psi[ n ][ k ], K_s * psi[ n ][ k ] ) );
                f.vectors[0] = psi[ n ][ k ];
                f.update_variables();
                f.call_after_solve();
            }
            
            if ( save_PGD_space_pvd or display_PGD_space_pvd )
                dp_space[ n ].add_mesh_iter( m, prefix_space + "_" + prefix_mode, lp_space, k );
            if ( save_PGD_parameter_pvd or display_PGD_parameter_pvd )
                dp_parameter[ n ].add_mesh_iter( m_param, prefix_parameter + "_" + prefix_mode, lp_parameter, k );
            
            while ( true ) {
                k++;
//                 cout << "Iteration k = " << k << endl << endl;
                string prefix_iteration = "iteration_" + to_string( k );
                psi[ n ][ k ].resize( f.vectors[0].size() );
                lambda[ n ][ k ].resize( f_unknown_param.vectors[0].size() );
                
                psi[ n ][ k ] = psi[ n ][ k-1 ];
                
                /// Construction et resolution du pb en parametre
                ///----------------------------------------------
                solve_PGD_param( m_param, f_unknown_param, n, k, nb_iterations, F_p, K_unk_p, K_k_p, F_s, K_unk_s, K_k_s, psi, lambda );
                
                /// Construction et resolution du pb en espace
                ///-------------------------------------------
                solve_PGD_space( m, f, n, k, nb_iterations, F_s, F_p, K_unk_p, K_k_p, want_iterative_solver_EF, criterium_iterative_solver_EF, list_elems_PGD_unknown_parameter, lambda, psi );
                
                /// Normalisation de la fonction psi en espace
                ///-------------------------------------------
                if ( want_normalization_space ) {
                    psi[ n ][ k ] /= sqrt( dot( psi[ n ][ k ], K_s * psi[ n ][ k ] ) );
                    f.vectors[ 0 ] = psi[ n ][ k ];
                    f.update_variables();
                    f.call_after_solve();
                }
                
//                 cout << "psi_" + to_string( k ) << endl;
//                 cout << psi[ n ][ k ] << endl;
//                 cout << "lambda_" + to_string( k ) << endl;
//                 cout << lambda[ n ][ k ] << endl << endl;
                
                if ( save_pvd_PGD_space or display_pvd_PGD_space )
                    dp_space[ n ].add_mesh_iter( m, prefix_space + "_" + prefix_mode + "_iteration", lp_space, k );
                if ( save_pvd_PGD_parameter or display_pvd_PGD_parameter )
                    dp_parameter[ n ].add_mesh_iter( m_param, prefix_parameter + "_" + prefix_mode + "_iteration", lp_parameter, k );
                
                if ( save_plot_PGD_parameter ) {
                    static GnuPlot gp;
                    gp.set_terminal( "postscript eps enhanced color" );
                    stringstream graph_name;
                    graph_name << "'" << prefix_parameter + "_" + prefix_mode + "_" + prefix_iteration << ".eps'";
                    gp.set_output(graph_name.str().c_str());
                    gp.set("key left bottom");
                    gp.set("format '%g'");
                    gp.set("font \"Times-Roman\"");
                    gp.set_xlabel("'{/Symbol q}'");
                    gp.set_ylabel("'{/Symbol l}'");
                    gp.hold_on();
                    gp.plot( vals_param, lambda[ n ][ k ], "title '{/Symbol l}({/Symbol q})' axes x1y1 w l lt 1 lw 1" );
                    gp.hold_off();
                }
                
                T alpha_p_unk = dot( psi[ n ][ k ], K_unk_s * psi[ n ][ k ] );
                T alpha_p_k = dot( psi[ n ][ k ], K_k_s * psi[ n ][ k ] );
                T alpha_p_unk_old = dot( psi[ n ][ k-1 ], K_unk_s * psi[ n ][ k-1 ] );
                T alpha_p_k_old = dot( psi[ n ][ k-1 ], K_k_s * psi[ n ][ k-1 ] );
                T beta_p_unk = dot( psi[ n ][ k ], K_unk_s * psi[ n ][ k-1 ] );
                T beta_p_k = dot( psi[ n ][ k ], K_k_s * psi[ n ][ k-1 ] );
                
                T alpha_s_unk = dot( lambda[ n ][ k ], K_unk_p * lambda[ n ][ k ] );
                T alpha_s_k = dot( lambda[ n ][ k ], K_k_p * lambda[ n ][ k ] );
                T alpha_s_unk_old = dot( lambda[ n ][ k-1 ], K_unk_p * lambda[ n ][ k-1 ] );
                T alpha_s_k_old = dot( lambda[ n ][ k-1 ], K_k_p * lambda[ n ][ k-1 ] );
                T beta_s_unk = dot( lambda[ n ][ k ], K_unk_p * lambda[ n ][ k-1 ] );
                T beta_s_k = dot( lambda[ n ][ k ], K_k_p * lambda[ n ][ k-1 ] );
                
                /// Stationnarite du produit fonction psi en espace * fonction lambda en parametre dans le processus iteratif associe au mode n
                ///----------------------------------------------------------------------------------------------------------------------------
                T dist = sqrt( fabs( alpha_s_unk * alpha_p_unk + alpha_s_k * alpha_p_k + alpha_s_unk_old * alpha_p_unk_old + alpha_s_k_old * alpha_p_k_old - 2 * ( beta_s_unk * beta_p_unk + beta_s_k * beta_p_k ) ) ) / ( 0.5 * sqrt( fabs( alpha_s_unk * alpha_p_unk + alpha_s_k * alpha_p_k + alpha_s_unk_old * alpha_p_unk_old + alpha_s_k_old * alpha_p_k_old + 2 * ( beta_s_unk * beta_p_unk + beta_s_k * beta_p_k ) ) ) );
//                 cout << "Distance a l'iteration " << k << " = " << dist << endl << endl;
                if ( dist < tol_local_convergence_criterium or k >= nb_iterations_max ) {
                    nb_iterations[ n ] = k;
                    cout << "Nb d'iterations associe au mode " << n << " = " << nb_iterations[ n ] << endl << endl;
                    break;
                }
            }
            
            if ( display_pvd_PGD_space )
                dp_space[ n ].exec( prefix_space + "_" + prefix_mode + ".pvd" );
            else if ( save_pvd_PGD_space )
                dp_space[ n ].make_pvd_file( prefix_space + "_" + prefix_mode + ".pvd" );
            if ( display_pvd_PGD_parameter )
                dp_parameter[ n ].exec( prefix_parameter + "_" + prefix_mode + ".pvd" );
            else if ( save_pvd_PGD_parameter )
                dp_parameter[ n ].make_pvd_file( prefix_parameter + "_" + prefix_mode + ".pvd" );
            
            /// Residu au sens faible associe a la solution u_n
            ///------------------------------------------------
            residual = 0.;
            T residual_sollicitation = 0.;
            for (unsigned i=0;i<n+1;++i) {
                T gamma_p = dot( F_s, psi[ i ][ nb_iterations[ i ] ] );
                T gamma_s = dot( F_p, lambda[ i ][ nb_iterations[ i ] ] );
                residual_sollicitation += gamma_p * gamma_s;
                residual -= gamma_p * gamma_s;
                for (unsigned j=0;j<n+1;++j) {
                    T alpha_p_unk = dot( psi[ i ][ nb_iterations[ i ] ], K_unk_s * psi[ j ][ nb_iterations[ j ] ] );
                    T alpha_p_k = dot( psi[ i ][ nb_iterations[ i ] ], K_k_s * psi[ j ][ nb_iterations[ j ] ] );
                    
                    T alpha_s_unk = dot( lambda[ i ][ nb_iterations[ i ] ], K_unk_p * lambda[ j ][ nb_iterations[ j ] ] );
                    T alpha_s_k = dot( lambda[ i ][ nb_iterations[ i ] ], K_k_p * lambda[ j ][ nb_iterations[ j ] ] );
                    
                    residual += alpha_p_unk * alpha_s_unk + alpha_p_k * alpha_s_k;
                }
            }
            
//             cout << "mode " << n << endl;
//             cout << "iteration " << k << endl;
//             cout << residual << endl;
//             cout << residual_sollicitation << endl;
            residual = fabs( residual ) / fabs( residual_sollicitation );
            residual_mode[ n ] = residual;
            cout << "Residu au sens faible associe a la solution u_" << n << " = " << residual_mode[ n ] << endl << endl;
            if ( /*fabs( residual_mode[ n ] ) < tol_global_convergence_criterium or*/ n >= nb_modes_max-1 ) {
                nb_modes = n+1;
                cout << "Nb de modes = " << nb_modes << endl << endl;
                break;
            }
            n++;
        }
        
        dep_psi.resize( nb_modes );
        dep_lambda.resize( nb_modes );
        for (unsigned n=0;n<nb_modes;++n) {
            dep_psi[ n ] = psi[ n ][ nb_iterations[ n ] ];
            dep_lambda[ n ] = lambda[ n ][ nb_iterations[ n ] ];
        }
        
        if ( want_verif_kinematic_PGD ) {
            /// Verification pour un jeu connu de parametres
            ///---------------------------------------------
            for (unsigned p=0;p<nb_vals_param_verif;++p) {
                unsigned ind_in_vals_param = rand() % nb_steps_param;
                for (unsigned i=0;i<m.elem_list.size();++i) {
                    if ( find( list_elems_PGD_unknown_parameter, _1 == i ) )
                        m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", 1. + vals_param[ ind_in_vals_param ] );
                    else
                        m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", 1. );
                }
                if ( want_iterative_solver_EF == 0 )
                    f.solve();
                else
                    f.solve( criterium_iterative_solver_EF );
                if ( save_PGD_space_verif_pvd or display_PGD_space_verif_pvd )
                    dp_space_verif[ p ].add_mesh_iter( m, prefix_space + "_verification_param_REF_" + to_string( vals_param[ ind_in_vals_param ] ), lp_space, 0 );
                f.vectors[0].set( 0. );
                for (unsigned n=0;n<nb_modes;++n) {
                    f.vectors[0] += dep_psi[ n ] * dep_lambda[ n ][ ind_in_vals_param ];
                }
                f.update_variables();
                f.call_after_solve();
                if ( save_PGD_space_verif_pvd or display_PGD_space_verif_pvd )
                    dp_space_verif[ p ].add_mesh_iter( m, prefix_space + "_verification_param_PGD_" + to_string( vals_param[ ind_in_vals_param ] ), lp_space, 1 );
                
                if ( display_pvd_PGD_space_verif )
                    dp_space_verif[ p ].exec( prefix_space + "_verification_param_" + to_string( vals_param[ ind_in_vals_param ] ) + ".pvd" );
                else if ( save_pvd_PGD_space_verif )
                    dp_space_verif[ p ].make_pvd_file( prefix_space + "_verification_param_" + to_string( vals_param[ ind_in_vals_param ] ) + ".pvd" );
            }
        }
        
    }
    t.stop();
    cout << "Temps de calcul du pb " << pb << " : " << t.res << endl << endl;
    
    Vec<T> dep_part;
    Vec<T> kappa;
    kappa.resize( f_known_param.vectors[0].size() );
    kappa.set( 1. );
    
    /// Construction d'un champ admissible
    ///-----------------------------------
    if ( want_PGD ) {
        for (unsigned i=0;i<m.elem_list.size();++i)
            m.elem_list[i]->set_field( "phi_elem_PGD_unknown_parameter", 1. );
        f.assemble();
        if ( want_iterative_solver_EF == 0 )
            f.solve_system();
        else
            f.solve_system( criterium_iterative_solver_EF );
        f.update_variables();
        f.call_after_solve();
        
        dep_part = f.vectors[0];
    }
    
    /// Verification de l'equilibre du pb direct
    ///-----------------------------------------
    check_equilibrium( f, "direct", verif_eq );
    
    if ( want_solve_ref and want_PGD == 0 ) {
        /// Resolution du pb de reference associe au pb direct
        ///---------------------------------------------------
        TicToc t_ref;
        t_ref.start();
        if ( want_iterative_solver_ref == 0 )
            f_ref.solve();
        else
            f_ref.solve( criterium_iterative_solver_ref );
        t_ref.stop();
        cout << "Temps de calcul du pb de reference associe au pb direct : " << t_ref.res << endl << endl;
        
        /// Verification de l'equilibre du pb de reference associe au pb direct
        ///--------------------------------------------------------------------
        check_equilibrium( f_ref, "de reference associe au pb direct", verif_eq );
    }
    
    ///---------------------------------------------------------------///
    /// Calcul et Affichage des informations relatives a la geometrie ///
    ///---------------------------------------------------------------///
    
    calcul_display_geometry( m, f, debug_geometry );
    
    ///---------------------------------------------------------------------------///
    /// Mesure de l'erreur de discretisation globale et locale associee pb direct ///
    ///---------------------------------------------------------------------------///
    
    if ( want_PGD == 0 )
        calcul_discretization_error( m, m_ref, f, f_ref, debug_discretization_error, want_global_discretization_error, want_local_discretization_error, want_solve_ref );
    
    T theta = 0.;
    Vec<T> theta_elem;
    Vec< Vec<T> > dep_part_hat;
    Vec< Vec< Vec<T> > > dep_psi_hat;
    dep_psi_hat.resize( nb_modes );
 /*   
    if ( want_global_estimation or ( want_local_estimation and want_handbook_only == 0 and want_interest_quantity_only == 0 ) ) {
        
        ///---------------------------------------------------------------------------------------------------------------///
        /// Construction d'un champ de contrainte admissible et Calcul d'un estimateur d'erreur globale associe pb direct ///
        ///---------------------------------------------------------------------------------------------------------------///
        
        calcul_global_error_estimation( f, m, "direct", method, cost_function, pen_N, solver, solver_minimisation, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, geometric_criterium, val_geometric_criterium, val_estimator_criterium, verif_compatibility_conditions, tol_compatibility_conditions, verif_eq_force_fluxes, tol_eq_force_fluxes, verif_solver, tol_solver, verif_solver_enhancement, tol_solver_enhancement, verif_solver_minimisation, tol_solver_minimisation, verif_solver_minimisation_enhancement, tol_solver_minimisation_enhancement, want_global_discretization_error, want_local_discretization_error, want_local_enrichment, theta, theta_elem, dep_hat, debug_geometry, debug_force_fluxes, debug_force_fluxes_enhancement, debug_criterium_enhancement, debug_error_estimate, debug_local_effectivity_index, debug_method, debug_method_enhancement );
        
    }
    
    TM m_adjoint, m_local_ref, m_lambda_min, m_lambda_max, m_lambda_opt, m_adjoint_lambda_min, m_adjoint_lambda_max, m_adjoint_lambda_opt, m_crown;
    static const bool want_compute_adjoint_ref = 0;
    static const bool want_global_discretization_error_adjoint = 0;
    static const bool want_local_discretization_error_adjoint = 0;
    static const bool save_adjoint_crown_vtu = 0;
    static const bool display_adjoint_crown_vtu = 0;
    
    if ( want_local_estimation ) {
        
        ///---------------------------------------///
        /// Construction de la quantite d'interet ///
        ///---------------------------------------///
        
        display_interest_quantity( interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_R1, radius_R2 );
        
        /// Definition de l'extracteur
        ///---------------------------
        if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" )
            create_structure_crown( m, m_crown, pos_crack_tip, radius_R1, radius_R2, spread_cut );
        TF f_crown( m_crown );
        define_extractor( m, m_crown, f, f_crown, interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_R1, radius_R2, want_local_enrichment );
        
        ///------------------------------------------------------///
        /// Calcul de la quantite d'interet locale approchee I_h ///
        ///------------------------------------------------------///
        
        T I_h;
        calcul_interest_quantity( m, m_crown, f, f_crown, "direct", interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_R1, radius_R2, I_h );
        
        ///------------------------------------------------------------///
        /// Calcul de la quantite d'interet locale (quasi-)exacte I_ex ///
        ///------------------------------------------------------------///
        
        T I_ex;
        if ( want_compute_local_ref ) {
            static const bool want_compute_local_ref_ref = 0;
            static const bool want_global_discretization_error_local_ref = 0;
            static const bool want_local_discretization_error_local_ref = 0;
            create_structure( m_local_ref, m_local_ref, "direct", structure, mesh_size, loading, deg_p, want_global_discretization_error_local_ref, want_local_discretization_error_local_ref, refinement_degree_ref, want_compute_local_ref_ref );
            
            Vec<unsigned> list_elems_local_ref_interest_quantity;
            unsigned node_local_ref_interest_quantity;
            create_structure_local_ref( m, m_local_ref, deg_p, refinement_degree_ref, interest_quantity, list_elems_interest_quantity, list_elems_local_ref_interest_quantity, node_interest_quantity, node_local_ref_interest_quantity, pos_crack_tip, radius_R1, radius_R2, spread_cut );
            
            /// Formulation du pb de reference local
            ///-------------------------------------
            TF f_local_ref( m_local_ref );
            
            /// Proprietes materiaux et Conditions aux limites du pb de reference local
            ///------------------------------------------------------------------------
            create_material_properties( f_local_ref, m_local_ref, structure, loading );
            create_boundary_conditions( f_local_ref, m_local_ref, boundary_condition_D, "direct", structure, loading, mesh_size );
            
            /// Resolution du pb de reference local
            ///------------------------------------
            TicToc t_local_ref;
            t_local_ref.start();
            if ( want_iterative_solver_ref == 0 )
                f_local_ref.solve();
            else
                f_local_ref.solve( criterium_iterative_solver_ref );
            t_local_ref.stop();
            cout << "Temps de calcul du pb de reference local associe au pb direct : " << t_local_ref.res << endl << endl;
            
            /// Verification de l'equilibre du pb de reference local associe au pb direct
            ///--------------------------------------------------------------------------
            check_equilibrium( f_local_ref, "de reference local associe au pb direct", verif_eq );
            
            /// Definition de l'extracteur du pb de reference local
            ///----------------------------------------------------
            TM m_crown_ref;
            if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" )
                create_structure_crown( m_local_ref, m_crown_ref, pos_crack_tip, radius_R1, radius_R2, spread_cut );
            TF f_crown_ref( m_crown_ref );
            define_extractor( m_local_ref, m_crown_ref, f_local_ref, f_crown_ref, interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_local_ref_interest_quantity, node_local_ref_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_R1, radius_R2, want_local_enrichment );
            
            /// Calcul de la quantite d'interet locale (quasi-)exacte I_ex
            ///-----------------------------------------------------------
            calcul_interest_quantity( m_local_ref, m_crown_ref, f_local_ref, f_crown_ref, "reference", interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_local_ref_interest_quantity, node_local_ref_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_R1, radius_R2, I_ex );
        }
        
        if ( want_interest_quantity_only == 0 ) {
            
            ///-------------------------------------///
            /// Construction de la solution adjoint ///
            ///-------------------------------------///
            
            Vec<unsigned> list_elems_adjoint_interest_quantity;
            Vec<unsigned> list_elems_adjoint_enrichment_zone_1;
            Vec<unsigned> list_elems_adjoint_enrichment_zone_2;
            Vec<unsigned> list_faces_adjoint_enrichment_zone_12;
            unsigned node_adjoint_interest_quantity;
            Vec<unsigned> list_nodes_adjoint_enrichment;
            
            /// Maillage du pb adjoint
            ///-----------------------
            create_structure( m_adjoint, m_local_ref, "adjoint", structure, mesh_size, loading, deg_p, want_global_discretization_error_adjoint, want_local_discretization_error_adjoint, refinement_degree_ref, want_compute_adjoint_ref );
            create_structure_adjoint( m, m_adjoint, deg_p, interest_quantity, direction_extractor, want_local_refinement_adjoint, local_refinement_adjoint_l_min, local_refinement_adjoint_k, pointwise_interest_quantity, list_elems_interest_quantity, list_elems_adjoint_interest_quantity, node_interest_quantity, node_adjoint_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_R1, radius_R2, spread_cut, want_local_enrichment, nb_layers_nodes_enrichment, list_elems_adjoint_enrichment_zone_1, list_elems_adjoint_enrichment_zone_2, list_faces_adjoint_enrichment_zone_12, list_nodes_adjoint_enrichment, debug_geometry, debug_geometry_adjoint );
            
            display_structure( m_adjoint, m_local_ref, "adjoint", structure, deg_p, want_compute_local_ref );
            display_params_adjoint( want_local_refinement_adjoint, local_refinement_adjoint_l_min, local_refinement_adjoint_k, spread_cut, want_local_enrichment, nb_layers_nodes_enrichment, list_elems_adjoint_enrichment_zone_1, list_elems_adjoint_enrichment_zone_2, list_faces_adjoint_enrichment_zone_12, list_nodes_adjoint_enrichment, want_local_improvement, local_improvement, shape, k_min, k_max, k_opt );
            
            /// Formulation du pb adjoint
            ///--------------------------
            TF f_adjoint( m_adjoint );
            
            /// Conditions aux limites du pb adjoint
            ///-------------------------------------
            create_material_properties( f_adjoint, m_adjoint, structure, loading );
            create_boundary_conditions( f_adjoint, m_adjoint, boundary_condition_D, "adjoint", structure, loading, mesh_size );
            create_load_conditions( m_adjoint, f_adjoint, m, m_crown, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, interest_quantity, direction_extractor, pointwise_interest_quantity, want_local_enrichment );
            
            /// Verification des contraintes cinematiques
            ///------------------------------------------
            check_constraints( f_adjoint, display_constraints );
            
            /// Resolution du pb adjoint
            ///-------------------------
            TicToc t_adjoint;
            t_adjoint.start();
            if ( want_iterative_solver_adjoint == 0 )
                f_adjoint.solve();
            else
                f_adjoint.solve( criterium_iterative_solver_adjoint );
            t_adjoint.stop();
            cout << "Temps de calcul du pb adjoint : " << t_adjoint.res << endl << endl;
            
            if ( want_local_enrichment )
                calcul_dep_tot_after_solve( m_adjoint );
            
            /// Verification de l'equilibre du pb adjoint
            ///------------------------------------------
            check_equilibrium( f_adjoint, "adjoint", verif_eq );
            
            ///-----------------------------------------------------------------------------------///
            /// Calcul et Affichage des informations relatives a la geometrie du maillage adjoint ///
            ///-----------------------------------------------------------------------------------///
            
            calcul_display_geometry( m_adjoint, f_adjoint, debug_geometry_adjoint );
            
            if ( want_handbook_only == 0 ) {
                
                ///-------------------------------------------------------------------------------------------------------------------///
                /// Construction d'un champ de contrainte admissible et Calcul d'un estimateur d'erreur globale associe au pb adjoint ///
                ///-------------------------------------------------------------------------------------------------------------------///
                
                T theta_adjoint = 0.;
                Vec<T> theta_adjoint_elem;
                Vec< Vec<T> > dep_adjoint_hat;
                calcul_global_error_estimation( f_adjoint, m_adjoint, "adjoint", method_adjoint, cost_function, pen_N, solver, solver_minimisation, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, geometric_criterium, val_geometric_criterium, val_estimator_criterium, verif_compatibility_conditions, tol_compatibility_conditions, verif_eq_force_fluxes, tol_eq_force_fluxes, verif_solver, tol_solver, verif_solver_enhancement, tol_solver_enhancement, verif_solver_minimisation, tol_solver_minimisation, verif_solver_minimisation_enhancement, tol_solver_minimisation_enhancement, want_global_discretization_error_adjoint, want_local_discretization_error_adjoint, want_local_enrichment, theta_adjoint, theta_adjoint_elem, dep_adjoint_hat, debug_geometry_adjoint, debug_force_fluxes_adjoint, debug_force_fluxes_enhancement_adjoint, debug_criterium_enhancement_adjoint, debug_error_estimate_adjoint, debug_local_effectivity_index_adjoint, debug_method_adjoint, debug_method_enhancement_adjoint );
                
                /// Construction de la correspondance entre maillages extraits et maillages initiaux direct/adjoint
                ///------------------------------------------------------------------------------------------------
                Vec<unsigned> correspondance_elem_m_adjoint_to_elem_m;
                correspondance_elem_m_adjoint_to_elem_m.resize( m_adjoint.elem_list.size() );
                
                Construct_Correspondance_Elem_Mesh_Extracted_To_Elem_Mesh construct_correspondance_elem_m_adjoint_to_elem_m;
                construct_correspondance_elem_m_adjoint_to_elem_m.correspondance_elem_m_extracted_to_elem_m = &correspondance_elem_m_adjoint_to_elem_m;
                apply_ij( m_adjoint.elem_list, m.elem_list, construct_correspondance_elem_m_adjoint_to_elem_m );
                
                ///-----------------------------------------------------------------------------------------------------///
                /// Calcul de la correction I_hh (avec ou sans introduction de sigma_hat_m) sur la quantite d'interet I ///
                ///-----------------------------------------------------------------------------------------------------///
                
                T I_hh = 0.;
                calcul_correction_interest_quantity( m, m_adjoint, f, f_adjoint, interest_quantity, method, method_adjoint, want_local_enrichment, theta, theta_adjoint, theta_adjoint_elem, correspondance_elem_m_adjoint_to_elem_m, dep_hat, dep_adjoint_hat, I_h, want_introduction_sigma_hat_m, I_hh );
                
                ///------------------------------------------------------------------------///
                /// Calcul standard des bornes d'erreur sur la quantite d'interet locale I ///
                ///------------------------------------------------------------------------///
                
                calcul_standard_local_error_bounds( m, m_adjoint, f, f_adjoint, method, theta, theta_adjoint, theta_adjoint_elem, correspondance_elem_m_adjoint_to_elem_m, dep_hat, I_h, I_hh, want_introduction_sigma_hat_m );
                
                ///------------------------------------------------------------------------///
                /// Calcul ameliore des bornes d'erreur sur la quantite d'interet locale I ///
                ///------------------------------------------------------------------------///
                if ( want_local_improvement ) {
                    calcul_enhanced_local_error_bounds( m, m_adjoint, f, f_adjoint, m_lambda_min, m_lambda_max, m_lambda_opt, m_adjoint_lambda_min, m_adjoint_lambda_opt, deg_p, method, method_adjoint, structure, loading, mesh_size, cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, val_geometric_criterium, val_estimator_criterium, geometric_criterium, deg_k, local_improvement, shape, k_min, k_max, k_opt, interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_R1, radius_R2, spread_cut, theta, theta_adjoint, theta_adjoint_elem, correspondance_elem_m_adjoint_to_elem_m, dep_hat, dep_adjoint_hat, I_h, I_hh, integration_k, integration_nb_steps, debug_method_adjoint, debug_method_enhancement_adjoint, debug_geometry_adjoint, debug_error_estimate_adjoint, want_introduction_sigma_hat_m, want_compute_eig_local_improvement, use_mask_eig_local_improvement );
                }
            }
        }
    }/*
    
    t_total.stop();
    cout << "Temps de calcul total : " << t_total.res << endl << endl;
    
    ///-----------///
    /// Affichage ///
    ///-----------///
    
//     display_vtu_pvd( m, m_ref, m_lambda_min, m_lambda_max, m_lambda_opt, m_crown, "direct", method, structure, loading, mesh_size, cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, val_geometric_criterium, val_estimator_criterium, geometric_criterium, deg_k, refinement_degree_ref, want_global_discretization_error, want_local_discretization_error, want_global_estimation, want_local_estimation, want_local_improvement, interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_R1, radius_R2, local_improvement, shape, k_min, k_max, k_opt, want_local_enrichment, nb_layers_nodes_enrichment, save_vtu, display_vtu, save_pvd, display_pvd, save_vtu_ref, display_vtu_ref, save_vtu_lambda, display_vtu_lambda, save_vtu_crown, display_vtu_crown );
//     if ( want_local_estimation and want_interest_quantity_only == 0 ) {
//         display_vtu_pvd( m_adjoint, m_local_ref, m_adjoint_lambda_min, m_adjoint_lambda_max, m_adjoint_lambda_opt, m_crown, "adjoint", method_adjoint, structure, loading, mesh_size, cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, val_geometric_criterium, val_estimator_criterium, geometric_criterium, deg_k, refinement_degree_ref, want_global_discretization_error_adjoint, want_local_discretization_error_adjoint, want_global_estimation, want_local_estimation, want_local_improvement, interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_R1, radius_R2, local_improvement, shape, k_min, k_max, k_opt, want_local_enrichment, nb_layers_nodes_enrichment, save_vtu_adjoint, display_vtu_adjoint, save_pvd_adjoint, display_adjoint_pvd, save_vtu_local_ref, display_vtu_local_ref, save_vtu_adjoint_lambda, display_vtu_adjoint_lambda, save_adjoint_crown_vtu, display_adjoint_crown_vtu );
//     }
    
}
