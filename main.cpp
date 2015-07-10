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
#include "build/problem_error_estimation/all_in_one.h" // sert a forcer le logiciel scons a generer le repertoire build et ses codes sources .h et .cpp correspondant a la formulation
#include "Structure.h"
#include "Material_properties.h"
#include "Boundary_conditions.h"
#include "GEOMETRY/Calcul_geometry.h"
#include "GEOMETRY/Geometry.h"
#include "DISCRETIZATION_ERROR/Calcul_discretization_error.h"
#include "Display.h"
#include "Calcul_global_error_estimation.h"
#include "Calcul_goal_oriented_error_estimation.h"

using namespace LMT;
using namespace std;

int main( int argc, char **argv ) {
    TicToc t_total;
    t_total.start();
    static const unsigned dim = 2;
    static const bool wont_add_nz = true; // pour completer plus vite la matrice de rigidite si ajout de termes non nuls. 
                                          // si wont_add_nz == true et " -DWITH_CHOLMOD " dans les parametres de compilation, pas de solveur iteratif
    typedef Mesh<Mesh_carac_error_estimation<double,dim> > TM; // definit un synonyme du type Mesh<Mesh_carac_error_estimation<double,dim> >
    typedef Formulation<TM,FormulationElasticity,DefaultBehavior,double,wont_add_nz> TF; // definit un synonyme du type Formulation<TM,FormulationElasticity>
    typedef TM::Pvec Pvec; // definit un synonyme du type TM::Pvec servant a representer la position d'un point
    typedef TM::TNode::T T;
    static const string structure = "square_32"; // structure 2D : plate_traction, plate_flexion, plate_hole, plate_crack, structure_crack, eprouvette, weight_sensor, circle
                                                     // structure 3D : beam_traction, beam_flexion, beam_hole, plate_hole, plate_hole_full, hub_rotor_helico, reactor_head, door_seal, spot_weld, blade, pipe, SAP, spherical_inclusions, sphere, sphere_center, sphere_hollow
    static const string mesh_size = "fine"; // mesh_size pour les structures plate_hole (2D ou  3D), plate_crack, structure_crack, test_specimen, weigth_sensor, spot_weld (3D), reactor_head (3D) : coarse, fine
    static const string loading = "pre_epsilon"; // loading pour la structure spot_weld (3D) : pull, shear, peeling et pour la structure plate_crack (2D) : pull, shear
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
    static const bool want_global_estimation = 1; // calcul de l'erreur globale au sens de la norme energetique usuelle
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
    
    ///--------------------------------------------///
    /// Construction de la solution elements finis ///
    ///--------------------------------------------///
    
    /// Maillage du pb direct
    ///----------------------
    TM m; // declaration d'un maillage de type TM
    TM m_ref;
    create_structure( m, m_ref, "direct", structure, mesh_size, loading, deg_p, want_global_discretization_error, want_local_discretization_error, refinement_degree_ref, want_solve_ref );
    display_structure( m, m_ref, "direct", structure, deg_p, want_solve_ref );
    
    /// Formulation du pb direct
    ///-------------------------
    TF f( m ); // creation d'une formulation du type TF avec le maillage m
    TF f_ref( m_ref );
    
    /// Proprietes materiaux et Conditions aux limites du pb direct
    ///------------------------------------------------------------
    create_material_properties( f, m, structure, loading );
    create_boundary_conditions( f, m, boundary_condition_D, "direct", structure, loading, mesh_size );
    if ( want_solve_ref ) {
        create_material_properties( f_ref, m_ref, structure, loading );
        create_boundary_conditions( f_ref, m_ref, boundary_condition_D, "direct", structure, loading, mesh_size );
    }
    
    /// Verification des contraintes cinematiques
    ///------------------------------------------
    check_constraints( f, display_constraints );
    
    /// Resolution du pb direct
    ///------------------------
    if ( structure.find("square") != string::npos ) {

    }
    else {
        TicToc t;
        t.start();
        if ( want_iterative_solver == 0 )
            f.solve();
        else
            f.solve( criterium_iterative_solver );
        t.stop();
        cout << "Temps de calcul du pb direct : " << t.res << endl << endl;
//    }
    
    /// Verification de l'equilibre du pb direct
    ///-----------------------------------------
    check_equilibrium( f, "direct", verif_eq );
    
    if ( want_solve_ref ) {
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
    
    calcul_discretization_error( m, m_ref, f, f_ref, debug_discretization_error, want_global_discretization_error, want_local_discretization_error, want_solve_ref );
    
    T theta = 0.;
    Vec<T> theta_elem;
    Vec< Vec<T> > dep_hat;
    
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
    }
    
    t_total.stop();
    cout << "Temps de calcul total : " << t_total.res << endl << endl;
    
    ///-----------///
    /// Affichage ///
    ///-----------///
    
    display_vtu_pvd( m, m_ref, m_lambda_min, m_lambda_max, m_lambda_opt, m_crown, "direct", method, structure, loading, mesh_size, cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, val_geometric_criterium, val_estimator_criterium, geometric_criterium, deg_k, refinement_degree_ref, want_global_discretization_error, want_local_discretization_error, want_global_estimation, want_local_estimation, want_local_improvement, interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_R1, radius_R2, local_improvement, shape, k_min, k_max, k_opt, want_local_enrichment, nb_layers_nodes_enrichment, save_vtu, display_vtu, save_pvd, display_pvd, save_vtu_ref, display_vtu_ref, save_vtu_lambda, display_vtu_lambda, save_vtu_crown, display_vtu_crown );
    if ( want_local_estimation and want_interest_quantity_only == 0 ) {
        display_vtu_pvd( m_adjoint, m_local_ref, m_adjoint_lambda_min, m_adjoint_lambda_max, m_adjoint_lambda_opt, m_crown, "adjoint", method_adjoint, structure, loading, mesh_size, cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, val_geometric_criterium, val_estimator_criterium, geometric_criterium, deg_k, refinement_degree_ref, want_global_discretization_error_adjoint, want_local_discretization_error_adjoint, want_global_estimation, want_local_estimation, want_local_improvement, interest_quantity, direction_extractor, pointwise_interest_quantity, list_elems_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_R1, radius_R2, local_improvement, shape, k_min, k_max, k_opt, want_local_enrichment, nb_layers_nodes_enrichment, save_vtu_adjoint, display_vtu_adjoint, save_pvd_adjoint, display_adjoint_pvd, save_vtu_local_ref, display_vtu_local_ref, save_vtu_adjoint_lambda, display_vtu_adjoint_lambda, save_adjoint_crown_vtu, display_adjoint_crown_vtu );
    }
    
}
