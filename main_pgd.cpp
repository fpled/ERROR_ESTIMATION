//
// C++ Implementation: main_pgd_cpp
//
// Description: Global/Goal-oriented error estimation methods for PGD
//
//
// Author: Pled Florent <florent.pled@univ-paris-est.fr>, (C) 2017
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "build/problem_space/all_in_one.h" // sert a forcer le logiciel scons a generer le repertoire build et ses codes sources .h et .cpp correspondant a la formulation
#include "build/problem_parameter/all_in_one.h" // sert a forcer le logiciel scons a generer le repertoire build et ses codes sources .h et .cpp correspondant a la formulation
#include "Mesh.h"
#include "Material_properties.h"
#include "Boundary_conditions.h"
#include "Display.h"
#include "CONNECTIVITY/Calcul_connectivity.h"
#include "VERIFICATION/Verification.h"
#include "DISCRETIZATION_ERROR/Calcul_discretization_error.h"
#include "EET/Construct_standard_force_fluxes_EET.h"
#include "EET/Construct_standard_force_fluxes_EET_PGD.h"
#include "EESPT/Construct_standard_force_fluxes_EESPT.h"
#include "EESPT/Construct_standard_force_fluxes_EESPT_PGD.h"
#include "ENHANCEMENT/Construct_enhanced_force_fluxes_EET_EESPT.h"
#include "CRITERIUM_ENHANCEMENT/Construct_criterium_enhancement.h"
#include "ECRE/Construct_K_hat.h"
#include "ECRE/Construct_F_hat.h"
#include "ECRE/Construct_dep_hat.h"
#include "ECRE/Calcul_error_estimate_prolongation_condition.h"
#include "ECRE/Calcul_error_estimate_prolongation_condition_PGD.h"
#include "SPET/Set_patch.h"
#include "SPET/Construct_K_hat.h"
#include "SPET/Construct_F_hat.h"
#include "SPET/Construct_F_hat_PGD.h"
#include "SPET/Construct_dep_hat.h"
#include "SPET/Calcul_error_estimate_partition_unity.h"
#include "SPET/Calcul_error_estimate_partition_unity_PGD.h"
#include "PGD/PGD.h"

#include "LMT/include/containers/gnuplot.h"
#include "LMT/include/containers/matlabplot.h"
#include "LMT/include/containers/matcholamd.h"
#include "LMT/include/containers/conjugate_gradient.h"
#include "LMT/include/containers/MatWithTinyBlocks.h"
#include "LMT/include/mesh/interpolation.h"

#include "LMT/include/util/MKL_direct_solver.h"
#include "LMT/include/util/MKL_iterative_solver.h"
#include "LMT/include/util/MUMPS_solver.h"

using namespace LMT;
using namespace std;

int main( int argc, char **argv ) {
    TicToc t_total;
    t_total.start();
    static const unsigned dim = 3;
    static const bool wont_add_nz = true;
    typedef Mesh<Mesh_carac_space<double,dim> > TM;
    typedef Formulation<TM,FormulationElasticity,DefaultBehavior,double,wont_add_nz> TF;
    typedef TM::Pvec Pvec;
    typedef TM::TNode::T T;
    typedef Mesh<Mesh_carac_parameter<double,1> > TM_param;
    typedef Formulation<TM_param,FormulationParam,DefaultBehavior,double,wont_add_nz> TF_param;
    typedef TM_param::Pvec Pvec_param;
    typedef TM::TElemListPtr TElemListPtr;
    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    static const string structure = "spherical_inclusions"; // structure
    // 2D : plate_traction, plate_flexion, plate_hole, plate_crack, structure_crack, test_specimen, weight_sensor, circular_inclusions, circular_holes,
    //      square_n (n=32,64,128,256,512,1024,2048,4096), square_init_n (n=32,64,128,256,512,1024,2048,4096)
    // 3D : beam_traction, beam_flexion, beam_hole, plate_hole, plate_hole_full, hub_rotor_helico, reactor_head, door_seal, spot_weld, blade, pipe, SAP, spherical_inclusions, spherical_holes,
    //      test_specimen_n (n=5,10,15,20,25,Q1_5,Q3_5,Q3_10,Q3_15,Q3_20,Q3_25,Q4_5,Q6_5,Q8_5), hashin_op_n (op=filtered,truncated,willot2015; n=32,64,128,256,512)
    static const string mesh_size = "fine"; // taille du maillage : coarse, fine
    // 2D : plate_hole, plate_crack, structure_crack, test_specimen, weigth_sensor
    // 3D : plate_hole, spot_weld, reactor_head
    static const string loading = "Step-2"; // chargement
    // spot_weld (3D) : pull, shear, peeling
    // plate_crack (2D) : pull, shear
    // test_specimen_n (3D) : Step-1, ..., Step-9
    static const unsigned deg_p = 1; // degre de l'analyse elements finis : 1, 2, ...
    static const string boundary_condition_D = "penalty"; // methode de prise en compte des conditions aux limites de Dirichlet (en deplacement) pour le pb direct : lagrange, penalty
    
    /// Global discretization error
    /// ---------------------------
    static const bool want_global_discretization_error = 0; // calcul de l'erreur de discretisation globale du pb direct
    static const bool want_local_discretization_error = 0; // calcul de l'erreur de discretisation locale du pb direct
    static const bool want_solve_ref = 0; // calcul d'une solution de reference sur un maillage de reference (tres fin)
    static const unsigned refinement_level_ref = 2; // degre du h-refinement pour la construction du maillage de reference du pb direct :
    // 1 -> sous-decoupage en 4/8 elements en 2D/3D
    // 2 -> sous-decoupage en 16/64 elements en 2D/3D
    // 3 -> sous-decoupage en 64/512 elements en 2D/3D
    // 4 -> sous-decoupage en 256/4096 elements en 2D/3D
    // 5 -> sous-decoupage en 1024/32768 elements en 2D/3D
    // 6 -> sous-decoupage en 4096/32768 elements en 2D/3D
    // 7 -> sous-decoupage en 16384/262144 elements en 2D/3D
    // 8 -> sous-decoupage en 65536/2097152 elements en 2D/3D
    
    /// Global error estimation method
    /// ------------------------------
    static const bool want_global_estimation = 1; // calcul d'un estimateur d'erreur globale (au sens de la norme energetique)
    static const string method = "EET"; //methode de construction de champs admissibles pour le pb direct : EET, SPET, EESPT
    static const string method_adjoint = "EET"; // methode de construction de champs admissibles pour le pb adjoint : EET, SPET, EESPT
    static const unsigned cost_function = 0; // fonction-cout pour les methodes EET, EESPT :
    // 0 : norme matricielle sans coefficient de ponderation (matrice identite)
    // 1 : norme matricielle avec coeff de ponderation (en 1/mes(face)^2)
    // 2 : norme energetique
    static const T penalty_val_N = 1e6; // coefficient de penalisation pour la prise en compte des conditions aux limites de Neumann (en effort) (methode EESPT)
    static const string solver = "LDL"; // solveur pour la resolution des pbs locaux avec blocage auto du noyau : CholMod (sym, def, pos), LDL (sym) // types de solveur sans blocage auto du noyau (-> ne marche pas!) : CholFactorize (sym, def, pos), LUFactorize, Inv, UMFPACK
    static const string solver_minimisation = "UMFPACK"; // solveur pour la resolution des pbs de minimisation : LDL (sym), UMFPACK, LUFactorize, Inv
    
    /// Enhanced technique
    /// ------------------
    static const bool enhancement_with_geometric_criterium = 0; // amelioration de la construction des densites d'effort (methodes EET, EESPT) basee sur un critere geometrique
    static const string geometric_criterium = "radius_ratio"; // critere d'amelioration geometrique : radius_ratio, edge_ratio
    static const T val_geometric_criterium = 0.34; // valeur du critere d'amelioration geometrique
    // critere radius_ratio : rapport entre rayon du cercle/sphere inscrit(e) et rayon du cercle/sphere circonscrit(e) à un élément
    // critere edge_ratio : rapport entre longueur/aire minimale et longueur/aire maximale des bords/faces d'un élément
    static const bool enhancement_with_estimator_criterium = 0; // amelioration de la construction des densites d'effort (methodes EET, EESPT) basee sur un critere sur l'estimateur d'erreur
    static const T val_estimator_criterium = 0.8; // valeur du critere d'amelioration sur l'estimateur d'erreur : rapport entre la contribution elementaire au carre a l'erreur estimee et la contribution elementaire maximale au carre
    
    /// Adaptive remeshing (mesh refinement)
    /// ------------------------------------
    static const bool want_remesh = 1; // remaillage adaptatif (raffinement du maillage)
    static const T tol_remesh = 1e-2; // tolerance pour le critère d'arrêt de l'algorithme de remaillage
    static const unsigned max_mesh = 6; // nb d'etapes de remaillage max
    static const unsigned max_iter_remesh = 10; // nb d'iterations max de l'algorithme de remaillage
    static const T k_remesh = 0.25; // rapport maximal entre la contribution élémentaire au carré à l'erreur estimée et la contribution élémentaire maximale au carré des barres qui ne seront pas divisées
    static const bool spread_cut_remesh = true; // propagation du raffinement au reste du maillage (étendue de la coupe si l'arête coupée n'est pas la plus longue de l'élément)
    
    /// Goal-oriented error estimation method
    /// -------------------------------------
    static const bool want_local_estimation = 0; // calcul de l'erreur locale sur une quantite d'interet
    static const bool want_introduction_sigma_hat_m = 1; // introduction de sigma_hat_m pour le calcul de l'erreur sur une quantite d'interet locale
    static const bool want_local_refinement = 1; // raffinement local du mailage adjoint
    static const bool want_local_enrichment = 0; // enrichissement local avec fonctions handbook
    static const bool want_local_improvement = 0; // amelioration des bornes pour le calcul de l'erreur locale sur une quantite d'interet
    static const bool want_solve_eig_local_improvement = 0; // resolution du pb aux valeurs propres generalisees pour calculer la constante dans l'amelioration des bornes d'erreur locale
    static const bool use_mask_eig_local_improvement = 0; // utilisation d'un masque (image) pour definir le maillage du pb aux valeurs propres generalisees
    static const bool want_solve_local_ref = 0; // calcul de la quantite d'interet (quasi-)exacte sur un maillage de reference
    static const string interest_quantity = "mean_sigma"; // quantite d'interet : mean_sigma, mean_epsilon, pointwise_dep, pointwise_sigma, pointwise_epsilon, SIF (stress intensity factor)
    static const string direction_extractor = "xx"; // direction de l'operateur d'extraction
    // quantites d'interet mean_sigma, mean_epsilon, pointwise_sigma, pointwise_epsilon : xx, yy, xy, zz, xz, yz
    // quantite d'interet pointwise_dep : x, y, z
    // quantite d'interet SIF (stress intensity factor) : I, II, III
    
    /// Zone of interest
    /// ----------------
    static const Vec<unsigned> elem_list_interest_quantity( 4886 ); // liste des elements definissant la zone d'interet (quantite d'interet mean_sigma, mean_epsilon)
    static const string pointwise_interest_quantity = "node"; // definition de la quantite d'interet ponctuelle : node, pos
    static const unsigned node_interest_quantity( 661 ); // noeud definissant la zone d'interet (quantite d'interet pointwise_dep, pointwise_sigma, pointwise_epsilon)
    static const Pvec pos_interest_quantity( 49.5, 135.5 ); // position definissant la zone d'interet (quantite d'interet pointwise_dep, pointwise_sigma, pointwise_epsilon)
    static const Pvec pos_crack_tip( 109., 105. ); // position de la pointe de fissure (quantite d'interet SIF) : ( 3.5, 0. ) pour plate_crack et ( 109., 105. ) pour structure_crack
    static const T angle_crack = atan2( -17, -3 ); // angle de la fissure (en rad) (quantite d'interet SIF) : 0. pour plate_crack et atan2( -17, -3 ) pour structure_crack
    static const T radius_Ri = 6; // rayon du cercle interieur a la couronne omega entourant la pointe de fissure (quantite d'interet SIF) : 1.6 pour plate_crack et 6 pour structure_crack
    static const T radius_Re = 8; // rayon du cercle exterieur a la couronne omega entourant la pointe de fissure (quantite d'interet SIF) : 3.4 pour plate_crack et 8 pour structure_crack
    
    /// Local refinement for adjoint problem
    /// ------------------------------------
    static const T l_min_refinement = 1.0; // longueur minimale des cotes des elements du maillage adjoint
    static const T k_refinement = 1.0; // coefficient d'augmentation de la longueur maximale des cotes des elements en fonction de la distance au point, au cercle, ... autour duquel on souhaite raffiner le maillage
    static const bool spread_cut = true; // propagation du raffinement au reste du maillage (étendue de la coupe si l'arête coupée n'est pas la plus longue de l'élément)
    
    /// Local enrichment with handbook functions
    /// ----------------------------------------
    static const unsigned nb_layers_nodes_enrichment = 2; // nombre de couches/rangées de noeuds enrichis par la PUM (pb direct)
    
    /// Improved goal-oriented error estimation methods based on Steklov/Rayleigh constants
    /// -----------------------------------------------------------------------------------
    static const string local_improvement = "rayleigh"; // amelioration du calcul de l'erreur locale sur une quantite d'interet, basee sur la constante de Steklov ou le quotient de Rayleigh : steklov, rayleigh
    static const string shape = "circle"; // forme geometrique des domaines homothetiques
    static const T k_min = 2.5; // parametre k_min du domaine homothetique (amelioration steklov) : facteur multiplicatif devant le rayon du cercle/sphere (shape circle/sphere)
    static const T k_max = 7.; // parametre k_max du domaines homothetique (amelioration steklov) : facteur multiplicatif devant le rayon du cercle/sphere (shape circle/sphere)
    static const T k_opt = 4.4; // parametre k_opt du domaine homothetique (amelioration rayleigh) : facteur multiplicatif devant le rayon du cercle/sphere (shape circle/sphere)
    static const string integration_k = "trapeze"; // type d'integration sur le parametre k (amelioration steklov) : gauss, trapeze, IPP
    static const unsigned integration_nb_points = 1000; // nb d'intervalles pour l'integration type trapeze sur le parametre k (amelioration steklov)
    
    /// Proper Generalized Decomposition - PGD
    /// --------------------------------------
    static const bool want_normalization = 1; // normalisation des fonctions en parametres
    static const unsigned max_mode = 6; // nb de modes max dans la decomposition
    static const unsigned max_iter = 4; // nb d'iterations max de l'algorithme de point fixe
    static const T tol_convergence_criterium_mode = 1e-4; // precision pour critere d'arret global (boucle sur les modes)
    static const T tol_convergence_criterium_iter = 1e-8; // precision pour critere d'arret local (processus iteratif)
    static const Vec<T,2> support_param( 1., 10. ); // support de l'espace des parametres
    static const unsigned nb_points_param = 100; // nb de points du maillage parametrique
    static const bool want_eval_PGD = 0; // evaluation de la decomposition PGD
    static const unsigned nb_vals = 3; // nb de valeurs des parametres pris aleatoirement pour l'evaluation de la decomposition PGD
    
    /// Verification equilibrium / solver
    /// ---------------------------------
    static const bool verif_eq = 0; // verification de l'equilibre global elements finis
    static const bool verif_compatibility_conditions = 1; // verification des conditions de compatibilite (equilibre elements finis) (methode EET)
    static const bool verif_eq_force_fluxes = 1; // verification de l'equilibre des densites d'effort (methodes EET, EESPT)
    static const T tol_compatibility_conditions = 1e-6; // tolerance pour la verification des conditions de compatibilite (equilibre elements finis) (methode EET)
    static const T tol_eq_force_fluxes = 1e-6; // tolerance pour la verification de l'equilibre des densites d'effort (methodes EET, EESPT)
    
    static const bool verif_solver = 1; // verification de la resolution des pbs locaux (methodes EET, SPET, EESPT)
    static const bool verif_solver_minimisation = 1; // verification de la resolution des pbs de minimisation (methodes EET, EESPT)
    static const bool verif_solver_enhancement = 1; // verification de la resolution des pbs locaux (amelioration des methodes EET, EESPT)
    static const bool verif_solver_minimisation_enhancement = 1; // verification de la resolution des pbs de minimisation (amelioration des methodes EET, EESPT)
    static const T tol_solver = 1e-6; // tolerance pour la verification de la resolution des pbs locaux (methodes EET, SPET, EESPT)
    static const T tol_solver_minimisation = 1e-6; // tolerance pour la verification de la resolution des pbs de minimisation (methodes EET, EESPT)
    static const T tol_solver_enhancement = 1e-6; // tolerance pour la verification de la resolution des pbs locaux (amelioration des methodes EET EESPT)
    static const T tol_solver_minimisation_enhancement = 1e-6; // tolerance pour la verification de la resolution des pbs de minimisation (amelioration des methodes EET, EESPT)
    
    /// Display outputs
    /// ---------------
    static const bool display_vtu = 0;
    static const bool display_pvd = 0;
    static const bool display_vtu_adjoint = 0;
    static const bool display_vtu_lambda = 0;
    static const bool display_vtu_adjoint_lambda = 0;
    static const bool display_vtu_crown = 0;
    
    static const bool display_pvd_space = 0;
    static const bool display_pvd_eval = 0;
    
    static const bool display_matlab = 0;
    
    /// -------------------------------------------- ///
    /// Construction de la solution PGD du pb direct ///
    /// -------------------------------------------- ///
    
    display_pb( dim, structure, deg_p  );
    
    /// Maillage en espace du pb direct
    /// -------------------------------
    TM m;
    set_mesh( m, structure, mesh_size, loading, deg_p, refinement_level_ref, want_global_discretization_error, want_local_discretization_error );
    string filename = set_filename( m, "direct", structure, loading, mesh_size, method, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, val_geometric_criterium, val_estimator_criterium, geometric_criterium, want_global_discretization_error, want_local_discretization_error, want_global_estimation, want_local_estimation, interest_quantity, direction_extractor, pointwise_interest_quantity, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_Ri, radius_Re, want_local_improvement, local_improvement, shape, k_min, k_max, k_opt, want_local_enrichment, nb_layers_nodes_enrichment );
    
    /// Formulation en espace du pb direct
    /// ----------------------------------
    TF f( m );
    
    /// Proprietes materiaux du pb direct
    /// ---------------------------------
    set_material_properties( f, m, structure );
    
    /// Conditions aux limites du pb direct
    /// -----------------------------------
    set_constraints( f, m, boundary_condition_D, "direct", structure, loading );
    set_load_conditions( m, structure, loading, mesh_size );
    
    /// Partition des elements du maillage en espace du pb direct
    /// ---------------------------------------------------------
    Vec< Vec<unsigned> > elem_group; // vecteur contenant une liste de pointeurs sur les elements definissant les zones avec parametres inconnus
    partition_elem_list( m, structure, elem_group );
    
    /// Maillage en parametre du pb direct
    /// ----------------------------------
    Vec<TM_param> m_param;
    m_param.resize( elem_group.size()-1 );
    for (unsigned p=0;p<elem_group.size()-1;++p)
        set_mesh_param( m_param[p], support_param[0], support_param[1], nb_points_param );
    
    /// Formulation en parametre du pb direct
    /// -------------------------------------
    Vec<TF_param> f_param;
    f_param.resize( elem_group.size()-1 );
    for (unsigned p=0;p<elem_group.size()-1;++p)
        f_param[p].set_mesh( &m_param[p] );
    
    /// Defintion des fonctions a variables separees
    /// --------------------------------------------
    Vec< Vec<T>, max_mode > dep_space;
    Vec< Vec< Vec<T>, max_mode > > dep_param;
    dep_param.resize( elem_group.size()-1 );
    Vec< Vec<T> > vals_param;
    vals_param.resize( elem_group.size()-1 );
    for (unsigned p=0;p<elem_group.size()-1;++p)
        vals_param[p] = generate( m_param[p].node_list, ExtractDMi<pos_DM>( 0 ) );
    
    Vec< DisplayParaview, max_mode > dp_space;
    Vec<string> lp_space;
    lp_space.push_back( "dep" );
    lp_space.push_back( "young_eff" );
    DisplayParaview dp;
    Vec<string> lp("all");
    
    /// Resolution du pb direct
    /// -----------------------
    TicToc t;
    t.start();
    
    /// Construction des opérateurs et du second membre en espace
    /// ---------------------------------------------------------
    Vec<T> F_space;
    Vec<TMatSymSparse> K_space;
    assemble_space( m, f, F_space, K_space, elem_group );
    
    /// Construction des opérateurs et seconds membres en parametre
    /// -----------------------------------------------------------
    Vec< Vec<T> > F_param;
    Vec< Vec<TMatSymSparse,2> > K_param;
    assemble_param( m_param, f_param, F_param, K_param, elem_group );
    
    /// Construction d'une solution elements finis en espace particuliere du pb direct
    /// ------------------------------------------------------------------------------
    f.solve();
    Vec<T> dep_space_FE_part = f.vectors[0];
    
    /// Verification de l'equilibre elements finis du pb direct
    /// -------------------------------------------------------
    if ( verif_eq )
        check_equilibrium( f, "direct" );
    
    /// Calcul de la norme du champ de deplacement approche du pb direct
    /// ----------------------------------------------------------------
    calcul_norm_dep( m, f, "direct" );
    
    // EET and EESPT methods
    Vec<bool> elem_flag_enh;
    Vec<bool> face_flag_enh;
    Vec<bool> elem_flag_bal;
    bool enhancement = 0;
    bool balancing = 0;
    
    // SPET method
    Vec<unsigned> connect_node_to_vertex_node;
    unsigned nb_vertex_nodes;
    Vec< Vec<unsigned> > elem_list_vertex_node;
    Vec< Vec<unsigned> > face_type;
    Vec<unsigned> nb_points_elem;
    Vec<unsigned> nb_points_patch;
    Vec< Vec< Vec<unsigned> > > patch_elem;
    Vec< Vec< Vec<unsigned> > > constrained_points_list_patch;
    Vec<unsigned> nb_constraints_patch;
    
    Vec< Mat<T, Sym<> > > K_hat;
    Vec< Vec<T> > F_hat;
    Vec< Vec<T> > dep_space_hat_part;
    
    if ( want_global_estimation or want_local_estimation ) {
        
        display_method( "direct", method, cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, solver, solver_minimisation );
        
        /// Construction d'un champ de contrainte admissible particulier
        /// ------------------------------------------------------------
        if ( method == "EET" or method == "EESPT" ) {
            Vec< Vec< Vec<T> > > force_fluxes_standard;
            if ( method == "EET" )
                construct_standard_force_fluxes_EET( m, f, "direct", cost_function, enhancement, face_flag_enh, solver_minimisation, force_fluxes_standard, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation, verif_compatibility_conditions, tol_compatibility_conditions );  
            else if ( method == "EESPT")
                construct_standard_force_fluxes_EESPT( m, f, "direct", cost_function, enhancement, face_flag_enh, solver_minimisation, penalty_val_N, force_fluxes_standard, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation );
            
            if ( verif_eq_force_fluxes )
                check_equilibrium_force_fluxes( m, f, "direct", force_fluxes_standard, tol_eq_force_fluxes, want_local_enrichment );
            
            construct_K_hat( m, f, K_hat );
            construct_F_hat( m, f, "direct", balancing, elem_flag_bal, elem_flag_enh, force_fluxes_standard, F_hat, want_local_enrichment );
            construct_dep_hat( m, f, solver, K_hat, F_hat, dep_space_hat_part, verif_solver, tol_solver );
        }
        else if ( method == "SPET" ) {
            set_patch( m, f, nb_vertex_nodes, face_type, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_elem, nb_points_patch, constrained_points_list_patch, nb_constraints_patch );
            
            construct_K_hat( m, f, nb_vertex_nodes, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, constrained_points_list_patch, K_hat );
            construct_F_hat( m, f, "direct", nb_vertex_nodes, face_type, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, constrained_points_list_patch, F_hat, want_local_enrichment );
            construct_dep_hat( m, f, solver, nb_vertex_nodes, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, nb_points_elem, K_hat, F_hat, dep_space_hat_part, verif_solver, tol_solver );
        }
    }
    
    Vec<T, max_mode+max_mesh > theta, theta_PGD, theta_dis;
    theta.set( 0. );
    theta_PGD.set( 0. );
    theta_dis.set( 0. );
    Vec<unsigned> modes;
    Vec< Vec<T>, max_mode+max_mesh > theta_elem, theta_elem_PGD, theta_elem_dis;
    Vec< Vec< Vec<T> >, max_mode+max_mesh > dep_space_hat, dep_space_FE;
    
    unsigned mode = 0;
    unsigned mesh = 0;
    while ( true ) {
        
        TicToc t_mode;
        t_mode.start();
        
        cout << "Mode #" << mode+1 << endl;
        cout << "-------" << endl << endl;
        
        /// Initialisation
        /// --------------
        unsigned iter = 0;
        
        for (unsigned p=0;p<elem_group.size()-1;++p) {
            dep_param[ p ][ mode ].resize( f_param[p].vectors[0].size() );
            dep_param[ p ][ mode ].set( 1. );
        }
        dep_space[ mode ].resize( f.vectors[0].size() );
        
        solve_space( m, f, mode, F_space, F_param, K_param, elem_group, dep_param, dep_space );
        
        string filename_mode;
        if ( want_remesh )
            filename_mode = filename + "_mode" + to_string(mode+1) + "_space_mesh" + to_string(mesh) + "_iter";
        else
            filename_mode = filename + "_mode" + to_string(mode+1) + "_space_iter";
        dp_space[ mode ].add_mesh_iter( m, filename_mode, lp_space, iter );
        
        /// Processus iteratif : Algorithme de point fixe
        /// ---------------------------------------------
        while ( true ) {
            ++iter;
            
            /// Construction et resolution des pbs en parametre
            /// -----------------------------------------------
            Vec< Vec<T> > dep_param_old;
            dep_param_old.resize( elem_group.size()-1 );
            for (unsigned p=0;p<elem_group.size()-1;++p) {
                dep_param_old[ p ] = dep_param[ p ][ mode ];
                solve_param( m_param[p], f_param[p], p, mode, F_space, F_param, K_space, K_param, elem_group, dep_space, dep_param, want_normalization );
            }
            
            /// Construction et resolution du pb en espace
            /// ------------------------------------------
            Vec<T> dep_space_old = dep_space[ mode ];
            solve_space( m, f, mode, F_space, F_param, K_param, elem_group, dep_param, dep_space );
            
//            cout << "Fonction en espace =" << endl;
//            cout << dep_space[ mode ] << endl << endl;
//            for (unsigned p=0;p<elem_group.size()-1;++p) {
//                cout << "Fonction en parametre " << p << " =" << endl;
//                cout << dep_param[ p ][ mode ] << endl << endl;
//            }
            
            dp_space[ mode ].add_mesh_iter( m, filename_mode, lp_space, iter );
            
            /// Stationnarite du produit mode en espace * modes en parametre dans le processus iteratif
            /// ---------------------------------------------------------------------------------------
            T stagnation_indicator = 0.;
            calc_stagnation_indicator( m, f, mode, K_param, elem_group, dep_space, dep_param, dep_space_old, dep_param_old, stagnation_indicator );
            
            cout << "Iteration #" << iter << " : stagnation = " << stagnation_indicator << endl;
            if ( iter >= max_iter ) { // if ( stagnation_indicator < tol_convergence_criterium_iter or iter >= max_iter )
                cout << endl;
                break;
            }
        }
        
        /// Residu (au sens faible) associe a la solution PGD
        /// -------------------------------------------------
        T residual = 0.;
        T error_indicator = 0.;
        calc_error_indicator( m, f, mode, F_space, F_param, K_param, elem_group, dep_space, dep_param, residual, error_indicator );
        
        cout << "Nb d'iterations = " << iter << " : residu = " << residual << ", erreur = " << error_indicator << endl << endl;
        
        t_mode.stop();
        cout << "temps de calcul du mode #" << mode+1 << " = " << t_mode.res << " s" << endl << endl;
        
        /// Estimation d'erreur globale
        /// ---------------------------
        if ( want_global_estimation or want_local_estimation ) {
            
            TicToc t_CRE;
            t_CRE.start();
            
            /// Construction d'un champ de contrainte admissible a zero
            /// -------------------------------------------------------
            reset_load_conditions( m );
            set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
            
            construct_dep_space_FE( m, F_param, K_param, elem_group, mode, dep_param, dep_space, dep_space_FE_part, dep_space_FE[ mode ] );
            
            if ( method == "EET" or method == "EESPT" ) {
                Vec< Vec< Vec<T> > > force_fluxes_standard;
                if ( method == "EET" )
                    construct_standard_force_fluxes_EET_PGD( m, f, "direct", cost_function, enhancement, face_flag_enh, solver_minimisation, force_fluxes_standard, dep_space_FE[ mode ], elem_group, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation, verif_compatibility_conditions, tol_compatibility_conditions );
                else if ( method == "EESPT" )
                    construct_standard_force_fluxes_EESPT_PGD( m, f, "direct", cost_function, enhancement, face_flag_enh, solver_minimisation, penalty_val_N, force_fluxes_standard, dep_space_FE[ mode ], elem_group, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation );
                
                if ( verif_eq_force_fluxes )
                    check_equilibrium_force_fluxes( m, f, "direct", force_fluxes_standard, tol_eq_force_fluxes, want_local_enrichment );
                
                construct_F_hat( m, f, "direct", balancing, elem_flag_bal, elem_flag_enh, force_fluxes_standard, F_hat, want_local_enrichment );
                construct_dep_hat( m, f, solver, K_hat, F_hat, dep_space_hat[ mode ], verif_solver, tol_solver );
            }
            else if ( method == "SPET" ) {
                construct_F_hat_PGD( m, f, "direct", nb_vertex_nodes, face_type, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, constrained_points_list_patch, F_hat, dep_space_FE[ mode ], elem_group, want_local_enrichment );
                construct_dep_hat( m, f, solver, nb_vertex_nodes, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, nb_points_elem, K_hat, F_hat, dep_space_hat[ mode ], verif_solver, tol_solver );
            }
            
            /// Calcul d'un estimateur d'erreur globale associe au pb direct
            /// ------------------------------------------------------------
            set_load_conditions( m, structure, loading, mesh_size );
            
            if ( method == "EET" or method == "EESPT" )
                calcul_error_estimate_prolongation_condition_PGD( m, f, "direct", theta[ mode+mesh ],  theta_PGD[ mode+mesh ], theta_dis[ mode+mesh ], theta_elem[ mode+mesh ], theta_elem_PGD[ mode+mesh ], theta_elem_dis[ mode+mesh ], dep_space, dep_param, dep_space_FE_part, dep_space_FE, dep_space_hat_part, dep_space_hat, elem_group, mode, want_global_discretization_error, want_local_discretization_error );
            else if ( method == "SPET" )
                calcul_error_estimate_partition_unity_PGD( m, f, "direct", theta[ mode+mesh ], theta_PGD[ mode+mesh ], theta_dis[ mode+mesh ], theta_elem[ mode+mesh ], theta_elem_PGD[ mode+mesh ], theta_elem_dis[ mode+mesh ], dep_space, dep_param, dep_space_FE_part, dep_space_FE, dep_space_hat_part, dep_space_hat, elem_group, mode, want_global_discretization_error, want_local_discretization_error );
            
            smoothing( m, ExtractDM< error_estimate_nodal_DM >(), ExtractDM< error_estimate_elem_DM >() );
            
            modes.push_back( mode+1 );
            
            t_CRE.stop();
            cout << "temps de calcul de la methode d'estimation d'erreur globale " << method << " = " << t_CRE.res << " s" << endl << endl;
            
        }
        
        /// Solution PGD
        /// ------------
        f.vectors[0].set( 0. );
        for (unsigned n=0;n<mode+1;++n) {
            Vec<T> dep_mode = dep_space[ n ];
            for (unsigned p=0;p<elem_group.size()-1;++p)
                dep_mode *= dep_param[ p ][ n ][ 0 ];
            f.vectors[0] +=  dep_mode;
        }
        f.update_variables();
        set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
        f.call_after_solve();
        
        if ( want_remesh )
            dp.add_mesh_iter( m, filename + "_solution_mode_" + to_string(mode+1) + "_mesh_" + to_string(mesh) + "_iter", lp, mode+1+mesh );
        else
            dp.add_mesh_iter( m, filename + "_solution_mode", lp, mode+1 );
        
        if ( want_remesh and theta_PGD[ mode+mesh ] < theta_dis[ mode+mesh ] ) {
            
            /// Adaptation du maillage associe au pb direct
            /// -------------------------------------------
            
            ++mesh;
            theta_PGD[ mode+mesh ] = theta_PGD[ mode+mesh-1 ];
            theta_dis[ mode+mesh ] = theta_dis[ mode+mesh-1 ];
            
            TM m_old = m;
            unsigned iter_remesh = 0;
            
            while ( theta_PGD[ mode+mesh ] < theta_dis[ mode+mesh ] and adapt_mesh( m, f, structure, method, ++iter_remesh, max_iter_remesh, k_remesh, spread_cut_remesh ) ) {
                
                f.set_mesh( &m );
                f.init();
                
                /// Proprietes materiaux du pb direct
                /// ---------------------------------
                //set_material_properties( f, m, structure );
                
                /// Conditions aux limites du pb direct
                /// -----------------------------------
                f.erase_constraints();
                set_constraints( f, m, boundary_condition_D, "direct", structure, loading, iter_remesh );
                reset_load_conditions( m );
                set_load_conditions( m, structure, loading, mesh_size );
                
                /// Partition des elements du maillage en espace du pb direct
                /// ---------------------------------------------------------
                partition_elem_list( m, structure, elem_group );
                
                /// Construction des opérateurs et du second membre en espace
                /// ---------------------------------------------------------
                assemble_space( m, f, F_space, K_space, elem_group );
                
                /// Construction d'une solution elements finis en espace particuliere du pb direct
                /// ------------------------------------------------------------------------------
                f.solve();
                dep_space_FE_part = f.vectors[0];
                
                /// Verification de l'equilibre elements finis du pb direct
                /// -------------------------------------------------------
                if ( verif_eq )
                    check_equilibrium( f, "direct" );
                
                /// Calcul de la norme du champ de deplacement approche du pb direct
                /// ----------------------------------------------------------------
                calcul_norm_dep( m, f, "direct" );
                
                /// Construction d'un champ de contrainte admissible particulier
                /// ------------------------------------------------------------
                if ( method == "EET" or method == "EESPT" ) {
                    Vec< Vec< Vec<T> > > force_fluxes_standard;
                    if ( method == "EET" )
                        construct_standard_force_fluxes_EET( m, f, "direct", cost_function, enhancement, face_flag_enh, solver_minimisation, force_fluxes_standard, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation, verif_compatibility_conditions, tol_compatibility_conditions );  
                    else if ( method == "EESPT")
                        construct_standard_force_fluxes_EESPT( m, f, "direct", cost_function, enhancement, face_flag_enh, solver_minimisation, penalty_val_N, force_fluxes_standard, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation );
                    
                    if ( verif_eq_force_fluxes )
                        check_equilibrium_force_fluxes( m, f, "direct", force_fluxes_standard, tol_eq_force_fluxes, want_local_enrichment );
                    
                    construct_K_hat( m, f, K_hat );
                    construct_F_hat( m, f, "direct", balancing, elem_flag_bal, elem_flag_enh, force_fluxes_standard, F_hat, want_local_enrichment );
                    construct_dep_hat( m, f, solver, K_hat, F_hat, dep_space_hat_part, verif_solver, tol_solver );
                }
                else if ( method == "SPET" ) {
                    set_patch( m, f, nb_vertex_nodes, face_type, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_elem, nb_points_patch, constrained_points_list_patch, nb_constraints_patch );
                    
                    construct_K_hat( m, f, nb_vertex_nodes, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, constrained_points_list_patch, K_hat );
                    construct_F_hat( m, f, "direct", nb_vertex_nodes, face_type, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, constrained_points_list_patch, F_hat, want_local_enrichment );
                    construct_dep_hat( m, f, solver, nb_vertex_nodes, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, nb_points_elem, K_hat, F_hat, dep_space_hat_part, verif_solver, tol_solver );
                }
                
                /// Projection des modes precedents sur le nouveau maillage
                /// -------------------------------------------------------
//                Mat<T, Gen<>, SparseLine<> > P;
//                interpolation_matrix( m_old, m, P );
//                for (unsigned n=0;n<mode;++n)
//                    dep_space[ n ]  = P*dep_space[ n ];
                
                /// Construction et resolution des pbs en espace associes aux differents modes
                /// --------------------------------------------------------------------------
                for (unsigned n=0;n<mode+1;++n)
                    solve_space( m, f, n, F_space, F_param, K_param, elem_group, dep_param, dep_space );
                
                /// Estimation d'erreur globale
                /// ---------------------------
                TicToc t_CRE;
                t_CRE.start();
                
                /// Construction des champs de contrainte admissibles a zero associes aux differents modes
                /// --------------------------------------------------------------------------------------
                reset_load_conditions( m );
                set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
                
                for (unsigned n=0;n<mode+1;++n) {
                    construct_dep_space_FE( m, F_param, K_param, elem_group, n, dep_param, dep_space, dep_space_FE_part, dep_space_FE[ n ] );
                    
                    if ( method == "EET" or method == "EESPT" ) {
                        Vec< Vec< Vec<T> > > force_fluxes_standard;
                        if ( method == "EET" )
                            construct_standard_force_fluxes_EET_PGD( m, f, "direct", cost_function, enhancement, face_flag_enh, solver_minimisation, force_fluxes_standard, dep_space_FE[ n ], elem_group, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation, verif_compatibility_conditions, tol_compatibility_conditions );
                        else if ( method == "EESPT" )
                            construct_standard_force_fluxes_EESPT_PGD( m, f, "direct", cost_function, enhancement, face_flag_enh, solver_minimisation, penalty_val_N, force_fluxes_standard, dep_space_FE[ n ], elem_group, want_local_enrichment, verif_solver_minimisation, tol_solver_minimisation );
                        
                        if ( verif_eq_force_fluxes )
                            check_equilibrium_force_fluxes( m, f, "direct", force_fluxes_standard, tol_eq_force_fluxes, want_local_enrichment );
                        
                        construct_F_hat( m, f, "direct", balancing, elem_flag_bal, elem_flag_enh, force_fluxes_standard, F_hat, want_local_enrichment );
                        construct_dep_hat( m, f, solver, K_hat, F_hat, dep_space_hat[ n ], verif_solver, tol_solver );
                    }
                    else if ( method == "SPET" ) {
                        construct_F_hat_PGD( m, f, "direct", nb_vertex_nodes, face_type, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, constrained_points_list_patch, F_hat, dep_space_FE[ n ], elem_group, want_local_enrichment );
                        construct_dep_hat( m, f, solver, nb_vertex_nodes, connect_node_to_vertex_node, elem_list_vertex_node, patch_elem, nb_points_patch, nb_points_elem, K_hat, F_hat, dep_space_hat[ n ], verif_solver, tol_solver );
                    }
                }
                
                /// Calcul d'un estimateur d'erreur globale associe au pb direct
                /// ------------------------------------------------------------
                set_load_conditions( m, structure, loading, mesh_size );
                
                if ( method == "EET" or method == "EESPT" )
                    calcul_error_estimate_prolongation_condition_PGD( m, f, "direct", theta[ mode+mesh ],  theta_PGD[ mode+mesh ], theta_dis[ mode+mesh ], theta_elem[ mode+mesh ], theta_elem_PGD[ mode+mesh ], theta_elem_dis[ mode+mesh ], dep_space, dep_param, dep_space_FE_part, dep_space_FE, dep_space_hat_part, dep_space_hat, elem_group, mode, want_global_discretization_error, want_local_discretization_error );
                else if ( method == "SPET" )
                    calcul_error_estimate_partition_unity_PGD( m, f, "direct", theta[ mode+mesh ], theta_PGD[ mode+mesh ], theta_dis[ mode+mesh ], theta_elem[ mode+mesh ], theta_elem_PGD[ mode+mesh ], theta_elem_dis[ mode+mesh ], dep_space, dep_param, dep_space_FE_part, dep_space_FE, dep_space_hat_part, dep_space_hat, elem_group, mode, want_global_discretization_error, want_local_discretization_error );
                
                smoothing( m, ExtractDM< error_estimate_nodal_DM >(), ExtractDM< error_estimate_elem_DM >() );
                
                t_CRE.stop();
                cout << "temps de calcul de la methode d'estimation d'erreur globale " << method << " = " << t_CRE.res << " s" << endl << endl;
            }
            
            modes.push_back( mode+1 );
            
            /// Modes PGD
            /// ---------
            for (unsigned n=0;n<mode+1;++n) {
                f.vectors[0] = dep_space[ n ];
                f.update_variables();
                
                filename_mode = filename + "_mode" + to_string(n+1) + "_space_mesh" + to_string(mesh) + "_iter";
                dp_space[ n ].add_mesh_iter( m, filename_mode, lp_space, iter+mesh );
            }
            
            /// Solution PGD
            /// ------------
            f.vectors[0].set( 0. );
            for (unsigned n=0;n<mode+1;++n) {
                Vec<T> dep_mode = dep_space[ n ];
                for (unsigned p=0;p<elem_group.size()-1;++p)
                    dep_mode *= dep_param[ p ][ n ][ 0 ];
                f.vectors[0] +=  dep_mode;
            }
            f.update_variables();
            set_field_alternativeontype( m, Number<1>(), 1., alpha_DM() );
            f.call_after_solve();
            
            if ( want_remesh )
                dp.add_mesh_iter( m, filename + "_solution_mode_" + to_string(mode+1) + "_mesh_" + to_string(mesh) + "_iter", lp, mode+1+mesh );
        }
        
        if ( mode+1 >= max_mode ) { // if ( error_indicator < tol_convergence_criterium_mode or mode+1 >= max_mode ) {
            cout << "Convergence de l'algorithme : nb de modes = " << mode+1 << ", residu = " << residual << ", erreur = " << error_indicator << endl << endl;
            break;
        }
        ++mode;
    }
    
    /// ---------------------- ///
    /// Sauvegarde / Affichage ///
    /// ---------------------- ///
    
    MatlabPlot mp(display_matlab);
    mp.cd_cwd();
    for (unsigned n=0;n<mode+1;++n) {
        // Paraview
        string filename_mode = filename + "_mode" + to_string(n+1) + "_space";
        if ( display_pvd_space )
            dp_space[ n ].exec( filename_mode );
        else
            dp_space[ n ].make_pvd_file( filename_mode );
        
        for (unsigned p=0;p<elem_group.size()-1;++p) {
            // Matlab
            string output = "'" + filename + "_mode" + to_string(n+1) + "_param" + to_string(p+1);
            string xlabel = "'$p_" + to_string(p+1) + "$'";
            string ylabel = "'$\\gamma_{" + to_string(p+1) + "," + to_string(n+1) + "}$'";
            string params = ",'LineStyle','-','Color',getfacecolor(" + to_string(p+4) + "),'LineWidth',1";
//            mp.save_plot( vals_param[p], dep_param[ p ][ n ], (output + ".fig'").c_str(), xlabel.c_str(), ylabel.c_str(), params.c_str() );
//            mp.save_plot( vals_param[p], dep_param[ p ][ n ], (output + ".epsc2'").c_str(), xlabel.c_str(), ylabel.c_str(), params.c_str() );
//            mp.save_plot( vals_param[p], dep_param[ p ][ n ], (output + ".tex'").c_str(), xlabel.c_str(), ylabel.c_str(), params.c_str() );
            mp.figure();
            mp.plot( vals_param[p], dep_param[ p ][ n ], params.c_str() );
            mp.grid_on();
            mp.box_on();
            mp.set_fontsize(16);
            mp.set_xlabel_interpreter(xlabel.c_str(),"'latex'");
            mp.set_ylabel_interpreter(ylabel.c_str(),"'latex'");
            mp.save_output((output + ".fig'").c_str());
            mp.save_output((output + ".epsc2'").c_str());
            mp.save_output((output + ".tex'").c_str());
            mp.close();
            
            // Gnuplot
//            output = "'gp_" + filename + "_mode" + to_string(n+1) + "_param" + to_string(p+1);
//            params = "notitle w l lt " + to_string(p+1) + " lw 1";
//            bool jump_lines = false;
////            save_plot( vals_param[p], dep_param[ p ][ n ], (output + ".tex'").c_str(), xlabel.c_str(), ylabel.c_str(), params.c_str() );
//            GnuPlot gp;
//            gp.set_output_terminal((output + ".tex'").c_str());
//            gp.set_xlabel(xlabel.c_str());
//            gp.set_ylabel(ylabel.c_str());
//            gp.plot( vals_param[p], dep_param[ p ][ n ], params.c_str(), jump_lines );
        }
    }
    
    if ( display_pvd )
        dp.exec( filename + "_solution" );
    else
        dp.make_pvd_file( filename + "_solution" );
    
    if ( want_global_estimation or want_local_estimation ) {
        Vec<T> theta_2, theta_2_PGD, theta_2_dis;
        theta_2.resize( modes.size(), 0. );
        theta_2_PGD.resize( modes.size(), 0. );
        theta_2_dis.resize( modes.size(), 0. );
        for (unsigned n=0;n<modes.size();++n) {
            theta_2[ n ] = pow( theta[ n ], 2 );
            theta_2_PGD[ n ] = pow( theta_PGD[ n ], 2 );
            theta_2_dis[ n ] = pow( theta_dis[ n ], 2 );
        }
        cout << modes << endl;
        cout << theta_2 << endl;
        cout << theta_2_PGD << endl;
        cout << theta_2_dis << endl;
        // Matlab
        string output = "'" + filename + "_estimates_global";
        string xlabel = "'$m$'";
        string ylabel = "'Error'";
        string legend = "'$E^2_{\\mathrm{CRE}}$','$\\eta^2_{\\mathrm{PGD}}$','$\\eta^2_{\\mathrm{dis}}$'";
        string params_CRE = ",'LineStyle','-','Color',getfacecolor(4),'LineWidth',1";
        string params_PGD = ",'LineStyle','-','Color',getfacecolor(5),'LineWidth',1";
        string params_dis = ",'LineStyle','-','Color',getfacecolor(6),'LineWidth',1";
        mp.figure();
        mp.semilogy( modes, theta_2, params_CRE.c_str() );
        mp.hold_on();
        mp.semilogy( modes, theta_2_PGD, params_PGD.c_str() );
        mp.semilogy( modes, theta_2_dis, params_dis.c_str() );
        mp.hold_off();
        mp.grid_on();
        mp.box_on();
        mp.set_fontsize(16);
        mp.set_xlabel_interpreter(xlabel.c_str(),"'latex'");
//        mp.set_ylabel_interpreter(ylabel.c_str(),"'latex'");
        mp.set_legend_interpreter(legend.c_str(),"'latex'");
        mp.save_output((output + ".fig'").c_str());
        mp.save_output((output + ".epsc2'").c_str());
        mp.save_output((output + ".tex'").c_str());
        mp.close();
    }
    
    if ( want_eval_PGD )
        eval_PGD( m_param, m, f, "direct", structure, boundary_condition_D, loading, mesh_size, elem_group, nb_vals, vals_param, mode, dep_space, dep_param, filename, display_pvd_eval );
    
    t.stop();
    cout << "temps de calcul de la resolution du pb direct = " << t.res << " s" << endl << endl;
    
    t_total.stop();
    cout << "temps de calcul total = " << t_total.res << " s" << endl << endl;
    
}
