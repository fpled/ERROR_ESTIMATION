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
    static const bool wont_add_nz = true;
    typedef Mesh<Mesh_carac_space<double,dim> > TM;
    typedef Formulation<TM,FormulationElasticity,DefaultBehavior,double,wont_add_nz> TF;
    typedef TM::Pvec Pvec;
    typedef TM::TNode::T T;
    typedef Mesh<Mesh_carac_parameter<double,1> > TM_param;
    typedef Formulation<TM_param,FormulationParam,DefaultBehavior,double,wont_add_nz> TF_param;
    typedef TM_param::Pvec Pvec_param;
    typedef TM::TElemListPtr TElemListPtr;
    static const string structure = "circular_inclusions"; // structure 2D : plate_traction, plate_flexion, plate_hole, plate_crack, structure_crack, eprouvette, weight_sensor, circle, circular_inclusions, circular_holes
                                                     // structure 3D : beam_traction, beam_flexion, beam_hole, plate_hole, plate_hole_full, hub_rotor_helico, reactor_head, door_seal, spot_weld, blade, pipe, SAP, spherical_inclusions, spherical_holes
    static const string mesh_size = "fine"; // mesh_size pour les structures plate_hole (2D ou 3D), plate_crack, structure_crack, test_specimen, weigth_sensor, spot_weld (3D), reactor_head (3D) : coarse, fine
    static const string loading = "pull"; // chargement pour la structure spot_weld (3D) : pull, shear, peeling et pour la structure plate_crack (2D) : pull, shear
    static const unsigned deg_p = 1; // degre de l'analyse elements finis : 1, 2, ...
    static const unsigned deg_k = 3; // degre supplementaire : 1, 2 , 3, ...
    static const string boundary_condition_D = "penalty"; // methode de prise en compte des conditions aux limites de Dirichlet (en deplacement) pour le pb direct : lagrange, penalty
    static const bool display_constraints = 0; // affichage des contraintes cinematiques
    
    /// Global discretization error
    ///----------------------------
    static const bool want_global_discretization_error = 0; // calcul de l'erreur de discretisation globale du pb direct
    static const bool want_local_discretization_error = 0; // calcul de l'erreur de discretisation locale du pb direct
    static const bool want_solve_ref = 0; // calcul d'une solution de reference sur un maillage de reference (tres fin)
    static const unsigned refinement_degree_ref = 2; // degre du h-refinement pour la construction du maillage de reference du pb de reference :
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
    static const bool want_global_estimation = 0; // calcul d'un estimateur d'erreur globale (au sens de la norme energetique)
    static const string method = "EET"; //methode de construction de champs admissibles pour le pb direct : EET, SPET, EESPT
    static const string method_adjoint = "EET"; // methode de construction de champs admissibles pour le pb adjoint : EET, SPET, EESPT

    static const unsigned cost_function = 0; // fonction-cout pour les methodes EET, EESPT :
                                             // 0 : norme matricielle sans coefficient de ponderation (matrice identite)
                                             // 1 : norme matricielle avec coeff de ponderation (en 1/mes(face)^2)
                                             // 2 : norme energetique
    static const T penalty_val_N = 1e6; // coefficient de penalisation pour la prise en compte des conditions aux limites de Neumann (en effort) (methode EESPT)
    static const string solver = "LDL"; // solveur pour la resolution des problemes locaux avec blocage auto du noyau : CholMod (sym, def, pos), LDL (sym) // types de solveur sans blocage auto du noyau (-> ne marche pas!) : CholFactorize (sym, def, pos), LUFactorize, Inv, UMFPACK
    static const string solver_minimisation = "UMFPACK"; // solveur pour la resolution des problemes de minimisation : LDL (sym), UMFPACK, LUFactorize, Inv
    
    /// Enhanced technique
    ///-------------------
    static const bool enhancement_with_geometric_criterium = 0; // amelioration de la construction des densites d'effort (methodes EET, EESPT) basee sur un critere geometrique
    static const string geometric_criterium = "radius_ratio"; // critere d'amelioration geometrique : radius_ratio, edge_ratio
    static const T val_geometric_criterium = 0.34; // valeur du critere d'amelioration geometrique pour le choix des elements dont les densites d'effort doivent etre ameliorees

    static const bool enhancement_with_estimator_criterium = 0; // amelioration de la construction des densites d'effort (methodes EET, EESPT) basee sur un critere sur l'estimateur d'erreur
    static const T val_estimator_criterium = 0.8; // valeur du critere d'amelioration sur l'estimateur d'erreur (% de la contribution maximale au carre a l'estimateur d'erreur) pour le choix des elements dont les densites d'effort doivent etre ameliorees
    
    /// Goal-oriented error estimation method
    ///--------------------------------------
    static const bool want_local_estimation = 0; // calcul de l'erreur locale sur une quantite d'interet
    static const bool want_interest_quantity_only = 0; // calcul de la quantite d'interet uniquement
    static const bool want_handbook_only = 0; // calcul de la solution handbook uniquement
    static const bool want_introduction_sigma_hat_m = 1; // introduction de sigma_hat_m pour le calcul de l'erreur sur une quantite d'interet locale
    static const bool want_local_refinement = 1; // raffinement local du mailage adjoint
    static const bool want_local_enrichment = 0; // enrichissement local avec fonctions handbook
    static const bool want_local_improvement = 0; // amelioration des bornes pour le calcul de l'erreur locale sur une quantite d'interet
    static const bool want_solve_eig_local_improvement = 0; // resolution du pb aux valeurs propres generalisees pour calculer la constante dans l'amelioration des bornes d'erreur locale
    static const bool use_mask_eig_local_improvement = 0; // utilisation d'un masque (image) pour definir le maillage du pb aux valeurs propres generalisees
    static const bool want_solve_local_ref = 0; // calcul de la quantite d'interet (quasi-)exacte sur un maillage de reference
    static const string interest_quantity = "mean_sigma"; // quantite d'interet : mean_sigma, mean_epsilon, pointwise_dep, pointwise_sigma, pointwise_epsilon, SIF (stress intensity factor)
    static const string direction_extractor = "xx"; // direction de l'operateur d'extraction pour quantite d'interet mean_sigma, mean_epsilon, pointwise_sigma, pointwise_epsilon : xx, yy, xy, zz, xz, yz
                                                   // direction de l'operateur d'extraction pour quantite d'interet pointwise_dep : x, y, z
                                                   // direction de l'operateur d'extraction pour quantite d'interet SIF (stress intensity factor) : I, II, III
    
    /// Zone of interest
    ///-----------------
    static const Vec<unsigned> elem_list_interest_quantity( 4886 ); // liste des elements definissant la zone d'interet (quantite d'interet mean_sigma, mean_epsilon)
    static const string pointwise_interest_quantity = "node"; // definition de la quantite d'interet ponctuelle : node, pos
    static const unsigned node_interest_quantity( 661 ); // noeud definissant la zone d'interet (quantite d'interet pointwise_dep, pointwise_sigma, pointwise_epsilon)
    static const Pvec pos_interest_quantity( 49.5, 135.5 ); // position definissant la zone d'interet (quantite d'interet pointwise_dep, pointwise_sigma, pointwise_epsilon)
    static const Pvec pos_crack_tip( 109, 105. ); // position de la pointe de fissure (quantite d'interet SIF) : ( 3.5, 0. ) pour plate_crack et ( 109., 105. ) pour structure_crack
    static const T angle_crack = atan2( -17, -3 ); // angle de la fissure (en rad) (quantite d'interet SIF) : 0. pour plate_crack et atan2( -17, -3 ) pour structure_crack
    static const T radius_Ri = 6; // rayon du cercle interieur a la couronne omega entourant la pointe de fissure (quantite d'interet SIF) : 1.6 pour plate_crack et 6 pour structure_crack
    static const T radius_Re = 8; // rayon du cercle exterieur a la couronne omega entourant la pointe de fissure (quantite d'interet SIF) : 3.4 pour plate_crack et 8 pour structure_crack
    
    /// Local refinement parameters for adjoint problem
    ///------------------------------------------------
    // Decoupe du cote d'un element (Bar) si sa longueur est superieure à d * k + l_min ou d est la distance entre le milieu du cote et le centre
    static const T l_min_refinement = 1.0; // longueur minimale des cotes des elements du maillage adjoint
    static const T k_refinement = 1.0; // coefficient d'augmentation de la longueur maximale des cotes des elements en fonction de la distance au point, au cercle, ... autour duquel on souhaite raffiner le maillage
    static const bool spread_cut = true; // propagation du raffinement au reste du maillage (étendue de la coupe si l'arête coupée n'est pas la plus longue de l'élément)
    
    /// Local enrichment with handbook functions
    ///-----------------------------------------
    static const unsigned nb_layers_nodes_enrichment = 2; // nombre de couches/rangées de noeuds enrichis par la PUM (pb direct)

    /// Improved goal-oriented error estimation methods based on Steklov/Rayleigh constants
    ///------------------------------------------------------------------------------------
    static const string local_improvement = "rayleigh"; // amelioration du calcul de l'erreur locale sur une quantite d'interet, basee sur la constante de Steklov ou le quotient de Rayleigh : steklov, rayleigh
    static const string shape = "circle"; // forme geometrique des domaines homothetiques
    static const T k_min = 2.5; // parametre k_min du domaine homothetique (amelioration steklov) : facteur multiplicatif devant le rayon du cercle/sphere (shape circle/sphere)
    static const T k_max = 7.; // parametre k_max du domaines homothetique (amelioration steklov) : facteur multiplicatif devant le rayon du cercle/sphere (shape circle/sphere)
    static const T k_opt = 4.4; // parametre k_opt du domaine homothetique (amelioration rayleigh) : facteur multiplicatif devant le rayon du cercle/sphere (shape circle/sphere)
    static const string integration_k = "trapeze"; // type d'integration sur le parametre k (amelioration steklov) : gauss, trapeze, IPP
    static const unsigned integration_nb_points = 100; // nb d'intervalles pour l'integration type trapeze sur le parametre k (amelioration steklov)
    
    /// Parameters PGD
    ///---------------
    static const bool want_PGD = 1; // methode de reduction de modele PGD
    static const bool want_normalization = 1; // normalisation
    static const unsigned max_mode = 5; // nb de modes max dans la decomposition
    static const unsigned max_iter = 5; // nb d'iterations max de l'algorithme de type point fixe
    static const T tol_convergence_criterium_mode = 1e-4; // precision pour critere d'arret global (boucle sur les modes)
    static const T tol_convergence_criterium_iter = 1e-8; // precision pour critere d'arret local (processus iteratif)
    static const T min_param = -2.; // valeur min sur l'intervalle des parametres
    static const T max_param = 2.; // valeur max sur l'intervalle des parametres
    static const unsigned nb_points_param = 100; // nb de points du maillage parametrique
    static const bool want_verif_PGD = 0; // verification de la decomposition PGD cinematique
    static const unsigned nb_vals_param_verif = 1; // nb de valeurs des parametres pris aleatoirement pour la verification de la decomposition PGD
    
    /// Parameters for iterative solver
    ///--------------------------------
    static const bool want_iterative_solver = 0; // solveur iteratif pour la resolution du pb direct
    static const bool want_iterative_solver_ref = 0; // solveur iteratif pour la resolution du pb de reference associe au pb direct
    static const bool want_iterative_solver_adjoint = 0; // solveur iteratif pour la resolution du pb adjoint
    static const T iterative_criterium = 1e-3; // tolerance du solveur iteratif pour la resolution du pb direct : residu en norme inf
    static const T iterative_criterium_ref = 1e-3; // tolerance du solveur iteratif pour la resolution du pb de reference associe au pb direct : residu en norme inf
    static const T iterative_criterium_adjoint = 1e-3; // tolerance du solveur itératif pour la resolution du pb adjoint : residu en norme inf
    
    /// Verification equilibre / solveur
    ///---------------------------------
    static const bool verif_eq = 1; // verification de l'equilibre global elements finis
    static const bool verif_compatibility_conditions = 1; // verification des conditions de compatibilite (equilibre elements finis) (methode EET)
    static const T tol_compatibility_conditions = 1e-6; // tolerance pour la verification des conditions de compatibilite (equilibre elements finis) (methode EET)
    static const bool verif_eq_force_fluxes = 1; // verification de l'equilibre des densites d'effort (methodes EET, EESPT)
    static const T tol_eq_force_fluxes = 1e-6; // tolerance pour la verification de l'equilibre des densites d'effort (methodes EET, EESPT)

    static const bool verif_solver = 1; // verification de la resolution des problemes locaux (methodes EET, SPET, EESPT)
    static const T tol_solver = 1e-6; // tolerance pour la verification de la resolution des problemes locaux (methodes EET, SPET, EESPT)
    static const bool verif_solver_enhancement = 0; // verification de la resolution des problemes locaux (amelioration des methodes EET, EESPT)
    static const T tol_solver_enhancement = 1e-6; // tolerance pour la verification de la resolution des problemes locaux (amelioration des methodes EET EESPT)
    static const bool verif_solver_minimisation = 1; // verification de la resolution des problemes de minimisation (methodes EET, EESPT)
    static const T tol_solver_minimisation = 1e-6; // tolerance pour la verification de la resolution des problemes de minimisation (methodes EET, EESPT)
    static const bool verif_solver_minimisation_enhancement = 0; // verification de la resolution des problemes de minimisation (amelioration des methodes EET, EESPT)
    static const T tol_solver_minimisation_enhancement = 1e-6; // tolerance pour la verification de la resolution des problemes de minimisation (amelioration des methodes EET, EESPT)
    
    /// Debug
    ///------
    static const bool debug_method = 0; // debug des methodes EET, SPET, EESPT (pb direct)
    static const bool debug_method_adjoint = 0; // debug des methodes EET, SPET, EESPT (pb adjoint)
    static const bool debug_method_enhancement = 0; // debug de l'amelioration des methodes EET, EESPT (pb direct)
    static const bool debug_method_enhancement_adjoint = 0; // debug de l'amelioration des methodes EET, EESPT (pb adjoint)
    static const bool debug_criterium_enhancement = 0; // debug du critere d'amelioration (pb direct)
    static const bool debug_criterium_enhancement_adjoint = 0; // debug du critere d'amelioration (pb adjoint)
    static const bool debug_geometry = 0; // debug de la geometrie (pb direct)
    static const bool debug_geometry_adjoint = 0; // debug de la geometrie (pb adjoint)
    static const bool debug_discretization_error = 0; // debug de l'erreur de discretisation (pb direct)
    static const bool debug_force_fluxes = 0; // debug des densites d'effort pour les methodes EET, EESPT (pb direct)
    static const bool debug_force_fluxes_adjoint = 0; // debug des densites d'effort pour les methodes EET, EESPT (pb adjoint)
    static const bool debug_force_fluxes_enhancement = 0; // debug de l'amelioration des densites d'effort pour les methodes EET, EESPT (pb direct)
    static const bool debug_force_fluxes_enhancement_adjoint = 0; // debug de l'amelioration des densites d'effort pour les methodes EET, EESPT (pb adjoint)
    static const bool debug_error_estimate = 0; // debug de l'estimateur d'erreur globale (pb direct)
    static const bool debug_error_estimate_adjoint = 0; // debug de l'estimateur d'erreur globale (pb adjoint)
    static const bool debug_local_effectivity_index = 0; // debug de l'indice d'efficacite local (pb direct)
    static const bool debug_local_effectivity_index_adjoint = 0; // debug de l'indice d'efficacite local (pb adjoint)
    
    /// Sauvegarde / Affichage
    ///-----------------------
    static const string pathname = "./RESULTS/";

    static const bool save_vtu = 1;
    static const bool display_vtu = 0;
    static const bool save_pvd = 0;
    static const bool display_pvd = 0;
    static const bool save_vtu_ref = 0;
    static const bool display_vtu_ref = 0;
    
    static const bool save_vtu_adjoint = 1;
    static const bool save_pvd_adjoint = 0;
    static const bool save_vtu_local_ref = 0;
    static const bool display_vtu_adjoint = 0;
    static const bool display_pvd_adjoint = 0;
    static const bool display_vtu_local_ref = 0;
    
    static const bool save_vtu_lambda = 1;
    static const bool save_vtu_adjoint_lambda = 1;
    static const bool display_vtu_lambda = 0;
    static const bool display_vtu_adjoint_lambda = 0;
    
    static const bool save_vtu_crown = 1;
    static const bool display_vtu_crown = 0;
    
    static const bool save_pvd_PGD_space = 1;
    static const bool save_pvd_PGD_param = 1;
    static const bool save_plot_PGD_param = 1;
    static const bool save_pvd_PGD_space_verif = 1;
    static const bool display_pvd_PGD_space = 0;
    static const bool display_pvd_PGD_param = 0;
    static const bool display_pvd_PGD_space_verif = 0;
    
    ///--------------------------------------------///
    /// Construction de la solution elements finis ///
    ///--------------------------------------------///
    
    /// Maillage en espace du pb direct
    ///--------------------------------
    TM m; // declaration d'un maillage de type TM
    TM m_ref;
    create_structure( m, m_ref, "direct", structure, mesh_size, loading, deg_p, want_global_discretization_error, want_local_discretization_error, refinement_degree_ref, want_solve_ref );
    display_structure( m, m_ref, "direct", structure, deg_p, want_solve_ref );
    
    /// Formulation en espace du pb direct
    ///-----------------------------------
    TF f( m ); // creation d'une formulation du type TF avec le maillage m
    TF f_ref( m_ref );
    
    /// Proprietes materiaux et Conditions aux limites du pb direct
    ///------------------------------------------------------------
    set_material_properties( f, m, structure );
    set_boundary_conditions( f, m, boundary_condition_D, "direct", structure, loading, mesh_size );
    Vec< Vec<unsigned> > elem_group; // vecteur contenant une liste de pointeurs sur les elements definissant les zones avec parametre inconnu
    partition_elem_list( m, structure, elem_group );
    if ( want_solve_ref ) {
        set_material_properties( f_ref, m_ref, structure );
        set_boundary_conditions( f_ref, m_ref, boundary_condition_D, "direct", structure, loading, mesh_size );
    }

    /// Maillage en parametre du pb direct
    ///-----------------------------------
    TM_param m_param;
    create_structure_param( m_param, min_param, max_param, nb_points_param );

    /// Formulation en parametre du pb direct
    ///--------------------------------------
    TF_param f_param( m_param );
    
    /// Defintion des fonctions a variables separees
    ///---------------------------------------------
    unsigned nb_modes; // nb de modes dans la decomposition
    Vec< Vec<T> > dep_space;
    Vec< Vec< Vec<T> > > dep_param;
    Vec< Vec< Vec<T>, max_iter+1 >, max_mode > psi;
    Vec< Vec< Vec< Vec<T>, max_iter+1 >, max_mode > > lambda;
    lambda.resize( elem_group.size()-1 );
    Vec<unsigned> nb_iterations; // nb d'iterations de l'algorithme de type point fixe pour chaque mode
    nb_iterations.resize( max_mode );
    nb_iterations.set( 1 );
    Vec<T> vals_param;
    for (unsigned i=0;i<m_param.node_list.size();++i)
        vals_param.push_back( m_param.node_list[i].pos[0] );
    Vec< DisplayParaview, max_mode > dp_space;
    Vec< Vec< DisplayParaview, max_mode > > dp_param;
    dp_param.resize( elem_group.size()-1 );
    Vec< DisplayParaview, nb_vals_param_verif > dp_space_verif;
    Vec<string> lp_space, lp_param;
    lp_space.push_back( "dep" );
    lp_space.push_back( "young_eff" );
    lp_param.push_back( "dep" );
    string prefix = define_prefix( m, pathname, "direct", structure, loading, mesh_size );
    
    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    Vec<TMatSymSparse> K_space;
    K_space.resize(elem_group.size());
    Vec< Vec<TMatSymSparse,2> > K_param;
    K_param.resize(elem_group.size()-1);
    Vec<T> F_space;
    Vec< Vec<T> > F_param;
    F_param.resize(elem_group.size()-1);
    
    /// Verification des conditions cinematiques
    ///-----------------------------------------
    if ( display_constraints )
        check_constraints( f );
    
    /// Resolution du pb direct
    ///------------------------
    TicToc t;
    t.start();
    if ( want_PGD == 0 ) {
        if ( want_iterative_solver == 0 )
            f.solve();
        else
            f.solve( iterative_criterium );
    }
    else {
        f.allocate_matrices();
        f.shift();
        f.assemble( false, true );
        F_space = f.sollicitation;
        for (unsigned g=0;g<elem_group.size();++g) {
            for (unsigned j=0;j<m.elem_list.size();++j) {
                if ( find( elem_group[g], _1 == j ) )
                    m.elem_list[j]->set_field( "alpha", 1. );
                else
                    m.elem_list[j]->set_field( "alpha", 0. );
            }
            f.assemble( true, false );
            K_space[g] = f.matrices(Number<0>());
        }

        f_param.allocate_matrices();
        f_param.shift();
        f_param.assemble( false, true );
        for (unsigned p=0;p<elem_group.size()-1;++p)
            F_param[p] = f_param.sollicitation;
        m_param.phi = 0;
        f_param.assemble( true, false );
        for (unsigned p=0;p<elem_group.size()-1;++p)
            K_param[p][0] = f_param.matrices(Number<0>());
        m_param.phi = 1;
        f_param.assemble( true, false );
        for (unsigned p=0;p<elem_group.size()-1;++p)
            K_param[p][1] = f_param.matrices(Number<0>());

        unsigned n = 0;
        while ( true ) {
            cout << "Mode " << n << endl << endl;

            /// Initialisation
            ///---------------
            unsigned k = 0;

            for (unsigned p=0;p<elem_group.size()-1;++p) {
                lambda[ p ][ n ][ k ].resize( f_param.vectors[0].size() );
                lambda[ p ][ n ][ k ].set( 1. );
            }
            psi[ n ][ k ].resize( f.vectors[0].size() );

            solve_space( m, f, n, k, nb_iterations, F_space, F_param, K_param, elem_group, lambda, psi, want_iterative_solver, iterative_criterium );

            if ( save_pvd_PGD_space or display_pvd_PGD_space )
                dp_space[ n ].add_mesh_iter( m, prefix + "_space_mode_" + to_string( n ) + "_iter", lp_space, k );
            if ( save_pvd_PGD_param or display_pvd_PGD_param ) {
                for (unsigned p=0;p<elem_group.size()-1;++p)
                    dp_param[ p ][ n ].add_mesh_iter( m_param, prefix + "_param_" + to_string( p ) + "_mode_" + to_string( n ) + "_iter", lp_param, k );
            }

            while ( true ) {
                k++;
                for (unsigned p=0;p<elem_group.size()-1;++p) {
                    lambda[ p ][ n ][ k ].resize( f_param.vectors[0].size() );
                    lambda[ p ][ n ][ k ] = lambda[ p ][ n ][ k-1 ];
                }
                psi[ n ][ k ].resize( f.vectors[0].size() );
                psi[ n ][ k ] = psi[ n ][ k-1 ];
                
                /// Construction et resolution des pbs en parametre
                ///------------------------------------------------
                for (unsigned p=0;p<elem_group.size()-1;++p)
                    solve_param( m_param, f_param, p, n, k, nb_iterations, F_space, F_param, K_space, K_param, elem_group, psi, lambda );
                
                /// Construction et resolution du pb en espace
                ///-------------------------------------------
                solve_space( m, f, n, k, nb_iterations, F_space, F_param, K_param, elem_group, lambda, psi, want_iterative_solver, iterative_criterium );
                
                /// Normalisation de la fonction psi en espace
                ///-------------------------------------------
//                if ( want_normalization ) {
//                    psi[ n ][ k ] /= sqrt( dot( psi[ n ][ k ], K_s * psi[ n ][ k ] ) );
//                    f.vectors[ 0 ] = psi[ n ][ k ];
//                    f.update_variables();
//                    f.call_after_solve();
//                }
                
//                cout << "Fonction en espace =" << endl;
//                cout << psi[ n ][ k ] << endl << endl;
//                for (unsigned p=0;p<elem_group.size()-1;++p) {
//                    cout << "Fonction en parametre " << p << " =" << endl;
//                    cout << lambda[ p ][ n ][ k ] << endl << endl;
//                }
                
                if ( save_pvd_PGD_space or display_pvd_PGD_space )
                    dp_space[ n ].add_mesh_iter( m, prefix + "_space_mode_" + to_string( n ) + "_iter", lp_space, k );
                if ( save_pvd_PGD_param or display_pvd_PGD_param ) {
                    for (unsigned p=0;p<elem_group.size()-1;++p)
                        dp_param[ p ][ n ].add_mesh_iter( m_param, prefix + "_param_" + to_string( p ) + "_mode_" + to_string( n ) + "_iter", lp_param, k );
                }
                
                if ( save_plot_PGD_param ) {
                    for (unsigned p=0;p<elem_group.size()-1;++p) {
                        static GnuPlot gp;
                        gp.set_terminal( "postscript eps enhanced color" );
                        stringstream graph_name;
                        graph_name << "'" << prefix + "_param_" + to_string( p ) + "_mode_" + to_string( n ) + "_iter_" + to_string( k ) << ".eps'";
                        gp.set_output(graph_name.str().c_str());
                        gp.set("key left bottom");
                        gp.set("format '%g'");
                        gp.set("font \"Times-Roman\"");
                        gp.set_xlabel("'{/Symbol q}'");
                        gp.set_ylabel("'{/Symbol l}'");
                        gp.hold_on();
                        gp.plot( vals_param, lambda[ p ][ n ][ k ], "title '{/Symbol l}({/Symbol q})' axes x1y1 w l lt 1 lw 1" );
                        gp.hold_off();
                    }
                }
                
                /// Stationnarite du produit fonction psi en espace * fonction lambda en parametre dans le processus iteratif associe au mode n
                ///----------------------------------------------------------------------------------------------------------------------------
                T stagnation_indicator = 0.;
                for (unsigned g=0;g<elem_group.size();++g) {
                    T alpha = 1.;
                    for (unsigned p=0;p<elem_group.size()-1;++p) {
                        if ( g == p )
                            alpha *= dot( lambda[ p ][ n ][ k ], K_param[p][1] * lambda[ p ][ n ][ k ] );
                        else
                            alpha *= dot( lambda[ p ][ n ][ k ], K_param[p][0] * lambda[ p ][ n ][ k ] );
                    }
                    for (unsigned j=0;j<m.elem_list.size();++j) {
                        if ( find( elem_group[g], _1 == j ) )
                            m.elem_list[j]->set_field( "alpha", alpha );
                    }
                }
                f.assemble( true, false );
                stagnation_indicator += dot( psi[ n ][ k ], f.matrices(Number<0>()) * psi[ n ][ k ] );
                for (unsigned g=0;g<elem_group.size();++g) {
                    T alpha = 1.;
                    for (unsigned p=0;p<elem_group.size()-1;++p) {
                        if ( g == p )
                            alpha *= dot( lambda[ p ][ n ][ k-1 ], K_param[p][1] * lambda[ p ][ n ][ k-1 ] );
                        else
                            alpha *= dot( lambda[ p ][ n ][ k-1 ], K_param[p][0] * lambda[ p ][ n ][ k-1 ] );
                    }
                    for (unsigned j=0;j<m.elem_list.size();++j) {
                        if ( find( elem_group[g], _1 == j ) )
                            m.elem_list[j]->set_field( "alpha", alpha );
                    }
                }
                f.assemble( true, false );
                stagnation_indicator += dot( psi[ n ][ k-1 ], f.matrices(Number<0>()) * psi[ n ][ k-1 ] );
                for (unsigned g=0;g<elem_group.size();++g) {
                    T alpha = 1.;
                    for (unsigned p=0;p<elem_group.size()-1;++p) {
                        if ( g == p )
                            alpha *= dot( lambda[ p ][ n ][ k ], K_param[p][1] * lambda[ p ][ n ][ k-1 ] );
                        else
                            alpha *= dot( lambda[ p ][ n ][ k ], K_param[p][0] * lambda[ p ][ n ][ k-1 ] );
                    }
                    for (unsigned j=0;j<m.elem_list.size();++j) {
                        if ( find( elem_group[g], _1 == j ) )
                            m.elem_list[j]->set_field( "alpha", alpha );
                    }
                }
                f.assemble( true, false );
                stagnation_indicator -= 2 * dot( psi[ n ][ k ], f.matrices(Number<0>()) * psi[ n ][ k-1 ] );
                stagnation_indicator = sqrt( fabs( stagnation_indicator ) );

                cout << "Iteration " << k << " : stagnation = " << stagnation_indicator << endl << endl;
                if ( k >= max_iter ) { //if ( stagnation_indicator < tol_convergence_criterium_iter or k >= max_iter ) {
                    nb_iterations[n] = k;
                    cout << "nb d'iterations associe au mode " << n << " = " << nb_iterations[n] << endl << endl;
                    break;
                }
            }
            
            if ( display_pvd_PGD_space )
                dp_space[ n ].exec( prefix + "_space_mode_" + to_string( n ) );
            else if ( save_pvd_PGD_space )
                dp_space[ n ].make_pvd_file( prefix + "_space_mode_" + to_string( n ) );
            if ( display_pvd_PGD_param ) {
                for (unsigned p=0;p<elem_group.size()-1;++p)
                    dp_param[ p ][ n ].exec( prefix + "_param_" + to_string( p )+ "_mode_" + to_string( n ) );
            }
            else if ( save_pvd_PGD_param ) {
                for (unsigned p=0;p<elem_group.size()-1;++p)
                    dp_param[ p ][ n ].make_pvd_file( prefix + "_param_" + to_string( p ) + "_mode_" + to_string( n ) );
            }
            
            /// Residu au sens faible associe a la solution au mode n
            ///------------------------------------------------------
            T residual = 0.;
            T sollicitation = 0.;
            for (unsigned i=0;i<n+1;++i) {
                T sollicitation_mode = dot( F_space, psi[ i ][ nb_iterations[i] ] );
                for (unsigned p=0;p<elem_group.size()-1;++p)
                    sollicitation_mode *= dot( F_param[p], lambda[ p ][ i ][ nb_iterations[i] ] );
                sollicitation += sollicitation_mode;
                residual -= sollicitation_mode;
                for (unsigned j=0;j<n+1;++j) {
                    for (unsigned g=0;g<elem_group.size();++g) {
                        T alpha = 1.;
                        for (unsigned p=0;p<elem_group.size()-1;++p) {
                            if ( g == p )
                                alpha *= dot( lambda[ p ][ i ][ nb_iterations[i] ], K_param[p][1] * lambda[ p ][ j ][ nb_iterations[j] ] );
                            else
                                alpha *= dot( lambda[ p ][ i ][ nb_iterations[i] ], K_param[p][0] * lambda[ p ][ j ][ nb_iterations[j] ] );
                        }
                        for (unsigned l=0;l<m.elem_list.size();++l) {
                            if ( find( elem_group[g], _1 == l ) )
                                m.elem_list[l]->set_field( "alpha", alpha );
                        }
                    }
                    f.assemble( true, false );
                    residual += dot( psi[ i ][ nb_iterations[i] ], f.matrices(Number<0>()) * psi[ j ][ nb_iterations[j] ] );
                }
            }
            
            T error_indicator = fabs( residual ) / fabs( sollicitation );
            cout << "Mode " << n << " : residu = " << residual << ", erreur = " << error_indicator << endl << endl;
            if ( n >= max_mode-1 ) { //if ( error_indicator < tol_convergence_criterium_mode or n >= max_mode-1 ) {
                nb_modes = n+1;
                cout << "nb de modes = " << nb_modes << endl << endl;
                break;
            }
            n++;
        }

        dep_space.resize( nb_modes );
        dep_param.resize( nb_modes );
        for (unsigned n=0;n<nb_modes;++n) {
            dep_space[ n ] = psi[ n ][ nb_iterations[n] ];
            dep_param[ n ] = lambda[ n ][ nb_iterations[n] ];
        }
        
        if ( want_verif_PGD ) {
            /// Verification pour un jeu connu de parametres
            ///---------------------------------------------
//            for (unsigned p=0;p<nb_vals_param_verif;++p) {
//                unsigned ind_in_vals_param = rand() % nb_points_param;
//                for (unsigned j=0;j<m.elem_list.size();++j) {
//                    if ( find( elem_group, _1 == j ) )
//                        m.elem_list[j]->set_field( "alpha", 1. + vals_param[ ind_in_vals_param ] );
//                    else
//                        m.elem_list[j]->set_field( "alpha", 1. );
//                }
//                if ( want_iterative_solver == 0 )
//                    f.solve();
//                else
//                    f.solve( iterative_criterium );
//                if ( save_pvd_PGD_space_verif or display_pvd_PGD_space_verif )
//                    dp_space_verif[ p ].add_mesh_iter( m, prefix + "_space_verif_param_REF_" + to_string( vals_param[ ind_in_vals_param ] ), lp_space, 0 );
//                f.vectors[0].set( 0. );
//                for (unsigned n=0;n<nb_modes;++n) {
//                    f.vectors[0] += dep_space[ n ] * dep_param[ n ][ ind_in_vals_param ];
//                }
//                f.update_variables();
//                f.call_after_solve();
//                if ( save_pvd_PGD_space_verif or display_pvd_PGD_space_verif )
//                    dp_space_verif[ p ].add_mesh_iter( m, prefix + "_space_verif_param_PGD_" + to_string( vals_param[ ind_in_vals_param ] ), lp_space, 1 );
                
//                if ( display_pvd_PGD_space_verif )
//                    dp_space_verif[ p ].exec( prefix + "_space_verif_param_" + to_string( vals_param[ ind_in_vals_param ] ) );
//                else if ( save_pvd_PGD_space_verif )
//                    dp_space_verif[ p ].make_pvd_file( prefix + "_space_verif_param_" + to_string( vals_param[ ind_in_vals_param ] ) );
//            }
        }

    }
    t.stop();
    cout << "Temps de calcul du pb direct : " << t.res << endl << endl;
    
    Vec<T> dep_part;
    Vec< Vec<T> > kappa;
    kappa.resize( elem_group.size()-1 );
    for (unsigned p=0;p<elem_group.size()-1;++p) {
        kappa[p].resize( f_param.vectors[0].size() );
        kappa[p].set( 1. );
    }
    
    /// Construction d'un champ admissible
    ///-----------------------------------
    if ( want_PGD ) {
        for (unsigned i=0;i<m.elem_list.size();++i)
            m.elem_list[i]->set_field( "alpha", 1. );
        f.assemble();
        if ( want_iterative_solver == 0 )
            f.solve_system();
        else
            f.solve_system( iterative_criterium );
        f.update_variables();
        f.call_after_solve();
        
        dep_part = f.vectors[0];
    }
    
    /// Verification de l'equilibre du pb direct
    ///-----------------------------------------
    if ( verif_eq )
        check_equilibrium( f, "direct" );
    
    if ( want_solve_ref and want_PGD == 0 ) {
        /// Resolution du pb de reference associe au pb direct
        ///---------------------------------------------------
        TicToc t_ref;
        t_ref.start();
        if ( want_iterative_solver_ref == 0 )
            f_ref.solve();
        else
            f_ref.solve( iterative_criterium_ref );
        t_ref.stop();
        cout << "Temps de calcul du pb de reference associe au pb direct : " << t_ref.res << endl << endl;
        
        /// Verification de l'equilibre du pb de reference associe au pb direct
        ///--------------------------------------------------------------------
        if ( verif_eq )
            check_equilibrium( f_ref, "de reference associe au pb direct" );
    }
    
    ///---------------------------------------------------------------///
    /// Calcul et Affichage des informations relatives a la geometrie ///
    ///---------------------------------------------------------------///
    
    calcul_display_geometry( m, f, debug_geometry );
    
    ///---------------------------------------------------------------------------///
    /// Mesure de l'erreur de discretisation globale et locale associee pb direct ///
    ///---------------------------------------------------------------------------///
    
    if ( want_PGD == 0 )
        calcul_discretization_error( m, m_ref, f, f_ref, want_global_discretization_error, want_local_discretization_error, want_solve_ref, debug_discretization_error );
    
    T theta = 0.;
    Vec<T> theta_elem;
    Vec< Vec<T> > dep_part_hat;
    Vec< Vec< Vec<T> > > dep_space_hat;
    dep_space_hat.resize( nb_modes );

//    if ( want_global_estimation or ( want_local_estimation and want_handbook_only == 0 and want_interest_quantity_only == 0 ) ) {
        
//        ///---------------------------------------------------------------------------------------------------------------///
//        /// Construction d'un champ de contrainte admissible et Calcul d'un estimateur d'erreur globale associe pb direct ///
//        ///---------------------------------------------------------------------------------------------------------------///
        
//        calcul_global_error_estimation( f, m, "direct", method, cost_function, penalty_val_N, solver, solver_minimisation, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, geometric_criterium, val_geometric_criterium, val_estimator_criterium, theta, theta_elem, dep_hat, verif_compatibility_conditions, tol_compatibility_conditions, verif_eq_force_fluxes, tol_eq_force_fluxes, verif_solver, tol_solver, verif_solver_enhancement, tol_solver_enhancement, verif_solver_minimisation, tol_solver_minimisation, verif_solver_minimisation_enhancement, tol_solver_minimisation_enhancement, want_global_discretization_error, want_local_discretization_error, want_local_enrichment, debug_geometry, debug_force_fluxes, debug_force_fluxes_enhancement, debug_criterium_enhancement, debug_error_estimate, debug_local_effectivity_index, debug_method, debug_method_enhancement );
        
//    }
    
//    TM m_adjoint, m_local_ref, m_lambda_min, m_lambda_max, m_lambda_opt, m_adjoint_lambda_min, m_adjoint_lambda_max, m_adjoint_lambda_opt, m_crown;
    
//    if ( want_local_estimation ) {
        
//        ///---------------------------------------///
//        /// Construction de la quantite d'interet ///
//        ///---------------------------------------///
        
//        display_interest_quantity( interest_quantity, direction_extractor, pointwise_interest_quantity, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_Ri, radius_Re );
        
//        /// Definition de l'extracteur
//        ///---------------------------
//        if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" )
//            create_structure_crown( m, m_crown, pos_crack_tip, radius_Ri, radius_Re, spread_cut );
//        TF f_crown( m_crown );
//        define_extractor( m, m_crown, f, f_crown, interest_quantity, direction_extractor, pointwise_interest_quantity, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_Ri, radius_Re, want_local_enrichment );
        
//        ///------------------------------------------------------///
//        /// Calcul de la quantite d'interet locale approchee I_h ///
//        ///------------------------------------------------------///
        
//        T I_h;
//        calcul_interest_quantity( m, m_crown, f, f_crown, "direct", interest_quantity, direction_extractor, pointwise_interest_quantity, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_Ri, radius_Re, I_h );
        
//        ///------------------------------------------------------------///
//        /// Calcul de la quantite d'interet locale (quasi-)exacte I_ex ///
//        ///------------------------------------------------------------///
        
//        T I_ex;
//        if ( want_solve_local_ref ) {
//            create_structure( m_local_ref, m_local_ref, "direct", structure, mesh_size, loading, deg_p );
            
//            Vec<unsigned> elem_list_local_ref_interest_quantity;
//            unsigned node_local_ref_interest_quantity;
//            create_structure_local_ref( m, m_local_ref, deg_p, refinement_degree_ref, interest_quantity, elem_list_interest_quantity, elem_list_local_ref_interest_quantity, node_interest_quantity, node_local_ref_interest_quantity, pos_crack_tip, radius_Ri, radius_Re, spread_cut );
            
//            /// Formulation du pb de reference local
//            ///-------------------------------------
//            TF f_local_ref( m_local_ref );
            
//            /// Proprietes materiaux et Conditions aux limites du pb de reference local
//            ///------------------------------------------------------------------------
//            set_material_properties( f_local_ref, m_local_ref, structure );
//            set_boundary_conditions( f_local_ref, m_local_ref, boundary_condition_D, "direct", structure, loading, mesh_size );
            
//            /// Resolution du pb de reference local
//            ///------------------------------------
//            TicToc t_local_ref;
//            t_local_ref.start();
//            if ( want_iterative_solver_ref == 0 )
//                f_local_ref.solve();
//            else
//                f_local_ref.solve( iterative_criterium_ref );
//            t_local_ref.stop();
//            cout << "Temps de calcul du pb de reference local associe au pb direct : " << t_local_ref.res << endl << endl;
            
//            /// Verification de l'equilibre du pb de reference local associe au pb direct
//            ///--------------------------------------------------------------------------
//              if ( verif_eq )
//                  check_equilibrium( f_local_ref, "de reference local associe au pb direct" );
            
//            /// Definition de l'extracteur du pb de reference local
//            ///----------------------------------------------------
//            TM m_crown_ref;
//            if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" )
//                create_structure_crown( m_local_ref, m_crown_ref, pos_crack_tip, radius_Ri, radius_Re, spread_cut );
//            TF f_crown_ref( m_crown_ref );
//            define_extractor( m_local_ref, m_crown_ref, f_local_ref, f_crown_ref, interest_quantity, direction_extractor, pointwise_interest_quantity, elem_list_local_ref_interest_quantity, node_local_ref_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_Ri, radius_Re, want_local_enrichment );
            
//            /// Calcul de la quantite d'interet locale (quasi-)exacte I_ex
//            ///-----------------------------------------------------------
//            calcul_interest_quantity( m_local_ref, m_crown_ref, f_local_ref, f_crown_ref, "reference", interest_quantity, direction_extractor, pointwise_interest_quantity, elem_list_local_ref_interest_quantity, node_local_ref_interest_quantity, pos_interest_quantity, pos_crack_tip, angle_crack, radius_Ri, radius_Re, I_ex );
//        }
        
//        if ( want_interest_quantity_only == 0 ) {
            
//            ///-------------------------------------///
//            /// Construction de la solution adjoint ///
//            ///-------------------------------------///
            
//            Vec<unsigned> elem_list_adjoint_interest_quantity;
//            Vec<unsigned> elem_list_adjoint_enrichment_zone_1;
//            Vec<unsigned> elem_list_adjoint_enrichment_zone_2;
//            Vec<unsigned> face_list_adjoint_enrichment_zone_12;
//            unsigned node_adjoint_interest_quantity;
//            Vec<unsigned> node_list_adjoint_enrichment;
            
//            /// Maillage du pb adjoint
//            ///-----------------------
//            create_structure( m_adjoint, m_local_ref, "adjoint", structure, mesh_size, loading, deg_p );
//            create_structure_adjoint( m, m_adjoint, deg_p, interest_quantity, direction_extractor, want_local_refinement, l_min_refinement, k_refinement, pointwise_interest_quantity, elem_list_interest_quantity, elem_list_adjoint_interest_quantity, node_interest_quantity, node_adjoint_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_Ri, radius_Re, spread_cut, want_local_enrichment, nb_layers_nodes_enrichment, elem_list_adjoint_enrichment_zone_1, elem_list_adjoint_enrichment_zone_2, face_list_adjoint_enrichment_zone_12, node_list_adjoint_enrichment, debug_geometry, debug_geometry_adjoint );
            
//            display_structure( m_adjoint, m_local_ref, "adjoint", structure, deg_p, want_solve_local_ref );
//            display_params_adjoint( want_local_refinement, l_min_refinement, k_refinement, spread_cut, want_local_enrichment, nb_layers_nodes_enrichment, elem_list_adjoint_enrichment_zone_1, elem_list_adjoint_enrichment_zone_2, face_list_adjoint_enrichment_zone_12, node_list_adjoint_enrichment, want_local_improvement, local_improvement, shape, k_min, k_max, k_opt );
            
//            /// Formulation du pb adjoint
//            ///--------------------------
//            TF f_adjoint( m_adjoint );
            
//            /// Proprietes materiaux et Conditions aux limites du pb adjoint
//            ///-------------------------------------------------------------
//            set_material_properties( f_adjoint, m_adjoint, structure );
//            set_boundary_conditions( f_adjoint, m_adjoint, boundary_condition_D, "adjoint", structure, loading, mesh_size );
//            set_load_conditions( m_adjoint, f_adjoint, m, m_crown, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, interest_quantity, direction_extractor, pointwise_interest_quantity, want_local_enrichment );
            
//            /// Verification des contraintes cinematiques
//            ///------------------------------------------
//            if ( display_constraints )
//                check_constraints( f_adjoint );
            
//            /// Resolution du pb adjoint
//            ///-------------------------
//            TicToc t_adjoint;
//            t_adjoint.start();
//            if ( want_iterative_solver_adjoint == 0 )
//                f_adjoint.solve();
//            else
//                f_adjoint.solve( iterative_criterium_adjoint );
//            t_adjoint.stop();
//            cout << "Temps de calcul du pb adjoint : " << t_adjoint.res << endl << endl;
            
//            if ( want_local_enrichment )
//                calcul_dep_tot_after_solve( m_adjoint );
            
//            /// Verification de l'equilibre du pb adjoint
//            ///------------------------------------------
//              if ( verif_eq )
//                  check_equilibrium( f_adjoint, "adjoint", verif_eq );
            
//            ///-----------------------------------------------------------------------------------///
//            /// Calcul et Affichage des informations relatives a la geometrie du maillage adjoint ///
//            ///-----------------------------------------------------------------------------------///
            
//            calcul_display_geometry( m_adjoint, f_adjoint, debug_geometry_adjoint );
            
//            if ( want_handbook_only == 0 ) {
                
//                ///-------------------------------------------------------------------------------------------------------------------///
//                /// Construction d'un champ de contrainte admissible et Calcul d'un estimateur d'erreur globale associe au pb adjoint ///
//                ///-------------------------------------------------------------------------------------------------------------------///
                
//                T theta_adjoint = 0.;
//                Vec<T> theta_adjoint_elem;
//                Vec< Vec<T> > dep_adjoint_hat;
//                calcul_global_error_estimation( f_adjoint, m_adjoint, "adjoint", method_adjoint, cost_function, penalty_val_N, solver, solver_minimisation, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, geometric_criterium, val_geometric_criterium, val_estimator_criterium, theta_adjoint, theta_adjoint_elem, dep_adjoint_hat, verif_compatibility_conditions, tol_compatibility_conditions, verif_eq_force_fluxes, tol_eq_force_fluxes, verif_solver, tol_solver, verif_solver_enhancement, tol_solver_enhancement, verif_solver_minimisation, tol_solver_minimisation, verif_solver_minimisation_enhancement, tol_solver_minimisation_enhancement, false, false, want_local_enrichment, debug_geometry_adjoint, debug_force_fluxes_adjoint, debug_force_fluxes_enhancement_adjoint, debug_criterium_enhancement_adjoint, debug_error_estimate_adjoint, debug_local_effectivity_index_adjoint, debug_method_adjoint, debug_method_enhancement_adjoint );
                
//                /// Construction de la correspondance entre maillages extraits et maillages initiaux direct/adjoint
//                ///------------------------------------------------------------------------------------------------
//                Vec<unsigned> correspondance_elem_m_adjoint_to_elem_m;
//                correspondance_elem_m_adjoint_to_elem_m.resize( m_adjoint.elem_list.size() );
                
//                Construct_Correspondance_Elem_Mesh_Extracted_To_Elem_Mesh construct_correspondance_elem_m_adjoint_to_elem_m;
//                construct_correspondance_elem_m_adjoint_to_elem_m.correspondance_elem_m_extracted_to_elem_m = &correspondance_elem_m_adjoint_to_elem_m;
//                apply_ij( m_adjoint.elem_list, m.elem_list, construct_correspondance_elem_m_adjoint_to_elem_m );
                
//                ///-----------------------------------------------------------------------------------------------------///
//                /// Calcul de la correction I_hh (avec ou sans introduction de sigma_hat_m) sur la quantite d'interet I ///
//                ///-----------------------------------------------------------------------------------------------------///
                
//                T I_hh = 0.;
//                calcul_correction_interest_quantity( m, m_adjoint, f, f_adjoint, interest_quantity, method, method_adjoint, theta, theta_adjoint, theta_adjoint_elem, correspondance_elem_m_adjoint_to_elem_m, dep_hat, dep_adjoint_hat, I_h, I_hh, want_local_enrichment, want_introduction_sigma_hat_m );
                
//                ///------------------------------------------------------------------------///
//                /// Calcul standard des bornes d'erreur sur la quantite d'interet locale I ///
//                ///------------------------------------------------------------------------///
                
//                calcul_standard_local_error_bounds( m, m_adjoint, f, f_adjoint, method, theta, theta_adjoint, theta_adjoint_elem, correspondance_elem_m_adjoint_to_elem_m, dep_hat, I_h, I_hh, want_introduction_sigma_hat_m );
                
//                ///------------------------------------------------------------------------///
//                /// Calcul ameliore des bornes d'erreur sur la quantite d'interet locale I ///
//                ///------------------------------------------------------------------------///
//                if ( want_local_improvement ) {
//                    calcul_enhanced_local_error_bounds( m, m_adjoint, f, f_adjoint, m_lambda_min, m_lambda_max, m_lambda_opt, m_adjoint_lambda_min, m_adjoint_lambda_opt, deg_p, method, method_adjoint, local_improvement, shape, k_min, k_max, k_opt, interest_quantity, pointwise_interest_quantity, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_Ri, radius_Re, spread_cut, theta, theta_adjoint, dep_hat, dep_adjoint_hat, I_h, I_hh, integration_k, integration_nb_points, debug_method_adjoint, debug_method_enhancement_adjoint, debug_error_estimate_adjoint, want_introduction_sigma_hat_m, want_solve_eig_local_improvement, use_mask_eig_local_improvement );
//                }
//            }
//        }
//    }

    t_total.stop();
    cout << "Temps de calcul total : " << t_total.res << endl << endl;
    
    ///-----------///
    /// Affichage ///
    ///-----------///
    
//     display_vtu_pvd( m, m_ref, m_lambda_min, m_lambda_max, m_lambda_opt, m_crown, pathname, "direct", method, structure, loading, mesh_size, cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, val_geometric_criterium, val_estimator_criterium, geometric_criterium, deg_k, refinement_degree_ref, want_global_discretization_error, want_local_discretization_error, want_global_estimation, want_local_estimation, want_local_improvement, interest_quantity, direction_extractor, pointwise_interest_quantity, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_Ri, radius_Re, local_improvement, shape, k_min, k_max, k_opt, want_local_enrichment, nb_layers_nodes_enrichment, save_vtu, display_vtu, save_pvd, display_pvd, save_vtu_ref, display_vtu_ref, save_vtu_lambda, display_vtu_lambda, save_vtu_crown, display_vtu_crown );
//     if ( want_local_estimation and want_interest_quantity_only == 0 ) {
//         display_vtu_pvd( m_adjoint, m_local_ref, m_adjoint_lambda_min, m_adjoint_lambda_max, m_adjoint_lambda_opt, m_crown, pathname, "adjoint", method_adjoint, structure, loading, mesh_size, cost_function, enhancement_with_geometric_criterium, enhancement_with_estimator_criterium, val_geometric_criterium, val_estimator_criterium, geometric_criterium, deg_k, refinement_degree_ref, want_global_discretization_error_adjoint, want_local_discretization_error_adjoint, want_global_estimation, want_local_estimation, want_local_improvement, interest_quantity, direction_extractor, pointwise_interest_quantity, elem_list_interest_quantity, node_interest_quantity, pos_interest_quantity, pos_crack_tip, radius_Ri, radius_Re, local_improvement, shape, k_min, k_max, k_opt, want_local_enrichment, nb_layers_nodes_enrichment, save_vtu_adjoint, display_vtu_adjoint, save_pvd_adjoint, display_pvd_adjoint, save_vtu_local_ref, display_vtu_local_ref, save_vtu_adjoint_lambda, display_vtu_adjoint_lambda );
//     }
    
}
