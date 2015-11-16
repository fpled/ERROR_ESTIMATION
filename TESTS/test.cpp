#include "../build/problem_error_estimation/all_in_one.h" // sert a forcer le logiciel scons a generer le repertoire build et ses codes sources .h et .cpp correspondant a la formulation
#include "../Structure.h"
#include "../Material_properties.h"
#include "../Boundary_conditions.h"
#include "../GEOMETRY/Calcul_geometry.h"
#include "../GEOMETRY/Geometry.h"
#include "../Display.h"

using namespace LMT;
using namespace std;

int main( int argc, char **argv ) {
    TicToc t_total;
    t_total.start();
    static const unsigned dim = 2;
    static const bool wont_add_nz = true;
    typedef Mesh<Mesh_carac_error_estimation<double,dim> > TM;
    typedef Formulation<TM,FormulationElasticity,DefaultBehavior,double,wont_add_nz> TF;
    typedef TM::Pvec Pvec;
    typedef TM::TNode::T T;
    static const string structure = "structure_crack"; // structure 2D : plate_traction, plate_flexion, plate_hole, plate_crack, structure_crack, eprouvette, weight_sensor, circular_inclusions, circular_holes
                                                     // structure 3D : beam_traction, beam_flexion, beam_hole, plate_hole, plate_hole_full, hub_rotor_helico, reactor_head, door_seal, spot_weld, blade, pipe, SAP, spherical_inclusions, spherical_holes
    static const string mesh_size = "fine"; // maillage pour les structures plate_hole (2D ou 3D), plate_crack, structure_crack, test_specimen, weigth_sensor, spot_weld (3D), reactor_head (3D) : coarse, fine
    static const string loading = "pull"; // chargement pour la structure spot_weld (3D) : pull, shear, peeling et pour la structure plate_crack (2D) : pull, shear
    static const unsigned deg_p = 1; // degre de l'analyse elements finis : 1, 2, ...
    static const unsigned deg_k = 3; // degre supplementaire : 1, 2 , 3, ...
    static const string boundary_condition_D = "penalty"; // methode de prise en compte des conditions aux limites de Dirichlet (en deplacement) pour le pb direct : lagrange, penalty
    static const bool display_constraints = 0; // affichage des contraintes cinematiques

    TM m; // declaration d'un maillage de type TM
    TM m_ref;
    create_structure( m, m_ref, "direct", structure, mesh_size, loading, deg_p );
    display_structure( m, m_ref, "direct", structure, deg_p );

    /// Formulation du pb direct
    /// ------------------------
    TF f( m ); // creation d'une formulation du type TF avec le maillage m

    /// Proprietes materiaux et Conditions aux limites du pb direct
    /// -----------------------------------------------------------
    set_material_properties( f, m, structure );
    set_boundary_conditions( f, m, boundary_condition_D, "direct", structure, loading, mesh_size );

    /// Verification des contraintes cinematiques
    /// -----------------------------------------
    if ( display_constraints )
        check_constraints( f );

    /// Resolution du pb direct
    /// -----------------------
    TicToc t;
    t.start();
    f.solve();
    t.stop();
    cout << "Temps de calcul du pb direct : " << t.res << endl << endl;

    cout << m.type_elements() << endl;

}

//#include "../LMT/include/mesh/meshcaracstd.h"

//using namespace std;
//using namespace LMT;

//int main() {
//    typedef Mesh<MeshCaracStd<2,2> > TM;
//    TM m;
//    m.add_node( TM::Pvec( 0, 1 ) );
//    PRINT( generate( m.node_list, ExtractDM<pos_DM>() ) );
//    PRINT( generate( m.node_list, ExtractDMi<pos_DM>( 0 ) ) );
//    PRINT( generate( m.node_list, ExtractDMi<pos_DM>( 1 ) ) );
//}

//int main () {
//    double random = rand();
//    cout << random << endl;
//    double random_max = RAND_MAX;
//    cout << random_max << endl;
//    cout << random/random_max << endl;
//}
