#include "../build/problem_error_estimation/all_in_one.h" // sert a forcer le logiciel scons a generer le repertoire build et ses codes sources .h et .cpp correspondant a la formulation
#include "../Mesh.h"
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
    static const string structure = "structure_crack";
    static const string mesh_size = "fine";
    static const string loading = "pull";
    static const unsigned deg_p = 1;
    static const unsigned deg_k = 3;
    static const string boundary_condition_D = "penalty";
    static const bool display_constraints = 0;

    display_pb( dim, structure, deg_p  );

    /// Maillage du pb direct
    /// ---------------------
    TM m; // declaration d'un maillage de type TM
    set_mesh( m, structure, mesh_size, loading, deg_p );

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
    cout << "Resolution du pb direct" << endl << endl;
    TicToc t;
    t.start();
    f.solve();
    t.stop();
    cout << "Temps de calcul de la resolution du pb direct = " << t.res << endl << endl;

//    f.allocate_matrices();
//    f.shift();
//    f.assemble();
////    f.solve_system();
//    f.get_initial_conditions();
//    f.update_variables();
//    f.call_after_solve();

    /// Verification de l'equilibre du pb direct
    /// ----------------------------------------
    check_equilibrium( f, "direct" );

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

//int main( int argc, char **argv, char** env ) {
//    /// Print all environment variables
//    while (*env)
//        printf("%s\n", *env++);
//    return 0;
//}
