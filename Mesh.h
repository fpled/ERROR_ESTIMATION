//
// C++ Interface: Maillage
//
// Description: creation du maillage
//
//
// Author: Pled Florent <pled@lmt.ens-cachan.fr>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Mesh_h
#define Mesh_h

#include "LMT/include/mesh/make_rect.h" // sert a definir la fonction make_rect() pour fabriquer un maillage rectangulaire/hypercubique
#include "LMT/include/mesh/read_msh_2.h" // sert a definir la fonction read_msh_2() pour charger un maillage a partir d'un fichier .msh
#include "LMT/include/mesh/read_avs.h" // sert a definir la fonction read_avs() pour charger un maillage a partir d'un fichier .avs
#include "LMT/include/mesh/read_inp.h" // sert a definir la fonction read_inp() pour charger un maillage a partir d'un fichier .inp
#include "LMT/include/mesh/ReaderINP.h" // sert a definir la fonction parse() pour charger un maillage et lire des donnees a partir d'un fichier .inp
#include "LMT/include/mesh/read_vtu.h" // sert a definir la fonction read_vtu() pour charger un maillage a partir d'un fichier .vtu
#include "LMT/include/mesh/read_vtk.h" // sert a definir la fonction read_vtk() pour charger un maillage a partir d'un fichier .vtk

#include "LMT/include/mesh/displayparaview.h" // sert a lire ou ecrire un maillage dans un fichier .vtk (.vti, .vtu, .vtr, .vts, .vtu) ou .pvd (collection de .vtu)
#include "LMT/include/mesh/write_avs.h" // sert a definir la fonction write_avs() pour ecrire un maillage dans un fichier .avs

#include "LMT/include/mesh/remove_lonely_nodes.h" // sert a retirer les noeuds seuls d'un maillage
#include "LMT/include/mesh/refinement.h" // sert a raffiner un maillage selon un critere donne
#include "LMT/include/mesh/sep_mesh.h" // sert a redefinir la peau d'un maillage selon un critere donne
#include "LMT/include/mesh/replace_Quad_by_Triangle.h" // sert a retirer tous les Quad d'un maillage et les remplacer par des Triangle
#include "LMT/include/mesh/replace_Hexa_by_Tetra.h" // sert a retirer tous les Hexa d'un maillage et les remplacer par des Tetra
#include "LMT/include/mesh/replace_Wedge_by_Tetra.h" // sert a retirer tous les Wedge d'un maillage et les remplacer par des Tetra

#include "LMT/include/util/Hdf.h" // sert a lire des donnees a partir d'un fichier .h5 ou .hdf5
#include "LMT/include/containers/Tens3.h"

#include "CONNECTIVITY/Calcul_connectivity.h"

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <string>
#include <hdf5.h>

using namespace LMT;
using namespace std;

/// Affichage des paramètre d'un maillage
/// -------------------------------------
template<class TM>
void display_mesh_carac( const TM &m ) {
    static const unsigned dim = TM::dim;
    cout << "nb de ddl     = " << m.node_list.size() * dim << endl;
    cout << "nb de noeuds  = " << m.node_list.size() << endl;
    cout << "nb d'elements = " << m.elem_list.size() << endl << endl;
}

/// Creation du maillage
/// --------------------
template<class TM>
void set_mesh( TM &m, const string &structure, const string &mesh_size, const string &loading, const unsigned &deg_p = 1, const unsigned &refinement_level_ref = 2, const bool want_global_discretization_error = false, const bool want_local_discretization_error = false, const bool disp = true ) {

    static const unsigned dim = TM::dim;
    typedef typename TM::Pvec Pvec;
    typedef typename TM::TNode::T T;

    /// Dimension 2
    /// -----------
    if ( dim == 2 ) {
        /// Plaque rectangulaire 2D en traction
        /// Plaque rectangulaire 2D en flexion
        /// -----------------------------------
        if ( structure == "plate_traction" or structure == "plate_flexion" ) {
            switch ( deg_p ) {
            case 1 :
                make_rect( m, Triangle(), Pvec( 0., 0. ), Pvec( 2., 1. ), Pvec( 3, 3 ) );
//                make_rect( m, Quad(), Pvec( 0., 0. ), Pvec( 2. , 1. ), Pvec( 3, 3 ) );
                break;
            case 2 :
                make_rect( m, Triangle_6(), Pvec( 0., 0. ), Pvec( 2., 1. ), Pvec( 3, 3 ) );
//                make_rect( m, Quad_8(), Pvec( 0., 0. ), Pvec( 2. , 1. ), Pvec( 3, 3 ) );
//                make_rect( m, Quad_9(), Pvec( 0., 0. ), Pvec( 2. , 1. ), Pvec( 3, 3 ) );
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Quart de plaque rectangulaire trouee 2D
        /// ---------------------------------------
        else if ( structure == "plate_hole" ) {
            switch ( deg_p ) {
            case 1 :
                if ( mesh_size == "coarse" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) ) {
//                        read_msh_2( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Triangle.msh" );
                        read_avs( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Quad.avs" );
                    }
                    else {
//                        read_vtu( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Triangle_direct_global_local_discretization_error.vtu" );
                        read_vtu( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Quad_direct_global_local_discretization_error.vtu" );
                    }
                }
                else if ( mesh_size == "fine" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) ) {
                        read_msh_2( m, "MESH/PLATE_HOLE_2D/plate_hole_fine_Triangle.msh" );
//                        read_avs( m, "MESH/PLATE_HOLE_2D/plate_hole_fine_Triangle_Quad.avs" );
                    }
                    else {
                        read_vtu( m, "MESH/PLATE_HOLE_2D/plate_hole_fine_Triangle_direct_global_local_discretization_error.vtu" );
//                        read_vtu( m, "MESH/PLATE_HOLE_2D/plate_hole_fine_Triangle_Quad_direct_global_local_discretization_error.vtu" );
                    }
                }
                else
                    cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
                break;
            case 2 :
                if ( mesh_size == "coarse" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) ) {
                        read_msh_2( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Triangle_6.msh" );
//                        read_avs( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Quad_8.avs" );
//                        read_avs( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Quad_9.avs" );
                    }
                    else {
                        read_vtu( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Triangle_6_direct_global_local_discretization_error.vtu" );
//                        read_vtu( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Quad_8_direct_global_local_discretization_error.vtu" );
//                        read_vtu( m, "MESH/PLATE_HOLE_2D/plate_hole_coarse_Quad_9_direct_global_local_discretization_error.vtu" );
                    }
                }
                else if ( mesh_size == "fine" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/PLATE_HOLE_2D/plate_hole_fine_Triangle_6.msh" );
                    else
                        read_vtu( m, "MESH/PLATE_HOLE_2D/plate_hole_fine_Triangle_6_direct_global_local_discretization_error.vtu" );
                }
                else
                    cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Plaque fissuree 2D
        /// ------------------
        else if ( structure == "plate_crack" ) {
            if ( mesh_size == "very_coarse" )
                read_msh_2( m, "MESH/PLATE_CRACK_2D/plate_crack_very_coarse_Triangle.msh" );
            else if ( mesh_size == "coarse" )
                read_msh_2( m, "MESH/PLATE_CRACK_2D/plate_crack_coarse_Triangle.msh" );
            else if ( mesh_size == "fine" )
                read_msh_2( m, "MESH/PLATE_CRACK_2D/plate_crack_fine_Triangle.msh" );
            else if ( mesh_size == "very_fine" )
                read_msh_2( m, "MESH/PLATE_CRACK_2D/plate_crack_very_fine_Triangle.msh" );
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Structure fissuree 2D
        /// ---------------------
        else if ( structure == "structure_crack" ) {
            switch ( deg_p ) {
            case 1 :
                if ( mesh_size == "coarse" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) ) {
//                        read_msh_2( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Triangle.msh" );
                        read_avs( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Quad.avs" );
                    }
                    else {
//                        read_vtu( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Triangle_direct_global_local_discretization_error.vtu" );
                        read_vtu( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Quad_direct_global_local_discretization_error.vtu" );
                    }
                }
                else if ( mesh_size == "fine" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_fine_Triangle.msh" );
                    else {
                        read_vtu( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_fine_Triangle_direct_global_local_discretization_error.vtu" );
//                        read_vtu( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Triangle_direct_global_local_discretization_error_old.vtu" );
                    }
                }
                else
                    cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
                break;
            case 2 :
                if ( mesh_size == "coarse" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) ) {
                        read_msh_2( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Triangle_6.msh" );
//                        read_avs( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Quad_8.avs" );
//                        read_avs( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Quad_9.avs" );
                    }
                    else {
                        read_vtu( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Triangle_6_direct_global_local_discretization_error.vtu" );
//                        read_vtu( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Quad_8_direct_global_local_discretization_error.vtu" );
//                        read_vtu( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_coarse_Quad_9_direct_global_local_discretization_error.vtu" );
                    }
                }
                else if ( mesh_size == "fine" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_fine_Triangle_6.msh" );
                    else
                        read_vtu( m, "MESH/STRUCTURE_CRACK_2D/structure_crack_fine_Triangle_6_direct_global_local_discretization_error.vtu" );
                }
                else
                    cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Eprouvette 2D
        /// -------------
        else if ( structure == "test_specimen" ) {
            switch ( deg_p ) {
            case 1 :
                if ( mesh_size == "coarse" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) ) {
//                        read_avs( m, "MESH/TEST_SPECIMEN_2D/test_specimen_coarse_Triangle.avs" );
                        read_msh_2( m, "MESH/TEST_SPECIMEN_2D/test_specimen_coarse_Triangle.msh" );
                    }
                    else {
//                        read_vtu( m, "MESH/TEST_SPECIMEN_2D/test_specimen_coarse_Triangle_direct_global_local_discretization_error_avs.vtu" );
                        read_vtu( m, "MESH/TEST_SPECIMEN_2D/test_specimen_coarse_Triangle_direct_global_local_discretization_error_msh.vtu" );
                    }
                }
                else if ( mesh_size == "fine" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/TEST_SPECIMEN_2D/test_specimen_fine_Triangle.msh" );
                    else
                        read_vtu( m, "MESH/TEST_SPECIMEN_2D/test_specimen_fine_Triangle_direct_global_local_discretization_error_msh.vtu" );
                }
                else
                    cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
                break;
            case 2 :
                if ( mesh_size == "coarse" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/TEST_SPECIMEN_2D/test_specimen_coarse_Triangle_6.msh" );
                    else
                        read_vtu( m, "MESH/TEST_SPECIMEN_2D/test_specimen_coarse_Triangle_6_direct_global_local_discretization_error.vtu" );
                }
                else if ( mesh_size == "fine" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/TEST_SPECIMEN_2D/test_specimen_fine_Triangle_6.msh" );
                    else
                        read_vtu( m, "MESH/TEST_SPECIMEN_2D/test_specimen_fine_Triangle_6_direct_global_local_discretization_error.vtu" );
                }
                else
                    cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Capteur d'effort 2D
        /// -------------------
        else if ( structure == "weight_sensor" ) {
            switch ( deg_p ) {
            case 1 :
                if ( mesh_size == "coarse" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/WEIGHT_SENSOR_2D/weight_sensor_coarse_Triangle.msh" );
                    else
                        read_vtu( m, "MESH/WEIGHT_SENSOR_2D/weight_sensor_coarse_Triangle_direct_global_local_discretization_error.vtu" );
                }
                else if ( mesh_size == "fine" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/WEIGHT_SENSOR_2D/weight_sensor_fine_Triangle.msh" );
                    else
                        read_vtu( m, "MESH/WEIGHT_SENSOR_2D/weight_sensor_fine_Triangle_direct_global_local_discretization_error.vtu" );
                }
                else
                    cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
                break;
            case 2 :
                if ( mesh_size == "coarse" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/WEIGHT_SENSOR_2D/weight_sensor_coarse_Triangle_6.msh" );
                    else
                        read_vtu( m, "MESH/WEIGHT_SENSOR_2D/weight_sensor_coarse_Triangle_6_direct_global_local_discretization_error.vtu" );
                }
                else if ( mesh_size == "fine" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_msh_2( m, "MESH/WEIGHT_SENSOR_2D/weight_sensor_fine_Triangle_6.msh" );
                    else
                        read_vtu( m, "MESH/WEIGHT_SENSOR_2D/weight_sensor_fine_Triangle_6_direct_global_local_discretization_error.vtu" );
                }
                else
                    cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Cercle 2D
        /// ---------
        else if ( structure == "circle" ) {
            if ( mesh_size == "very_coarse" )
                read_msh_2( m, "MESH/CIRCLE_2D/circle_very_coarse_Triangle.msh" );
            else if ( mesh_size == "coarse" )
                read_msh_2( m, "MESH/CIRCLE_2D/circle_coarse_Triangle.msh" );
            else if ( mesh_size == "fine" )
                read_msh_2( m, "MESH/CIRCLE_2D/circle_fine_Triangle.msh" );
            else if ( mesh_size == "very_fine" )
                read_msh_2( m, "MESH/CIRCLE_2D/circle_very_fine_Triangle.msh" );
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Inclusions circulaires 2D
        /// -------------------------
        else if ( structure == "circular_inclusions" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                read_msh_2( m, "MESH/CIRCULAR_INCLUSIONS_2D/circular_inclusions_Triangle.msh" );
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/CIRCULAR_INCLUSIONS_2D/circular_inclusions_Triangle_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/CIRCULAR_INCLUSIONS_2D/circular_inclusions_Triangle_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
                    read_msh_2( m, "MESH/CIRCULAR_INCLUSIONS_2D/circular_inclusions_Triangle.msh" );
                    break;
                }
            }
        }
        /// Trous circulaires 2D
        /// --------------------
        else if ( structure == "circular_holes" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                read_msh_2( m, "MESH/CIRCULAR_HOLES_2D/circular_holes_Triangle.msh" );
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/CIRCULAR_HOLES_2D/circular_holes_Triangle_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/CIRCULAR_HOLES_2D/circular_holes_Triangle_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
                    read_msh_2( m, "MESH/CIRCULAR_HOLES_2D/circular_holes_Triangle.msh" );
                    break;
                }
            }
        }
        /// Carre 2D
        /// --------
        else if ( structure.find("square") != string::npos ) {
            size_t off = structure.rfind( "_" );
            string str = structure.substr( off+1 );
            istringstream buffer(str);
            int N;
            buffer >> N;
            make_rect( m, Quad(), Pvec( 0., 0. ), Pvec( 1. , 1. ), Pvec( N+1, N+1 ) );
        }
        else
            cerr << "structure " << structure << " not implemented..." << endl << endl;
    }
    /// Dimension 3
    /// -----------
    else if ( dim == 3 ) {
        /// Barre rectangulaire 3D en traction
        /// Barre rectangulaire 3D en flexion
        /// ----------------------------------
        if ( structure == "beam_traction" or structure == "beam_flexion" ) {
            switch ( deg_p ) {
            case 1 :
                make_rect( m, Tetra(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 10, 5, 5 ) );
//                make_rect( m, Hexa(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 3, 3, 3 ) );
                break;
            case 2 :
                make_rect( m, Tetra_10(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 2, 2, 2 ) );
//                make_rect( m, Hexa_20(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 2, 2, 2 ) );
//                make_rect( m, Hexa_27(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 2, 2, 2 ) );
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Barre rectangulaire trouee 3D
        /// -----------------------------
        else if ( structure == "beam_hole" ) {
            if ( mesh_size == "coarse" ) {
                if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                    read_msh_2( m, "MESH/BEAM_HOLE_3D/beam_hole_coarse_Tetra.msh" );
                else {
                    switch ( refinement_level_ref ) {
                    case 1 :
                        read_vtu( m, "MESH/BEAM_HOLE_3D/beam_hole_coarse_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                        break;
                    case 2 :
                        read_vtu( m, "MESH/BEAM_HOLE_3D/beam_hole_coarse_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                        break;
                    default :
                        read_msh_2( m, "MESH/BEAM_HOLE_3D/beam_hole_coarse_Tetra.msh" );
                        break;
                    }
                }
            }
            else if ( mesh_size == "fine" ) {
                if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                    read_msh_2( m, "MESH/BEAM_HOLE_3D/beam_hole_fine_Tetra.msh" );
                else {
                    switch ( refinement_level_ref ) {
                    case 1 :
                        read_vtu( m, "MESH/BEAM_HOLE_3D/beam_hole_fine_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                        break;
                    case 2 :
                        read_vtu( m, "MESH/BEAM_HOLE_3D/beam_hole_fine_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                        break;
                    default :
                        read_msh_2( m, "MESH/BEAM_HOLE_3D/beam_hole_fine_Tetra.msh" );
                        break;
                    }
                }
            }
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Quart de plaque rectangulaire trouee 3D
        /// ---------------------------------------
        else if ( structure == "plate_hole" ) {
            if ( mesh_size == "coarse" ) {
                if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                    read_msh_2( m, "MESH/PLATE_HOLE_3D/plate_hole_coarse_Tetra.msh" );
                else {
                    switch ( refinement_level_ref ) {
                    case 1 :
                        read_vtu( m, "MESH/PLATE_HOLE_3D/plate_hole_coarse_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                        break;
                    case 2 :
                        read_vtu( m, "MESH/PLATE_HOLE_3D/plate_hole_coarse_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                        break;
                    case 3 :
                        read_vtu( m, "MESH/PLATE_HOLE_3D/plate_hole_coarse_Tetra_direct_global_local_discretization_error_ref_3.vtu" );
                        break;
                    default :
                        read_msh_2( m, "MESH/PLATE_HOLE_3D/plate_hole_coarse_Tetra.msh" );
                        break;
                    }
                }
            }
            else if ( mesh_size == "fine" ) {
                if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                    read_msh_2( m, "MESH/PLATE_HOLE_3D/plate_hole_fine_Tetra.msh" );
                else {
                    switch ( refinement_level_ref ) {
                    case 1 :
                        read_vtu( m, "MESH/PLATE_HOLE_3D/plate_hole_fine_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                        break;
                    case 2 :
                        read_vtu( m, "MESH/PLATE_HOLE_3D/plate_hole_fine_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                        break;
                    default :
                        read_msh_2( m, "MESH/PLATE_HOLE_3D/plate_hole_fine_Tetra.msh" );
                        break;
                    }
                }
            }
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Plaque rectangulaire trouee complete 3D
        /// ---------------------------------------
        else if ( structure == "plate_hole_full" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                read_avs( m, "MESH/PLATE_HOLE_3D/plate_hole_full_Hexa.avs" );
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/PLATE_HOLE_3D/plate_hole_full_Hexa_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/PLATE_HOLE_3D/plate_hole_full_Hexa_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
                    read_avs( m, "MESH/PLATE_HOLE_3D/plate_hole_full_Hexa.avs" );
                    break;
                }
            }
        }
        /// Moyeu-rotor 3D de l'helicoptere NH90B
        /// -------------------------------------
        else if ( structure == "hub_rotor_helico" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                read_avs( m, "MESH/HUB_ROTOR_HELICO_3D/hub_rotor_helico_Tetra.avs" );
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/HUB_ROTOR_HELICO_3D/hub_rotor_helico_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/HUB_ROTOR_HELICO_3D/hub_rotor_helico_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
                    read_avs( m, "MESH/HUB_ROTOR_HELICO_3D/hub_rotor_helico_Tetra.avs" );
                    break;
                }
            }
        }
        /// Quart de tete de reacteur nucleaire 3D
        /// --------------------------------------
        else if ( structure == "reactor_head" ) {
            if ( mesh_size == "coarse" ) {
                if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                    read_msh_2( m, "MESH/REACTOR_HEAD_3D/reactor_head_coarse_Tetra.msh" );
                else {
                    switch ( refinement_level_ref ) {
                    case 1 :
                        read_vtu( m, "MESH/REACTOR_HEAD_3D/reactor_head_coarse_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                        break;
                    case 2 :
                        read_vtu( m, "MESH/REACTOR_HEAD_3D/reactor_head_coarse_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                        break;
                    default :
                        read_msh_2( m, "MESH/REACTOR_HEAD_3D/reactor_head_coarse_Tetra.msh" );
                        break;
                    }
                }
            }
            else if ( mesh_size == "fine" ) {
                if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                    read_avs( m, "MESH/REACTOR_HEAD_3D/reactor_head_fine_Tetra.avs" );
                else {
                    switch ( refinement_level_ref ) {
                    case 1 :
                        read_vtu( m, "MESH/REACTOR_HEAD_3D/reactor_head_fine_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                        break;
                    case 2 :
                        read_vtu( m, "MESH/REACTOR_HEAD_3D/reactor_head_fine_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                        break;
                    default :
                        read_avs( m, "MESH/REACTOR_HEAD_3D/reactor_head_fine_Tetra.avs" );
                        break;
                    }
                }
            }
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Joint de porte 3D de vehicule auto
        /// ----------------------------------
        else if ( structure == "door_seal" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                read_avs( m, "MESH/DOOR_SEAL_3D/door_seal_Hexa.avs" );
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/DOOR_SEAL_3D/door_seal_Hexa_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/DOOR_SEAL_3D/door_seal_Hexa_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
                    read_avs( m, "MESH/DOOR_SEAL_3D/door_seal_Hexa.avs" );
                    break;
                }
            }
        }
        /// Point de soudure 3D
        /// -------------------
        else if ( structure == "spot_weld" ) {
            if ( mesh_size == "coarse" ) {
                if ( loading == "pull" or loading == "shear" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_avs( m, "MESH/SPOT_WELD_3D/spot_weld_pull_shear_coarse_Hexa.avs" );
                    else {
                        switch ( refinement_level_ref ) {
                        case 1 :
                            if ( loading == "pull" )
                                read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_pull_coarse_Hexa_direct_global_local_discretization_error_ref_1.vtu" );
                            else if ( loading == "shear" )
                                read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_shear_coarse_Hexa_direct_global_local_discretization_error_ref_1.vtu" );
                            break;
                        case 2 :
                            if ( loading == "pull" )
                                read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_pull_coarse_Hexa_direct_global_local_discretization_error_ref_2.vtu" );
                            else if ( loading == "shear" )
                                read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_shear_coarse_Hexa_direct_global_local_discretization_error_ref_2.vtu" );
                            break;
                        default :
                            read_avs( m, "MESH/SPOT_WELD_3D/spot_weld_pull_shear_coarse_Hexa.avs" );
                            break;
                        }
                    }
                }
                else if ( loading == "peeling" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_avs( m, "MESH/SPOT_WELD_3D/spot_weld_peeling_coarse_Hexa.avs" );
                    else {
                        switch ( refinement_level_ref ) {
                        case 1 :
                            read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_peeling_coarse_Hexa_direct_global_local_discretization_error_ref_1.vtu" );
                            break;
                        case 2 :
                            read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_peeling_coarse_Hexa_direct_global_local_discretization_error_ref_2.vtu" );
                            break;
                        default :
                            read_avs( m, "MESH/SPOT_WELD_3D/spot_weld_peeling_coarse_Hexa.avs" );
                            break;
                        }
                    }
                }
                else
                    cerr << "loading " << loading << " not implemented..." << endl << endl;
            }
            else if ( mesh_size == "fine" ) {
                if ( loading == "pull" or loading == "shear" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_avs( m, "MESH/SPOT_WELD_3D/spot_weld_pull_shear_fine_Hexa.avs" );
                    else {
                        switch ( refinement_level_ref ) {
                        case 1 :
                            if ( loading == "pull" )
                                read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_pull_fine_Hexa_direct_global_local_discretization_error_ref_1.vtu" );
                            else if ( loading == "shear" )
                                read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_shear_fine_Hexa_direct_global_local_discretization_error_ref_1.vtu" );
                            break;
                        case 2 :
                            if ( loading == "pull" )
                                read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_pull_fine_Hexa_direct_global_local_discretization_error_ref_2.vtu" );
                            else if ( loading == "shear" )
                                read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_shear_fine_Hexa_direct_global_local_discretization_error_ref_2.vtu" );
                            break;
                        default :
                            read_avs( m, "MESH/SPOT_WELD_3D/spot_weld_pull_shear_fine_Hexa.avs" );
                            break;
                        }
                    }
                }
                else if ( loading == "peeling" ) {
                    if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                        read_avs( m, "MESH/SPOT_WELD_3D/spot_weld_peeling_fine_Hexa.avs" );
                    else {
                        switch ( refinement_level_ref ) {
                        case 1 :
                            read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_peeling_fine_Hexa_direct_global_local_discretization_error_ref_1.vtu" );
                            break;
                        case 2 :
                            read_vtu( m, "MESH/SPOT_WELD_3D/spot_weld_peeling_fine_Hexa_direct_global_local_discretization_error_ref_2.vtu" );
                            break;
                        default :
                            read_avs( m, "MESH/SPOT_WELD_3D/spot_weld_peeling_fine_Hexa.avs" );
                            break;
                        }
                    }
                }
                else
                    cerr << "loading " << loading << " not implemented..." << endl << endl;
            }
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Ailette 3D
        /// ----------
        else if ( structure == "blade" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) ) {
//                read_avs( m, "MESH/BLADE_3D/blade_Tetra.avs" );
//                read_inp( m, "MESH/BLADE_3D/blade_Tetra.inp" );
                ReaderINP<TM> RI( m, "MESH/BLADE_3D/blade_Tetra.inp" );
            }
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/BLADE_3D/blade_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/BLADE_3D/blade_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
//                    read_avs( m, "MESH/BLADE_3D/blade_Tetra.avs" );
//                    read_inp( m, "MESH/BLADE_3D/blade_Tetra.inp" );
                    ReaderINP<TM> RI( m, "MESH/BLADE_3D/blade_Tetra.inp" );
                    break;
                }
            }
        }
        /// Quart de conduite 3D
        /// --------------------
        else if ( structure == "pipe" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) ) {
//                read_avs( m, "MESH/PIPE_3D/pipe_Tetra.avs" );
//                read_inp( m, "MESH/PIPE_3D/pipe_Tetra.inp" );
                ReaderINP<TM> RI( m, "MESH/PIPE_3D/pipe_Tetra.inp" );
            }
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/PIPE_3D/pipe_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/PIPE_3D/pipe_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
                    //                read_avs( m, "MESH/PIPE_3D/pipe_Tetra.avs" );
                    //                read_inp( m, "MESH/PIPE_3D/pipe_Tetra.inp" );
                    ReaderINP<TM> RI( m, "MESH/PIPE_3D/pipe_Tetra.inp" );
                    break;
                }
            }
        }
        /// SAP 3D
        /// ------
        else if ( structure == "SAP" ) {
            if ( mesh_size == "coarse" ) {
                if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                    read_msh_2( m, "MESH/SAP_3D/SAP_coarse_Tetra.msh" );
                else {
                    switch ( refinement_level_ref ) {
                    case 1 :
                        read_vtu( m, "MESH/SAP_3D/SAP_coarse_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                        break;
                    case 2 :
                        read_vtu( m, "MESH/SAP_3D/SAP_coarse_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                        break;
                    default :
                        read_msh_2( m, "MESH/SAP_3D/SAP_coarse_Tetra.msh" );
                        break;
                    }
                }
            }
            else if ( mesh_size == "fine" ) {
                if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                    read_msh_2( m, "MESH/SAP_3D/SAP_fine_Tetra.msh" );
                else {
                    switch ( refinement_level_ref ) {
                    case 1 :
                        read_vtu( m, "MESH/SAP_3D/SAP_fine_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                        break;
                    case 2 :
                        read_vtu( m, "MESH/SAP_3D/SAP_fine_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                        break;
                    default :
                        read_msh_2( m, "MESH/SAP_3D/SAP_fine_Tetra.msh" );
                        break;
                    }
                }
            }
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Inclusions sphériques 3D
        /// ------------------------
        else if ( structure == "spherical_inclusions" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                read_msh_2( m, "MESH/SPHERICAL_INCLUSIONS_3D/spherical_inclusions_Tetra.msh" );
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/SPHERICAL_INCLUSIONS_3D/spherical_inclusions_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/SPHERICAL_INCLUSIONS_3D/spherical_inclusions_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
                    read_msh_2( m, "MESH/SPHERICAL_INCLUSIONS_3D/spherical_inclusions_Tetra.msh" );
                    break;
                }
            }
        }
        /// Trous sphériques 3D
        /// -------------------
        else if ( structure == "spherical_holes" ) {
            if ( not( want_global_discretization_error ) and not( want_local_discretization_error ) )
                read_msh_2( m, "MESH/SPHERICAL_HOLES_3D/spherical_holes_Tetra.msh" );
            else {
                switch ( refinement_level_ref ) {
                case 1 :
                    read_vtu( m, "MESH/SPHERICAL_HOLES_3D/spherical_holes_Tetra_direct_global_local_discretization_error_ref_1.vtu" );
                    break;
                case 2 :
                    read_vtu( m, "MESH/SPHERICAL_HOLES_3D/spherical_holes_Tetra_direct_global_local_discretization_error_ref_2.vtu" );
                    break;
                default :
                    read_msh_2( m, "MESH/SPHERICAL_HOLES_3D/spherical_holes_Tetra.msh" );
                    break;
                }
            }
        }
        /// Eprouvette 3D
        /// -------------
        else if ( structure.find("test_specimen") != string::npos ) {
//            string filename = "MESH/TEST_SPECIMEN_3D/" + structure + "_Hexa.inp";
            string filename = "MESH/TEST_SPECIMEN_3D/" + structure + "_Tetra.inp";
//            read_inp( m, filename );
            ReaderINP<TM> RI( m, filename.c_str() );
//            replace_Hexa_by_Tetra( m );
        }
        /// Cube 3D
        /// -------
        else if ( structure.find("cube") != string::npos ) {
            size_t off = structure.rfind( "_" );
            string str = structure.substr( off+1 );
            istringstream buffer(str);
            int N;
            buffer >> N;
            make_rect( m, Hexa(), Pvec( 0., 0., 0. ), Pvec( 1. , 1., 1. ), Pvec( N+1, N+1, N+1 ) );
        }
        /// Sphere 3D avec noeud au centre
        /// ------------------------------
        else if ( structure == "sphere_center" ) {
            if ( mesh_size == "very_coarse" )
                read_msh_2( m, "MESH/SPHERE_3D/sphere_center_very_coarse_Tetra.msh" );
            else if ( mesh_size == "coarse" )
                read_msh_2( m, "MESH/SPHERE_3D/sphere_center_coarse_Tetra.msh" );
            else if ( mesh_size == "fine" )
                read_msh_2( m, "MESH/SPHERE_3D/sphere_center_fine_Tetra.msh" );
            else if ( mesh_size == "very_fine" )
                read_msh_2( m, "MESH/SPHERE_3D/sphere_center_very_fine_Tetra.msh" );
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Sphere 3D sans noeud au centre
        /// ------------------------------
        else if ( structure == "sphere" ) {
            if ( mesh_size == "very_coarse" )
                read_msh_2( m, "MESH/SPHERE_3D/sphere_very_coarse_Tetra.msh" );
            else if ( mesh_size == "coarse" )
                read_msh_2( m, "MESH/SPHERE_3D/sphere_coarse_Tetra.msh" );
            else if ( mesh_size == "fine" )
                read_msh_2( m, "MESH/SPHERE_3D/sphere_fine_Tetra.msh" );
            else if ( mesh_size == "very_fine" )
                read_msh_2( m, "MESH/SPHERE_3D/sphere_very_fine_Tetra.msh" );
            else
                cerr << "mesh_size " << mesh_size << " not implemented..." << endl << endl;
        }
        /// Sphere 3D creuse avec noeud au centre
        /// -------------------------------------
        else if ( structure == "sphere_hollow" )
            read_msh_2( m, "MESH/SPHERE_3D/sphere_hollow_very_fine_Tetra.msh" );
        else
            cerr << "structure " << structure << " not implemented..." << endl << endl;
    }
    else
        cerr << "dim " << dim << " not implemented..." << endl << endl;

    if ( m.node_list.size() ) {
        if ( remove_lonely_nodes( m ) )
            cerr << "Des noeuds seuls ont ete retires du maillage du pb direct..." << endl << endl;
        if ( not m.elem_list.size() ) {
            cerr << "Arret brutal, car il n'y a aucun element dans le maillage du pb direct..." << endl << endl;
            throw "Baleinou sous caillou...";
        }
        if ( disp ) {
            cout << "maillage du pb direct :" << endl;
            display_mesh_carac( m );
        }
    }
    else {
        cerr << "Arret brutal, car il n'y a aucun noeud dans le maillage du pb direct..." << endl << endl;
        throw "Ourson sous gravillon...";
    }
}

/// Creation du maillage de reference pour le calcul de la mesure de l'erreur de discretisation
/// -------------------------------------------------------------------------------------------
template<class TM>
void set_mesh_ref( TM &m_ref, const TM &m, const string &structure, const unsigned &deg_p = 1, const unsigned &refinement_level_ref = 2, const bool disp = true ) {

    static const unsigned dim = TM::dim;
    typedef typename TM::Pvec Pvec;
    typedef typename TM::TNode::T T;

    /// Dimension 2
    /// -----------
    if ( dim == 2 ) {
        /// Plaque rectangulaire 2D en traction
        /// Plaque rectangulaire 2D en flexion
        /// -----------------------------------
        if ( structure == "plate_traction" or structure == "plate_flexion" ) {
            switch ( deg_p ) {
            case 1 :
                make_rect( m_ref, Triangle(), Pvec( 0., 0. ), Pvec( 2. , 1. ), Pvec( 101, 101 ) );
//                make_rect( m_ref, Quad(), Pvec( 0., 0. ), Pvec( 2. , 1. ), Pvec( 101, 101 ) );
                break;
            case 2 :
                make_rect( m_ref, Triangle_6(), Pvec( 0., 0. ), Pvec( 2. , 1. ), Pvec( 101, 101 ) );
//                make_rect( m_ref, Quad_8(), Pvec( 0., 0. ), Pvec( 2. , 1. ), Pvec( 101, 101 ) );
//                make_rect( m_ref, Quad_9(), Pvec( 0., 0. ), Pvec( 2. , 1. ), Pvec( 101, 101 ) );
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Quart de plaque rectangulaire trouee 2D
        /// ---------------------------------------
        else if ( structure == "plate_hole" ) {
            switch ( deg_p ) {
            case 1 :
                read_msh_2( m_ref, "MESH/PLATE_HOLE_2D/plate_hole_very_fine_Triangle.msh" );
                break;
            case 2 :
                read_msh_2( m_ref, "MESH/PLATE_HOLE_2D/plate_hole_very_fine_Triangle_6.msh" );
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Plaque fissuree 2D
        /// ------------------
        else if ( structure == "plate_crack" ) {
            read_msh_2( m_ref, "MESH/PLATE_CRACK_2D/plate_crack_very_fine_Triangle.msh" );
        }
        /// Structure fissuree 2D
        /// ---------------------
        else if ( structure == "structure_crack" ) {
            switch ( deg_p ) {
            case 1 :
//                read_msh_2( m_ref, "MESH/STRUCTURE_CRACK_2D/structure_crack_very_fine_Triangle.msh" );
                read_msh_2( m_ref, "MESH/STRUCTURE_CRACK_2D/structure_crack_very_fine_Triangle_6.msh" );
                break;
            case 2 :
                read_msh_2( m_ref, "MESH/STRUCTURE_CRACK_2D/structure_crack_very_fine_Triangle_6.msh" );
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Eprouvette 2D
        /// -------------
        else if ( structure == "test_specimen" ) {
            switch ( deg_p ) {
            case 1 :
                read_msh_2( m_ref, "MESH/TEST_SPECIMEN_2D/test_specimen_very_fine_Triangle.msh" );
                break;
            case 2 :
                read_msh_2( m_ref, "MESH/TEST_SPECIMEN_2D/test_specimen_very_fine_Triangle_6.msh" );
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Capteur d'effort 2D
        /// -------------------
        else if ( structure == "weight_sensor" ) {
            switch ( deg_p ) {
            case 1 :
                read_msh_2( m_ref, "MESH/WEIGHT_SENSOR_2D/weight_sensor_very_fine_Triangle.msh" );
                break;
            case 2 :
                read_msh_2( m_ref, "MESH/WEIGHT_SENSOR_2D/weight_sensor_very_fine_Triangle_6.msh" );
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Cercle 2D
        /// ---------
        else if ( structure == "circle" ) {
            read_msh_2( m_ref, "MESH/CIRCLE_2D/circle_very_fine_Triangle.msh" );
        }
        /// Inclusions circulaires 2D
        /// -------------------------
        else if ( structure == "circular_inclusions" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Trous circulaires 2D
        /// --------------------
        else if ( structure == "circular_holes" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Carre 2D
        /// --------
        else if ( structure.find("square") != string::npos ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        else
            cerr << "structure " << structure << " not implemented..." << endl << endl;
    }
    /// Dimension 3
    /// -----------
    else if ( dim == 3 ) {
        /// Barre rectangulaire 3D en traction
        /// Barre rectangulaire 3D en flexion
        /// ----------------------------------
        if ( structure == "beam_traction" or structure == "beam_flexion" ) {
            switch ( deg_p ) {
            case 1 :
                make_rect( m_ref, Tetra(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 101, 101, 101 ) );
//                make_rect( m_ref, Hexa(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 101, 101, 101 ) );
                break;
            case 2 :
                make_rect( m_ref, Tetra_10(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 51, 51, 51 ) );
//                make_rect( m_ref, Hexa_20(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 51, 51, 51 ) );
//                make_rect( m_ref, Hexa_27(), Pvec( 0., 0., 0. ), Pvec( 2., 1., 1. ), Pvec( 51, 51, 51 ) );
                break;
            default :
                cerr << "deg " << deg_p << " not implemented..." << endl << endl;
                break;
            }
        }
        /// Barre rectangulaire trouee 3D
        /// -----------------------------
        else if ( structure == "beam_hole" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Quart de plaque rectangulaire trouee 3D
        /// ---------------------------------------
        else if ( structure == "plate_hole" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Plaque rectangulaire trouee complete 3D
        /// ---------------------------------------
        else if ( structure == "plate_hole_full" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Moyeu-rotor 3D de l'helicoptere NH90B
        /// -------------------------------------
        else if ( structure == "hub_rotor_helico" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Quart de tete de reacteur nucleaire 3D
        /// --------------------------------------
        else if ( structure == "reactor_head" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Joint de porte 3D de vehicule auto
        /// ----------------------------------
        else if ( structure == "door_seal" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Point de soudure 3D
        /// -------------------
        else if ( structure == "spot_weld" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Ailette 3D
        /// ----------
        else if ( structure == "blade" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Quart de conduite 3D
        /// --------------------
        else if ( structure == "pipe" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// SAP 3D
        /// ------
        else if ( structure == "SAP" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Inclusions sphériques 3D
        /// ------------------------
        else if ( structure == "spherical_inclusions" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Trous sphériques 3D
        /// -------------------
        else if ( structure == "spherical_holes" ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Eprouvette 3D
        /// -------------
        else if ( structure.find("test_specimen") != string::npos ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Cube 3D
        /// -------
        else if ( structure.find("cube") != string::npos ) {
            m_ref = m;
            for (unsigned n=0;n<refinement_level_ref;++n)
                divide( m_ref );
        }
        /// Sphere 3D avec noeud au centre
        /// ------------------------------
        else if ( structure == "sphere_center" ) {
            read_msh_2( m_ref, "MESH/SPHERE_3D/sphere_center_very_fine_Tetra.msh" );
        }
        /// Sphere 3D sans noeud au centre
        /// ------------------------------
        else if ( structure == "sphere" ) {
            read_msh_2( m_ref, "MESH/SPHERE_3D/sphere_very_fine_Tetra.msh" );
        }
        /// Sphere 3D creuse avec noeud au centre
        /// -------------------------------------
        else if ( structure == "sphere_hollow" ) {
            read_msh_2( m_ref, "MESH/SPHERE_3D/sphere_hollow_very_fine_Tetra.msh" );
        }
        else
            cerr << "structure " << structure << " not implemented..." << endl << endl;
    }
    else
        cerr << "dim " << dim << " not implemented..." << endl << endl;

    if ( m_ref.node_list.size() ) {
        if ( remove_lonely_nodes( m_ref ) )
            cerr << "Des noeuds seuls ont ete retires du maillage de reference du pb direct..." << endl << endl;
        if ( not m_ref.elem_list.size() ) {
            cerr << "Arret brutal, car il n'y a aucun element dans le maillage de reference du pb direct..." << endl << endl;
            throw "Baleinou sous caillou...";
        }
        if ( disp ) {
            cout << "maillage de reference du pb direct :" << endl;
            display_mesh_carac( m_ref );
        }
    }
    else {
        cerr << "Arret brutal, car il n'y a aucun noeud dans le maillage de reference du pb direct..." << endl << endl;
        throw "Ourson sous gravillon...";
    }

}

/// Adaptation du maillage
/// ----------------------
template<class TM,class TF,class T>
bool adapt_mesh( TM &m, const TF &f, const string &structure, const string &method = "EET", unsigned &iter = 1, const unsigned &max_iter = 10, const T &k = 0.75, const bool spread_cut = false, const bool disp = true ) {

    static const unsigned dim = TM::dim;

    if ( iter > max_iter )
        return 1;

    bool res = 0;
    if ( dim == 3 and structure.find("test_specimen") != string::npos ) {
        if ( iter <= max_iter ) {
            if ( method.find("EET") != string::npos )
                res = refinement_if_constraints_or_elem_field_sup( m, f, theta_elem_EET_DM(), k, spread_cut );
            else if ( method.find("SPET") != string::npos )
                res = refinement_if_constraints_or_elem_field_sup( m, f, theta_elem_SPET_DM(), k, spread_cut );
            else if ( method.find("EESPT") != string::npos )
                res = refinement_if_constraints_or_elem_field_sup( m, f, theta_elem_EESPT_DM(), k, spread_cut );
            if ( disp ) {
                cout << "Adaptation de maillage #" << iter << endl;
                cout << "--------------------------" << endl << endl;
                cout << "raffinement du maillage autour des noeuds contraints (appartenant a un bord de Dirichlet) :" << endl;
                cout << "nb de noeuds constraints = " << f.nb_constrained_nodes() << endl;
                cout << "raffinement du maillage autour des elements contribuant majoritairement a l'erreur estimee :" << endl;
                cout << "rapport maximal entre la contribution elementaire au carre a l'erreur estimee et la contribution elementaire maximale au carre : k = " << k << endl;
                cout << "propagation du raffinement au reste du maillage = " << spread_cut << endl << endl;
            }
        }
        else {
            if ( method.find("EET") != string::npos )
                res = refinement_if_elem_field_sup( m, theta_elem_EET_DM(), k, spread_cut );
            else if ( method.find("SPET") != string::npos )
                res = refinement_if_elem_field_sup( m, theta_elem_SPET_DM(), k, spread_cut );
            else if ( method.find("EESPT") != string::npos )
                res = refinement_if_elem_field_sup( m, theta_elem_EESPT_DM(), k, spread_cut );
            if ( disp ) {
                cout << "Adaptation de maillage #" << iter << endl;
                cout << "--------------------------" << endl << endl;
                cout << "raffinement du maillage autour des elements contribuant majoritairement a l'erreur estimee :" << endl;
                cout << "rapport maximal entre la contribution elementaire au carre a l'erreur estimee et la contribution elementaire maximale au carre : k = " << k << endl;
                cout << "propagation du raffinement au reste du maillage = " << spread_cut << endl << endl;
            }
        }

//        res = refinement_if_constraints_or_nodal_field_sup( m, f, ExtractDM< theta_nodal_DM >(), k, spread_cut );
//        if ( disp ) {
//            cout << "Adaptation de maillage #" << iter << endl;
//            cout << "--------------------------" << endl << endl;
//            cout << "raffinement du maillage autour des noeuds contraints (appartenant a un bord de Dirichlet) :" << endl;
//            cout << "nb de noeuds constraints = " << f.nb_constrained_nodes() << endl;
//            cout << "raffinement du maillage autour des noeuds contribuant majoritairement a l'erreur estimee lissee :" << endl;
//            cout << "rapport maximal entre la contribution nodale au carre a l'erreur estimee et la contribution nodale maximale au carre : k = " << k << endl;
//            cout << "propagation du raffinement au reste du maillage = " << spread_cut << endl << endl;
//        }
    }
    else {
        if ( method.find("EET") != string::npos )
            res = refinement_if_elem_field_sup( m, theta_elem_EET_DM(), k, spread_cut );
        else if ( method.find("SPET") != string::npos )
            res = refinement_if_elem_field_sup( m, theta_elem_SPET_DM(), k, spread_cut );
        else if ( method.find("EESPT") != string::npos )
            res = refinement_if_elem_field_sup( m, theta_elem_EESPT_DM(), k, spread_cut );
        if ( disp ) {
            cout << "Adaptation de maillage #" << iter << endl;
            cout << "--------------------------" << endl << endl;
            cout << "raffinement du maillage autour des elements contribuant majoritairement a l'erreur estimee :" << endl;
            cout << "rapport maximal entre la contribution elementaire au carre a l'erreur estimee et la contribution elementaire maximale au carre : k = " << k << endl;
            cout << "propagation du raffinement au reste du maillage = " << spread_cut << endl << endl;
        }

//        res = refinement_if_nodal_field_sup( m, ExtractDM< theta_nodal_DM >(), k, spread_cut );
//        if ( disp ) {
//            cout << "Adaptation de maillage #" << iter << endl;
//            cout << "--------------------------" << endl << endl;
//            cout << "raffinement du maillage autour des noeuds contribuant majoritairement a l'erreur estimee lissee :" << endl;
//            cout << "rapport maximal entre la contribution nodale au carre a l'erreur estimee et la contribution nodale maximale au carre : k = " << k << endl;
//            cout << "propagation du raffinement au reste du maillage = " << spread_cut << endl << endl;
//        }
    }

    if ( not res )
        cerr << "Pas d'adaptation de maillage..." << endl << endl;
    if ( disp )
        display_mesh_carac( m );

    return res;
}

/// Creation du maillage adjoint
/// ----------------------------
template<class TM, class T, class Pvec>
void set_mesh_adjoint( TM &m_adjoint, TM &m, const string &interest_quantity, const string &direction_extractor, const bool &want_local_refinement, const T &l_min, const T &k, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list, Vec<unsigned> &elem_list_adjoint, const unsigned &node, unsigned &node_adjoint, const Pvec &pos_interest_quantity, const Pvec &pos_crack_tip, const T &radius_Ri, const T &radius_Re, const bool &spread_cut, const bool &want_local_enrichment, const unsigned &nb_layers_nodes_enrichment, Vec<unsigned> &elem_list_adjoint_enrichment_zone_1, Vec<unsigned> &elem_list_adjoint_enrichment_zone_2, Vec<unsigned> &face_list_adjoint_enrichment_zone_12, Vec<unsigned> &node_list_adjoint_enrichment, const bool disp = true ) {

    static const unsigned dim = TM::dim;

    m_adjoint = m;

    if ( interest_quantity.find("mean") != string::npos ) {
        for (unsigned n=0;n<elem_list.size();++n) {
            if ( elem_list[ n ] > m.elem_list.size() ) {
                cerr << "Arret brutal, car l'element " << elem_list[ n ] << " definissant la zone d'interet ne fait pas partie de la liste des elements du maillage direct..." << endl << endl;
                throw "Baleinou sous caillou...";
            }
        }
    }
    else if ( interest_quantity.find("pointwise") != string::npos ) {
        if ( pointwise_interest_quantity == "node" and node > m.node_list.size() ) {
            cerr << "Arret brutal, car le noeud " << node << " definissant la zone d'interet ne fait pas partie de la liste des noeuds du maillage direct..." << endl << endl;
            throw "Ourson sous gravillon...";
        }
    }
    
    /// Raffinement local du maillage adjoint
    /// -------------------------------------
    if ( want_local_refinement ) {
        if ( interest_quantity.find("mean") != string::npos ) {
            for (unsigned n=0;n<elem_list.size();++n) {
                while( refinement_point( m_adjoint, l_min, k, center( *m.elem_list[ elem_list[ n ] ] ), spread_cut ) );
                for (unsigned i=0;i<(m.elem_list[ elem_list[ n ] ]->nb_nodes_virtual());++i) {
                    while( refinement_point( m_adjoint, l_min, k, m.elem_list[ elem_list[ n ] ]->node_virtual(i)->pos, spread_cut ) );
                }
            }
        }
        else if ( interest_quantity.find("pointwise") != string::npos ) {
            Pvec pos;
            if ( pointwise_interest_quantity == "node" )
                pos = m.node_list[ node ].pos;
            else if ( pointwise_interest_quantity == "pos" )
                pos = pos_interest_quantity;
            else {
                cerr << "Arret brutal, car la definition de la quantite d'interet ponctuelle " << interest_quantity << " n'est pas implementee..." << endl << endl;
                throw "Ourson sous gravillon...";
            }
            while( refinement_point( m_adjoint, l_min, k, pos, spread_cut ) );
        }
        else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
            while( refinement_point( m_adjoint, l_min, k, pos_crack_tip, spread_cut ) or refinement_circle( m_adjoint, l_min, k, pos_crack_tip, radius_Ri, spread_cut ) or refinement_circle( m_adjoint, l_min, k, pos_crack_tip, radius_Re, spread_cut ) );
        }
        else
            divide( m_adjoint );

        if ( disp ) {
            cout << "raffinement local du maillage adjoint :" << endl;
            cout << "longueur minimale des cotes des elements du maillage adjoint : l_min = " << l_min << endl;
            cout << "coefficient d'augmentation de la longueur maximale des cotes des elements : k = " << k << endl;
            cout << "propagation du raffinement au reste du maillage = " << spread_cut << endl << endl;
        }
    }
    
    /// Zone d'interet du maillage adjoint
    /// ----------------------------------
    if ( interest_quantity.find("mean") != string::npos )
        apply( m_adjoint.elem_list, Construct_Elem_List_Ref(), m, elem_list, elem_list_adjoint );
    else if ( interest_quantity.find("pointwise") != string::npos )
        apply( m_adjoint.node_list, Construct_Node_Ref(), m, node, node_adjoint );
    
    /// Enrichissement local du maillage adjoint : definition des zones d'enrichissement
    /// --------------------------------------------------------------------------------
    if ( want_local_enrichment ) {
        m.update_node_neighbours();
        m.update_node_parents();
        Vec<unsigned> node_list_enrichment;
        Vec<unsigned> elem_list_enrichment_zone_1;
        /// Definition de la zone d'enrichissement zone_1 du pb direct
        /// ----------------------------------------------------------
        if ( interest_quantity.find("mean") != string::npos ) {
            for (unsigned n=0;n<elem_list.size();++n) {
                unsigned node_layer_cpt_enrichment = 1;
                elem_list_enrichment_zone_1.push_back( elem_list[ n ] );
                for (unsigned i=0;i<(m.elem_list[ elem_list[ n ] ]->nb_nodes_virtual());++i)
                    node_list_enrichment.push_back( m.elem_list[ elem_list[ n ] ]->node_virtual(i)->number );
                while ( node_layer_cpt_enrichment < nb_layers_nodes_enrichment ) {
                    Vec<unsigned> node_list_enrichment_tmp = node_list_enrichment;
                    for (unsigned i=0;i<node_list_enrichment_tmp.size();++i) {
                        for (unsigned j=0;j<m.get_node_neighbours( node_list_enrichment_tmp[ i ] ).size();++j) {
                            if ( not find( node_list_enrichment, _1 == m.get_node_neighbours( node_list_enrichment_tmp[ i ] )[ j ]->number ) )
                                node_list_enrichment.push_back( m.get_node_neighbours( node_list_enrichment[ i ] )[ j ]->number );
                        }
                        for (unsigned p=0;p<m.get_node_parents( node_list_enrichment_tmp[ i ] ).size();++p) {
                            if ( not find( elem_list_enrichment_zone_1, _1 == m.get_node_parents( node_list_enrichment_tmp[ i ] )[ p ]->number ) )
                                elem_list_enrichment_zone_1.push_back( m.get_node_parents( node_list_enrichment_tmp[ i ] )[ p ]->number );
                        }
                    }
                    node_layer_cpt_enrichment++;
                }
            }
        }
        else if ( interest_quantity.find("pointwise") != string::npos ) {
            if ( pointwise_interest_quantity == "node" ) {
                unsigned node_layer_cpt_enrichment = 1;
                for (unsigned n=0;n<m.get_node_parents( node ).size();++n)
                    elem_list_enrichment_zone_1.push_back( m.get_node_parents( node )[ n ]->number );
                node_list_enrichment.push_back( node );
                for (unsigned j=0;j<m.get_node_neighbours( node ).size();++j)
                    node_list_enrichment.push_back( m.get_node_neighbours( node )[ j ]->number );
                while ( node_layer_cpt_enrichment < nb_layers_nodes_enrichment ) {
                    Vec<unsigned> node_list_enrichment_tmp = node_list_enrichment;
                    for (unsigned j=0;j<node_list_enrichment_tmp.size();++j) {
                        for (unsigned k=0;k<m.get_node_neighbours( node_list_enrichment_tmp[ j ] ).size();++k) {
                            if ( not find( node_list_enrichment, _1 == m.get_node_neighbours( node_list_enrichment_tmp[ j ] )[ k ]->number ) )
                                node_list_enrichment.push_back( m.get_node_neighbours( node_list_enrichment[ j ] )[ k ]->number );
                        }
                        for (unsigned p=0;p<m.get_node_parents( node_list_enrichment_tmp[ j ] ).size();++p) {
                            if ( not find( elem_list_enrichment_zone_1, _1 == m.get_node_parents( node_list_enrichment_tmp[ j ] )[ p ]->number ) )
                                elem_list_enrichment_zone_1.push_back( m.get_node_parents( node_list_enrichment_tmp[ j ] )[ p ]->number );
                        }
                    }
                    node_layer_cpt_enrichment++;
                }
            }
            else if ( pointwise_interest_quantity == "pos" ) {
                unsigned num_elem = 0;
                for (unsigned n=0;n<m.elem_list.size();++n) {
                    if ( length( center( *m.elem_list[ n ] ) - pos_interest_quantity ) < length( center( *m.elem_list[ num_elem ] ) - pos_interest_quantity ) )
                        num_elem = n;
                }
                unsigned node_layer_cpt_enrichment = 1;
                elem_list_enrichment_zone_1.push_back( num_elem );
                for (unsigned i=0;i<(m.elem_list[ num_elem ]->nb_nodes_virtual());++i)
                    node_list_enrichment.push_back( m.elem_list[ num_elem ]->node_virtual(i)->number );
                while ( node_layer_cpt_enrichment < nb_layers_nodes_enrichment ) {
                    Vec<unsigned> node_list_enrichment_tmp = node_list_enrichment;
                    for (unsigned i=0;i<node_list_enrichment_tmp.size();++i) {
                        for (unsigned j=0;j<m.get_node_neighbours( node_list_enrichment_tmp[ i ] ).size();++j) {
                            if ( not find( node_list_enrichment, _1 == m.get_node_neighbours( node_list_enrichment_tmp[ i ] )[ j ]->number ) )
                                node_list_enrichment.push_back( m.get_node_neighbours( node_list_enrichment[ i ] )[ j ]->number );
                        }
                        for (unsigned p=0;p<m.get_node_parents( node_list_enrichment_tmp[ i ] ).size();++p) {
                            if ( not find( elem_list_enrichment_zone_1, _1 == m.get_node_parents( node_list_enrichment_tmp[ i ] )[ p ]->number ) )
                                elem_list_enrichment_zone_1.push_back( m.get_node_parents( node_list_enrichment_tmp[ i ] )[ p ]->number );
                        }
                    }
                    node_layer_cpt_enrichment++;
                }
            }
            else
                cerr << "pointwise quantity of interest " << interest_quantity << " defined by " << pointwise_interest_quantity << " not implemented..." << endl << endl;
        }
        /// Definition de la zone d'enrichissement zone_2 du pb direct
        /// ----------------------------------------------------------
        Vec<unsigned> elem_list_enrichment_zone_2;
        for (unsigned i=0;i<node_list_enrichment.size();++i) {
            for (unsigned n=0;n<m.get_node_parents( node_list_enrichment[ i ] ).size();++n) {
                if ( not find( elem_list_enrichment_zone_1, _1 == m.get_node_parents( node_list_enrichment[ i ] )[ n ]->number ) and not find( elem_list_enrichment_zone_2, _1 == m.get_node_parents( node_list_enrichment[ i ] )[ n ]->number ) )
                    elem_list_enrichment_zone_2.push_back( m.get_node_parents( node_list_enrichment[ i ] )[ n ]->number );
            }
        }
        /// Definition de la zone d'enrichissement zone_12 du pb direct
        /// -----------------------------------------------------------
        Vec<unsigned> child_cpt;
        Vec< Vec<unsigned> > child_list;
        construct_child( m, child_cpt, child_list );
        Vec<unsigned> face_list_enrichment_zone_1;
        for (unsigned n=0;n<elem_list_enrichment_zone_1.size();++n) {
            for (unsigned k=0;k<child_cpt[ elem_list_enrichment_zone_1[ n ] ];++k) {
                if ( not find( face_list_enrichment_zone_1, _1 == child_list[ elem_list_enrichment_zone_1[ n ] ][ k ] ) )
                    face_list_enrichment_zone_1.push_back( child_list[ elem_list_enrichment_zone_1[ n ] ][ k ] );
            }
        }
        Vec<unsigned> face_list_enrichment_zone_12;
        for (unsigned n=0;n<elem_list_enrichment_zone_2.size();++n) {
            for (unsigned k=0;k<child_cpt[ elem_list_enrichment_zone_2[ n ] ];++k) {
                if ( find( face_list_enrichment_zone_1, _1 == child_list[ elem_list_enrichment_zone_2[ n ] ][ k ] ) )
                    face_list_enrichment_zone_12.push_back( child_list[ elem_list_enrichment_zone_2[ n ] ][ k ] );
            }
        }
        /// Definition des level set du pb direct
        /// -------------------------------------
        for (unsigned i=0;i<node_list_enrichment.size();++i)
            m.node_list[ node_list_enrichment[ i ] ].phi_nodal_handbook = 1.;
        for (unsigned n=0;n<elem_list_enrichment_zone_1.size();++n)
            m.elem_list[ elem_list_enrichment_zone_1[ n ] ]->set_field( "phi_elem_handbook_zone_1", 1 );
        for (unsigned n=0;n<elem_list_enrichment_zone_2.size();++n)
            m.elem_list[ elem_list_enrichment_zone_2[ n ] ]->set_field( "phi_elem_handbook_zone_2", 1 );
        for (unsigned k=0;k<face_list_enrichment_zone_12.size();++k)
            m.sub_mesh(Number<1>()).elem_list[ face_list_enrichment_zone_12[ k ] ]->set_field( "phi_surf_handbook_zone_12", 1 );
        
        /// Definition de la zone d'enrichissement zone_1 du pb adjoint
        /// -----------------------------------------------------------
        apply( m_adjoint.elem_list, Construct_Elem_List_Ref(), m, elem_list_enrichment_zone_1, elem_list_adjoint_enrichment_zone_1 );
        for (unsigned n=0;n<elem_list_adjoint_enrichment_zone_1.size();++n) {
            for (unsigned i=0;i<(m_adjoint.elem_list[ elem_list_adjoint_enrichment_zone_1[ n ] ]->nb_nodes_virtual());++i) {
                if ( not find( node_list_adjoint_enrichment, _1 == m_adjoint.elem_list[ elem_list_adjoint_enrichment_zone_1[ n ] ]->node_virtual(i)->number ) )
                    node_list_adjoint_enrichment.push_back( m_adjoint.elem_list[ elem_list_adjoint_enrichment_zone_1[ n ] ]->node_virtual(i)->number );
            }
        }
        /// Definition de la zone d'enrichissement zone_2 du pb adjoint
        /// -----------------------------------------------------------
        m_adjoint.update_node_parents();
        for (unsigned i=0;i<node_list_adjoint_enrichment.size();++i) {
            for (unsigned n=0;n<m_adjoint.get_node_parents( node_list_adjoint_enrichment[ i ] ).size();++n) {
                if ( not find( elem_list_adjoint_enrichment_zone_1, _1 == m_adjoint.get_node_parents( node_list_adjoint_enrichment[ i ] )[ n ]->number ) and not find( elem_list_adjoint_enrichment_zone_2, _1 == m_adjoint.get_node_parents( node_list_adjoint_enrichment[ i ] )[ n ]->number ) )
                    elem_list_adjoint_enrichment_zone_2.push_back( m_adjoint.get_node_parents( node_list_adjoint_enrichment[ i ] )[ n ]->number );
            }
        }
        /// Definition de la zone d'enrichissement zone_12 du pb adjoint
        /// ------------------------------------------------------------
        Vec<unsigned> child_cpt_adjoint;
        Vec< Vec<unsigned> > child_list_adjoint;
        construct_child( m_adjoint, child_cpt_adjoint, child_list_adjoint );
        Vec<unsigned> face_list_adjoint_enrichment_zone_1;
        for (unsigned n=0;n<elem_list_adjoint_enrichment_zone_1.size();++n) {
            for (unsigned k=0;k<child_cpt_adjoint[ elem_list_adjoint_enrichment_zone_1[ n ] ];++k) {
                if ( not find( face_list_adjoint_enrichment_zone_1, _1 == child_list_adjoint[ elem_list_adjoint_enrichment_zone_1[ n ] ][ k ] ) )
                    face_list_adjoint_enrichment_zone_1.push_back( child_list_adjoint[ elem_list_adjoint_enrichment_zone_1[ n ] ][ k ] );
            }
        }
        for (unsigned n=0;n<elem_list_adjoint_enrichment_zone_2.size();++n) {
            for (unsigned k=0;k<child_cpt_adjoint[ elem_list_adjoint_enrichment_zone_2[ n ] ];++k) {
                if ( find( face_list_adjoint_enrichment_zone_1, _1 == child_list_adjoint[ elem_list_adjoint_enrichment_zone_2[ n ] ][ k ] ) )
                    face_list_adjoint_enrichment_zone_12.push_back( child_list_adjoint[ elem_list_adjoint_enrichment_zone_2[ n ] ][ k ] );
            }
        }
        /// Definition des level sets du pb adjoint
        /// ---------------------------------------
        if ( interest_quantity == "pointwise_dep" ) {
            switch ( dim ) {
            case 1 :
                if ( direction_extractor == "x" )
                    m_adjoint.phi_handbook_pointwise_force_in_infinite_domain[ 0 ] = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            case 2 :
                if ( direction_extractor == "x" )
                    m_adjoint.phi_handbook_pointwise_force_in_infinite_domain[ 0 ] = 1;
                else if ( direction_extractor == "y" )
                    m_adjoint.phi_handbook_pointwise_force_in_infinite_domain[ 1 ] = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            case 3 :
                if ( direction_extractor == "x" )
                    m_adjoint.phi_handbook_pointwise_force_in_infinite_domain[ 0 ] = 1;
                else if ( direction_extractor == "y" )
                    m_adjoint.phi_handbook_pointwise_force_in_infinite_domain[ 1 ] = 1;
                else if ( direction_extractor == "z" )
                    m_adjoint.phi_handbook_pointwise_force_in_infinite_domain[ 2 ] = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            }
        }
        else if ( interest_quantity == "pointwise_epsilon" ) {
            switch ( dim ) {
            case 1 :
                if ( direction_extractor == "xx" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 0, 0 ) = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            case 2 :
                if ( direction_extractor == "xx" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 0, 0 ) = 1;
                else if ( direction_extractor == "xy" or direction_extractor == "yx" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 0, 1 ) = 1;
                else if ( direction_extractor == "yy" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 1, 1 ) = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            case 3 :
                if ( direction_extractor == "xx" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 0, 0 ) = 1;
                else if ( direction_extractor == "xy" or direction_extractor == "yx" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 0, 1 ) = 1;
                else if ( direction_extractor == "yy" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 1, 1 ) = 1;
                else if ( direction_extractor == "xz" or direction_extractor == "zx" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 0, 2 ) = 1;
                else if ( direction_extractor == "yz" or direction_extractor == "zy" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 1, 2 ) = 1;
                else if ( direction_extractor == "zz" )
                    m_adjoint.phi_handbook_pointwise_pre_sigma_in_infinite_domain( 2, 2 ) = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            }
        }
        else if ( interest_quantity == "pointwise_sigma" ) {
            switch ( dim ) {
            case 1 :
                if ( direction_extractor == "xx" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 0, 0 ) = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            case 2 :
                if ( direction_extractor == "xx" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 0, 0 ) = 1;
                else if ( direction_extractor == "xy" or direction_extractor == "yx" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 0, 1 ) = 1;
                else if ( direction_extractor == "yy" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 1, 1 ) = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            case 3 :
                if ( direction_extractor == "xx" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 0, 0 ) = 1;
                else if ( direction_extractor == "xy" or direction_extractor == "yx" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 0, 1 ) = 1;
                else if ( direction_extractor == "yy" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 1, 1 ) = 1;
                else if ( direction_extractor == "xz" or direction_extractor == "zx" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 0, 2 ) = 1;
                else if ( direction_extractor == "yz" or direction_extractor == "zy" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 1, 2 ) = 1;
                else if ( direction_extractor == "zz" )
                    m_adjoint.phi_handbook_pointwise_pre_epsilon_in_infinite_domain( 2, 2 ) = 1;
                else {
                    cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << interest_quantity << " en dimension " << dim << " n'est pas implementee..." << endl << endl;
                    throw "Anguille sous coquille...";
                }
                break;
            }
        }
        for (unsigned i=0;i<node_list_adjoint_enrichment.size();++i)
            m_adjoint.node_list[ node_list_adjoint_enrichment[ i ] ].phi_nodal_handbook = 1.;
        for (unsigned n=0;n<elem_list_adjoint_enrichment_zone_1.size();++n)
            m_adjoint.elem_list[ elem_list_adjoint_enrichment_zone_1[ n ] ]->set_field( "phi_elem_handbook_zone_1", 1 );
        for (unsigned n=0;n<elem_list_adjoint_enrichment_zone_2.size();++n)
            m_adjoint.elem_list[ elem_list_adjoint_enrichment_zone_2[ n ] ]->set_field( "phi_elem_handbook_zone_2", 1 );
        for (unsigned k=0;k<face_list_adjoint_enrichment_zone_12.size();++k)
            m_adjoint.sub_mesh(Number<1>()).elem_list[ face_list_adjoint_enrichment_zone_12[ k ] ]->set_field( "phi_surf_handbook_zone_12", 1 );
        LMT::sort( node_list_adjoint_enrichment );
        LMT::sort( elem_list_adjoint_enrichment_zone_1 );
        LMT::sort( elem_list_adjoint_enrichment_zone_2 );
        LMT::sort( face_list_adjoint_enrichment_zone_12 );

        if ( disp ) {
            cout << "enrichissement local avec fonctions handbook : " << endl;
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
    }

    if ( m_adjoint.node_list.size() ) {
        if ( remove_lonely_nodes( m_adjoint ) )
            cerr << "Des noeuds seuls ont ete retires du maillage du pb adjoint..." << endl << endl;
        if ( not m_adjoint.elem_list.size() ) {
            cerr << "Arret brutal, car il n'y a aucun element dans le maillage du pb adjoint..." << endl << endl;
            throw "Baleinou sous caillou...";
        }
        if ( disp ) {
            cout << "maillage du pb adjoint :" << endl;
            display_mesh_carac( m_adjoint );
        }
    }
    else {
        cerr << "Arret brutal, car il n'y a aucun noeud dans le maillage du pb adjoint..." << endl << endl;
        throw "Ourson sous gravillon...";
    }

}

/// Creation du maillage de reference pour le calcul de la quantite d'interet quasi-exacte
/// --------------------------------------------------------------------------------------
template<class TM, class T, class Pvec>
void set_mesh_local_ref( TM &m_ref, TM &m, const unsigned &refinement_level_ref, const string &interest_quantity, const Vec<unsigned> &elem_list, Vec<unsigned> &elem_list_ref, const unsigned &node, unsigned &node_ref, const Pvec &pos_crack_tip, const T &radius_Ri, const T &radius_Re, const bool spread_cut = false, const bool disp = true ) {

    static const unsigned dim = TM::dim;

    m_ref = m;

    for (unsigned n=0;n<refinement_level_ref;++n)
        divide( m_ref );

//    if ( interest_quantity.find("mean") != string::npos ) {
//        for (unsigned n=0;n<elem_list.size();++n) {
//            while( refinement_point( m_ref, 0.05, 0.05, center( *m.elem_list[ elem_list[ n ] ] ), spread_cut ) );
//            for (unsigned i=0;i<(m.elem_list[ elem_list[ n ] ]->nb_nodes_virtual());++i) {
//                while( refinement_point( m_ref, 0.05, 0.05, m.elem_list[ elem_list[ n ] ]->node_virtual(i)->pos, spread_cut ) );
//            }
//        }
//    }
//    else if ( interest_quantity.find("pointwise") != string::npos ) {
//        while( refinement_point( m_ref, 0.5, 0.5, m.node_list[ node ].pos, spread_cut ) );
//    }
//    else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
//        while( refinement_point( m_ref, 0.2, 0.2, pos_crack_tip, spread_cut ) or refinement_circle( m_ref, 0.2, 0.2, pos_crack_tip, radius_Ri, spread_cut ) or refinement_circle( m_ref, 0.2, 0.2, pos_crack_tip, radius_Re, spread_cut )  );
//    }

    if ( m_ref.node_list.size() ) {
        if ( remove_lonely_nodes( m_ref ) )
            cerr << "Des noeuds seuls ont ete retires du maillage de reference du pb direct..." << endl << endl;
        if ( not m_ref.elem_list.size() ) {
            cerr << "Arret brutal, car il n'y a aucun element dans le maillage de reference du pb direct..." << endl << endl;
            throw "Baleinou sous caillou...";
        }
        if ( disp ) {
            cout << "maillage de reference du pb direct :" << endl;
            display_mesh_carac( m_ref );
        }
    }
    else {
        cerr << "Arret brutal, car il n'y a aucun noeud dans le maillage de reference du pb direct..." << endl << endl;
        throw "Ourson sous gravillon...";
    }
    
    if ( interest_quantity.find("mean") != string::npos )
        apply( m_ref.elem_list, Construct_Elem_List_Ref(), m, elem_list, elem_list_ref );
    else if ( interest_quantity.find("pointwise") != string::npos )
        apply( m_ref.node_list, Construct_Node_Ref(), m, node, node_ref );
}

/// Decoupe du maillage autour d'un domaine repere par son centre domain_center et sa taille domain_length
/// ------------------------------------------------------------------------------------------------------
template<class TM, class T, class Pvec>
void set_mesh_domain( TM &m_domain, TM &m, const string &shape, const T &k, const Vec<T> &domain_length, const Pvec &domain_center, const bool spread_cut = false, const bool disp = false ) {
    
    if ( shape.find("circle") != string::npos or shape.find("sphere") != string::npos ) {
        for (unsigned i=0;i<m.node_list.size();++i)
            m.node_list[i].phi_domain = k * domain_length[ 0 ] - length( m.node_list[ i ].pos - domain_center );
        m_domain = m;
        level_set_cut( m_domain, ExtractDM< phi_domain_DM >(), spread_cut );
    }
    if ( disp ) {
        cout << "maillage du sous-domaine :" << endl;
        display_mesh_carac( m_domain );
    }
}

/// Creation du maillage de la couronne entourant la pointe de fissure pour la quantite d'interet SIF
/// -------------------------------------------------------------------------------------------------
template<class TM, class T, class Pvec>
void set_mesh_crown( TM &m_crown, TM &m, const Pvec &pos_crack_tip, const T &radius_Ri, const T &radius_Re, const bool spread_cut = false, const bool disp = false ) {

    for (unsigned i=0;i<m.node_list.size();++i) {
        m.node_list[i].phi_crown_int = length( m.node_list[ i ].pos - pos_crack_tip ) - radius_Ri;
        m.node_list[i].phi_crown_ext = radius_Re - length( m.node_list[ i ].pos - pos_crack_tip );
    }
    m_crown = m;
    level_set_cut( m_crown, ExtractDM< phi_crown_int_DM >(), spread_cut );
    level_set_cut( m_crown, ExtractDM< phi_crown_ext_DM >(), spread_cut );
    if ( disp ) {
        cout << "maillage de la couronne entourant la pointe de fissure :" << endl;
        display_mesh_carac( m_crown );
    }
}

/// Creation du maillage parametrique
/// ---------------------------------
template<class TM_param, class T, class TT>
void set_mesh_param( TM_param &m_param, const T &min_param, const T &max_param, const TT &nb_points, const bool disp = false ) {
    typedef typename TM_param::Pvec Pvec_param;
    make_rect( m_param, Bar(), Pvec_param( min_param ), Pvec_param( max_param ), Pvec_param( nb_points ) );
    if ( disp ) {
        cout << "maillage du pb en paramètre :" << endl;
        display_mesh_carac( m_param );
    }
}

#endif // Mesh_h
