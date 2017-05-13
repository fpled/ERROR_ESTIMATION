//
// C++ Interface: Construct_standard_force_fluxes_EESPT_PGD
//
// Description: construction standard des densites d'effort par la methode EESPT dans le cadre des methodes PGD
//
//
// Author: Pled Florent <florent.pled@univ-paris-est.fr>, (C) 2017
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_standard_force_fluxes_EESPT_PGD_h
#define Construct_standard_force_fluxes_EESPT_PGD_h

#include "EESPT.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"
#include "../LMT/include/containers/matumfpack.h"

using namespace LMT;
using namespace std;

/// Construction standard des densites d'effort par la methode EESPT
/// ----------------------------------------------------------------
template<class TM, class TF, class T, class TVV, class TVVV>
void construct_standard_force_fluxes_EESPT_PGD( TM &m, const TF &f, const string &pb, const unsigned &cost_function, const bool &enhancement, const Vec<bool> &flag_face_enh, const string &solver_minimisation, const T &penalty_val_N, TVVV &force_fluxes, const TVV &dep, const Vec< Vec<unsigned> > &elem_group, const bool want_local_enrichment = false, const bool verif_solver_minimisation = false, const T tol_solver_minimisation = 1e-6, const bool disp = false ) {
    
    static const unsigned dim = TM::dim;
    
    typedef Mat<T, Gen<>, SparseLine<> > TMatGenSparse;
    
    TicToc t_force_fluxes_std;
    t_force_fluxes_std.start();
    
    Vec<unsigned> node_cpt_face;
    Vec< Vec<unsigned> > node_list_face;
    construct_nodes_connected_to_face( m, node_cpt_face, node_list_face );
    
    Vec<bool> correspondance_node_to_vertex_node;
    Vec<unsigned> connect_node_to_vertex_node;
    unsigned nb_vertex_nodes = match_node_to_vertex_node( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node );
    
    Vec<unsigned> face_cpt_vertex_node;
    Vec< Vec<unsigned> > face_list_vertex_node;
    construct_faces_connected_to_vertex_node( m, nb_vertex_nodes, correspondance_node_to_vertex_node, connect_node_to_vertex_node, face_cpt_vertex_node, face_list_vertex_node );
    
    Vec<unsigned> vertex_node_cpt_elem;
    Vec< Vec<unsigned> > vertex_node_list_elem;
    construct_vertex_nodes_connected_to_elem( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, vertex_node_cpt_elem, vertex_node_list_elem );
    
    vertex_node_cpt_elem.free();
    
    Vec<unsigned> vertex_node_cpt_face;
    Vec< Vec<unsigned> > vertex_node_list_face;
    construct_vertex_nodes_connected_to_face( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node, vertex_node_cpt_face, vertex_node_list_face );
    
    correspondance_node_to_vertex_node.free();
    
    Vec<unsigned> elem_cpt_node;
    Vec< Vec<unsigned> > elem_list_node;
    construct_elems_connected_to_node( m, elem_cpt_node, elem_list_node );
    
    elem_list_node.free();
    
    Vec< Vec<unsigned> > face_type;
    construct_face_type( m, f, face_type );
    
    if ( disp ) {
        cout << "Construction standard des densites d'effort" << endl;
        cout << "-------------------------------------------" << endl << endl;
    }
    
    /// Reperage pour chaque face k et chaque direction d de l'indice de debut de ligne dans les matrices A[ j ][ d ] : face_ind[ k ][ d ]
    /// Calcul du nb de lignes de la matrice A[ j ][ d ] : nb_unk[ j ][ d ]
    /// ----------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Calcul du vecteur nb_unk" << endl << endl;
    
    Vec< Vec<unsigned> > nb_unk;
    nb_unk.resize( nb_vertex_nodes );
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        nb_unk[ j ].resize( dim );
        nb_unk[ j ].set( 0 );
    }
    
    Vec< Vec< Vec<unsigned> > > face_ind;
    face_ind.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        face_ind[ k ].resize( dim );
    }
    
    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, Calcul_Face_Ind_EESPT(), face_ind, nb_unk, connect_node_to_vertex_node ); // nb_unk[ j ][ d ] contient le nb de lignes de la matrice A[ j ][ d ] associee au j eme noeud sommet dans la direction d // face_ind[ k ][ d ][ numero_noeud_sommet_dans_face (0 ou 1) ] = indice de debut de ligne dans la matrice A[ j ][ d ] pour chaque face k du maillage et chaque la direction d
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "nb d'inconnues lambda_F_hat associees au noeud sommet " << j << " dans la direction " << d << " = " << nb_unk[ j ][ d ] << endl;
            }
        }
        cout << endl;
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "indice de debut de ligne de la face " << k << " dans les matrices A[ noeud sommet connecte a la face " << k << " ] dans la direction " << d << " = " << face_ind[ k ][ d ] << endl;
            }
        }
        cout << endl;
    }
    
    /// Reperage pour chaque element n et chaque direction d de l'indice de debut de colonne dans les matrices A[ j ][ d ] et de debut de ligne dans les vecteurs R[ j ][ d ] : vertex_nodal_ind[ n ][ d ]
    /// Calcul du nb de colonnes de la matrice A[ j ][ d ] : nb_eq[ j ][ d ]
    /// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Calcul du vecteur nb_eq" << endl << endl;
    
    Vec< Vec<unsigned> > nb_eq;
    nb_eq.resize( nb_vertex_nodes );
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        nb_eq[ j ].resize( dim );
        nb_eq[ j ].set( 0 );
    }
    
    Vec< Vec< Vec<unsigned> > > vertex_nodal_ind;
    vertex_nodal_ind.resize( m.elem_list.size() );
    for (unsigned n=0;n<m.elem_list.size();++n) {
        vertex_nodal_ind[ n ].resize( dim );
    }
    
    apply( m.elem_list, Calcul_Vertex_Nodal_Ind_EESPT(), vertex_nodal_ind, nb_eq, connect_node_to_vertex_node ); // nb_eq[ j ][ d ] contient le nb de colonnes de la matrice A[ j ][ d ] associee au j eme noeud sommet dans la direction d // vertex_nodal_ind[ n ][ d ][ numero_noeud_sommet_dans_elem (0 ou 1 ou 2) ] = indice de debut de colonne dans la matrice A[ j ][ d ] pour chaque element n du maillage et chaque direction d
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "nb d'equations associees au noeud sommet " << j << " dans la direction " << d << " = " << nb_eq[ j ][ d ] << endl << endl;
            }
        }
        for (unsigned n=0;n<m.elem_list.size();++n) {
            for (unsigned d=0;d<dim;++d) {
                cout << "indice de debut de colonne de l'element " << n << " dans les matrices A[ noeud sommet connecte a l'element " << n << " dans la direction " << d << " ] = " << vertex_nodal_ind[ n ][ d ] << endl << endl;
            }
        }
    }
    
    /// Construction des matrices A[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// --------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices A" << endl << endl;
    
    Vec< Vec< TMatGenSparse > > A;
    A.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        A[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            A[ j ][ d ].resize( nb_unk[ j ][ d ], nb_eq[ j ][ d ] );
        }
    }
    
    Calcul_Vertex_Nodal_Matrix_A calcul_vertex_nodal_matrix_A;
    calcul_vertex_nodal_matrix_A.face_ind = &face_ind;
    calcul_vertex_nodal_matrix_A.vertex_nodal_ind = &vertex_nodal_ind;
    calcul_vertex_nodal_matrix_A.vertex_node_list_elem = &vertex_node_list_elem;
    calcul_vertex_nodal_matrix_A.node_list_face = &node_list_face;
    
    apply( m.elem_list, calcul_vertex_nodal_matrix_A, m, connect_node_to_vertex_node, A );
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension de la matrice A associee au noeud sommet " << j << " dans la direction " << d << " = ( " << nb_unk[ j ][ d ] << ", " << nb_eq[ j ][ d ] << " )" << endl;
                cout << "matrice A associee au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << A[ j ][ d ] << endl;
            }
        }
    }
    
    
    /// Construction des vecteurs R[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// --------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs R" << endl << endl;
    
    Vec< Vec< Vec<T> > > R;
    R.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        R[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            R[ j ][ d ].resize( nb_eq[ j ][ d ] );
            R[ j ][ d ].set( 0. );
        }
    }
    
    Calcul_Vertex_Nodal_Vector_R_PGD<TVV> calcul_vertex_nodal_vector_R_PGD;
    calcul_vertex_nodal_vector_R_PGD.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_vertex_nodal_vector_R_PGD.vertex_nodal_ind = &vertex_nodal_ind;
    calcul_vertex_nodal_vector_R_PGD.vertex_node_list_elem = &vertex_node_list_elem;
    calcul_vertex_nodal_vector_R_PGD.node_list_face = &node_list_face;
    calcul_vertex_nodal_vector_R_PGD.elem_cpt_node = &elem_cpt_node;
    calcul_vertex_nodal_vector_R_PGD.pb = &pb;
    calcul_vertex_nodal_vector_R_PGD.want_local_enrichment = &want_local_enrichment;
    calcul_vertex_nodal_vector_R_PGD.dep = &dep;
    calcul_vertex_nodal_vector_R_PGD.elem_group = &elem_group;
    
    apply( m.elem_list, calcul_vertex_nodal_vector_R_PGD, m, f, R );
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur R associe au noeud sommet " << j << " dans la direction " << d << " = " << nb_eq[ j ][ d ] << endl;
                cout << "vecteur R associe au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << R[ j ][ d ] << endl << endl;
            }
        }
    }
    
    vertex_node_list_elem.free();
    elem_cpt_node.free();
    
    /// Construction de la liste des noeuds connectes a un noeud sommet : node_list_vertex_node[ j ]
    /// Construction de la liste des noeuds connectes a un noeud sommet et appartenant au bord du patch : edge_node_list_vertex_node[ j ]
    /// ---------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )  {
        cout << "Construction de la liste des noeuds connectes a un noeud sommet" << endl;
        cout << "Construction de la liste des noeuds connectes a un noeud sommet et appartenant au bord du patch" << endl << endl;
    }
    
    Vec< Vec<unsigned> > node_list_vertex_node;
    node_list_vertex_node.resize( nb_vertex_nodes );
    
    Vec< Vec<unsigned> > edge_node_list_vertex_node;
    edge_node_list_vertex_node.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned k=0;k<face_cpt_vertex_node[ j ];++k) {
            for (unsigned i=0;i<node_cpt_face[ face_list_vertex_node[ j ][ k ] ];++i) {
                node_list_vertex_node[ j ].push_back( node_list_face[ face_list_vertex_node[ j ][ k ] ][ i ] );
            }
            if ( face_type[ face_list_vertex_node[ j ][ k ] ][ 0 ] != 0 ) { // si la face consideree est sur le bord delta_Omega
                for (unsigned i=0;i<node_cpt_face[ face_list_vertex_node[ j ][ k ] ];++i) {
                    edge_node_list_vertex_node[ j ].push_back( node_list_face[ face_list_vertex_node[ j ][ k ] ][ i ] );
                }
            }
        }
        remove_doubles( node_list_vertex_node[ j ] );
        remove_doubles( edge_node_list_vertex_node[ j ] );
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            cout << "liste des noeuds connectes au noeud sommet " << j << " = " << node_list_vertex_node[ j ] << endl ;
            cout << "liste des noeuds connectes au noeud sommet " << j << " et appartenant au bord du patch = " << edge_node_list_vertex_node[ j ] << endl << endl;
        }
    }
    
    node_cpt_face.free();
    face_cpt_vertex_node.free();
    face_list_vertex_node.free();
    
    /// Suppression du noyau des matrices A[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// Stockage des inconnues non bloquees : eq_indep[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// Calcul du nb d'equations independantes : nb_eq_indep[ i ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// ---------------------------------------------------------------------------------------------------------------------------
    
    if ( disp ) {
        cout << "Suppression du noyau des matrices A" << endl;
        cout << "Calcul du vecteur nb_eq_indep" << endl << endl;
    }
    
    Vec< Vec<unsigned> > nb_eq_indep;
    nb_eq_indep.resize( nb_vertex_nodes );
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        nb_eq_indep[ j ].resize( dim );
        nb_eq_indep[ j ].set( 0 );
    }
    
    Vec< Vec< Vec<unsigned> > > eq_indep; // liste des inconnues non bloquees avec redondance de certaines inconnues
    eq_indep.resize( nb_vertex_nodes );
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        eq_indep[ j ].resize( dim );
    }
    
    Vec< Vec<bool> > node_flag;
    Vec< Vec<unsigned> > elem_flag;
    node_flag.resize( nb_vertex_nodes );
    elem_flag.resize( nb_vertex_nodes );
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        node_flag[ j ].resize( m.node_list.size() );
        node_flag[ j ].set( 0 );
        elem_flag[ j ].resize( m.node_list.size() );
        elem_flag[ j ].set( 0 );
    }
    
    Remove_Kernel remove_kernel;
    remove_kernel.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    remove_kernel.vertex_nodal_ind = &vertex_nodal_ind;
    remove_kernel.node_list_vertex_node = &node_list_vertex_node;
    remove_kernel.edge_node_list_vertex_node = &edge_node_list_vertex_node;
    remove_kernel.node_flag = &node_flag;
    remove_kernel.elem_flag = &elem_flag;
    
    apply( m.elem_list, remove_kernel, m, eq_indep );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            remove_doubles( eq_indep[ j ][ d ] );
            nb_eq_indep[ j ][ d ] = eq_indep[ j ][ d ].size();
        }
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "nb d'equations associees au noeud sommet " << j << " dans la direction " << d << " avant suppression du noyau = " << nb_eq[ j ][ d ] << endl;
                cout << "nb d'equations supprimees associees au noeud sommet " << j << " dans la direction " << d << " = " << ( nb_eq[ j ][ d ] - nb_eq_indep[ j ][ d ] ) << endl;
                cout << "nb d'equations independantes associees au noeud sommet " << j << " dans la direction " << d << " apres suppression du noyau = " << nb_eq_indep[ j ][ d ] << endl;
                cout << "indices des colonnes dans la matrice A_tilde et des lignes dans le vecteur R_tilde associes au noeud sommet " << j << " dans la direction " << d << " = " << eq_indep[ j ][ d ] << endl << endl;
            }
        }
    }
    
    nb_eq.free();
    vertex_nodal_ind.free();
    node_flag.free();
    elem_flag.free();
    node_list_vertex_node.free();
    edge_node_list_vertex_node.free();
    
    /// Construction des matrices A_tilde[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// --------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices A_tilde" << endl << endl;
    
    Vec< Vec< TMatGenSparse > > A_tilde;
    A_tilde.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        A_tilde[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            A_tilde[ j ][ d ].resize( nb_unk[ j ][ d ], nb_eq_indep[ j ][ d ] );
        }
    }
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            for (unsigned k=0;k<nb_eq_indep[ j ][ d ];++k) {
                A_tilde[ j ][ d ].col( k ) = A[ j ][ d ].col( eq_indep[ j ][ d ][ k ] );
            }
        }
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension de la matrice A_tilde associee au noeud sommet " << j << " dans la direction " << d << " = ( " << nb_unk[ j ][ d ] << ", " << nb_eq_indep[ j ][ d ] << " )" << endl;
                cout << "matrice A_tilde associee au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << A_tilde[ j ][ d ] << endl;
            }
        }
    }
    
    if ( not disp )
        A.free();
    
    /// Construction des vecteurs R_tilde[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// --------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs R_tilde" << endl << endl;
    
    Vec< Vec< Vec<T> > > R_tilde;
    R_tilde.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        R_tilde[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            R_tilde[ j ][ d ].resize( nb_eq_indep[ j ][ d ] );
            R_tilde[ j ][ d ].set( 0. );
        }
    }
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            for (unsigned k=0;k<nb_eq_indep[ j ][ d ];++k) {
                R_tilde[ j ][ d ][ k ] += R[ j ][ d ][ eq_indep[ j ][ d ][ k ] ];
            }
        }
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur R_tilde associe au noeud sommet " << j << " dans la direction " << d << " = " << nb_eq_indep[ j ][ d ] << endl;
                cout << "vecteur R_tilde associe au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << R_tilde[ j ][ d ] << endl << endl;
            }
        }
    }
    
    eq_indep.free();
    if ( not disp )
        R.free();
    
    /// Construction des vecteurs lambda_F[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// ---------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs lambda_F" << endl << endl;
    
    Vec< Vec< Vec<T> > > lambda_F;
    lambda_F.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        lambda_F[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            lambda_F[ j ][ d ].resize( nb_unk[ j ][ d ] );
            lambda_F[ j ][ d ].set( 0. );
        }
    }
    
    Vec<unsigned> degre_p;
    degre_p.resize( m.elem_list.size() );
    degre_p.set( 0 );
    apply( m.elem_list, Get_Elem_Degree(), degre_p );
    
    if ( degre_p[ 0 ] == 1 ) {
        
        /// Cas p = 1
        /// ---------
        
        if ( disp )
            cout << "Cas : degre p = 1" << endl << endl;
        
        Vec< Vec< Vec<T> > > lambda_F_face;
        lambda_F_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );
        
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            lambda_F_face[ k ].resize( dim );
            for (unsigned d=0;d<dim;++d) {
                lambda_F_face[ k ][ d ].resize( vertex_node_cpt_face[ k ] * m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
                lambda_F_face[ k ][ d ].set( 0. );
            }
        }
        
        /// Construction des matrices B[ k ][ d ] pour chaque face k du maillage et chaque direction d
        /// ------------------------------------------------------------------------------------------
        
        if ( disp )
            cout << "Construction des matrices B" << endl << endl;
        
        Vec< Vec< Mat<T, Gen<>, SparseUMFPACK > > > B;
        B.resize( m.sub_mesh(Number<1>()).elem_list.size() );
        
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            B[ k ].resize( dim );
            for (unsigned d=0;d<dim;++d) {
                B[ k ][ d ].resize( vertex_node_cpt_face[ k ] * m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
            }
        }
        
        apply( m.sub_mesh(Number<1>()).elem_list, Calcul_Skin_Elem_Matrix_B_p_1(), m, B );
        
        if ( disp ) {
            for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
                for (unsigned d=0;d<dim;++d) {
                    cout << "dimension de la matrice B associe a la face " << k << " dans la direction " << d << " = ( " << vertex_node_cpt_face[ k ] * m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << ", " << vertex_node_cpt_face[ k ] * m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << " )" << endl;
                    cout << "matrice B associe a la face " << k << " dans la direction " << d << " =" << endl;
                    cout << B[ k ][ d ] << endl;
                }
            }
        }
        
        /// Construction des vecteurs Q[ k ][ d ] pour chaque face k du maillage et chaque direction d
        /// ------------------------------------------------------------------------------------------
        
        if ( disp )
            cout << "Construction des vecteurs Q" << endl << endl;
        
        Vec< Vec< Vec<T> > > Q;
        Q.resize( m.sub_mesh(Number<1>()).elem_list.size() );
        
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            Q[ k ].resize( dim );
            for (unsigned d=0;d<dim;++d) {
                Q[ k ][ d ].resize( vertex_node_cpt_face[ k ] * m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
                Q[ k ][ d ].set( 0. );
            }
        }
        
        Calcul_Skin_Elem_Vector_Q_p_1_PGD<TVV> calcul_skin_elem_vector_Q_p_1_PGD;
        calcul_skin_elem_vector_Q_p_1_PGD.connect_node_to_vertex_node = &connect_node_to_vertex_node;
        calcul_skin_elem_vector_Q_p_1_PGD.face_type = &face_type;
        calcul_skin_elem_vector_Q_p_1_PGD.face_ind = &face_ind;
        calcul_skin_elem_vector_Q_p_1_PGD.node_list_face = &node_list_face;
        calcul_skin_elem_vector_Q_p_1_PGD.pb = &pb;
        calcul_skin_elem_vector_Q_p_1_PGD.want_local_enrichment = &want_local_enrichment;
        calcul_skin_elem_vector_Q_p_1_PGD.dep = &dep;
        calcul_skin_elem_vector_Q_p_1_PGD.elem_group = &elem_group;
        
        apply( m.elem_list, calcul_skin_elem_vector_Q_p_1_PGD, m, f, Q );
        
        if ( disp ) {
            for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
                for (unsigned d=0;d<dim;++d) {
                    cout << "dimension du vecteur Q associe a la face " << k << " dans la direction " << d << " = " << vertex_node_cpt_face[ k ] * m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << endl;
                    cout << "vecteur Q associe a la face " << k << " dans la direction " << d << " =" << endl;
                    cout << Q[ k ][ d ] << endl << endl;
                }
            }
        }
        
        /// Resolution des systemes lineaires B[ k ][ d ] * lambda_F_face[ k ][ d ] = Q[ k ][ d ] sur chaque face k du maillage et dans chaque direction d
        /// Construction des vecteurs lambda_F_face[ k ][ d ] pour chaque face k du maillage et chaque direction d
        /// ----------------------------------------------------------------------------------------------------------------------------------------------
        
        if ( disp ) {
            cout << "Resolution des systemes lineaires B * lambda_F_face = Q" << endl;
            cout << "Construction des vecteurs lambda_F_face" << endl << endl;
        }
        
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                B[ k ][ d ].get_factorization();
                lambda_F_face[ k ][ d ] = B[ k ][ d ].solve( Q[ k ][ d ] );
            }
        }
        
        if ( disp ) {
            for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
                for (unsigned d=0;d<dim;++d) {
                    cout << "dimension du vecteur lambda_F_face associe a la face " << k << " dans la direction " << d << " = " << vertex_node_cpt_face[ k ] * m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << endl;
                    cout << "vecteur lambda_F_face associe a la face " << k << " dans la direction " << d << " =" << endl;
                    cout << lambda_F_face[ k ][ d ] << endl << endl;
                }
            }
        }
        
        B.free();
        Q.free();
        vertex_node_cpt_face.free();
        
        /// Construction des vecteurs lambda_F[ j ][ d ] pour chaque sommet j du maillage et chaque direction d
        /// ---------------------------------------------------------------------------------------------------
        
        if ( disp )
            cout << "Construction des vecteurs lambda_F" << endl << endl;
        
        Calcul_Vertex_Nodal_Vector_lambda_F_p_1 calcul_vertex_nodal_vector_lambda_F_p_1;
        calcul_vertex_nodal_vector_lambda_F_p_1.face_ind = &face_ind;
        calcul_vertex_nodal_vector_lambda_F_p_1.connect_node_to_vertex_node = &connect_node_to_vertex_node;
        
        apply( m.sub_mesh(Number<1>()).elem_list, calcul_vertex_nodal_vector_lambda_F_p_1, lambda_F_face, lambda_F );
        
        lambda_F_face.free();
        
    }
    else {
        
        /// Cas p >= 2
        /// ----------
        
        if ( disp )
            cout << "Cas : degre p >= 2" << endl << endl;
        
        if ( disp )
            cout << "Construction des vecteurs lambda_F" << endl << endl;
        
        Calcul_Vertex_Nodal_Vector_lambda_F_p_2_PGD<TVV> calcul_vertex_nodal_vector_lambda_F_p_2_PGD;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.connect_node_to_vertex_node = &connect_node_to_vertex_node;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.face_type = &face_type;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.face_ind = &face_ind;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.node_list_face = &node_list_face;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.vertex_node_list_face = &vertex_node_list_face;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.pb = &pb;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.want_local_enrichment = &want_local_enrichment;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.dep = &dep;
        calcul_vertex_nodal_vector_lambda_F_p_2_PGD.elem_group = &elem_group;
        
        apply( m.elem_list, calcul_vertex_nodal_vector_lambda_F_p_2_PGD, m, f, lambda_F );
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur lambda_F associe au noeud sommet " << j << " dans la direction " << d << " = " << nb_unk[ j ][ d ] << endl;
                cout << "vecteur lambda_F associe au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << lambda_F[ j ][ d ] << endl << endl;
            }
        }
    }
    
    degre_p.free();
    vertex_node_list_face.free();
    node_list_face.free();
    
    /// Condition de minimsation
    /// si minimisation[ j ][ d ] = 0, pas d'etape de minimisation d'une fonction-cout pour le noeud sommet j du maillage dans la direction d
    /// si minimisation[ j ][ d ] = 1, etape de minimisation d'une fonction-cout pour le noeud sommet j du maillage dans la direction d
    /// -------------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Condition de minimsation" << endl << endl;
    
    Vec< Vec<bool> > minimisation;
    minimisation.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        minimisation[ j ].resize( dim );
        minimisation[ j ].set( 0 );
    }
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            if ( nb_eq_indep[ j ][ d ] < nb_unk[ j ][ d ] ) {
                minimisation[ j ][ d ] = 1;
            }
        }
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "etape de minimisation associee au noeud sommet " << j << " dans la direction " << d << " = " << minimisation[ j ][ d ] << endl << endl;
            }
        }
    }
    
    /// Construction des matrices de minimisation M[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// ------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices M" << endl << endl;
    
    Vec< Vec< Mat< T, Diag<> > > > M;
    M.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        M[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            if ( minimisation[ j ][ d ] ) {
                M[ j ][ d ].resize( nb_unk[ j ][ d ] );
                M[ j ][ d ].set( 0. );
            }
        }
    }
    
    Calcul_Vertex_Nodal_Matrix_M<T> calcul_vertex_nodal_matrix_M;
    calcul_vertex_nodal_matrix_M.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_vertex_nodal_matrix_M.face_type = &face_type;
    calcul_vertex_nodal_matrix_M.face_ind = &face_ind;
    calcul_vertex_nodal_matrix_M.cost_function = &cost_function;
    calcul_vertex_nodal_matrix_M.penalty_val_N = &penalty_val_N;
    calcul_vertex_nodal_matrix_M.minimisation = &minimisation;
    
    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, calcul_vertex_nodal_matrix_M, m, f, M );
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                if ( minimisation[ j ][ d ] ) {
                    cout << "dimension de la matrice de minimisation M associee au noeud sommet " << j << " dans la direction " << d << " = ( " << nb_unk[ j ][ d ] << ", " << nb_unk[ j ][ d ] << " )" << endl;
                    cout << "matrice M associe au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                    cout << M[ j ][ d ] << endl;
                }
            }
        }
    }
    
    face_type.free();
    
    /// Construction des matrices K[ j ][ d ] et des vecteurs F[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// ------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices K et des vecteurs F" << endl << endl;
    
    Vec< Vec< Mat<T> > > K;
    K.resize( nb_vertex_nodes );
    
    Vec< Vec< Vec<T> > > F;
    F.resize( nb_vertex_nodes );
    
    Vec< Vec< Vec<T> > > U;
    U.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        K[ j ].resize( dim );
        F[ j ].resize( dim );
        U[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            if ( not minimisation[ j ][ d ] ) {
                K[ j ][ d ].resize( nb_unk[ j ][ d ] );
                F[ j ][ d ].resize( nb_unk[ j ][ d ] );
                U[ j ][ d ].resize( nb_unk[ j ][ d ] );
            }
            else {
                K[ j ][ d ].resize( nb_unk[ j ][ d ] + nb_eq_indep[ j ][ d ] );
                K[ j ][ d ].set( 0. );
                F[ j ][ d ].resize( nb_unk[ j ][ d ] + nb_eq_indep[ j ][ d ] );
                U[ j ][ d ].resize( nb_unk[ j ][ d ] + nb_eq_indep[ j ][ d ] );
            }
            K[ j ][ d ].set( 0. );
            F[ j ][ d ].set( 0. );
            U[ j ][ d ].set( 0. );
        }
    }
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            if ( not minimisation[ j ][ d ] ) {
                Vec<unsigned> vec_unk = range( nb_unk[ j ][ d ] );
                Vec<unsigned> vec_eq_indep = range( nb_eq_indep[ j ][ d ] );
                TMatGenSparse trans_A_tilde = trans( A_tilde[ j ][ d ] );
                K[ j ][ d ]( vec_eq_indep, vec_unk ) = trans_A_tilde( vec_eq_indep, vec_unk ) * 1.;
                F[ j ][ d ][ vec_eq_indep ] = R_tilde[ j ][ d ][ vec_eq_indep ] * 1.;
            }
            else {
                Vec<unsigned> vec_unk = range( nb_unk[ j ][ d ] );
                Vec<unsigned> vec_eq_indep = range( nb_eq_indep[ j ][ d ] );
                Vec<unsigned> vec_unk_to_unk_plus_eq_indep = range( nb_unk[ j ][ d ], nb_unk[ j ][ d ] + nb_eq_indep[ j ][ d ] );
                TMatGenSparse trans_A_tilde = trans( A_tilde[ j ][ d ] );
                K[ j ][ d ]( vec_unk, vec_unk ) = M[ j ][ d ][ vec_unk ] * 1.;
                K[ j ][ d ]( vec_unk_to_unk_plus_eq_indep, vec_unk ) = trans_A_tilde( vec_eq_indep, vec_unk ) * 1.;
                K[ j ][ d ]( vec_unk, vec_unk_to_unk_plus_eq_indep ) = A_tilde[ j ][ d ]( vec_unk, vec_eq_indep ) * 1.;
                F[ j ][ d ][ vec_unk ] = ( M[ j ][ d ] * lambda_F[ j ][ d ] )[ vec_unk ] * 1.;
                F[ j ][ d ][ vec_unk_to_unk_plus_eq_indep ] = R_tilde[ j ][ d ][ vec_eq_indep ] * 1.;
            }
        }
    }
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension de la matrice K associee au noeud sommet " << j << " dans la direction " << d << " = ( " << K[ j ][ d ].nb_rows() << ", " << K[ j ][ d ].nb_cols() << " ) "<< endl;
                cout << "matrice K associe au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << K[ j ][ d ] << endl;
            }
        }
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur F associe au noeud sommet " << j << " dans la direction " << d << " = " << F[ j ][ d ].size() << endl;
                cout << "vecteur F associe au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << F[ j ][ d ] << endl << endl;
            }
        }
    }
    
    M.free();
    lambda_F.free();
    nb_eq_indep.free();
    if ( not disp ) {
        A_tilde.free();
        R_tilde.free();
    }
    
    /// Resolution des problemes de minimsation K[ j ][ d ] * U[ j ][ d ] = F[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// Construction des vecteurs U[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// --------------------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp ) {
        cout << "Resolution des pbs de minimisation K * U = F" << endl;
        cout << "Construction des vecteurs U" << endl << endl;
    }
    
    TicToc t_solve_minimization;
    t_solve_minimization.start();
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            if ( not minimisation[ j ][ d ] ) {
                #ifdef WITH_UMFPACK
                Mat<T, Gen<>, SparseUMFPACK > K_UMFPACK = K[ j ][ d ];
                K_UMFPACK.get_factorization();
                U[ j ][ d ] = K_UMFPACK.solve( F[ j ][ d ] );
                #endif
            }
            else {
                if ( solver_minimisation == "LDL" ) {
                    #ifdef WITH_LDL
                    Mat<T, Sym<>, SparseLine<> > K_LDL = K[ j ][ d ];
                    U[ j ][ d ] = F[ j ][ d ];
                    LDL_solver ls;
                    Vec< Vec<T> > Ker;
                    Vec<int> Pivots;
                    ls.get_factorization( K_LDL, Ker, Pivots );
                    ls.solve( U[ j ][ d ] );
                    #endif
                }
                else if ( solver_minimisation == "UMFPACK" ) {
                    #ifdef WITH_UMFPACK
                    Mat<T, Gen<>, SparseUMFPACK > K_UMFPACK = K[ j ][ d ];
                    K_UMFPACK.get_factorization();
                    U[ j ][ d ] = K_UMFPACK.solve( F[ j ][ d ] );
                    #endif
                }
                else if ( solver_minimisation == "Inv" ) {
                    Mat<T, Sym<>, SparseLine<> > K_Inv = K[ j ][ d ];
                    U[ j ][ d ] = inv( K_Inv ) * F[ j ][ d ];
                }
                else if ( solver_minimisation == "LUFactorize" ) {
                    Mat<T, Sym<>, SparseLine<> > K_sym = K[ j ][ d ];
                    Mat<T, Gen<>, SparseLU > K_LU = K_sym;
//                    Vec<int> vector_permutation;
                    lu_factorize( K_LU/*, vector_permutation*/ );
                    solve_using_lu_factorize( K_LU, /*vector_permutation,*/ F[ j ][ d ], U[ j ][ d ] );
                }
                else if ( solver_minimisation == "CholFactorize" ) {
                    Mat<T , Sym<>, SparseLine<> > K_sym = K[ j ][ d ];
                    Mat<T , Sym<>, SparseLine<> > K_fact = K_sym;
                    incomplete_chol_factorize( K_fact );
                    solve_using_incomplete_chol_factorize( K_fact, K_sym, F[ j ][ d ], U[ j ][ d ] );
                }
                else {
                    cerr << "Bing. Error : solveur " << solver_minimisation << " pour la minimisation non implemente" << endl << endl;
                }
            }
        }
    }
    
    t_solve_minimization.stop();
    if ( disp )
        cout << "temps de calcul de la resolution des pbs de minimisation = " << t_solve_minimization.res << endl << endl;
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur U associe au noeud sommet " << j << " dans la direction " << d << " = " << U[ j ][ d ].size() << endl;
                cout << "vecteur U associe au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << U[ j ][ d ] << endl << endl;
            }
        }
    }
    if ( verif_solver_minimisation ) {
        if ( disp )
            cout << "Verification de la resolution des pbs de minimisation : tolerance = " << tol_solver_minimisation << endl << endl;
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                T residual = norm_2( K[ j ][ d ] * U[ j ][ d ] - F[ j ][ d ] );
                T b = norm_2( F[ j ][ d ] );
                if ( residual / b > tol_solver_minimisation ) {
                    cout << "residu associe au noeud sommet " << j << " dans la direction " << d << " :" << endl;
//                    cout << "K * U - F =" << endl;
//                    cout << K[ j ][ d ] * U[ j ][ d ] - F[ j ][ d ] << endl;
                    cout << "norme du residu = " << residual << endl;
                    cout << "norme du residu relatif = " << residual / b << endl << endl;
                }
            }
        }
    }
    
    K.free();
    F.free();
    minimisation.free();
    
    /// Construction des vecteurs lambda_F_hat[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// -------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs lambda_F_hat" << endl << endl;
    
    Vec< Vec< Vec<T> > > lambda_F_hat;
    lambda_F_hat.resize( nb_vertex_nodes );
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        lambda_F_hat[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            lambda_F_hat[ j ][ d ].resize( nb_unk[ j ][ d ] );
            lambda_F_hat[ j ][ d ].set( 0. );
        }
    }
    
    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        for (unsigned d=0;d<dim;++d) {
            Vec<unsigned> vec_unk_lambda_F_hat = range( nb_unk[ j ][ d ] );
            lambda_F_hat[ j ][ d ][ vec_unk_lambda_F_hat ] = U[ j ][ d ][ vec_unk_lambda_F_hat ] * 1.;
        }
    }
    
    Reset_Vertex_Nodal_Vector_lambda_F_hat reset_vertex_nodal_vector_lambda_F_hat;
    reset_vertex_nodal_vector_lambda_F_hat.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    reset_vertex_nodal_vector_lambda_F_hat.enhancement = &enhancement;
    reset_vertex_nodal_vector_lambda_F_hat.flag_face_enh = &flag_face_enh;
    reset_vertex_nodal_vector_lambda_F_hat.face_ind = &face_ind;
    
    apply( m.sub_mesh(Number<1>()).elem_list, reset_vertex_nodal_vector_lambda_F_hat, m, lambda_F_hat );
    
    if ( disp ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur lambda_F_hat associee au noeud sommet " << j << " dans la direction " << d << " = " << nb_unk[ j ][ d ] << endl;
                cout << "vecteur lambda_F_hat associe au noeud sommet " << j << " dans la direction " << d << " =" << endl;
                cout << lambda_F_hat[ j ][ d ] << endl << endl;
                if ( norm_2( trans( A[ j ][ d ] ) * lambda_F_hat[ j ][ d ] - R[ j ][ d ] ) / norm_2( R[ j ][ d ] ) > tol_solver_minimisation ) {
                    cout << "A^T * lambda_F_hat - R =" << endl;
                    cout << trans( A[ j ][ d ] ) * lambda_F_hat[ j ][ d ] - R[ j ][ d ] << endl << endl;
                    cout << "A_tilde^T * lambda_F_hat - R_tilde =" << endl;
                    cout << trans( A_tilde[ j ][ d ] ) * lambda_F_hat[ j ][ d ] - R_tilde[ j ][ d ]<< endl << endl;
                }
            }
        }
    }
    
    nb_unk.free();
    U.free();
    
    /// Construction des vecteurs force_fluxes[ k ][ d ] pour chaque face k du maillage et chaque direction d
    /// -----------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs force_fluxes_std" << endl << endl;
    
    force_fluxes.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        force_fluxes[ k ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            force_fluxes[ k ][ d ].resize( m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
            force_fluxes[ k ][ d ].set( 0. );
        }
    }
    
    Calcul_Skin_Elem_Force_Fluxes calcul_skin_elem_force_fluxes;
    calcul_skin_elem_force_fluxes.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    calcul_skin_elem_force_fluxes.face_ind = &face_ind;
    
    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, calcul_skin_elem_force_fluxes, lambda_F_hat, force_fluxes );
    
    if ( disp ) {
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur des densites d'effort associe a la face " << k << " dans la direction " << d << " = " << m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << endl;
                cout << "vecteur des densites d'effort associe a la face " << k << " dans la direction " << d << " =" << endl;
                cout << force_fluxes[ k ][ d ] << endl << endl;
            }
        }
    }
    
    face_ind.free();
    lambda_F_hat.free();
    connect_node_to_vertex_node.free();
    
    t_force_fluxes_std.stop();
    cout << "temps de calcul de la construction standard des densites d'effort = " << t_force_fluxes_std.res << endl << endl;
    
}

#endif // Construct_standard_force_fluxes_EESPT_PGD_h
