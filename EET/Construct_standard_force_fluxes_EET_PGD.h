//
// C++ Interface: Construct_standard_force_fluxes_EET_PGD
//
// Description: construction standard des densites d'effort par la methode EET dans le cadre des methodes PGD
//
//
// Author: Pled Florent <florent.pled@univ-paris-est.fr>, (C) 2017
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef Construct_standard_force_fluxes_EET_PGD_h
#define Construct_standard_force_fluxes_EET_PGD_h

#include "EET.h"
#include "../CONNECTIVITY/Calcul_connectivity.h"
#include "../LMT/include/containers/matumfpack.h"

using namespace LMT;
using namespace std;

/// Construction standard des densites d'effort par la methode EET
/// --------------------------------------------------------------
template<class TM, class TF, class T, class TVV, class TVVV>
void construct_standard_force_fluxes_EET_PGD( TM &m, const TF &f, const string &pb, const unsigned &cost_function, const bool enhancement, const Vec<bool> &flag_face_enh, const string &solver_minimisation, TVVV &force_fluxes, const TVV &dep, const Vec< Vec<unsigned> > &elem_group, const bool want_local_enrichment = false, const bool verif_solver_minimisation = false, const T tol_solver_minimisation = 1e-6, const bool verif_compatibility_conditions = false, const T tol_compatibility_conditions = 1e-6, const bool disp = false ) {
    
    static const unsigned dim = TM::dim;
    
    typedef Mat<T, Gen<>, SparseLine<> > TMatGenSparse;
    typedef Mat<T, Sym<>, SparseLine<> > TMatSymSparse;
    typedef Mat< T, Diag<> > TMatDiag;
    typedef Mat<T, Gen<>, SparseUMFPACK > TMatGenSparseUMFPACK;
    typedef Mat<T, Gen<>, SparseLU > TMatGenSparseLU;
    
    TicToc t_force_fluxes_std;
    t_force_fluxes_std.start();
    
    Vec<unsigned> elem_cpt_node;
    Vec< Vec<unsigned> > elem_list_node;
    construct_elems_connected_to_node( m, elem_cpt_node, elem_list_node );
    
    elem_list_node.free();
    
    Vec<unsigned> face_cpt_node;
    Vec< Vec<unsigned> > face_list_node;
    construct_faces_connected_to_node( m, face_cpt_node, face_list_node );
    
    Vec<unsigned> node_cpt_face;
    Vec< Vec<unsigned> > node_list_face;
    construct_nodes_connected_to_face( m, node_cpt_face, node_list_face );
    
    node_cpt_face.free();
    
    Vec<bool> correspondance_node_to_vertex_node;
    Vec<unsigned> connect_node_to_vertex_node;
    unsigned nb_vertex_nodes = match_node_to_vertex_node( m, correspondance_node_to_vertex_node, connect_node_to_vertex_node );
    
    connect_node_to_vertex_node.free();
    
    Vec< Vec<unsigned> > face_type;
    construct_face_type( m, f, face_type );
    
    Vec< Vec<unsigned> > node_type;
    construct_node_type( m, f, face_type, node_type );
    
    /// -------------------------------------------------- ///
    /// Construction des projections des densites d'effort ///
    /// -------------------------------------------------- ///
    
    if ( disp )
        cout << "Construction des projections des densites d'effort standard" << endl << endl;
    
    /// Reperage pour chaque face k et chaque direction d de l'indice de debut de ligne dans les vecteurs b[ i ][ d ] et de debut de colonne dans les matrices B[ i ][ d ] : face_ind[ k ][ d ]
    /// Calcul du nb de lignes du vecteur b[ i ][ d ] et de colonnes de la matrice B[ i ][ d ] : nb_unk[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Calcul des vecteurs nb_unk" << endl << endl;
    
    Vec< Vec<unsigned> > nb_unk;
    nb_unk.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        nb_unk[ i ].resize( dim );
        nb_unk[ i ].set( 0 );
    }
    
    Vec< Vec< Vec<unsigned> > > face_ind;
    face_ind.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        face_ind[ k ].resize( dim );
    }
    
    apply( m.sub_mesh(Number<1>()).elem_list, Calcul_Face_Ind_EET(), face_ind, nb_unk ); // nb_unk[ i ][ d ] contient le nb de lignes du vecteur b[ i ][ d ] et le nb de colonnes de la matrice B[ i ][ d ] associee au i eme noeud dans la direction d // face_ind[ k ][ d ][ numero_noeud_dans_face (0 ou 1) ] = indice de debut de ligne dans le vecteur b[ i ][ d ] et de debut de colonne dans la matrice B[ i ][ d ] pour chaque face k du maillage et chaque direction d
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "nb d'inconnues associees au noeud " << i << " dans la direction " << d << " = " << nb_unk[ i ][ d ] << endl << endl;
            }
        }
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "indice de debut de ligne de la face " << k << " dans la direction " << d << " dans les vecteurs b[ noeud connecte a la face " << k << " ][ " << d << " ] et de debut de colonne dans les matrices B[ noeud connecte a la face " << k << " ][ " << d << " ] = " << face_ind[ k ][ d ] << endl << endl;
            }
        }
    }
    
    /// Reperage pour chaque element n et chaque direction d de l'indice de debut de ligne dans les vecteurs r[ i ][ d ] et dans les matrices B[ i ][ d ] : elem_ind[ n ][ d ]
    /// Calcul du nb de lignes du vecteur r[ i ][ d ] et de la matrice B[ i ][ d ] : nb_eq[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Calcul des vecteurs nb_eq" << endl << endl;
    
    Vec< Vec<unsigned> > nb_eq;
    nb_eq.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        nb_eq[ i ].resize( dim );
        nb_eq[ i ].set( 0 );
    }
    
    Vec< Vec< Vec<unsigned> > > elem_ind;
    elem_ind.resize( m.elem_list.size() );
    
    for (unsigned n=0;n<m.elem_list.size();++n) {
        elem_ind[ n ].resize( dim );
    }
    
    apply( m.elem_list, Calcul_Elem_Ind_EET(), elem_ind, nb_eq ); // nb_eq[ i ][ d ] contient le nb de lignes du vecteur r[ i ][ d ] et de la matrice B[ i ][ d ] associee au i eme noeud dans la direction d // elem_ind[ n ][ d ][ numero_noeud_dans_elem (0 ou 1 ou 2) ] = indice de debut de ligne dans le vecteur r[ i ][ d ] et dans la matrices B[ i ][ d ] pour chaque element n du maillage et chaque direction d
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "nb d'equations associees au noeud " << i << " dans la direction " << d << " = " << nb_eq[ i ][ d ] << endl << endl;
            }
        }
        for (unsigned n=0;n<m.elem_list.size();++n) {
            for (unsigned d=0;d<dim;++d) {
                cout << "indice de debut de ligne de l'element " << n << " dans la direction " << d << " dans les vecteurs r[ noeud connecte a l'element " << n << " ][ " << d << " ] et dans les matrices B[ noeud connecte a l'element " << n << " ][ " << d << " ] = " << elem_ind[ n ][ d ] << endl << endl;
            }
        }
    }
    
    /// Construction des matrices B[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// -------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices B" << endl << endl;
    
    Vec< Vec< TMatGenSparse > > B;
    B.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        B[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            B[ i ][ d ].resize( nb_eq[ i ][ d ], nb_unk[ i ][ d ] );
        }
    }
    
    Calcul_Nodal_Matrix_B calcul_nodal_matrix_B;
    calcul_nodal_matrix_B.elem_ind = &elem_ind;
    calcul_nodal_matrix_B.face_ind = &face_ind;
    calcul_nodal_matrix_B.node_list_face = &node_list_face;
    
    apply( m.elem_list, calcul_nodal_matrix_B, m, B );
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension de la matrice B associe au noeud " << i << " dans la direction " << d << " = ( " << nb_eq[ i ][ d ] << ", " << nb_unk[ i ][ d ] << " )" << endl;
                cout << "matrice B associe au noeud " << i << " dans la direction " << d << " =" << endl;
                cout << B[ i ][ d ] << endl;
            }
        }
    }
    
    /// Construction des vecteurs r[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// -------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs r" << endl << endl;
    
    Vec< Vec< Vec<T> > > r;
    r.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        r[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            r[ i ][ d ].resize( nb_eq[ i ][ d ] );
            r[ i ][ d ].set( 0. );
        }
    }
    
    Calcul_Nodal_Vector_r_PGD<TVV> calcul_nodal_vector_r_PGD;
    calcul_nodal_vector_r_PGD.elem_ind = &elem_ind;
    calcul_nodal_vector_r_PGD.node_list_face = &node_list_face;
    calcul_nodal_vector_r_PGD.elem_cpt_node = &elem_cpt_node;
    calcul_nodal_vector_r_PGD.pb = &pb;
    calcul_nodal_vector_r_PGD.want_local_enrichment = &want_local_enrichment;
    calcul_nodal_vector_r_PGD.dep = &dep;
    calcul_nodal_vector_r_PGD.elem_group = &elem_group;
    
    apply( m.elem_list, calcul_nodal_vector_r_PGD, m, f, r );
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur r associe au noeud " << i << " dans la direction " << d << " = " << nb_eq[ i ][ d ] << endl;
                cout << "vecteur r associe au noeud " << i << " dans la direction " << d << " =" << endl;
                cout << r[ i ][ d ] << endl << endl;
            }
        }
    }
    
    elem_ind.free();
    
    /// Verification des conditions de compatibilite (equilibre elements finis) pour les noeuds interieurs (type 0)
    /// -----------------------------------------------------------------------------------------------------------
    
    if ( verif_compatibility_conditions ) {
        
        if ( disp )
            cout << "Verification des conditions de compatibilite (equilibre elements finis) pour les noeuds interieurs (type 0) : tolerance = " << tol_compatibility_conditions << endl << endl;
        
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                if ( node_type[ i ][ d ] == 0 ) {
                    T residual = 0;
                    for (unsigned n=0;n<elem_cpt_node[ i ];n++)
                        residual += r[ i ][ d ][ n ];
                    if ( fabs( residual ) > tol_compatibility_conditions ) {
                        cout << "verification des conditions de compatibilite au noeud " << i << " dans la direction " << d << " :" << endl;
                        cout << residual << " != 0" << endl << endl;
                    }
                }
            }
        }
    }
    
    /// Reperage pour chaque face k connectee au noeud i et chaque direction d de l'indice de debut de ligne dans les matrices C[ i ][ d ] et de debut de ligne dans les vecteurs q[ i ][ d ] : nodal_ind[ i ][ d ][ k ]
    /// Calcul du nb de lignes de la matrice C[ i ][ d ] et du vecteur q[ i ][ d ] : nb_eq_imp[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Calcul des vecteurs nb_eq_imp" << endl << endl;
    
    Vec< Vec<unsigned> > nb_eq_imp;
    nb_eq_imp.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        nb_eq_imp[ i ].resize( dim );
        nb_eq_imp[ i ].set( 0 );
    }
    
    Vec< Vec< Vec<unsigned> > > nodal_ind;
    nodal_ind.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        nodal_ind[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            nodal_ind[ i ][ d ].resize( face_cpt_node[ i ] );
            nodal_ind[ i ][ d ].set( 0 );
        }
    }
    
    Calcul_Nodal_Ind calcul_nodal_ind;
    calcul_nodal_ind.face_cpt_node = &face_cpt_node;
    calcul_nodal_ind.face_list_node = &face_list_node;
    calcul_nodal_ind.node_type = &node_type;
    calcul_nodal_ind.face_type = &face_type;
    
    apply( m.node_list, calcul_nodal_ind, nb_eq_imp, nodal_ind ); // nb_eq_imp[ i ][ d ] contient le nb de lignes de la matrice C[ i ][ d ] et du vecteur q[ i ][ d ] associee au i eme noeud dans la direction d // nodal_ind[ i ][ d ][ k ] = indice de debut de ligne dans la matrice C[ i ][ d ] et dans le vecteur q[ i ][ d ] associes au ieme noeud pour chaque face k du maillage et chaque direction d
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "nb d'equations imposees associees au noeud " << i << " dans la direction " << d << " = " << nb_eq_imp[ i ][ d ] << endl;
                if ( node_type[ i ][ d ] == 2 or node_type[ i ][ d ] == 12 ) {
                    for (unsigned k=0;k<face_cpt_node[ i ];++k) {
                        if ( face_type[ face_list_node[ i ][ k ] ][ d ] == 2 ) {
                            cout << "indice de debut de ligne de la face " << face_list_node[ i ][ k ] << " dans la direction " << d << " dans la matrice C[ " << i << " ][ " << d << " ] = " << nodal_ind[ i ][ d ][ k ] << endl;
                        }
                    }
                }
                cout << endl;
            }
        }
    }
    
    /// Construction des matrices C[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// -------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices C" << endl << endl;
    
    Vec< Vec< TMatGenSparse > > C;
    C.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        C[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            C[ i ][ d ].resize( nb_eq_imp[ i ][ d ], nb_unk[ i ][ d ] );
        }
    }
    
    Calcul_Nodal_Matrix_C calcul_nodal_matrix_C;
    calcul_nodal_matrix_C.face_type = &face_type;
    calcul_nodal_matrix_C.nodal_ind = &nodal_ind;
    calcul_nodal_matrix_C.face_ind = &face_ind;
    calcul_nodal_matrix_C.face_list_node = &face_list_node;
    
    apply( m.sub_mesh(Number<1>()).elem_list, calcul_nodal_matrix_C, m, C );
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                if ( node_type[ i ][ d ] == 2 or node_type[ i ][ d ] == 12 ) {
                    cout << "dimension de la matrice C associe au noeud " << i << " dans la direction " << d << " = ( " << nb_eq_imp[ i ][ d ] << ", " << nb_unk[ i ][ d ] << " )" << endl;
                    cout << "matrice C associe au noeud " << i << " dans la direction " << d << " =" << endl;
                    cout << C[ i ][ d ] << endl;
                }
            }
        }
    }
    
    /// Construction des vecteurs q[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// -------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs q" << endl << endl;
    
    Vec< Vec< Vec<T> > > q;
    q.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        q[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            q[ i ][ d ].resize( nb_eq_imp[ i ][ d ] );
            q[ i ][ d ].set( 0. );
        }
    }
    
    Calcul_Nodal_Vector_q_PGD<TVV> calcul_nodal_vector_q_PGD;
    calcul_nodal_vector_q_PGD.face_type = &face_type;
    calcul_nodal_vector_q_PGD.nodal_ind = &nodal_ind;
    calcul_nodal_vector_q_PGD.node_list_face = &node_list_face;
    calcul_nodal_vector_q_PGD.face_list_node = &face_list_node;
    calcul_nodal_vector_q_PGD.dep = &dep;
    calcul_nodal_vector_q_PGD.elem_group = &elem_group;
    
    apply( m.elem_list, calcul_nodal_vector_q_PGD, m, f, q );
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                if ( node_type[ i ][ d ] == 2 or node_type[ i ][ d ] == 12 ) {
                    cout << "dimension du vecteur q associe au noeud " << i << " dans la direction " << d << " = " << nb_eq_imp[ i ][ d ] << endl;
                    cout << "vecteur q associe au noeud " << i << " dans la direction " << d << " =" << endl;
                    cout << q[ i ][ d ] << endl << endl;
                }
            }
        }
    }
    
    /// Verification des conditions de compatibilite (equilibre elements finis) pour les noeuds sur la frontière de Neumann (type 2)
    /// ----------------------------------------------------------------------------------------------------------------------------
    
    if ( verif_compatibility_conditions ) {
        
        if ( disp )
            cout << "Verification des conditions de compatibilite (equilibre elements finis) pour les noeuds sur la frontière de Neumann (type 2) : tolerance = " << tol_compatibility_conditions << endl << endl;
        
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                if ( node_type[ i ][ d ] == 2 ) {
                    T sum_r = 0.;
                    for (unsigned n=0;n<elem_cpt_node[ i ];n++)
                        sum_r += r[ i ][ d ][ n ];
                    T sum_q = 0.;
                    for (unsigned k=0;k<face_cpt_node[ i ];k++) {
                        if ( face_type[ face_list_node[ i ][ k ] ][ d ] == 2 )
                            sum_q += q[ i ][ d ][ nodal_ind[ i ][ d ][ k ] ];
                    }
                    if ( fabs( sum_r - sum_q ) > tol_compatibility_conditions ) {
                        cout << "verification des conditions de compatibilite au noeud " << i << " dans la direction " << d << " :" << endl;
                        cout << sum_r << " != " << sum_q << endl << endl;
                    }
                }
            }
        }
    }
    
    nodal_ind.free();
    elem_cpt_node.free();
    face_cpt_node.free();
    face_list_node.free();
    
    /// Calcul du nb d'equations independantes : nb_eq_indep[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// --------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Calcul des vecteurs nb_eq_indep" << endl << endl;
    
    Vec< Vec<unsigned> > nb_eq_indep;
    nb_eq_indep.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        nb_eq_indep[ i ].resize( dim );
        nb_eq_indep[ i ].set( 0 );
    }
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        for (unsigned d=0;d<dim;++d) {
            if ( node_type[ i ][ d ] == 0 or node_type[ i ][ d ] == 2 ) { // si node_type[ i ][ d ] = 0 ou 2, l'ecriture de l'equilibre elements finis dans la direction d supprime 1 equation : nb_eq_indep[ i ][ d ] = nb_eq[ i ][ d ] - 1
                nb_eq_indep[ i ][ d ] = nb_eq[ i ][ d ] - 1;
            }
            else {
                nb_eq_indep[ i ][ d ] = nb_eq[ i ][ d ];
            }
        }
    }
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "nb d'equations independantes associees au noeud " << i << " dans la direction " << d << " = " << nb_eq_indep[ i ][ d ] << endl << endl;
            }
        }
    }
    
    /// Condition de minimsation
    /// si minimisation[ i ][ d ] = 0, pas d'etape de minimisation d'une fonction-cout pour le noeud i du maillage dans la direction d
    /// si minimisation[ i ][ d ] = 1, etape de minimisation d'une fonction-cout pour le noeud i du maillage dans la direction d
    /// ------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Condition de minimsation" << endl << endl;
    
    Vec< Vec<bool> > minimisation;
    minimisation.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        minimisation[ i ].resize( dim );
        minimisation[ i ].set( 0 );
    }
    
    /// Cas node_type[ i ][ d ] = 0 ou 1
    /// --------------------------------
    // si i est un noeud sommet (en 2D ou 3D) ou i est un noeud interieur a une arete (en 3D)
    //     etape de minimsation : minimisation[ i ][ d ] = 1
    // si i est un noeud non sommet (en 2D) ou i est un noeud interieur a une face (en 3D)
    //     pas d'etape de minimsation : minimisation[ i ][ d ] = 0
    // autrement dit :
    // si nb_eq_indep[ i ][ d ] < nb_unk[ i ][ d ]
    //     etape de minimsation : minimisation[ i ][ d ] = 1
    // sinon ( nb_eq_indep[ i ][ d ] = nb_unk[ i ][ d ] )
    //     pas d'etape de minimsation : minimisation[ i ][ d ] = 0
    /// Cas node_type[ i ][ d ] = 12 ou 2
    /// ---------------------------------
    // si nb_eq_imp[ i ][ d ] = nb_unk[ i ][ d ] ( toutes les faces connectees au noeud i dans la direction d sont de type 2, possible uniquement pour node_type[ i ][ d ] = 2 )
    //     pas d'etape de minimsation : minimisation[ i ][ d ] = 0
    // sinon ( nb_eq_imp[ i ][ d ] < nb_unk[ i ][ d ] )
    //     si nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] >= nb_unk[ i ][ d ]
    //         pas d'etape de minimsation : minimisation[ i ][ d ] = 0
    //         si nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] = nb_unk[ i ][ d ]
    //             resolution explicite
    //         si nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] > nb_unk[ i ][ d ]
    //             troncature de B[ i ][ d ] et r[ i ][ d ] puis resolution explicite
    //     si nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] < nb_unk[ i ][ d ]
    //         etape de minimsation : minimisation[ i ][ d ] = 1
    for (unsigned i=0;i<m.node_list.size();++i) {
        for (unsigned d=0;d<dim;++d) {
//            if ( node_type[ i ][ d ] == 0 or node_type[ i ][ d ] == 1 ) {
//                if ( nb_eq_indep[ i ][ d ] < nb_unk[ i ][ d ] ) {
//                    minimisation[ i ][ d ] = 1;
//                }
//            }
//            else if ( node_type[ i ][ d ] == 2 or node_type[ i ][ d ] == 12 ) {
//                if ( nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] < nb_unk[ i ][ d ] ) {
//                    minimisation[ i ][ d ] = 1;
//                }
//            }
            if ( nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] < nb_unk[ i ][ d ] ) {
                minimisation[ i ][ d ] = 1;
            }
        }
    }
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "etape de minimisation associee au noeud " << i << " dans la direction " << d << " = " << minimisation[ i ][ d ] << endl << endl;
            }
        }
    }
    
    /// Construction des matrices de minimisation M[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// -----------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices M" << endl << endl;
    
    Vec< Vec< TMatDiag > > M;
    M.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        M[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            if ( minimisation[ i ][ d ] ) {
                M[ i ][ d ].resize( nb_unk[ i ][ d ] );
                M[ i ][ d ].set( 0. );
            }
        }
    }
    
    Calcul_Nodal_Matrix_M calcul_nodal_matrix_M;
    calcul_nodal_matrix_M.face_ind = &face_ind;
    calcul_nodal_matrix_M.cost_function = &cost_function;
    calcul_nodal_matrix_M.minimisation = &minimisation;
    
    m.update_skin();
    apply( m.sub_mesh(Number<1>()).elem_list, calcul_nodal_matrix_M, m, f, M );
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                if ( minimisation[ i ][ d ] ) {
                    cout << "dimension de la matrice de minimisation M associee au noeud " << i << " dans la direction " << d << " = ( " << nb_unk[ i ][ d ] << ", " << nb_unk[ i ][ d ] << " )" << endl;
                    cout << "matrice M associe au noeud " << i << " dans la direction " << d << " =" << endl;
                    cout << M[ i ][ d ] << endl;
                }
            }
        }
    }
    
    /// Construction des vecteurs de minimisation b[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// -----------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs b" << endl << endl;
    
    Vec< Vec< Vec<T> > > b;
    b.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        b[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            if ( minimisation[ i ][ d ] ) {
                b[ i ][ d ].resize( nb_unk[ i ][ d ] );
                b[ i ][ d ].set( 0. );
            }
        }
    }
    
    Calcul_Nodal_Vector_b_PGD<TVV> calcul_nodal_vector_b_PGD;
    calcul_nodal_vector_b_PGD.minimisation = &minimisation;
    calcul_nodal_vector_b_PGD.face_type = &face_type;
    calcul_nodal_vector_b_PGD.face_ind = &face_ind;
    calcul_nodal_vector_b_PGD.node_list_face = &node_list_face;
    calcul_nodal_vector_b_PGD.pb = &pb;
    calcul_nodal_vector_b_PGD.want_local_enrichment = &want_local_enrichment;
    calcul_nodal_vector_b_PGD.dep = &dep;
    calcul_nodal_vector_b_PGD.elem_group = &elem_group;
    
    apply( m.elem_list, calcul_nodal_vector_b_PGD, m, f, b );
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                if ( minimisation[ i ][ d ] ) {
                    cout << "dimension du vecteur de minimisation b associee au noeud " << i << " dans la direction " << d << " = " << nb_unk[ i ][ d ] << endl;
                    cout << "vecteur b associe au noeud " << i << " dans la direction " << d << " =" << endl;
                    cout << b[ i ][ d ] << endl << endl;
                }
            }
        }
    }
    
    node_list_face.free();
    face_type.free();
    
    /// Construction des matrices K[ i ][ d ] et des vecteurs F[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// -----------------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices K et des vecteurs F" << endl << endl;
    
    Vec< Vec< Mat<T> > > K;
    K.resize( m.node_list.size() );
    
    Vec< Vec< Vec<T> > > F;
    F.resize( m.node_list.size() );
    
    Vec< Vec< Vec<T> > > U;
    U.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        K[ i ].resize( dim );
        F[ i ].resize( dim );
        U[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
//            if ( node_type[ i ][ d ] == 0 or node_type[ i ][ d ] == 1 ) {
//                if ( not minimisation[ i ][ d ] ) {
//                    K[ i ][ d ].resize( nb_unk[ i ][ d ] );
//                    F[ i ][ d ].resize( nb_unk[ i ][ d ] );
//                    U[ i ][ d ].resize( nb_unk[ i ][ d ] );
//                }
//                else {
//                    K[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_indep[ i ][ d ] );
//                    F[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_indep[ i ][ d ] );
//                    U[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_indep[ i ][ d ] );
//                }
//                else if ( node_type[ i ][ d ] == 2 or node_type[ i ][ d ] == 12 ) {
//                    if ( not minimisation[ i ][ d ] ) {
//                        K[ i ][ d ].resize( nb_unk[ i ][ d ] );
//                        F[ i ][ d ].resize( nb_unk[ i ][ d ] );
//                        U[ i ][ d ].resize( nb_unk[ i ][ d ] );
//                    }
//                    else {
//                        K[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] );
//                        F[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] );
//                        U[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] );
//                    }
//                }
            if ( not minimisation[ i ][ d ] ) {
                K[ i ][ d ].resize( nb_unk[ i ][ d ] );
                F[ i ][ d ].resize( nb_unk[ i ][ d ] );
                U[ i ][ d ].resize( nb_unk[ i ][ d ] );
                
            }
            else {
                K[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] );
                F[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] );
                U[ i ][ d ].resize( nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] );
            }
            K[ i ][ d ].set( 0. );
            F[ i ][ d ].set( 0. );
            U[ i ][ d ].set( 0. );
        }
    }
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        for (unsigned d=0;d<dim;++d) {
//            if ( node_type[ i ][ d ] == 0 or node_type[ i ][ d ] == 1 ) {
//                if ( not minimisation[ i ][ d ] ) {
//                    Vec<unsigned> vec_unk = range( nb_unk[ i ][ d ] );
//                    K[ i ][ d ]( vec_unk, vec_unk ) = B[ i ][ d ]( vec_unk, vec_unk ) * 1.;
//                    F[ i ][ d ][ vec_unk ] = r[ i ][ d ][ vec_unk ] * 1.;
//                }
//                else {
//                    Vec<unsigned> vec_unk = range( nb_unk[ i ][ d ] );
//                    Vec<unsigned> vec_eq_indep = range( nb_eq_indep[ i ][ d ] );
//                    Vec<unsigned> vec_unk_to_unk_plus_eq_indep = range( nb_unk[ i ][ d ], nb_unk[ i ][ d ] + nb_eq_indep[ i ][ d ] );
//                    TMatGenSparse trans_B = trans( B[ i ][ d ] );
//                    K[ i ][ d ]( vec_unk, vec_unk ) = M[ i ][ d ][ vec_unk ] * 1.;
//                    K[ i ][ d ]( vec_unk_to_unk_plus_eq_indep, vec_unk ) = B[ i ][ d ]( vec_eq_indep, vec_unk ) * 1.;
//                    K[ i ][ d ]( vec_unk, vec_unk_to_unk_plus_eq_indep ) = trans_B( vec_unk, vec_eq_indep ) * 1.;
//                    F[ i ][ d ][ vec_unk ] = ( M[ i ][ d ] * b[ i ][ d ] )[ vec_unk ] * 1.;
//                    F[ i ][ d ][ vec_unk_to_unk_plus_eq_indep ] = r[ i ][ d ][ vec_eq_indep ] * 1.;
//                }
//            }
//            else if ( node_type[ i ][ d ] == 2 or node_type[ i ][ d ] == 12 ) {
//                if ( nb_eq_imp[ i ][ d ] == nb_unk[ i ][ d ] ) { // possible uniquement si node_type[ i ][ d ] = 2
//                    Vec<unsigned> vec_unk = range( nb_unk[ i ][ d ] );
//                    K[ i ][ d ]( vec_unk, vec_unk ) = C[ i ][ d ]( vec_unk, vec_unk ) * 1.;
//                    F[ i ][ d ][ vec_unk ] = q[ i ][ d ][ vec_unk ] * 1.;
//                }
//                else if ( nb_eq_imp[ i ][ d ][ d ] + nb_eq_indep[ i ][ d ] >= nb_unk[ i ][ d ] ) {
//                    Vec<unsigned> vec_unk = range( nb_unk[ i ][ d ] );
//                    Vec<unsigned> vec_eq_imp = range( nb_eq_imp[ i ][ d ] );
//                    Vec<unsigned> vec_eq_indep = range( nb_eq_indep[ i ][ d ] );
//                    Vec<unsigned> vec_eq_imp_to_unk = range( nb_eq_imp[ i ][ d ], nb_unk[ i ][ d ] );
//                    Vec<unsigned> vec_unk_minus_eq_imp = range( nb_unk[ i ][ d ] - nb_eq_imp[ i ][ d ] );
//                    K[ i ][ d ]( vec_eq_imp, vec_unk ) = C[ i ][ d ]( vec_eq_imp, vec_unk ) * 1.;
//                    K[ i ][ d ]( vec_eq_imp_to_unk, vec_unk ) = B[ i ][ d ]( vec_unk_minus_eq_imp, vec_unk ) * 1.;
//                    F[ i ][ d ][ vec_eq_imp ] = q[ i ][ d ][ vec_eq_imp ] * 1.;
//                    F[ i ][ d ][ vec_eq_imp_to_unk ] = r[ i ][ d ][ vec_unk_minus_eq_imp ] * 1.;
//                }
//                else if ( nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] < nb_unk[ i ][ d ] ) {
//                    Vec<unsigned> vec_unk = range( nb_unk[ i ][ d ] );
//                    Vec<unsigned> vec_eq_imp = range( nb_eq_imp[ i ][ d ] );
//                    Vec<unsigned> vec_eq_indep = range( nb_eq_indep[ i ][ d ] );
//                    Vec<unsigned> vec_unk_to_unk_plus_eq_imp = range( nb_unk[ i ][ d ], nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] );
//                    Vec<unsigned> vec_unk_plus_eq_imp_to_unk_plus_eq_imp_plus_eq_indep = range( nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ], nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] );
//                    TMatGenSparse trans_B = trans( B[ i ][ d ] );
//                    TMatGenSparse trans_C = trans( C[ i ][ d ] );
//                    K[ i ][ d ]( vec_unk, vec_unk ) = M[ i ][ d ][ vec_unk ] * 1.;
//                    K[ i ][ d ]( vec_unk_to_unk_plus_eq_imp, vec_unk ) = C[ i ][ d ]( vec_eq_imp, vec_unk ) * 1.;
//                    K[ i ][ d ]( vec_unk_plus_eq_imp_to_unk_plus_eq_imp_plus_eq_indep, vec_unk ) = B[ i ][ d ]( vec_eq_indep, vec_unk ) * 1.;
//                    K[ i ][ d ]( vec_unk, vec_unk_to_unk_plus_eq_imp ) = trans_C( vec_unk, vec_eq_imp ) * 1.;
//                    K[ i ][ d ]( vec_unk, vec_unk_plus_eq_imp_to_unk_plus_eq_imp_plus_eq_indep ) = trans_B( vec_unk, vec_eq_indep ) * 1.;
//                    F[ i ][ d ][ vec_unk ] = ( M[ i ][ d ] * b[ i ][ d ] )[ vec_unk ] * 1.;
//                    F[ i ][ d ][ vec_unk_to_unk_plus_eq_imp ] = q[ i ][ d ][ vec_eq_imp ] * 1.;
//                    F[ i ][ d ][ vec_unk_plus_eq_imp_to_unk_plus_eq_imp_plus_eq_indep ] = r[ i ][ d ][ vec_eq_indep ] * 1.;
//                }
//            }
            if ( not minimisation[ i ][ d ] ) {
                Vec<unsigned> vec_unk = range( nb_unk[ i ][ d ] );
                Vec<unsigned> vec_eq_imp = range( nb_eq_imp[ i ][ d ] );
                Vec<unsigned> vec_eq_indep = range( nb_eq_indep[ i ][ d ] );
                Vec<unsigned> vec_eq_imp_to_unk = range( nb_eq_imp[ i ][ d ], nb_unk[ i ][ d ] );
                Vec<unsigned> vec_unk_minus_eq_imp = range( nb_unk[ i ][ d ] - nb_eq_imp[ i ][ d ] );
                K[ i ][ d ]( vec_eq_imp, vec_unk ) = C[ i ][ d ]( vec_eq_imp, vec_unk ) * 1.;
                K[ i ][ d ]( vec_eq_imp_to_unk, vec_unk ) = B[ i ][ d ]( vec_unk_minus_eq_imp, vec_unk ) * 1.;
                F[ i ][ d ][ vec_eq_imp ] = q[ i ][ d ][ vec_eq_imp ] * 1.;
                F[ i ][ d ][ vec_eq_imp_to_unk ] = r[ i ][ d ][ vec_unk_minus_eq_imp ] * 1.;
            }
            else {
                Vec<unsigned> vec_unk = range( nb_unk[ i ][ d ] );
                Vec<unsigned> vec_eq_imp = range( nb_eq_imp[ i ][ d ] );
                Vec<unsigned> vec_eq_indep = range( nb_eq_indep[ i ][ d ] );
                Vec<unsigned> vec_unk_to_unk_plus_eq_imp = range( nb_unk[ i ][ d ], nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] );
                Vec<unsigned> vec_unk_plus_eq_imp_to_unk_plus_eq_imp_plus_eq_indep = range( nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ], nb_unk[ i ][ d ] + nb_eq_imp[ i ][ d ] + nb_eq_indep[ i ][ d ] );
                TMatGenSparse trans_B = trans( B[ i ][ d ] );
                TMatGenSparse trans_C = trans( C[ i ][ d ] );
                K[ i ][ d ]( vec_unk, vec_unk ) = M[ i ][ d ][ vec_unk ] * 1.;
                K[ i ][ d ]( vec_unk_to_unk_plus_eq_imp, vec_unk ) = C[ i ][ d ]( vec_eq_imp, vec_unk ) * 1.;
                K[ i ][ d ]( vec_unk_plus_eq_imp_to_unk_plus_eq_imp_plus_eq_indep, vec_unk ) = B[ i ][ d ]( vec_eq_indep, vec_unk ) * 1.;
                K[ i ][ d ]( vec_unk, vec_unk_to_unk_plus_eq_imp ) = trans_C( vec_unk, vec_eq_imp ) * 1.;
                K[ i ][ d ]( vec_unk, vec_unk_plus_eq_imp_to_unk_plus_eq_imp_plus_eq_indep ) = trans_B( vec_unk, vec_eq_indep ) * 1.;
                F[ i ][ d ][ vec_unk ] = ( M[ i ][ d ] * b[ i ][ d ] )[ vec_unk ] * 1.;
                F[ i ][ d ][ vec_unk_to_unk_plus_eq_imp ] = q[ i ][ d ][ vec_eq_imp ] * 1.;
                F[ i ][ d ][ vec_unk_plus_eq_imp_to_unk_plus_eq_imp_plus_eq_indep ] = r[ i ][ d ][ vec_eq_indep ] * 1.;
            }
        }
    }
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension de la matrice K associee au noeud " << i << " dans la direction " << d << " = ( " << K[ i ][ d ].nb_rows() << ", " << K[ i ][ d ].nb_cols() << " ) "<< endl;
                cout << "matrice K associe au noeud " << i << " dans la direction " << d << " =" << endl;
                cout << K[ i ][ d ] << endl;
            }
        }
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur F associe au noeud " << i << " dans la direction " << d << " = " << F[ i ][ d ].size() << endl;
                cout << "vecteur F associe au noeud " << i << " dans la direction " << d << " =" << endl;
                cout << F[ i ][ d ] << endl << endl;
            }
        }
    }
    
    nb_eq.free();
    r.free();
    B.free();
    nb_eq_imp.free();
    C.free();
    q.free();
    nb_eq_indep.free();
    M.free();
    b.free();
    node_type.free();
    
    /// Resolution des problemes de minimisation K[ i ][ d ] * U[ i ][ d ] = F[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// Construction des vecteurs U[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// --------------------------------------------------------------------------------------------------------------------------------------
    
    if ( disp ) {
        cout << "Resolution des pbs de minimisation K * U = F" << endl;
        cout << "Construction des vecteurs U" << endl << endl;
    }
    
    TicToc t_solve_minimization;
    t_solve_minimization.start();
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        for (unsigned d=0;d<dim;++d) {
            if ( K[ i ][ d ].nb_rows() != 0 ) {
                if ( not minimisation[ i ][ d ] ) {
                    #ifdef WITH_UMFPACK
                    TMatGenSparseUMFPACK K_UMFPACK = K[ i ][ d ];
                    K_UMFPACK.get_factorization();
                    U[ i ][ d ] = K_UMFPACK.solve( F[ i ][ d ] );
                    #endif
                }
                else {
                    if ( solver_minimisation == "LDL" ) {
                        #ifdef WITH_LDL
                        TMatSymSparse K_LDL = K[ i ][ d ];
                        U[ i ][ d ] = F[ i ][ d ];
                        LDL_solver ls;
                        Vec< Vec<T> > Ker;
                        Vec<int> Pivots;
                        ls.get_factorization( K_LDL, Ker, Pivots );
                        ls.solve( U[ i ][ d ] );
                        #endif
                    }
                    else if ( solver_minimisation == "UMFPACK" ) {
                        #ifdef WITH_UMFPACK
                        TMatGenSparseUMFPACK K_UMFPACK = K[ i ][ d ];
                        K_UMFPACK.get_factorization();
                        U[ i ][ d ] = K_UMFPACK.solve( F[ i ][ d ] );
                        #endif
                    }
                    else if ( solver_minimisation == "Inv" ) {
                        TMatSymSparse K_Inv = K[ i ][ d ];
                        U[ i ][ d ] = inv( K_Inv ) * F[ i ][ d ];
                    }
                    else if ( solver_minimisation == "LUFactorize" ) {
                        TMatSymSparse K_sym = K[ i ][ d ];
                        TMatGenSparseLU K_LU = K_sym;
//                        Vec<int> vector_permutation;
                        lu_factorize( K_LU/*, vector_permutation*/ );
                        solve_using_lu_factorize( K_LU, /*vector_permutation,*/ F[ i ][ d ], U[ i ][ d ] );
                    }
                    else {
                        cerr << "Bing. Error : solveur " << solver_minimisation << " pour la minimisation non implemente" << endl << endl;
                    }
                }
            }
        }
    }
    
    t_solve_minimization.stop();
    if ( disp )
        cout << "temps de calcul de la resolution des pbs de minimisation = " << t_solve_minimization.res << endl << endl;
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur U associe au noeud " << i << " dans la direction " << d << " = " << U[ i ][ d ].size() << endl;
                cout << "vecteur U associe au noeud " << i << " dans la direction " << d << " =" << endl;
                cout << U[ i ][ d ] << endl << endl;
            }
        }
    }
    
    if ( verif_solver_minimisation ) {
        if ( disp )
            cout << "Verification de la resolution des pbs de minimisation : tolerance = " << tol_solver_minimisation << endl << endl;
        
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                T residual = norm_2( K[ i ][ d ] * U[ i ][ d ] - F[ i ][ d ] );
                T b = norm_2( F[ i ][ d ] );
                if ( residual / b > tol_solver_minimisation ) {
                    cout << "residu associe au noeud " << i << " dans la direction " << d << " :" << endl;
//                    cout << "K * U - F =" << endl;
//                    cout << K[ i ][ d ] * U[ i ][ d ] - F[ i ][ d ] << endl;
                    cout << "norme du residu = " << residual << endl;
                    cout << "norme du residu relatif = " << residual / b << endl << endl;
                }
            }
        }
    }
    
    K.free();
    F.free();
    minimisation.free();
    
    /// Construction des vecteurs de projection b_hat[ i ][ d ] pour chaque noeud i du maillage et chaque direction d
    /// -------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs b_hat" << endl << endl;
    
    Vec< Vec< Vec<T> > > b_hat;
    b_hat.resize( m.node_list.size() );
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        b_hat[ i ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            b_hat[ i ][ d ].resize( nb_unk[ i ][ d ] );
            b_hat[ i ][ d ].set( 0. );
        }
    }
    
    for (unsigned i=0;i<m.node_list.size();++i) {
        for (unsigned d=0;d<dim;++d) {
            Vec<unsigned> vec_unk = range( nb_unk[ i ][ d ] );
            b_hat[ i ][ d ][ vec_unk ] = U[ i ][ d ][ vec_unk ] * 1.;
        }
    }
    
    Reset_Nodal_Vector_b_hat reset_nodal_vector_b_hat;
    reset_nodal_vector_b_hat.enhancement = &enhancement;
    reset_nodal_vector_b_hat.flag_face_enh = &flag_face_enh;
    reset_nodal_vector_b_hat.face_ind = &face_ind;
    
    apply( m.sub_mesh(Number<1>()).elem_list, reset_nodal_vector_b_hat, m, b_hat );
    
    if ( disp ) {
        for (unsigned i=0;i<m.node_list.size();++i) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur b_hat associe au noeud " << i << " dans la direction " << d << " = " << nb_unk[ i ][ d ] << endl;
                cout << "vecteur b_hat associe au noeud " << i << " dans la direction " << d << " =" << endl;
                cout << b_hat[ i ][ d ] << endl << endl;
            }
        }
    }
    
    nb_unk.free();
    U.free();
    correspondance_node_to_vertex_node.free();
    
    /// Construction des vecteurs de projection b_face[ k ][ d ] pour chaque face k du maillage et chaque direction d
    /// -------------------------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des vecteurs b_face" << endl << endl;
    
    Vec< Vec< Vec<T> > > b_face;
    b_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        b_face[ k ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            b_face[ k ][ d ].resize( m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
            b_face[ k ][ d ].set( 0. );
        }
    }
    
    apply( m.sub_mesh(Number<1>()).elem_list, Calcul_Skin_Elem_Vector_b_face(), face_ind, b_hat, b_face );
    
    if ( disp ) {
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur b_face associe a la face " << k << " dans la direction " << d << " = " << m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << endl;
                cout << "vecteur b_face associe a la face " << k << " dans la direction " << d << " =" << endl;
                cout << b_face[ k ][ d ] << endl << endl;
            }
        }
    }
    
    face_ind.free();
    b_hat.free();
    
    /// ----------------------------------------------------------------------------------- ///
    /// Construction standard (de la partie standard si amelioration) des densites d'effort ///
    /// ----------------------------------------------------------------------------------- ///
    
    if ( disp ) {
        cout << "Construction standard des densites d'effort" << endl;
        cout << "-------------------------------------------" << endl << endl;
    }
    
    /// Construction des matrices K_face[ k ][ d ] pour chaque face k du maillage et chaque direction d
    /// -----------------------------------------------------------------------------------------------
    
    if ( disp )
        cout << "Construction des matrices K_face" << endl << endl;
    
    Vec< Vec< TMatGenSparseUMFPACK > > K_face;
    K_face.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        K_face[ k ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            K_face[ k ][ d ].resize( m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
        }
    }
    
    apply( m.sub_mesh(Number<1>()).elem_list, Calcul_Skin_Elem_Matrix_K_face(), K_face );
    
    if ( disp ) {
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension de la matrice K_face associe a la face " << k << " dans la direction " << d << " = ( " << m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << ", " << m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << " )" << endl;
                cout << "matrice K_face associe a la face " << k << " dans la direction " << d << " =" << endl;
                cout << K_face[ k ][ d ] << endl;
            }
        }
    }
    
    /// Resolution des systemes lineaires K_face[ k ][ d ] * force_fluxes[ k ][ d ] = b_face[ k ][ d ]
    /// Construction des vecteurs force_fluxes[ k ][ d ] pour chaque face k du maillage et chaque direction d
    /// -----------------------------------------------------------------------------------------------------
    
    if ( disp ) {
        cout << "Resolution des systemes lineaires K_face * force_fluxes_std = b_face" << endl;
        cout << "Construction des vecteurs force_fluxes_std" << endl << endl;
    }
    
    force_fluxes.resize( m.sub_mesh(Number<1>()).elem_list.size() );
    
    for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
        force_fluxes[ k ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            force_fluxes[ k ][ d ].resize( m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() );
//            force_fluxes[ k ][ d ].set( 0. );
            K_face[ k ][ d ].get_factorization();
            force_fluxes[ k ][ d ] = K_face[ k ][ d ].solve( b_face[ k ][ d ] );
        }
    }
    
    if ( disp ) {
        for (unsigned k=0;k<m.sub_mesh(Number<1>()).elem_list.size();++k) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur des densites d'effort standard associe a la face " << k << " dans la direction " << d << " = " << m.sub_mesh(Number<1>()).elem_list[k]->nb_nodes_virtual() << endl;
                cout << "vecteur des densites d'effort standard associe a la face " << k << " dans la direction " << d << " =" << endl;
                cout << force_fluxes[ k ][ d ] << endl << endl;
            }
        }
    }
    
    t_force_fluxes_std.stop();
    cout << "temps de calcul de la construction standard des densites d'effort = " << t_force_fluxes_std.res << endl << endl;
    
}

#endif // Construct_standard_force_fluxes_EET_PGD_h
