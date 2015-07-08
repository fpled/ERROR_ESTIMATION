#ifndef Test_cpp_h
#define Test_cpp_h

// Eta
//----
struct Eta {
    template<class TE,class TM> void operator()( const TE &elem, TM &m ) const {
        for (unsigned i=0;i<3;++i) {
            PRINT( elem );
            PRINT( i ); // 0, 1 ou 2 pour un Triangle
            m.update_elem_children();
            PRINT( *m.get_children_of( elem, Number<1>() )[i] ); // *m.get_children_of( elem, Number<1>() )[0] renvoie Bar 1 2 pour elem : Triangle 1 2 0
                                                                 // *m.get_children_of( elem, Number<1>() )[1] renvoie Bar 2 0 pour elem : Triangle 1 2 0
                                                                 // *m.get_children_of( elem, Number<1>() )[2] renvoie Bar 0 1 pour elem : Triangle 1 2 0
                                                                 
                                                                 // *m.get_children_of( elem, Number<1>() )[0] renvoie Bar 1 3 pour elem : Triangle 1 3 2
                                                                 // *m.get_children_of( elem, Number<1>() )[1] renvoie Bar 2 0 pour elem : Triangle 1 3 2
                                                                 // *m.get_children_of( elem, Number<1>() )[2] renvoie Bar 0 1 pour elem : Triangle 1 3 2
//             PRINT( m.get_children_of( elem, Number<1>() )[i]->number ); // renvoie le numero de la i ème face (0, 1 ou 2) de elem dans le maillage
//             PRINT( *m.sub_mesh(Number<1>()).get_parents_of_EA( m.get_children_of( elem, Number<1>() )[i] )[0] ); // *m.sub_mesh(Number<1>()).get_parents_of_EA( m.get_children_of( elem, Number<1>() )[0] )[0] renvoie le 1er element ([0]) auquel est connecte le 1er enfant ([0]) de elem, c'est-a-dire Triangle 1 2 0
            int eta = ( m.sub_mesh(Number<1>()).get_parents_of_EA( m.get_children_of( elem, Number<1>() )[i] )[0] == &elem ? 1 : -1 ); // m.get_children_of( elem, Number<1>() )[i] renvoie l'adresse d'un pointeur vers le i ème enfant de elem (il s'agit de l'adresse d'un pointeur vers un ElementAncestor) //   m.sub_mesh(Number<1>()).get_parents_of_EA( truc )[0] renvoie le 1er element auquel est connecte truc (qui est un ElementAncestor); s'il s'agit de elem , eta = 1, sinon eta = -1
            PRINT( eta );
//             m.sub_mesh(Number<1>()).get_parents_of( child_elem )[i] renvoie le ieme parent de l'interface child_elem (definie comme un Element)
//             cout << *m.sub_mesh(Number<1>()).get_parents_of_EA( m.sub_mesh(Number<1>()).elem_list[i] )[0] << endl << endl; renvoie le 1er parent de la ieme face du maillage
        }
    }
};

// Vec<T> eig_val;
// Mat<T> eig_vec;
// 
// get_eig_sym( A, eig_val, eig_vec );
// 
// cout << "nb valeurs propres de A = " << eig_val.size() << endl;
// cout << "valeurs propres de A : " << eig_val << endl << endl;
// 
// get_eig_sym( B, eig_val, eig_vec );
// 
// cout << "nb valeurs propres de B = " << eig_val.size() << endl;
// cout << "valeurs propres de B : " << eig_val << endl << endl;

// cout << m.get_node_neighbours( 100 ).size() << endl;
// cout << *m.get_node_neighbours( 100 )[ 0 ] << endl;
// cout << m.get_node_neighbours( 100 )[ 0 ]->number << endl;
// cout << *m.get_node_neighbours( 100 )[ 1 ] << endl;
// cout << m.get_node_neighbours( 100 )[ 1 ]->number << endl;


// m.sub_mesh(Number<1>()).update_node_parents();
// for (unsigned k=0;k<m.sub_mesh(Number<1>()).get_node_parents( node.number ).size();++k) {
//     if ( (*type_face)[ m.sub_mesh(Number<1>()).get_node_parents( node.number )[ k ]->number ][ d ] == 2 ) {
//         C_i_ind[ node.number ][ d ][ k ] += nb_eq_imp[ node.number ][ d ];
//         nb_eq_imp[ node.number ][ d ]++;
//     }
// }


// for(unsigned n=0;n<m.elem_list.size();++n) {
//     Mat<T,Sym<dim> > epsilon_mat, epsilon_mat_new, sigma_mat, sigma_mat_new;
//     //m.elem_list[n]->set_field( "epsilon", cos( 3.14159 * i / 7 ) );
//     epsilon_mat  = m.elem_list[n]->get_field( "epsilon", StructForType<Mat<T,Sym<dim> > >() );
//     epsilon_mat_new  = m_new.elem_list[n]->get_field( "epsilon", StructForType<Mat<T,Sym<dim> > >() );
//     sigma_mat = m.elem_list[n]->get_field( "sigma", StructForType<Mat<T,Sym<dim> > >() );
//     sigma_mat_new = m_new.elem_list[n]->get_field( "sigma", StructForType<Mat<T,Sym<dim> > >() );
//     PRINT( epsilon_mat );
//     PRINT( epsilon_mat_new );
//     PRINT( sigma_mat );
//     PRINT( sigma_mat_new );
//     
//     Vec<T,unsigned(dim*(dim+1)/2) > epsilon_vec, epsilon_vec_new, sigma_vec, sigma_vec_new;
//     epsilon_vec  = m.elem_list[n]->get_field( "epsilon", StructForType<Vec<T,unsigned(dim*(dim+1)/2)> >() );
//     epsilon_vec_new  = m_new.elem_list[n]->get_field( "epsilon", StructForType<Vec<T,unsigned(dim*(dim+1)/2)> >() );
//     sigma_vec = m.elem_list[n]->get_field( "sigma", StructForType<Vec<T,unsigned(dim*(dim+1)/2)> >() );
//     sigma_vec_new = m_new.elem_list[n]->get_field( "sigma", StructForType<Vec<T,unsigned(dim*(dim+1)/2)> >() );
//     PRINT( epsilon_vec );
//     PRINT( epsilon_vec_new );
//     PRINT( sigma_vec );
//     PRINT( sigma_vec_new );
// }

// TNode *n[] = {
//     ban.get_node( e.node(0)->pos ),
//     ban.get_node( e.node(1)->pos ),
//     ban.get_node( e.node(2)->pos ),
//     ban.get_node( e.node(3)->pos ),
//     ban.get_node( e.node(4)->pos ),
//     ban.get_node( e.node(5)->pos ),
//     ban.get_node( e.node(6)->pos ),
//     ban.get_node( e.node(7)->pos )
// };
// m_ref.add_element( Hexa(), DefaultBehavior(), node_Hexa[0], node_Hexa[1], node_Hexa[2], node_Hexa[3], node_Hexa[4], node_Hexa[5], node_Hexa[6], node_Hexa[7] );
// m_ref.add_element( Hexa(), DefaultBehavior(), n );

// m.update_elem_children();
// m.update_elem_children( Number<2>() );
// for (unsigned i=0;i<6;++i) {
//     PRINT( *m.get_children_of( elem, Number<1>() )[i] );
// }
// for (unsigned i=0;i<12;++i) {
//     PRINT( *m.get_children_of( elem, Number<2>() )[i] );
// }

// const double* data = gauss_point_for_order( 5, Tetra() );
// unsigned cpt = 1;
// for (unsigned i=0; i<100; i+=(1+dim)) {
//     if ( data[i] != 0 ) {
//         cout << "point de gauss : " << cpt << endl;
//         cout << "poids : " << data[i] << endl;
//         for (unsigned d=0; d<dim; ++d) {
//             cout << "valeur["+to_string(d)+"] : " << data[i+d+1] << endl;
//         }
//         cout << endl ;
//     }
//     else break;
//     cpt++;
// }



//     typedef ImgInterp<double,2> TI;
//     TI bin( "/home/leclerc/Data/Croix/masque_0.png" );
//     TI cut = img_dist_from_front( bin, 100, 128.0 );
//     TI stp; stp.resize( cut.sizes, -1 );
//     // premier raffinement
//     Raf raf;
//     raf.cut = cut;
//     raf.v = 50;
//     refinement( m, raf, true );
//     raf.v = 25;
//     refinement( m, raf, true );
//     // coupe
//     LevelSetImageRefinement<TI> lr( cut, stp );
//     refinement( m, lr, true );
//     display( m );


//  struct Raf {
//     template<class TE>
//     double operator()( const TE &e ) const {
//         return abs( cut( center( e ) ) ) < v;
//     }
//     ImgInterp<double,2> cut;
//     double v;
// };


    /// Construction de la correspondance entre noeuds sommets et noeuds sommets appartenant a delta_Omega : correspondance_vertex_node_to_skin_vertex_node[ j ]
    ///---------------------------------------------------------------------------------------------------------------------------------------------------------

    cout << "Construction de la correspondance entre noeuds sommets et noeuds sommets sur delta_Omega" << endl << endl;

    Vec<bool> correspondance_vertex_node_to_skin_vertex_node;
    correspondance_vertex_node_to_skin_vertex_node.resize( nb_vertex_nodes, 1 );

    m.update_skin();
    apply( m.skin.node_list, Construct_Correspondance_Vertex_Node_To_Skin_Vertex_Node(), correspondance_node_to_vertex_node, connect_node_to_vertex_node, correspondance_vertex_node_to_skin_vertex_node ); // si correspondance_vertex_node_to_skin_vertex_node[ j ] != 0, alors le noeud sommet j est un noeud interne ; sinon le noeud sommet j appartient a delta_Omega

    if ( debug_method ) {
        cout << "correspondance entre noeuds sommets et noeuds sommets du maillage appartenant a delta Omega : " << correspondance_vertex_node_to_skin_vertex_node << endl << endl;
    }

    /// Reperage pour chaque face k et chaque direction d de l'indice de debut de ligne dans les matrices C[ j ][ d ] et dans les vecteurs q[ j ][ d ] : C_ind[ j ][ d ][ k ]
    /// Calcul du nb de lignes de la matrice C[ j ][ d ] et du vecteur q[ j ][ d ] : nb_eq_imp[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if ( boundary_condition_force_EESPT == "lagrange" ) {
        cout << "Calcul du vecteur nb_eq_imp" << endl << endl;
    }

    Vec< Vec<unsigned> > nb_eq_imp;
    nb_eq_imp.resize( nb_vertex_nodes );

    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        nb_eq_imp[ j ].resize( dim, 0 );
    }

    Vec< Vec< Vec<unsigned> > > C_ind; 
    C_ind.resize( nb_vertex_nodes );

    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        C_ind[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            C_ind[ j ][ d ].resize( m.sub_mesh(Number<1>()).elem_list.size(), 0 );
        }
    }

    Construction_C_Ind constr_C_ind;
    constr_C_ind.cpt_nodes_face = &cpt_nodes_face;
    constr_C_ind.cpt_faces_node = &cpt_faces_node;
    constr_C_ind.list_faces_node = &list_faces_node;
    constr_C_ind.correspondance_node_to_vertex_node = &correspondance_node_to_vertex_node;
    constr_C_ind.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    constr_C_ind.boundary_condition_force_EESPT = &boundary_condition_force_EESPT;

    apply( m.node_list, constr_C_ind, type_face, nb_eq_imp, C_ind ); // nb_eq_imp[ j ][ d ] contient le nb de lignes de la matrice C[ j ][ d ] et du vecteur q[ j ][ d ] associee au j eme noeud sommet dans la direction d // C_ind[ j ][ d ][ k ] = indice de debut de ligne dans la matrice C[ j ][ d ] et dans le vecteur q[ j ][ d ] associes au jeme noeud sommet pour chaque face k du maillage et chaque direction d

    if ( debug_method and boundary_condition_force_EESPT == "lagrange" ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "nb d'equations imposees C associees au noeud sommet " << j << " dans la direction " << d << " : " << nb_eq_imp[ j ][ d ] << endl;
                for (unsigned k=0;k<cpt_faces_vertex_node[ j ];++k) {
                    if ( type_face[ list_faces_vertex_node[ j ][ k ] ][ d ] == 2 ) {
                        cout << "indice de debut de ligne de la face " << list_faces_vertex_node[ j ][ k ] << " dans la direction " << d << " dans la matrice C[ " << j << " ][ " << d <<  " ] : " << C_ind[ j ][ d ][ list_faces_vertex_node[ j ][ k ] ] << endl;
                    }
                }
                cout << endl << endl;
            }
        }
    }

    /// Construction des matrices C[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    ///---------------------------------------------------------------------------------------------------

    if ( boundary_condition_force_EESPT == "lagrange" ) {
        cout << "Construction des matrices C" << endl << endl;
    }

    Vec< Vec< Mat<double, Gen<>, SparseLine<> > > > C;
    C.resize( nb_vertex_nodes );

    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        C[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            C[ j ][ d ].resize( nb_eq_imp[ j ][ d ], nb_unk[ j ][ d ] );
        }
    }

    Construction_C constr_C;
    constr_C.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    constr_C.type_face  = &type_face;
    constr_C.C_ind  = &C_ind;
    constr_C.lambda_i_Fh_hat_ind  = &lambda_i_Fh_hat_ind;
    constr_C.boundary_condition_force_EESPT = &boundary_condition_force_EESPT;

    apply( m.sub_mesh(Number<1>()).elem_list, constr_C, m, C );

    if ( debug_method and boundary_condition_force_EESPT == "lagrange" ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension de la matrices C associe au noeud sommet " << j << " dans la direction " << d <<  " : ( " << nb_eq_imp[ j ][ d ] << ", " << nb_unk[ j ][ d ] << " )" << endl;
                cout << "matrice C associe au noeud sommet " << j << " dans la direction " << d << " :" << endl;
                cout << C[ j ][ d ] << endl << endl;
            }
        }
    }

    /// Construction des vecteurs q[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    ///---------------------------------------------------------------------------------------------------

    if ( boundary_condition_force_EESPT == "lagrange" ) {
        cout << "Construction des vecteurs q" << endl << endl;
    }

    Vec< Vec< Vec<double> > > q;
    q.resize( nb_vertex_nodes );

    for (unsigned j=0;j<nb_vertex_nodes;++j) {
        q[ j ].resize( dim );
        for (unsigned d=0;d<dim;++d) {
            q[ j ][ d ].resize( nb_eq_imp[ j ][ d ], 0. );
        }
    }

    Construction_q constr_q;
    constr_q.connect_node_to_vertex_node = &connect_node_to_vertex_node;
    constr_q.type_face  = &type_face;
    constr_q.lambda_i_Fh_hat_ind = &lambda_i_Fh_hat_ind;
    constr_q.C_ind  = &C_ind;
    constr_q.boundary_condition_force_EESPT = &boundary_condition_force_EESPT;

    apply( m.sub_mesh(Number<1>()).elem_list, constr_q, lambda_i_Fh, q );

    if ( debug_method and boundary_condition_force_EESPT == "lagrange" ) {
        for (unsigned j=0;j<nb_vertex_nodes;++j) {
            for (unsigned d=0;d<dim;++d) {
                cout << "dimension du vecteur q associe au noeud sommet " << j << " dans la direction " << d << " : " << nb_eq_imp[ j ][ d ] << endl;
                cout << "vecteur q associe au noeud sommet " << j << " dans la direction " << d << " :" << endl;
                cout << q[ j ][ d ] << endl << endl;
            }
        }
    }

/// Reperage pour chaque face k et chaque direction d de l'indice de debut de ligne dans les matrices C[ j ][ d ] et dans les vecteurs q[ j ][ d ] : C_ind[ j ][ d ][ k ]
/// Calcul du nb de lignes de la matrice C[ j ][ d ] et du vecteur q[ j ][ d ] : nb_eq_imp[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Construction_C_Ind {
    const Vec<unsigned>* cpt_nodes_face;
    const Vec<unsigned>* cpt_faces_node;
    const Vec< Vec<unsigned> >* list_faces_node;
    const Vec<bool>* correspondance_node_to_vertex_node;
    const Vec<unsigned>* connect_node_to_vertex_node;
    const string* boundary_condition_force_EESPT;
    template<class TN> void operator()( const TN &node, const Vec< Vec<unsigned> > &type_face, Vec< Vec<unsigned> > &nb_eq_imp, Vec< Vec< Vec<unsigned> > > &C_ind ) const {
        if ( *boundary_condition_force_EESPT == "lagrange" ) {
            if ( (*correspondance_node_to_vertex_node)[ node.number ] ) {
                for (unsigned d=0;d<TN::dim;++d) {
                    for (unsigned k=0;k<(*cpt_faces_node)[ node.number ];++k) {
                        if ( type_face[ (*list_faces_node)[ node.number ][ k ] ][ d ] == 2 ) {
                            C_ind[ (*connect_node_to_vertex_node)[ node.number ] ][ d ][ (*list_faces_node)[ node.number ][ k ] ] += nb_eq_imp[ (*connect_node_to_vertex_node)[ node.number ] ][ d ];
                            nb_eq_imp[ (*connect_node_to_vertex_node)[ node.number ] ][ d ] += (*cpt_nodes_face)[ (*list_faces_node)[ node.number ][ k ] ];
                        }
                    }
                }
            }
        }
    }
};

/// Construction des matrices C[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
///---------------------------------------------------------------------------------------------------
template<class TE, class TM, class TV, class TVV, class TVVV, class S, class TTMVV> 
void construction_C( const TE &child_elem , const TM &m, const TV &connect_node_to_vertex_node, const TVV &type_face, const TVVV &C_ind, const TVVV &lambda_i_Fh_hat_ind, const S &boundary_condition_force_EESPT, TTMVV &C ) {}

struct Construction_C {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* type_face;
    const Vec< Vec< Vec<unsigned> > >* C_ind;
    const Vec< Vec< Vec<unsigned> > >* lambda_i_Fh_hat_ind;
    const string* boundary_condition_force_EESPT;
    template<class TE, class TM> void operator()( const TE &child_elem, const TM &m, Vec< Vec< Mat<double, Gen<>, SparseLine<> > > > &C ) const {
        construction_C( child_elem, m, *connect_node_to_vertex_node, *type_face, *C_ind, *lambda_i_Fh_hat_ind, *boundary_condition_force_EESPT, C );
    }
};

/// Construction des vecteurs q[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
///---------------------------------------------------------------------------------------------------
template<class TE, class TV, class TVV, class TVVV, class S, class TTVVV> 
void construction_q( const TE &child_elem , const TV &connect_node_to_vertex_node, const TVV &type_face, const TVVV &lambda_i_Fh_hat_ind, const TVVV &C_ind, const TTVVV &lambda_i_Fh, const S &boundary_condition_force_EESPT, TTVVV &q ) {}

struct Construction_q {
    const Vec<unsigned>* connect_node_to_vertex_node;
    const Vec< Vec<unsigned> >* type_face;
    const Vec< Vec< Vec<unsigned> > >* lambda_i_Fh_hat_ind;
    const Vec< Vec< Vec<unsigned> > >* C_ind;
    const string* boundary_condition_force_EESPT;
    template<class TE> void operator()( const TE &child_elem, const Vec< Vec< Vec<double> > > &lambda_i_Fh, Vec< Vec< Vec<double> > > &q ) const {
        construction_q( child_elem, *connect_node_to_vertex_node, *type_face, *lambda_i_Fh_hat_ind, *C_ind, lambda_i_Fh, *boundary_condition_force_EESPT, q );
    }
};

    /// Construction de la liste des noeuds soumis a des contraintes cinematiques sur chaque patch : constrained_nodes_patch[ j ] pour chaque noeud sommet j du maillage
    ///-----------------------------------------------------------------------------------------------------------------------------------------------------------------
// struct Construct_Constrained_Nodes_Patch {
//     Vec<unsigned>* correspondance_node_to_vertex_node;
//     Vec<unsigned>* connect_node_to_vertex_node;
//     template<class TE, class TF> void operator()( const TE &elem, const TF &f, Vec< Vec<unsigned> > &constrained_nodes_patch ) const {
//         unsigned elem_nb_nodes = TE::nb_nodes; // nb de noeuds de l'element elem
//         for (unsigned n=0;n<elem_nb_nodes;++n) {
//             if ( (*correspondance_node_to_vertex_node)[ elem.node( n )->number ] ) {
//                 for (unsigned i=0;i<elem_nb_nodes;++i) {
//                     if ( f.constrained_nodes()[ elem.node( i )->number ] ) {
//                         constrained_nodes_patch[ (*connect_node_to_vertex_node)[ elem.node( n )->number ] ].push_back( elem.node( i )->number );
//                     }
//                 }
//             }
//         }
//     }
// };

    /// Construction base lineaire sur le patch associe au noeud sommet j
    ///------------------------------------------------------------------
//         Vec< Mat<double> > R;
//         R.resize( nb_vertex_nodes );
// 
//         Vec< Mat<unsigned, Gen<dim,dim>, SparseLine<> > > F_linear;
//         F_linear.resize( dim * ( dim ) );
// 
//         for (unsigned d=0; d<dim; d++) {
//             for (unsigned e=0; e<dim; e++) {
//                 F_linear[ d * dim + e ]( d, e ) = 1;
//             }
//         }
// 
//         if ( debug_method ) {
//             for (unsigned d=0;d<(dim * ( dim ));++d) {
//                 cout << "dimension de la matrice F_linear " << d << " : ( " << dim << " , " << dim << " )" << endl;
//                 cout << "matrice F_linear " << d << " :" << endl;
//                 cout << F_linear[ d ] << endl << endl;
//             }
//         }
// 
//         for (unsigned j=0;j<nb_vertex_nodes;++j) {
//             R[ j ].resize( nb_points_patch[ j ] * dim, dim * ( dim + 1 ) );
//             R[ j ].set( 0 );
//             for (unsigned p=0;p<nb_points_patch[ j ];++p) {
//                 for (unsigned d=0;d<dim;++d) {
//                     R[ j ].col( d )[ p * dim + d ] = 1;
//                     for (unsigned e=0;e<dim;++e) {
//                         for (unsigned f=0;f<dim;++f) {
//                             R[ j ].col( dim + d * dim + e )[ p * dim + f ] = ( F_linear[ d * dim + e ] * pos_patch[ j ][ p ] )[ f ];
//                         }
//                     }
//                 }
//             }
//             orthonormalisation_schmidt( R[ j ] );
//         }
// 
//         if ( debug_method ) {
//             for (unsigned j=0;j<nb_vertex_nodes;++j) {
//                 cout << "dimension de la base des mouvements lineaires R associe au noeud sommet " << j << " : " << nb_points_patch[ j ] * dim << " , " << dim * ( dim + 1 ) << endl;
//                 cout << "matrice R associe au noeud sommet " << j << " :" << endl;
//                 cout << R[ j ] << endl << endl;
//             }
//         }
// 
//         for (unsigned j=0;j<nb_vertex_nodes;++j) {
//             F[ j ] -= R[ j ] * trans(R[ j ]) * F[ j ];
//         }

//     cout << f.matrices(Number<0>()) << endl << endl;
//     cout << f.sollicitation << endl << endl;
//     cout << f.vectors[0] << endl << endl;
//     cout << f.matrices(Number<0>()).size() << endl << endl;
//     cout << f.sollicitation.size() << endl << endl;
//     cout << f.vectors[0].size() << endl << endl;
//     cout << sqrt( dot( f.matrices(Number<0>()) * f.vectors[0] - f.sollicitation, f.matrices(Number<0>()) * f.vectors[0] - f.sollicitation ) ) << endl << endl;
//     m.f_vol[1] = -0.01;
//     cout << m.node_list[0].dep << endl; // valeur du deplacement du noeud 0 du maillage



//     TicToc t_EESPT_solve_0;
//     t_EESPT_solve_0.start();

    /// Resolution de A_tilde[ j ][ d ]^T * P[ j ][ d ]^-1 * A_tilde[ j ][ d ] * vh_star[ j ][ d ] = R_tilde[ j ][ d ] - A_tilde[ j ][ d ]^T * lambda_i_Fh[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    /// Construction des vecteurs vh_star[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    ///-----------------------------------------------------------------------------------------------------------------------------------------------------------------
//     cout << "Resolution des systemes A_tilde^T * P^-1 * A_tilde * vh_star = R_tilde - A_tilde^T * lambda_i_Fh" << endl;
//     cout << "Construction des vecteurs vh_star" << endl << endl;
// 
//     Vec< Vec< Vec<double> > > vh_star; 
//     vh_star.resize( nb_vertex_nodes );
// 
//     for (unsigned j=0;j<nb_vertex_nodes;++j) {
//         vh_star[ j ].resize( dim );
//         for (unsigned d=0;d<dim;++d) {
//             vh_star[ j ][ d ].resize( unk_vh_star_free[ j ][ d ].size(), 0. );
//         }
//     }
// 
//     Vec< Vec< Mat<double, Sym<>, SparseLine<> > > > inv_P;
//     inv_P.resize( nb_vertex_nodes );
// 
//     for (unsigned j=0;j<nb_vertex_nodes;++j) {
//         inv_P[ j ].resize( dim );
//         for (unsigned d=0;d<dim;++d) {
//             inv_P[ j ][ d ].resize( nb_unk_lambda_i_Fh_hat[ j ][ d ] );
//         }
//     }
//     for (unsigned j=0;j<nb_vertex_nodes;++j) {
//         for (unsigned d=0;d<dim;++d) {
//             for (unsigned i=0;i<nb_unk_lambda_i_Fh_hat[ j ][ d ];++i) {
//                 inv_P[ j ][ d ]( i, i ) = 1/( P[ j ][ d ][ i ] );
//             }
//         }
//     }
// 
//     P.free();
// 
//     Vec< Vec< Mat<double, Sym<>, SparseLine<> > > > M; // SparseLine<> : rapide a completer
//     M.resize( nb_vertex_nodes );
// 
//     for (unsigned j=0;j<nb_vertex_nodes;++j) {
//         M[ j ].resize( dim );
//         for (unsigned d=0;d<dim;++d) {
//             M[ j ][ d ].resize( unk_vh_star_free[ j ][ d ].size() );
//         }
//     }
// 
//     for (unsigned j=0;j<nb_vertex_nodes;++j) {
//         for (unsigned d=0;d<dim;++d) {
//             M[ j ][ d ] = trans(A_tilde[ j ][ d ]) * inv_P[ j ][ d ] * A_tilde[ j ][ d ];
//         }
//     }
// 
//     for (unsigned j=0;j<nb_vertex_nodes;++j) {
//         for (unsigned d=0;d<dim;++d) {
//             if ( solver == "CholMod" ) {
//                 Mat<double, Sym<>, SparseCholMod > N = M[ j ][ d ];
//                 N.get_factorization();
//                 vh_star[ j ][ d ] = N.solve( R_tilde[ j ][ d ] - trans(A_tilde[ j ][ d ]) * lambda_i_Fh[ j ][ d ] );
//                 N.clear();
//             }
//             else if ( solver == "LDL" ) {
//                 LDL_solver ls;
//                 ls.get_factorization( M[ j ][ d ] );
//                 vh_star[ j ][ d ] = R_tilde[ j ][ d ] - trans(A_tilde[ j ][ d ]) * lambda_i_Fh[ j ][ d ];
//                 ls.solve( vh_star[ j ][ d ] );
//             }
//             M[ j ][ d ].clear();
//         }
//         M[ j ].free();
//     }
// 
//     if ( debug_method ) {
//         for (unsigned j=0;j<nb_vertex_nodes;++j) {
//             for (unsigned d=0;d<dim;++d) {
//                 cout << "dimension du vecteur vh_star associee au noeud sommet " << j << " dans la direction " << d << " : " << unk_vh_star_free[ j ][ d ].size() << endl;
//                 cout << "vecteur vh_star associe au noeud sommet " << j << " dans la direction " << d << " :" << endl;
//                 cout << vh_star[ j ][ d ] << endl << endl;
//             }
//         }
//     }
// 
//     unk_vh_star_free.free();
//     R_tilde.free();
//     M.free();
    /// Construction des vecteurs lambda_i_Fh_hat[ j ][ d ] pour chaque noeud sommet j du maillage et chaque direction d
    ///-----------------------------------------------------------------------------------------------------------------
//     cout << "Construction des vecteurs lambda_i_Fh_hat" << endl << endl;
// 
//     Vec< Vec< Vec<double> > > lambda_i_Fh_hat; 
//     lambda_i_Fh_hat.resize( nb_vertex_nodes );
// 
//     for (unsigned j=0;j<nb_vertex_nodes;++j) {
//         lambda_i_Fh_hat[ j ].resize( dim );
//         for (unsigned d=0;d<dim;++d) {
//             lambda_i_Fh_hat[ j ][ d ].resize( nb_unk_lambda_i_Fh_hat[ j ][ d ], 0. );
//         }
//     }
// 
//     for (unsigned j=0;j<nb_vertex_nodes;++j) {
//         for (unsigned d=0;d<dim;++d) {
//             lambda_i_Fh_hat[ j ][ d ] = lambda_i_Fh[ j ][ d ] + inv_P[ j ][ d ] * A_tilde[ j ][ d ] * vh_star[ j ][ d ];
//         }
//     }
// 
//     if ( debug_method ) {
//         for (unsigned j=0;j<nb_vertex_nodes;++j) {
//             for (unsigned d=0;d<dim;++d) {
//                 cout << "dimension du vecteur lambda_i_Fh_hat associee au noeud sommet " << j << " dans la direction " << d << " : " << nb_unk_lambda_i_Fh_hat[ j ][ d ] << endl;
//                 cout << "vecteur lambda_i_Fh_hat associe au noeud sommet " << j << " dans la direction " << d << " :" << endl;
//                 cout << lambda_i_Fh_hat[ j ][ d ] << endl << endl;
//             }
//         }
//     }
// 
//     nb_unk_lambda_i_Fh_hat.free();
//     A_tilde.free();
//     lambda_i_Fh.free();
//     vh_star.free();
//     inv_P.free();
// 
//     t_EESPT_solve_0.stop();
//     cout << "Temps de calcul de la resolution 0 des systemes pour la technique EESPT : " << t_EESPT_solve_0.res << endl << endl;

#endif // Test_cpp_h
