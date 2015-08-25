#ifndef Test_cpp_h
#define Test_cpp_h

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
//            PRINT( m.get_children_of( elem, Number<1>() )[i]->number ); // renvoie le numero de la i ème face (0, 1 ou 2) de elem dans le maillage
//            PRINT( *m.sub_mesh(Number<1>()).get_parents_of_EA( m.get_children_of( elem, Number<1>() )[i] )[0] ); // *m.sub_mesh(Number<1>()).get_parents_of_EA( m.get_children_of( elem, Number<1>() )[0] )[0] renvoie le 1er element ([0]) auquel est connecte le 1er enfant ([0]) de elem, c'est-a-dire Triangle 1 2 0
            int eta = ( m.sub_mesh(Number<1>()).get_parents_of_EA( m.get_children_of( elem, Number<1>() )[i] )[0] == &elem ? 1 : -1 ); // m.get_children_of( elem, Number<1>() )[i] renvoie l'adresse d'un pointeur vers le i ème enfant de elem (il s'agit de l'adresse d'un pointeur vers un ElementAncestor) //   m.sub_mesh(Number<1>()).get_parents_of_EA( truc )[0] renvoie le 1er element auquel est connecte truc (qui est un ElementAncestor); s'il s'agit de elem , eta = 1, sinon eta = -1
            PRINT( eta );
//            m.sub_mesh(Number<1>()).get_parents_of( child_elem )[i] renvoie le ieme parent de l'interface child_elem (definie comme un Element)
//            cout << *m.sub_mesh(Number<1>()).get_parents_of_EA( m.sub_mesh(Number<1>()).elem_list[i] )[0] << endl << endl; renvoie le 1er parent de la ieme face du maillage
        }
    }
};

//Vec<T> eig_val;
//Mat<T> eig_vec;

//get_eig_sym( A, eig_val, eig_vec );

//cout << "nb valeurs propres de A = " << eig_val.size() << endl;
//cout << "valeurs propres de A : " << eig_val << endl << endl;

//cout << m.get_node_neighbours( 100 ).size() << endl;
//cout << *m.get_node_neighbours( 100 )[ 0 ] << endl;
//cout << m.get_node_neighbours( 100 )[ 0 ]->number << endl;
//cout << *m.get_node_neighbours( 100 )[ 1 ] << endl;
//cout << m.get_node_neighbours( 100 )[ 1 ]->number << endl;


//m.sub_mesh(Number<1>()).update_node_parents();
//for (unsigned k=0;k<m.sub_mesh(Number<1>()).get_node_parents( node.number ).size();++k) {
//    if ( (*face_type)[ m.sub_mesh(Number<1>()).get_node_parents( node.number )[ k ]->number ][ d ] == 2 ) {
//        C_i_ind[ node.number ][ d ][ k ] += nb_eq_imp[ node.number ][ d ];
//        nb_eq_imp[ node.number ][ d ]++;
//    }
//}


//for (unsigned n=0;n<m.elem_list.size();++n) {
//    Mat<T,Sym<dim> > epsilon_mat, epsilon_mat_new, sigma_mat, sigma_mat_new;
//    //m.elem_list[n]->set_field( "epsilon", cos( 3.14159 * i / 7 ) );
//    epsilon_mat  = m.elem_list[n]->get_field( "epsilon", StructForType<Mat<T,Sym<dim> > >() );
//    epsilon_mat_new  = m_new.elem_list[n]->get_field( "epsilon", StructForType<Mat<T,Sym<dim> > >() );
//    sigma_mat = m.elem_list[n]->get_field( "sigma", StructForType<Mat<T,Sym<dim> > >() );
//    sigma_mat_new = m_new.elem_list[n]->get_field( "sigma", StructForType<Mat<T,Sym<dim> > >() );
//    PRINT( epsilon_mat );
//    PRINT( epsilon_mat_new );
//    PRINT( sigma_mat );
//    PRINT( sigma_mat_new );

//    Vec<T,unsigned(dim*(dim+1)/2) > epsilon_vec, epsilon_vec_new, sigma_vec, sigma_vec_new;
//    epsilon_vec  = m.elem_list[n]->get_field( "epsilon", StructForType<Vec<T,unsigned(dim*(dim+1)/2)> >() );
//    epsilon_vec_new  = m_new.elem_list[n]->get_field( "epsilon", StructForType<Vec<T,unsigned(dim*(dim+1)/2)> >() );
//    sigma_vec = m.elem_list[n]->get_field( "sigma", StructForType<Vec<T,unsigned(dim*(dim+1)/2)> >() );
//    sigma_vec_new = m_new.elem_list[n]->get_field( "sigma", StructForType<Vec<T,unsigned(dim*(dim+1)/2)> >() );
//    PRINT( epsilon_vec );
//    PRINT( epsilon_vec_new );
//    PRINT( sigma_vec );
//    PRINT( sigma_vec_new );
//}

//TNode *n[] = {
//    ban.get_node( e.node(0)->pos ),
//    ban.get_node( e.node(1)->pos ),
//    ban.get_node( e.node(2)->pos ),
//    ban.get_node( e.node(3)->pos ),
//    ban.get_node( e.node(4)->pos ),
//    ban.get_node( e.node(5)->pos ),
//    ban.get_node( e.node(6)->pos ),
//    ban.get_node( e.node(7)->pos )
//};
//m_ref.add_element( Hexa(), DefaultBehavior(), node_Hexa[0], node_Hexa[1], node_Hexa[2], node_Hexa[3], node_Hexa[4], node_Hexa[5], node_Hexa[6], node_Hexa[7] );
//m_ref.add_element( Hexa(), DefaultBehavior(), n );

//m.update_elem_children();
//m.update_elem_children( Number<2>() );
//for (unsigned i=0;i<6;++i) {
//    PRINT( *m.get_children_of( elem, Number<1>() )[i] );
//}
//for (unsigned i=0;i<12;++i) {
//    PRINT( *m.get_children_of( elem, Number<2>() )[i] );
//}

//const double* data = gauss_point_for_order( 5, Tetra() );
//unsigned cpt = 1;
//for (unsigned i=0; i<100; i+=(1+dim)) {
//    if ( data[i] != 0 ) {
//        cout << "point de gauss : " << cpt << endl;
//        cout << "poids : " << data[i] << endl;
//        for (unsigned d=0; d<dim; ++d) {
//            cout << "valeur["+to_string(d)+"] : " << data[i+d+1] << endl;
//        }
//        cout << endl ;
//    }
//    else break;
//    cpt++;
//}


//struct Raf {
//    template<class TE>
//    double operator()( const TE &e ) const {
//        return abs( cut( center( e ) ) ) < v;
//    }
//    ImgInterp<double,2> cut;
//    double v;
//};

//typedef ImgInterp<double,2> TI;
//TI bin( "/home/leclerc/Data/Croix/masque_0.png" );
//TI cut = img_dist_from_front( bin, 100, 128.0 );
//TI stp; stp.resize( cut.sizes, -1 );
//Raf raf; // premier raffinement
//raf.cut = cut;
//raf.v = 50;
//refinement( m, raf, true );
//raf.v = 25;
//refinement( m, raf, true );
//LevelSetImageRefinement<TI> lr( cut, stp ); // coupe
//refinement( m, lr, true );
//display( m );


/// Construction de la liste des noeuds soumis a des contraintes cinematiques sur chaque patch : constrained_nodes_patch[ j ] pour chaque noeud sommet j du maillage
/// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
//struct Construct_Constrained_Nodes_Patch {
//    Vec<unsigned>* correspondance_node_to_vertex_node;
//    Vec<unsigned>* connect_node_to_vertex_node;
//    template<class TE, class TF> void operator()( const TE &elem, const TF &f, Vec< Vec<unsigned> > &constrained_nodes_patch ) const {
//        unsigned elem_nb_nodes = TE::nb_nodes; // nb de noeuds de l'element elem
//        for (unsigned n=0;n<elem_nb_nodes;++n) {
//            if ( (*correspondance_node_to_vertex_node)[ elem.node( n )->number ] ) {
//                for (unsigned i=0;i<elem_nb_nodes;++i) {
//                    if ( f.constrained_nodes()[ elem.node( i )->number ] ) {
//                        constrained_nodes_patch[ (*connect_node_to_vertex_node)[ elem.node( n )->number ] ].push_back( elem.node( i )->number );
//                    }
//                }
//            }
//        }
//    }
//};

//     cout << f.matrices(Number<0>()) << endl << endl;
//     cout << f.sollicitation << endl << endl;
//     cout << f.vectors[0] << endl << endl;
//     cout << f.matrices(Number<0>()).size() << endl << endl;
//     cout << f.sollicitation.size() << endl << endl;
//     cout << f.vectors[0].size() << endl << endl;
//     cout << sqrt( dot( f.matrices(Number<0>()) * f.vectors[0] - f.sollicitation, f.matrices(Number<0>()) * f.vectors[0] - f.sollicitation ) ) << endl << endl;
//     m.f_vol[1] = -0.01;
//     cout << m.node_list[0].dep << endl; // valeur du deplacement du noeud 0 du maillage

#endif // Test_cpp_h
