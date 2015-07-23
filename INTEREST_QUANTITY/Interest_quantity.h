//
// C++ Interface: Interest_quantity
//
// Description: calcul de la quantite d'interet
//
//
// Author: Pled Florent These 2011 <pled@lmt.ens-cachan.fr>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef Interest_quantity_h
#define Interest_quantity_h

#include "../GEOMETRY/Geometry.h"
#include "../CRITERIUM_ENHANCEMENT/Criterium_enhancement.h"

using namespace LMT;
using namespace std;

/// Definition de l'extracteur associe a une quantite d'interet = facteur d'intensite de contrainte (SIF)
///------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class S, class TTV, class TT>
void def_extractor_SIF( TE &elem, const TM &m, const TF &f, const S &interest_quantity, const S &direction_extractor, const TTV &pos_crack_tip, const TT &angle_crack, const TT &radius_Ri, const TT &radius_Re ) {}

template<class T, class Pvec>
struct Define_Extractor_SIF {
    const string* interest_quantity;
    const string* direction_extractor;
    const Pvec* pos_crack_tip;
    const T* angle_crack;
    const T* radius_Ri;
    const T* radius_Re;
    template<class TE, class TM, class TF> void operator()( TE &elem, const TM &m, const TF &f ) const {
        def_extractor_SIF( elem, m, f, *interest_quantity, *direction_extractor, *pos_crack_tip, *angle_crack, *radius_Ri, *radius_Re );
    }
};

/// Calcul des fonctions d'enrichissement (handbook)
///-------------------------------------------------
template<class TE, class TM, class TF>
void calc_dep_handbook_in_infinite_domain( TE &elem_adjoint, const TM &m_adjoint, const TF &f_adjoint ) {}

struct Calcul_Dep_Handbook_In_Infinite_Domain {
    template<class TE, class TM, class TF> void operator()( TE &elem_adjoint, const TM &m_adjoint, const TF &f_adjoint ) const {
        calc_dep_handbook_in_infinite_domain( elem_adjoint, m_adjoint, f_adjoint );
    }
};

/// Construction du chargement du pb adjoint a partir de l'extracteur
///------------------------------------------------------------------
struct Construct_Extractor_Mean_Epsilon {
template<class TE, class TM> void operator()( TE &elem_adjoint, const TM &m, const Vec<unsigned> &elem_list_interest_quantity ) const {
    typedef typename TE::T T;
    for (unsigned n=0;n<elem_list_interest_quantity.size();++n) {
        Vec<Vec<T,TE::dim>, TE::nb_nodes > pos_nodes;
        for (unsigned i=0;i<(m.elem_list[ elem_list_interest_quantity[ n ] ]->nb_nodes_virtual());++i)
            pos_nodes[i] = m.elem_list[ elem_list_interest_quantity[ n ] ]->node_virtual(i)->pos;
        if ( is_inside_linear( typename TE::NE(), pos_nodes, center( elem_adjoint ) ) ) {
            Vec<T,unsigned(TE::dim*(TE::dim+1)/2) > extractor = m.elem_list[ elem_list_interest_quantity[ n ] ]->get_field( "pre_sigma", StructForType<Vec<T,unsigned(TE::dim*(TE::dim+1)/2)> >() );
            elem_adjoint.set_field( "pre_sigma", extractor );
//                cout << "pre sigma = " << endl;
//                cout << extractor << endl;
//                Mat<T,Sym<TE::dim > > pre_sig = elem_adjoint.get_field( "pre_sigma", StructForType<Mat<T,Sym<TE::dim > > >() );
//                cout << pre_sig << endl;
        }
    }
}
};

struct Construct_Extractor_Mean_Sigma {
    template<class TE, class TM> void operator()( TE &elem_adjoint, const TM &m, const Vec<unsigned> &elem_list_interest_quantity ) const {
        typedef typename TE::T T;
        for (unsigned n=0;n<elem_list_interest_quantity.size();++n) {
            Vec<Vec<T,TE::dim>, TE::nb_nodes > pos_nodes;
            for (unsigned i=0;i<(m.elem_list[ elem_list_interest_quantity[ n ] ]->nb_nodes_virtual());++i)
                pos_nodes[i] = m.elem_list[ elem_list_interest_quantity[ n ] ]->node_virtual(i)->pos;
            if ( is_inside_linear( typename TE::NE(), pos_nodes, center( elem_adjoint ) ) ) {
                Vec<T,unsigned(TE::dim*(TE::dim+1)/2) > extractor = m.elem_list[ elem_list_interest_quantity[ n ] ]->get_field( "pre_epsilon", StructForType<Vec<T,unsigned(TE::dim*(TE::dim+1)/2)> >() );
                elem_adjoint.set_field( "pre_epsilon", extractor );
//                cout << "pre epsilon = " << endl;
//                cout << extractor << endl;
//                Mat<T,Sym<TE::dim > > pre_eps = elem_adjoint.get_field( "pre_epsilon", StructForType<Mat<T,Sym<TE::dim > > >() );
//                cout << pre_eps << endl;
            }
        }
    }
};

struct Construct_Extractor_Pointwise_Dep {
    template<class TN, class TM> void operator()( TN &node_adjoint, const TM &m, const unsigned &node ) const {
        if ( norm_2( node_adjoint.pos - m.node_list[ node ].pos ) / norm_2( m.node_list[ node ].pos ) < 1e-6 ) {
            Vec<typename TN::T,TN::dim> extractor = m.node_list[ node ].pre_f_nodal;
            node_adjoint.pre_f_nodal = extractor;
        }
    }
};

struct Construct_Extractor_SIF {
    template<class TE, class TM> void operator()( TE &elem_adjoint, const TM &m_crown ) const {
        typedef typename TE::T T;
        for (unsigned n=0;n<m_crown.elem_list.size();++n) {
            Vec<Vec<T,TE::dim>, TE::nb_nodes > pos_nodes;
            for (unsigned i=0;i<(m_crown.elem_list[ n ]->nb_nodes_virtual());++i)
                pos_nodes[i] = m_crown.elem_list[ n ]->node_virtual(i)->pos;
            if ( is_inside_linear( typename TE::NE(), pos_nodes, center( elem_adjoint ) ) ) {
                Vec<T,unsigned(TE::dim*(TE::dim+1)/2) > extractor_pre_sigma = m_crown.elem_list[ n ]->get_field( "pre_sigma", StructForType<Vec<T,unsigned(TE::dim*(TE::dim+1)/2)> >() );
                elem_adjoint.set_field( "pre_sigma", extractor_pre_sigma );
                Vec<T,unsigned(TE::dim) > extractor_pre_f_vol = m_crown.elem_list[ n ]->get_field( "pre_f_vol", StructForType<Vec<T,unsigned(TE::dim)> >() );
                elem_adjoint.set_field( "pre_f_vol", extractor_pre_f_vol );
//                cout << "pre sigma = " << endl;
//                cout << extractor_pre_sigma << endl;
//                Mat<T,Sym<TE::dim > > pre_sigma = elem_adjoint.get_field( "pre_sigma", StructForType<Mat<T,Sym<TE::dim > > >() );
//                cout << pre_sigma << endl;
            }
        }
    }
};

/// Calcul d'une quantite d'interet = valeur moyenne d'une composante du champ de contrainte (mean_sigma) ou de deformation (mean_epsilon)
///---------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TTWW, class TT>
void calc_interest_quantity_mean_sigma_epsilon( const TE &elem, const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, TT &I ) {}

struct Calcul_Interest_Quantity_Mean_Sigma_Epsilon {
    const Vec<unsigned>* elem_list_interest_quantity;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, T &I ) const {
        if ( find( *elem_list_interest_quantity, _1 == elem.number ) ) {
            Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
            calc_interest_quantity_mean_sigma_epsilon( elem, m, f, f.vectors, ind, I );
        }
    }
};

/// Calcul d'une quantite d'interet = valeur ponctuelle d'une composante du champ de deplacement (pointwise_dep) definie a partir du numero d'un noeud du maillage
///---------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Calcul_Interest_Quantity_Pointwise_Dep_Node {
    const string* interest_quantity;
    template<class TN, class T> void operator()( TN &node, const string &direction_extractor, const unsigned &node_interest_quantity, T &I ) const {
        if ( node.number == node_interest_quantity ) {
            if ( TN::dim == 1 ) {
                if ( direction_extractor == "x" )
                    I += node.dep[0];
                else {
                    std::cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << *interest_quantity << " en dimension " << TN::dim << " n'est pas implementee..." << std::endl << std::endl;
                    throw "Anguille sous coquille...";
                }
            }
            else if ( TN::dim == 2 ) {
                if ( direction_extractor == "x" )
                    I += node.dep[0];
                else if ( direction_extractor == "y" )
                    I += node.dep[1];
                else {
                    std::cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << *interest_quantity << " en dimension " << TN::dim << " n'est pas implementee..." << std::endl << std::endl;
                    throw "Anguille sous coquille...";
                }
            }
            else if ( TN::dim == 3 ) {
                if ( direction_extractor == "x" )
                    I += node.dep[0];
                else if ( direction_extractor == "y" )
                    I += node.dep[1];
                else if ( direction_extractor == "z" )
                    I += node.dep[2];
                else {
                    std::cerr << "Arret brutal, car la direction " << direction_extractor << " pour la quantite d'interet " << *interest_quantity << " en dimension " << TN::dim << " n'est pas implementee..." << std::endl << std::endl;
                    throw "Anguille sous coquille...";
                }
            }
        }
    }
};

/// Calcul d'une quantite d'interet = valeur ponctuelle d'une composante du champ de deplacement (pointwise_dep) ou de contrainte (pointwise_sigma) ou de deformation (pointwise_epsilon) definie a partir de la position d'un point
///---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TM, class TF, class TTWW, class S, class Pvec, class TT>
void calc_interest_quantity_pointwise_dep_sigma_epsilon( const TE &elem, const TM &m, const TF &f, const TTWW &vectors, const Vec<unsigned> &indices, const S &interest_quantity, const S &direction_extractor, const Pvec &pos, TT &I ) {}

template<class Pvec>
struct Calcul_Interest_Quantity_Pointwise_Dep_Sigma_Epsilon {
    const string* interest_quantity;
    const string* direction_extractor;
    const Pvec* pos_interest_quantity;
    template<class TE, class TM, class TF, class T> void operator()( const TE &elem, const TM &m, const TF &f, T &I ) const {
        Vec<Vec<T,TE::dim>, TE::nb_nodes > pos_nodes;
        for (unsigned i=0;i<(m.elem_list[ elem.number ]->nb_nodes_virtual());++i)
            pos_nodes[i] = m.elem_list[ elem.number ]->node_virtual(i)->pos;
        if ( is_inside_linear( typename TE::NE(), pos_nodes, *pos_interest_quantity ) ) {
            Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f.indices_for_element( elem );
            calc_interest_quantity_pointwise_dep_sigma_epsilon( elem, m, f, f.vectors, ind, *interest_quantity, *direction_extractor, *pos_interest_quantity, I );
        }
    }
};

/// Calcul d'une quantite d'interet = facteur d'intensite de contrainte (SIF)
///--------------------------------------------------------------------------
template<class TE, class TM, class TF, class S, class TTV, class TT, class TTWW>
void calc_interest_quantity_SIF( const TE &elem, const TM &m_crown, const TF &f_crown, const S &direction_extractor, const TTV &pos_crack_tip, const TT &angle_crack, const TT &radius_Ri, const TT &radius_Re, const TTWW &vectors, const Vec<unsigned> &indices, TT &I ) {}

template<class T, class Pvec>
struct Calcul_Interest_Quantity_SIF {
    const string* direction_extractor;
    const Pvec* pos_crack_tip;
    const T* angle_crack;
    const T* radius_Ri;
    const T* radius_Re;
    template<class TE, class TM, class TF> void operator()( const TE &elem, const TM &m_crown, const TF &f_crown, T &I ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = f_crown.indices_for_element( elem );
        calc_interest_quantity_SIF( elem, m_crown, f_crown, *direction_extractor, *pos_crack_tip, *angle_crack, *radius_Ri, *radius_Re, f_crown.vectors, ind, I );
    }
};

/// Calcul de la correction I_hh sur la quantite d'interet (sans introduction de sigma_hat_m)
///------------------------------------------------------------------------------------------
template<class TE, class TE_adjoint, class TM, class TF, class TTWW, class S, class B, class TTVV, class TT>
void calc_elem_correction_interest_quantity_wo_sigma_hat_m( const TE &elem, const TE_adjoint &elem_adjoint, const TM &m, const TM &m_adjoint, const TF &f, const TF &f_adjoint, const TTWW &vectors, const TTWW &adjoint_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &adjoint_indices, const S &method, const B &want_local_enrichment, const TTVV &dep_hat, TT &I_hh  ) {}

template<class TM, class TF, class T>
struct Calcul_Elem_Correction_Interest_Quantity_wo_sigma_hat_m {
    const TM* m;
    const TM* m_adjoint;
    const TF* f;
    const TF* f_adjoint;
    const string* method;
    const bool* want_local_enrichment;
    const Vec< Vec<T> >* dep_hat;
    T* I_hh;
    template<class TE, class TE_adjoint> void operator()( const TE &elem, const TE_adjoint &elem_adjoint ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = (*f).indices_for_element( elem );
        Vec<unsigned,TE_adjoint::nb_nodes+1+TF::nb_global_unknowns> ind_adjoint = (*f_adjoint).indices_for_element( elem_adjoint );
        calc_elem_correction_interest_quantity_wo_sigma_hat_m( elem, elem_adjoint, *m, *m_adjoint, *f, *f_adjoint, (*f).vectors, (*f_adjoint).vectors, ind, ind_adjoint, *method, *want_local_enrichment, *dep_hat, *I_hh );
    }
};

template<class TM, class TF, class T>
struct Calcul_Correction_Interest_Quantity_wo_sigma_hat_m {
    const TM* m;
    const TM* m_adjoint;
    const TF* f;
    const TF* f_adjoint;
    const string* interest_quantity;
    const string* method;
    const bool* want_local_enrichment;
    const Vec<unsigned>* correspondance_elem_m_adjoint_to_elem_m;
    const Vec< Vec<T> >* dep_hat;
    T* I_hh;
    template<class TE_adjoint> void operator()( const TE_adjoint &elem_adjoint ) const {
        unsigned num_elem = (*correspondance_elem_m_adjoint_to_elem_m)[ elem_adjoint.number ];
        
        Calcul_Elem_Correction_Interest_Quantity_wo_sigma_hat_m<TM, TF, T> calcul_elem_correction_interest_quantity_wo_sigma_hat_m;
        calcul_elem_correction_interest_quantity_wo_sigma_hat_m.m = m;
        calcul_elem_correction_interest_quantity_wo_sigma_hat_m.m_adjoint = m_adjoint;
        calcul_elem_correction_interest_quantity_wo_sigma_hat_m.f = f;
        calcul_elem_correction_interest_quantity_wo_sigma_hat_m.f_adjoint = f_adjoint;
        calcul_elem_correction_interest_quantity_wo_sigma_hat_m.method = method;
        calcul_elem_correction_interest_quantity_wo_sigma_hat_m.want_local_enrichment = want_local_enrichment;
        calcul_elem_correction_interest_quantity_wo_sigma_hat_m.dep_hat = dep_hat;
        calcul_elem_correction_interest_quantity_wo_sigma_hat_m.I_hh = I_hh;
        
        apply_on_number( m->elem_list, num_elem, calcul_elem_correction_interest_quantity_wo_sigma_hat_m, elem_adjoint );
    }
};

/// Calcul de la correction I_hh sur la quantite d'interet (avec introduction de sigma_hat_m)
///------------------------------------------------------------------------------------------
template<class TE, class TE_adjoint, class TM, class TF, class TTWW, class S, class B, class TTVV, class TT>
void calc_elem_correction_interest_quantity_w_sigma_hat_m( const TE &elem, const TE_adjoint &elem_adjoint, const TM &m, const TM &m_adjoint, const TF &f, const TF &f_adjoint, const TTWW &vectors, const TTWW &adjoint_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &adjoint_indices, const S &method, const S &method_adjoint, const B &want_local_enrichment, const TTVV &dep_hat, const TTVV &dep_adjoint_hat, TT &I_hh ) {}

template<class TM, class TF,class T>
struct Calcul_Elem_Correction_Interest_Quantity_w_sigma_hat_m {
    const TM* m;
    const TM* m_adjoint;
    const TF* f;
    const TF* f_adjoint;
    const string* interest_quantity;
    const string* method;
    const string* method_adjoint;
    const bool* want_local_enrichment;
    const Vec< Vec<T> >* dep_hat;
    const Vec< Vec<T> >* dep_adjoint_hat;
    T* I_hh;
    template<class TE, class TE_adjoint> void operator()( const TE &elem, const TE_adjoint &elem_adjoint ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = (*f).indices_for_element( elem );
        Vec<unsigned,TE_adjoint::nb_nodes+1+TF::nb_global_unknowns> ind_adjoint = (*f_adjoint).indices_for_element( elem_adjoint );
        calc_elem_correction_interest_quantity_w_sigma_hat_m( elem, elem_adjoint, *m, *m_adjoint, *f, *f_adjoint, (*f).vectors, (*f_adjoint).vectors, ind, ind_adjoint, *method, *method_adjoint, *want_local_enrichment, *dep_hat, *dep_adjoint_hat, *I_hh );
    }
};

template<class TM, class TF, class T>
struct Calcul_Correction_Interest_Quantity_w_sigma_hat_m {
    const TM* m;
    const TM* m_adjoint;
    const TF* f;
    const TF* f_adjoint;
    const string* interest_quantity;
    const string* method;
    const string* method_adjoint;
    const bool* want_local_enrichment;
    const Vec<unsigned>* correspondance_elem_m_adjoint_to_elem_m;
    const Vec< Vec<T> >* dep_hat;
    const Vec< Vec<T> >* dep_adjoint_hat;
    T* I_hh;
    template<class TE_adjoint> void operator()( const TE_adjoint &elem_adjoint ) const {
        unsigned num_elem = (*correspondance_elem_m_adjoint_to_elem_m)[ elem_adjoint.number ];
        
        Calcul_Elem_Correction_Interest_Quantity_w_sigma_hat_m<TM, TF, T> calcul_elem_correction_interest_quantity_w_sigma_hat_m;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.m = m;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.m_adjoint = m_adjoint;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.f = f;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.f_adjoint = f_adjoint;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.method = method;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.method_adjoint = method_adjoint;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.want_local_enrichment = want_local_enrichment;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.dep_hat = dep_hat;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.dep_adjoint_hat = dep_adjoint_hat;
        calcul_elem_correction_interest_quantity_w_sigma_hat_m.I_hh = I_hh;
        
        apply_on_number( m->elem_list, num_elem, calcul_elem_correction_interest_quantity_w_sigma_hat_m, elem_adjoint );
    }
};

/// Construction du centre et de la taille du domaine Omega_lambda selon la forme et la quantite d'interet
///-------------------------------------------------------------------------------------------------------
template<class TM, class T, class Pvec>
void construct_center_length_domain( TM &m, const unsigned &deg_p, const string &shape, const string &interest_quantity, const string &pointwise_interest_quantity, const Vec<unsigned> &elem_list, const unsigned &node, const Pvec &pos, const Pvec &pos_crack_tip, const T &radius_Re, Pvec &domain_center, Vec<T> &domain_length ) {
    
    static const unsigned dim = TM::dim;
    
    m.update_node_neighbours();
    
    if (deg_p == 1) {
        if ( shape.find("circle") != string::npos or shape.find("sphere") != string::npos )
            domain_length.resize( 1 );
        else if ( shape == "square" or shape == "rectangle" or shape == "cube" or shape == "cuboid" )
            domain_length.resize( dim );
        else
            cerr << "forme " << shape << " non implementee pour le calcul de la constante dans l'amelioration..." << endl << endl;
        domain_length.set( 0. );
        T domain_size = 0.;
        if ( interest_quantity.find("mean") != string::npos ) {
            /// Construction du vecteur de vecteurs circum_center et du vecteur circum_radius
            /// circum_center[ n ] : position du centre du cercle/sphere circonscrit au n ieme element de la liste elem_list (elem_list[ n ])
            /// circum_radius[ n ] : rayon du cercle/sphere circonscrit au n ieme element de la liste elem_list (elem_list[ n ])
            ///--------------------------------------------------------------------------------------------------------------------------------
            Vec< Vec<T> > circum_center;
            circum_center.resize( elem_list.size() );

            for (unsigned n=0;n<elem_list.size();++n) {
                circum_center[ n ].resize( dim );
                circum_center[ n ].set( 0. );
            }

            Vec<T> circum_radius;
            circum_radius.resize( elem_list.size() );
            circum_radius.set( 0. );

            Construct_Circum_Center_Radius_Elem_List construct_circum_center_radius_elem_list;
            construct_circum_center_radius_elem_list.elem_list = &elem_list;
            
            apply( m.elem_list, construct_circum_center_radius_elem_list, circum_center, circum_radius );
            
//             cout << "Construction du vecteur de vecteurs circum_center et du vecteur circum_radius" << endl << endl;
//             for (unsigned n=0;n<elem_list.size();++n) {
//                 if ( dim == 2 ) {
//                     cout << "position du centre du cercle circonscrit a l'element " << elem_list[ n ] << " : " << circum_center[ n ] << endl;
//                     cout << "rayon du cercle circonscrit a l'element " << elem_list[ n ] << " : " << circum_radius[ n ] << endl;
//                 }
//                 else if ( dim == 3 ) {
//                     cout << "position du centre de la sphere circonscrite a l'element " << elem_list[ n ] << " : " << circum_center[ n ] << endl;
//                     cout << "rayon de la sphere circonscrite a l'element " << elem_list[ n ] << " : " << circum_radius[ n ] << endl;
//                 }
//                 cout << endl << endl;
//             }
            
            for (unsigned d=0;d<dim;++d) {
                Vec<T> pos_circum_centers;
                pos_circum_centers.reserve( elem_list.size() );
                for (unsigned n=0;n<elem_list.size();++n)
                    pos_circum_centers.push_back( circum_center[ n ][ d ] );
                domain_center[ d ] = mean( pos_circum_centers );
            }
            for (unsigned n=0;n<elem_list.size();++n)
                domain_size = max( domain_size, length( circum_center[ n ] - domain_center ) + circum_radius[ n ] );
        }
        else if ( interest_quantity.find("pointwise") != string::npos ) {
            if ( pointwise_interest_quantity == "node" ) {
                domain_center = m.node_list[ node ].pos;
                for (unsigned i=0;i<m.get_node_neighbours( node ).size();++i)
                    domain_size = max( domain_size, length( m.get_node_neighbours( node )[ i ]->pos - domain_center ) );
            }
            else if ( pointwise_interest_quantity == "pos" ) {
                domain_center = pos;
                Vec<unsigned> elem_list_pos;
                apply( m.elem_list, Construct_Elem_List_Pos(), pos, elem_list_pos );
                for (unsigned n=0;n<elem_list_pos.size();++n) {
                    for (unsigned i=0;i<(m.elem_list[ elem_list_pos[ n ] ]->nb_nodes_virtual());++i)
                        domain_size = max( domain_size, length( m.elem_list[ elem_list_pos[ n ] ]->node_virtual(i)->pos - domain_center ) );
                }
            }
        }
        else if ( interest_quantity == "SIF" or interest_quantity == "stress_intensity_factor" ) {
            domain_center = pos_crack_tip;
            domain_size = radius_Re;
        }
        if ( shape.find("circle") != string::npos or shape.find("sphere") != string::npos )
            domain_length[ 0 ] = domain_size;
        else if ( shape == "square" or shape == "rectangle" or shape == "cube" or shape == "cuboid" ) {
            for (unsigned d=0;d<dim;++d)
                domain_length[ d ] = domain_size;
        }
    }
}

/// Calcul des champs de contrainte sigma_lambda et sigma_hat_lambda, Transfert du champ de deplacement dep_hat a dep_hat_lambda, Calcul d'un estimateur d'erreur globale sur le domaine Omega_lambda
///--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TE_lambda, class TM, class TF, class TTWW, class TTVV, class S, class TTV, class TT>
void calc_elem_error_estimate_lambda( const TE &elem, TE_lambda &elem_lambda, const TM &m, const TM &m_lambda, const TF &f, const TF &f_lambda, const TTWW &vectors, const TTWW &lambda_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &lambda_indices, const TTVV &dep_hat, const TTVV &dep_hat_lambda, const S &method, TTV &theta_elem, TT &theta ) {}

template<class TM, class TF, class T>
struct Calcul_Elem_Error_Estimate_Lambda {
    const TM* m;
    const TM* m_lambda;
    const TF* f;
    const TF* f_lambda;
    const string* method;
    const Vec< Vec<T> >* dep_hat;
    Vec< Vec<T> >* dep_hat_lambda;
    Vec<T>* theta_elem;
    T* theta;
    template<class TE, class TE_lambda> void operator()( const TE &elem, TE_lambda &elem_lambda ) const {
        Vec<unsigned,TE_lambda::nb_nodes+1+TF::nb_global_unknowns> ind_lambda = (*f_lambda).indices_for_element( elem_lambda );
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = (*f).indices_for_element( elem );
        calc_elem_error_estimate_lambda( elem, elem_lambda, *m, *m_lambda, *f, *f_lambda, (*f).vectors, (*f_lambda).vectors, ind, ind_lambda, *dep_hat, *dep_hat_lambda, *method, *theta_elem, *theta );
    }
};

template<class TM, class TF, class T>
struct Calcul_Error_Estimate_Lambda {
    const TM* m;
    const TM* m_lambda;
    const TF* f;
    const TF* f_lambda;
    const string* method;
    const Vec<unsigned>* correspondance_elem_m_lambda_to_elem_m;
    const Vec< Vec<T> >* dep_hat;
    Vec< Vec<T> >* dep_hat_lambda;
    Vec<T>* theta_elem;
    T* theta;
    template<class TE_lambda> void operator()( TE_lambda &elem_lambda ) const {
        unsigned num_elem = (*correspondance_elem_m_lambda_to_elem_m)[ elem_lambda.number ];
        
        Calcul_Elem_Error_Estimate_Lambda<TM, TF, T> calcul_elem_error_estimate_lambda;
        calcul_elem_error_estimate_lambda.m = m;
        calcul_elem_error_estimate_lambda.m_lambda = m_lambda;
        calcul_elem_error_estimate_lambda.f = f;
        calcul_elem_error_estimate_lambda.f_lambda = f_lambda;
        calcul_elem_error_estimate_lambda.method = method;
        calcul_elem_error_estimate_lambda.dep_hat = dep_hat;
        calcul_elem_error_estimate_lambda.dep_hat_lambda = dep_hat_lambda;
        calcul_elem_error_estimate_lambda.theta_elem = theta_elem;
        calcul_elem_error_estimate_lambda.theta = theta;
        
        apply_on_number( m->elem_list, num_elem, calcul_elem_error_estimate_lambda, elem_lambda );
    }
};

/// Calcul d'un estimateur pondere d'erreur globale sur le domaine Omega_lambda
///----------------------------------------------------------------------------
template<class TE, class TE_lambda, class TM, class TF, class TTWW, class TTVV, class S, class Pvec, class TTV, class TT>
void calc_elem_weighted_error_estimate_lambda( const TE &elem, const TE_lambda &elem_lambda, const TM &m, const TM &m_lambda, const TF &f, const TF &f_lambda, const TTWW &vectors, const TTWW &lambda_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &lambda_indices, const TTVV &dep_hat, S &method, const TT &h, const Pvec &domain_center, const TTV &domain_length, const TT &k_min, TT &weighted_theta ) {}

template<class TM, class TF, class T, class Pvec>
struct Calcul_Elem_Weighted_Error_Estimate_Lambda {
    const TM* m;
    const TM* m_lambda;
    const TF* f;
    const TF* f_lambda;
    const string* method;
    const Vec< Vec<T> >* dep_hat;
    const T* h;
    const Pvec* domain_center;
    const Vec<T>* domain_length;
    const T* k_min;
    T* weighted_theta;
    template<class TE, class TE_lambda> void operator()( const TE &elem, const TE_lambda &elem_lambda ) const {
        Vec<unsigned,TE_lambda::nb_nodes+1+TF::nb_global_unknowns> ind_lambda = (*f_lambda).indices_for_element( elem_lambda );
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = (*f).indices_for_element( elem );
        calc_elem_weighted_error_estimate_lambda( elem, elem_lambda, *m, *m_lambda, *f, *f_lambda, (*f).vectors, (*f_lambda).vectors, ind, ind_lambda, *dep_hat, *method, *h, *domain_center, *domain_length, *k_min, *weighted_theta );
    }
};

template<class TM, class TF, class T, class Pvec>
struct Calcul_Weighted_Error_Estimate_Lambda {
    const TM* m;
    const TM* m_lambda;
    const TF* f;
    const TF* f_lambda;
    const string* method;
    const Vec<unsigned>* correspondance_elem_m_lambda_to_elem_m;
    const Vec< Vec<T> >* dep_hat;
    const T* h;
    const Pvec* domain_center;
    const Vec<T>* domain_length;
    const T* k_min;
    T* weighted_theta;
    template<class TE_lambda> void operator()( const TE_lambda &elem_lambda ) const {
        unsigned num_elem = (*correspondance_elem_m_lambda_to_elem_m)[ elem_lambda.number ];
        
        Calcul_Elem_Weighted_Error_Estimate_Lambda<TM, TF, T, Pvec> calcul_elem_weighted_error_estimate_lambda;
        calcul_elem_weighted_error_estimate_lambda.m = m;
        calcul_elem_weighted_error_estimate_lambda.m_lambda = m_lambda;
        calcul_elem_weighted_error_estimate_lambda.f = f;
        calcul_elem_weighted_error_estimate_lambda.f_lambda = f_lambda;
        calcul_elem_weighted_error_estimate_lambda.method = method;
        calcul_elem_weighted_error_estimate_lambda.dep_hat = dep_hat;
        calcul_elem_weighted_error_estimate_lambda.h = h;
        calcul_elem_weighted_error_estimate_lambda.domain_center = domain_center;
        calcul_elem_weighted_error_estimate_lambda.domain_length = domain_length;
        calcul_elem_weighted_error_estimate_lambda.k_min = k_min;
        calcul_elem_weighted_error_estimate_lambda.weighted_theta = weighted_theta;
        
        apply_on_number( m->elem_list, num_elem, calcul_elem_weighted_error_estimate_lambda, elem_lambda );
    }
};

/// Calcul d'un estimateur d'erreur globale sur le bord du domaine Omega_lambda
///----------------------------------------------------------------------------
template<class TE, class TE_lambda, class TM, class TF, class TTWW, class TTVV, class S, class Pvec, class TTV, class TT>
void calc_skin_elem_error_estimate_lambda_boundary( const TE &elem, const TE_lambda &elem_lambda, const TM &m, const TM &m_lambda, const TF &f, const TF &f_lambda, const TTWW &vectors, const TTWW &lambda_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &lambda_indices, const TTVV &dep_hat, S &method, const Pvec &domain_center, TTV &theta_boundary_face, TT &theta_boundary ) {}

template<class TM, class TF, class T, class Pvec>
struct Calcul_Skin_Elem_Error_Estimate_Lambda_Boundary {
    const TM* m;
    const TM* m_lambda;
    const TF* f;
    const TF* f_lambda;
    const string* method;
    const Vec< Vec<T> >* dep_hat;
    const Pvec* domain_center;
    Vec<T>* theta_boundary_face;
    T* theta_boundary;
    template<class TE, class TE_lambda> void operator()( const TE &elem, const TE_lambda &elem_lambda ) const {
        Vec<unsigned,TE_lambda::nb_nodes+1+TF::nb_global_unknowns> ind_lambda = (*f_lambda).indices_for_element( elem_lambda );
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = (*f).indices_for_element( elem );
        calc_skin_elem_error_estimate_lambda_boundary( elem, elem_lambda, *m, *m_lambda, *f, *f_lambda, (*f).vectors, (*f_lambda).vectors, ind, ind_lambda, *dep_hat, *method, *domain_center, *theta_boundary_face, *theta_boundary );
    }
};

template<class TM, class TF, class T, class Pvec>
struct Calcul_Error_Estimate_Lambda_Boundary {
    const TM* m;
    const TM* m_lambda;
    const TF* f;
    const TF* f_lambda;
    const string* method;
    const Vec<unsigned>* correspondance_elem_m_lambda_to_elem_m;
    const Vec< Vec<T> >* dep_hat;
    const Pvec* domain_center;
    Vec<T>* theta_boundary_face;
    T* theta_boundary;
    template<class TE_lambda> void operator()( const TE_lambda &elem_lambda ) const {
        unsigned num_elem = (*correspondance_elem_m_lambda_to_elem_m)[ elem_lambda.number ];
        
        Calcul_Skin_Elem_Error_Estimate_Lambda_Boundary<TM, TF, T, Pvec> calcul_skin_elem_error_estimate_lambda_boundary;
        calcul_skin_elem_error_estimate_lambda_boundary.m = m;
        calcul_skin_elem_error_estimate_lambda_boundary.m_lambda = m_lambda;
        calcul_skin_elem_error_estimate_lambda_boundary.f = f;
        calcul_skin_elem_error_estimate_lambda_boundary.f_lambda = f_lambda;
        calcul_skin_elem_error_estimate_lambda_boundary.method = method;
        calcul_skin_elem_error_estimate_lambda_boundary.dep_hat = dep_hat;
        calcul_skin_elem_error_estimate_lambda_boundary.domain_center = domain_center;
        calcul_skin_elem_error_estimate_lambda_boundary.theta_boundary_face = theta_boundary_face;
        calcul_skin_elem_error_estimate_lambda_boundary.theta_boundary = theta_boundary;
        
        apply_on_number( m->elem_list, num_elem, calcul_skin_elem_error_estimate_lambda_boundary, elem_lambda );
    }
};

/// Calcul de la correction I_hhh sur la quantite d'interet locale I
///-----------------------------------------------------------------
template<class TE, class TE_adjoint, class TM, class TF, class TTWW, class S, class TTVV, class TT>
void calc_elem_correction_interest_quantity_lambda( const TE &elem, const TE_adjoint &elem_adjoint, const TM &m, const TM &m_adjoint, const TF &f, const TF &f_adjoint, const TTWW &vectors, const TTWW &adjoint_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &adjoint_indices, const S &method, const S &method_adjoint, const TTVV &dep_hat, const TTVV &dep_adjoint_hat, TT &I_hhh ) {}

template<class TM, class TF, class T>
struct Calcul_Correction_Interest_Quantity_Lambda {
    const TM* m;
    const TM* m_adjoint;
    const TF* f;
    const TF* f_adjoint;
    const string* interest_quantity;
    const string* method;
    const string* method_adjoint;
    const Vec< Vec<T> >* dep_hat;
    const Vec< Vec<T> >* dep_adjoint_hat;
    T* I_hhh;
    template<class TE, class TE_adjoint> void operator()( const TE &elem, unsigned i, const TE_adjoint &elem_adjoint, unsigned j ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = (*f).indices_for_element( elem );
        Vec<unsigned,TE_adjoint::nb_nodes+1+TF::nb_global_unknowns> ind_adjoint = (*f_adjoint).indices_for_element( elem_adjoint );
        calc_elem_correction_interest_quantity_lambda( elem, elem_adjoint, *m, *m_adjoint, *f, *f_adjoint, (*f).vectors, (*f_adjoint).vectors, ind, ind_adjoint, *method, *method_adjoint, *dep_hat, *dep_adjoint_hat, *I_hhh );
    }
};

/// Calcul de la projection sur les elements du maillage adjoint de la contribution elementaire a l'estimateur d'erreur globale du maillage direct
///-----------------------------------------------------------------------------------------------------------------------------------------------
template<class TE, class TE_adjoint, class TM, class TF, class TTWW, class TTVV, class S, class TTV>
void calc_elem_error_estimate_proj_on_adjoint( const TE &elem, const TE_adjoint &elem_adjoint, const TM &m, const TM &m_adjoint, const TF &f, const TF &f_adjoint, const TTWW &vectors, const TTWW &adjoint_vectors, const Vec<unsigned> &indices, const Vec<unsigned> &adjoint_indices, const S &method, const TTVV &dep_hat, TTV &theta_elem_proj_on_adjoint ) {}

template<class TM, class TF, class T>
struct Calcul_Elem_Error_Estimate_Proj_On_Adjoint {
    const TM* m;
    const TM* m_adjoint;
    const TF* f;
    const TF* f_adjoint;
    const string* method;
    const Vec< Vec<T> >* dep_hat;
    Vec<T>* theta_elem_proj_on_adjoint;
    template<class TE, class TE_adjoint> void operator()( const TE &elem, const TE_adjoint &elem_adjoint ) const {
        Vec<unsigned,TE::nb_nodes+1+TF::nb_global_unknowns> ind = (*f).indices_for_element( elem );
        Vec<unsigned,TE_adjoint::nb_nodes+1+TF::nb_global_unknowns> ind_adjoint = (*f_adjoint).indices_for_element( elem_adjoint );
        calc_elem_error_estimate_proj_on_adjoint( elem, elem_adjoint, *m, *m_adjoint, *f, *f_adjoint, (*f).vectors, (*f_adjoint).vectors, ind, ind_adjoint, *method, *dep_hat, *theta_elem_proj_on_adjoint );
    }
};

template<class TM, class TF, class T>
struct Calcul_Error_Estimate_Proj_On_Adjoint {
    const TM* m;
    const TM* m_adjoint;
    const TF* f;
    const TF* f_adjoint;
    const string* method;
    const Vec<unsigned>* correspondance_elem_m_adjoint_to_elem_m;
    const Vec< Vec<T> >* dep_hat;
    Vec<T>* theta_elem_proj_on_adjoint;
    template<class TE_adjoint> void operator()( const TE_adjoint &elem_adjoint ) const {
        unsigned num_elem = (*correspondance_elem_m_adjoint_to_elem_m)[ elem_adjoint.number ];
        
        Calcul_Elem_Error_Estimate_Proj_On_Adjoint<TM, TF, T> calcul_elem_error_estimate_proj_on_adjoint;
        calcul_elem_error_estimate_proj_on_adjoint.m = m;
        calcul_elem_error_estimate_proj_on_adjoint.m_adjoint = m_adjoint;
        calcul_elem_error_estimate_proj_on_adjoint.f = f;
        calcul_elem_error_estimate_proj_on_adjoint.f_adjoint = f_adjoint;
        calcul_elem_error_estimate_proj_on_adjoint.method = method;
        calcul_elem_error_estimate_proj_on_adjoint.dep_hat = dep_hat;
        calcul_elem_error_estimate_proj_on_adjoint.theta_elem_proj_on_adjoint = theta_elem_proj_on_adjoint;
        
        apply_on_number( m->elem_list, num_elem, calcul_elem_error_estimate_proj_on_adjoint, elem_adjoint );
    }
};

#endif // Interest_quantity_h
