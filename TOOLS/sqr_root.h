#ifndef sqr_root_h
#define sqr_root_h

#include "eigen_problem_using_lapack.h"

/// partie positive d'une matrice diagonale
template<class TM> TM pos_part_mat_diag(const TM &M)
{
    unsigned dim = M.nb_cols();

    TM res;
    res.resize(dim);
    res.set(0.);
    for(unsigned k=0; k<dim; ++k)
    {
        res(k, k) = LMT::max(M(k, k), 0.);
    }
    return res;
}

/// partie négative d'une matrice diagonale
template<class TM> TM neg_part_mat_diag(const TM &M)
{
    unsigned dim = M.nb_cols();

    TM res;
    res.resize(dim);
    res.set(0.);
    for(unsigned k=0; k<dim; ++k)
    {
        res(k, k) = LMT::min(M(k, k), 0.);
    }
    return res;
}

/// partie positive d'une matrice symetrique definie positive using_lapack
template<class TM> TM pos_part_mat_sym_using_lapack(const TM &M)
{
    unsigned dim = M.nb_cols();

    LMT::Vec<double> V;
    V.resize(dim);
    LMT::Mat<double, LMT::Gen<> > P;
    P.resize(dim);
    eigen_problem_using_lapack(M, V, P);

    LMT::Mat<double, LMT::Diag<> > D = LMT::trans(P) * M * P;
    return P * pos_part_mat_diag(D) * LMT::trans(P);
}

/// partie négative d'une matrice symetrique definie positive using_lapack
template<class TM> TM neg_part_mat_sym_using_lapack(const TM &M)
{
    unsigned dim = M.nb_cols();

    LMT::Vec<double> V;
    V.resize(dim);
    LMT::Mat<double, LMT::Gen<> > P;
    P.resize(dim);
    eigen_problem_using_lapack(M, V, P);

    LMT::Mat<double, LMT::Diag<> > D = LMT::trans(P) * M * P;
    return P * neg_part_mat_diag(D) * LMT::trans(P);
}

/// partie positive d'un vecteur colonne using_lapack
template<class TV> TV pos_part_vec_col_using_lapack(const TV &V)
{
    return mat_sym_to_vec_col(pos_part_mat_sym_using_lapack(vec_col_to_mat_sym(V)));
}

/// partie negative d'un vecteur colonne using_lapack
template<class TV> TV neg_part_vec_col_using_lapack(const TV &V)
{
    return mat_sym_to_vec_col(neg_part_mat_sym_using_lapack(vec_col_to_mat_sym(V)));
}

/// partie positive generalisee using_lapack
template<class T, int n> LMT::Vec<T, n> gen_pos_part_using_lapack(const LMT::Vec<T, n> &sigma, const LMT::Mat<T, LMT::Gen<n> > &H)
{
    LMT::Vec<T, n> tmp = H * sigma;
    return inv(H) * pos_part_vec_col_using_lapack(tmp);
}

/// partie negative generalisee using_lapack
template<class T, int n> LMT::Vec<T, n> gen_neg_part_using_lapack(const LMT::Vec<T, n> &sigma, const LMT::Mat<T, LMT::Gen<n> > &H)
{
    LMT::Vec<T, n> tmp = H * sigma;
    return inv(H) * neg_part_vec_col_using_lapack(tmp);
}

/// inverse d'une matrice diagonale definie positive
template<class TM> TM inv_mat_diag(const TM &M)
{
    unsigned dim = M.nb_cols();

    TM res;
    res.resize(dim);
    res.set(0.);
    for(unsigned k=0; k<dim; ++k)
    {
        res(k, k) = 1. / M(k, k);
    }
    return res;
}

/// racine carrée d'une matrice diagonale definie positive
template<class TM> TM sqr_root_mat_diag(const TM &M)
{
    unsigned dim = M.nb_cols();

    TM res;
    res.resize(dim);
    res.set(0.);
    for(unsigned k=0; k<dim; ++k)
    {
        res(k, k) = sqrt(M(k, k));
    }
    return res;
}

/// racine carrée d'une matrice symetrique definie positive using_lapack
template<class TM> TM sqr_root_mat_sym_using_lapack(const TM &M)
{
    unsigned dim = M.nb_cols();

    LMT::Vec<double> V;
    V.resize(dim);
    LMT::Mat<double> P;
    P.resize(dim);
    eigen_problem_using_lapack(M, V, P);

    LMT::Mat<double, LMT::Diag<> > D = LMT::trans(P) * M * P;
    return P * sqr_root_mat_diag(D) * LMT::trans(P);
}

#endif
