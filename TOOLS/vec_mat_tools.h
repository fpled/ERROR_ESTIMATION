#ifndef vec_mat_tools_h
#define vec_mat_tools_h

#include <containers/mat.h>

/// partie symetrique d'une matrice
template<class TM> TM sym_part_mat(const TM &M)
{
    return (M + LMT::trans(M)) / 2;
}

/// prod_contr
template<class TM> typename TM::T prod_contr (const TM &M1, const TM &M2)
{
    unsigned dim = M1.nb_rows();

    typename TM::T res = 0.;
    for (unsigned i=0; i<dim; ++i)
    {
        for (unsigned j=0; j<dim; ++j)
        {
            res += M1(i,j) * M2(i,j);
        }
    }
    return res;
}

/// self_prod_contr
template<class TM> typename TM::T self_prod_contr (const TM &M)
{
    return prod_contr(M, M);
}

/// trace_vec_col
template<class TV> typename LMT::TypeReduction<LMT::Multiplies,TV>::T trace_vec_col (const TV &V)
{
    unsigned dim_vec = V.size();
    unsigned dim = std::ceil(double(dim_vec)/2);

    typename LMT::TypeReduction<LMT::Multiplies,TV>::T res = 0.;
    for (unsigned i=0; i<dim; ++i)
    {
        res += V[i];
    }
    return res;
}

/// trace_sym_col
template<class TV> typename LMT::TypeReduction<LMT::Multiplies,TV>::T trace_sym_col (const TV &V1, const TV &V2)
{
    unsigned dim_vec = V1.size();
    unsigned dim = std::ceil(double(dim_vec)/2);

    typename LMT::TypeReduction<LMT::Multiplies,TV>::T res = 0.;
    for (unsigned i=0; i<dim; ++i)
    {
        res += V1[i] * V2[i];
    }
    for (unsigned i=dim; i<dim_vec; ++i)
    {
        res += 2 * V1[i] * V2[i];
    }
    return res;
}

/// self_trace_sym_col
template<class TV> typename LMT::TypeReduction<LMT::Multiplies,TV>::T self_trace_sym_col (const TV &V)
{
    return trace_sym_col(V, V);
}

/// norm_inf_mat
template<class TM> typename TM::T norm_inf_mat (const TM &M)
{
    typename TM::T res = 0.;
    for (unsigned i=0; i<M.nb_rows(); ++i)
    {
        for (unsigned j=0; j<M.nb_cols(); ++j)
        {
            res = LMT::max(res, LMT::abs(M(i, j)));
        }
    }
    return res;
}

/// pow_mat_sym
template<class TM> TM pow_mat_sym(const TM &M, const unsigned &p)
{
    LMT::Vec<typename TM::T> ones;
    ones.resize(M.nb_cols());
    ones.set(1.);
    LMT::Mat<typename TM::T> res = diag(ones);
    for (unsigned k=0; k<p; ++k)
    {
        res *= M;
    }
    return res;
}

/// tens_prod
template<class TV> LMT::Mat<typename LMT::TypeReduction<LMT::Multiplies,TV>::T> tens_prod(const TV &V1, const TV &V2) // pour le type on a pas de typename TV::T car les vecteurs peuvent être hétérogènes
{
    unsigned dim = V1.size();

    LMT::Mat<typename LMT::TypeReduction<LMT::Multiplies,TV>::T> res( dim );
    for (unsigned i=0; i<dim; ++i)
    {
        for (unsigned j=0; j<dim; ++j)
        {
            res(i, j) = V1[i] * V2[j];
        }
    }
    return res;
}

/// self_tens_prod
template<class TV> LMT::Mat<typename LMT::TypeReduction<LMT::Multiplies,TV>::T> self_tens_prod(const TV &V) // pour le type on a pas de typename TV::T car les vecteurs peuvent être hétérogènes
{
    return tens_prod(V, V);
}

/// mat_sym_to_vec_col
template<class TM> LMT::Vec<typename TM::T> mat_sym_to_vec_col(const TM &M)
{
    unsigned dim = M.nb_cols();

    if ( dim == 2 )
    {
        LMT::Vec<typename TM::T> V;
        V.resize(3);
        V[0] = M(0, 0);
        V[1] = M(1, 1);
        V[2] = M(0, 1);
        return V;
    }

    else if ( dim == 3 )
    {
        LMT::Vec<typename TM::T> V;
        V.resize(6);
        V[0] = M(0, 0);
        V[1] = M(1, 1);
        V[2] = M(2, 2);
        V[3] = M(0, 1);
        V[4] = M(0, 2);
        V[5] = M(1, 2);
        return V;
    }
}

/// vec_col_to_mat_sym
template<class TV> LMT::Mat<typename LMT::TypeReduction<LMT::Multiplies,TV>::T, LMT::Sym<> > vec_col_to_mat_sym(const TV &V)
{
    unsigned dim_vec = V.size();

    if (dim_vec == 3)
    {
        LMT::Mat<typename LMT::TypeReduction<LMT::Multiplies,TV>::T, LMT::Sym<> > M;
        M.resize(2);
        M.set(0.);
        M(0, 0) = V[0];
        M(1, 1) = V[1];
        M(0, 1) = V[2];
        return M;
    }

    else if (dim_vec == 6)
    {
        LMT::Mat<typename LMT::TypeReduction<LMT::Multiplies,TV>::T, LMT::Sym<> > M;
        M.resize(3);
        M.set(0.);
        M(0, 0) = V[0];
        M(1, 1) = V[1];
        M(2, 2) = V[2];
        M(0, 1) = V[3];
        M(0, 2) = V[4];
        M(1, 2) = V[5];
        return M;
    }
}

#endif
