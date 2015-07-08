#ifndef eigen_problem_using_lapack_h
#define eigen_problem_using_lapack_h

#include <containers/mat.h>

extern "C" {
    void dsyev_(const char &JOBZ, const char &UPLO, const int &N, double* A, const int &LDA, double* W, double* WORK, const int &LWORK, int &INFO);
    void dsygv_(const int &ITYPE, const char &JOBZ, const char &UPLO, const int &N, double* A, const int &LDA, double* B, const int &LDB, double* W, double* WORK, const int &LWORK, int &INFO);
};

template<class MatType, class VecType> void eigen_values_using_lapack(const MatType &A, VecType &V)
{
    unsigned dim = A.nb_cols();
//     PRINT(dim);

    LMT::Vec<double> A_data;
    A_data.resize(dim*dim);
    unsigned A_cpt = 0;
    for(unsigned j=0; j<dim; j++)
    {
        for(unsigned i=0; i<dim; i++)
        {
            A_data[A_cpt++] = A(i, j);
        }
    }
//     PRINT(A);
//     PRINT(A_data);

    unsigned lwork = 3*dim-1;
//     PRINT(lwork);

    LMT::Vec<double> work;
    work.resize(lwork);
//     PRINT(work);

    int info = 0;
//     PRINT(info);

    dsyev_('N', 'U', dim, A_data.ptr(), dim, V.ptr(), work.ptr(), lwork, info);
//     PRINT(info);
//     PRINT(A_data);
//     PRINT(V);
}

template<class MatType1, class VecType, class MatType2> void eigen_problem_using_lapack(const MatType1 &A, VecType &V, MatType2 &P)
{
    unsigned dim = A.nb_cols();
//     PRINT(dim);

    LMT::Vec<double> A_data;
    A_data.resize(dim*dim);
    unsigned A_cpt = 0;
    for(unsigned j=0; j<dim; j++)
    {
        for(unsigned i=0; i<dim; i++)
        {
            A_data[A_cpt++] = A(i, j);
        }
    }
//     PRINT(A);
//     PRINT(A_data);

//     PRINT(V);

    unsigned lwork = 3*dim-1;
//     PRINT(lwork);

    LMT::Vec<double> work;
    work.resize(lwork);
    work.set(0.);
//     PRINT(work);

    int info = 0;
//     PRINT(info);

    dsyev_('V', 'U', dim, A_data.ptr(), dim, V.ptr(), work.ptr(), lwork, info);
//     PRINT(info);
//     PRINT(A_data);
//     PRINT(V);

    unsigned P_cpt = 0;
    for(unsigned j=0; j<dim; j++)
    {
        for(unsigned i=0; i<dim; i++)
        {
            P(i, j) = A_data[P_cpt++];
        }
    }
//     PRINT(P);
}

template<class MatType1, class MatType2, class VecType> void generalized_eigen_values_using_lapack(const MatType1 &A, const MatType2 &B, VecType &V)
{
    unsigned dim = A.nb_cols();
//     PRINT(dim);

    LMT::Vec<double> A_data;
    A_data.resize(dim*dim);
    unsigned A_cpt = 0;
    for(unsigned j=0; j<dim; j++)
    {
        for(unsigned i=0; i<dim; i++)
        {
            A_data[A_cpt++] = A(i, j);
        }
    }
//     PRINT(A);
//     PRINT(A_data);

    LMT::Vec<double> B_data;
    B_data.resize(dim*dim);
    unsigned B_cpt = 0;
    for(unsigned j=0; j<dim; j++)
    {
        for(unsigned i=0; i<dim; i++)
        {
            B_data[B_cpt++] = B(i, j);
        }
    }
//     PRINT(B);
//     PRINT(B_data);

    unsigned lwork = 3*dim-1;
//     PRINT(lwork);

    LMT::Vec<double> work;
    work.resize(lwork);
//     PRINT(work);

    int info = 0;
//     PRINT(info);

    dsygv_(1, 'N', 'U', dim, A_data.ptr(), dim, B_data.ptr(), dim, V.ptr(), work.ptr(), lwork, info);
//     PRINT(info);
//     PRINT(A_data);
//     PRINT(B_data);
//     PRINT(V);
}

template<class MatType1, class MatType2, class VecType, class MatType3> void generalized_eigen_problem_using_lapack(const MatType1 &A, const MatType2 &B, VecType &V, MatType3 &P)
{
    unsigned dim = A.nb_cols();
//     PRINT(dim);

    LMT::Vec<double> A_data;
    A_data.resize(dim*dim);
    unsigned A_cpt = 0;
    for(unsigned j=0; j<dim; j++)
    {
        for(unsigned i=0; i<dim; i++)
        {
            A_data[A_cpt++] = A(i, j);
        }
    }
//     PRINT(A);
//     PRINT(A_data);

    LMT::Vec<double> B_data;
    B_data.resize(dim*dim);
    unsigned B_cpt = 0;
    for(unsigned j=0; j<dim; j++)
    {
        for(unsigned i=0; i<dim; i++)
        {
            B_data[B_cpt++] = B(i, j);
        }
    }
//     PRINT(B);
//     PRINT(B_data);

    unsigned lwork = 3*dim-1;
//     PRINT(lwork);

    LMT::Vec<double> work;
    work.resize(lwork);
//     PRINT(work);

    int info = 0;
//     PRINT(info);

    dsygv_(1, 'V', 'U', dim, A_data.ptr(), dim, B_data.ptr(), dim, V.ptr(), work.ptr(), lwork, info);
//     PRINT(info);
//     PRINT(A_data);
//     PRINT(B_data);
//     PRINT(V);

    unsigned P_cpt = 0;
    for(unsigned j=0; j<dim; j++)
    {
        for(unsigned i=0; i<dim; i++)
        {
            P(i, j) = A_data[P_cpt++];
        }
    }
//     PRINT(P);
}

#endif
