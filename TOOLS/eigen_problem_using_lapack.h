#ifndef EIGEN_PROBLEM_USING_LAPACK_H
#define EIGEN_PROBLEM_USING_LAPACK_H

#include <containers/mat.h>

extern "C" {
    int dsyev_(
    char *jobz, char *uplo, int *n, double *a, int *lda, 
    double *w, double *work, int *lwork,
    int *info );
    int dsygv_(
    int *itype, char *jobz, char *uplo, int *n, double *a, int *lda, double *b, int *ldb, 
    double *w, double *work, int *lwork,
    int *info );
}

namespace LMT {

template<class TM, class TV> void eigen_values_using_lapack( const TM &A, TV &eig_val ) {
    char jobz = 'N'; // Compute eigenvalues only.
    char uplo = 'L'; // Lower triangle of A is stored.
    int n = A.nb_rows();
    Vec<double> A_data; A_data.resize( n * n );
    for(unsigned j=0,c=0;j<A.nb_cols();++j)
        for(unsigned i=0;i<A.nb_rows();++i,++c)
            A_data[ c ] = A( i, j );
    eig_val.resize( n );
    int lwork = 3*n-1;
    Vec<double> work; work.resize( lwork );
    int info;
    dsyev_( &jobz, &uplo, &n, A_data.ptr(), &n, eig_val.ptr(), work.ptr(), &lwork, &info );
}

template<class TM1, class TV, class TM2> void eigen_problem_using_lapack( const TM1 &A, TV &eig_val, TM2 &eig_vec ) {
    char jobz = 'V'; // Compute eigenvalues and eigenvectors.
    char uplo = 'L'; // Lower triangle of A is stored.
    int n = A.nb_rows();
    Vec<double> A_data; A_data.resize( n * n );
    for(unsigned j=0,c=0;j<A.nb_cols();++j)
        for(unsigned i=0;i<A.nb_rows();++i,++c)
            A_data[ c ] = A( i, j );
    eig_val.resize( n );
    int lwork = 3*n-1;
    Vec<double> work; work.resize( lwork );
    int info;
    dsyev_( &jobz, &uplo, &n, A_data.ptr(), &n, eig_val.ptr(), work.ptr(), &lwork, &info );
    if ( info )
        throw "pb with eig with lapack";
    //
    eig_vec.resize( n, n );
    for(int i=0,c=0;i<n;++i)
        for(int j=0;j<n;++j,++c)
            eig_vec( i, j ) = A_data[ c ];
}

template<class TM1, class TM2, class TV> void generalized_eigen_values_using_lapack( const TM1 &A, const TM2 &B, TV &eig_val ) {
    char jobz = 'N'; // Compute eigenvalues only.
    char uplo = 'L'; // Lower triangles of A and B are stored.
    int n = A.nb_rows();
    Vec<double> A_data; A_data.resize( n * n );
    for(unsigned j=0,c=0;j<A.nb_cols();++j)
        for(unsigned i=0;i<A.nb_rows();++i,++c)
            A_data[ c ] = A( i, j );
    Vec<double> B_data; B_data.resize( n * n );
    for(unsigned j=0,c=0;j<B.nb_cols();++j)
        for(unsigned i=0;i<B.nb_rows();++i,++c)
            B_data[ c ] = B( i, j );
    eig_val.resize( n );
    int lwork = 3*n-1;
    Vec<double> work; work.resize(lwork);
    int info;
    int itype = 1;
    dsygv_( &itype, &jobz, &uplo, &n, A_data.ptr(), &n, B_data.ptr(), &n, eig_val.ptr(), work.ptr(), &lwork, &info );
    if ( info )
        throw "pb with eig with lapack";
}

template<class TM1, class TM2, class TV, class TM3> void generalized_eigen_problem_using_lapack( const TM1 &A, const TM2 &B, TV &eig_val, TM3 &eig_vec ) {
    char jobz = 'V'; // Compute eigenvalues and eigenvectors.
    char uplo = 'L'; // Lower triangles of A and B are stored.
    int n = A.nb_rows();
    Vec<double> A_data; A_data.resize( n * n );
    for(unsigned j=0,c=0;j<A.nb_cols();++j)
        for(unsigned i=0;i<A.nb_rows();++i,++c)
            A_data[ c ] = A( i, j );
    Vec<double> B_data; B_data.resize( n * n );
    for(unsigned j=0,c=0;j<B.nb_cols();++j)
        for(unsigned i=0;i<B.nb_rows();++i,++c)
            B_data[ c ] = B( i, j );
    eig_val.resize( n );
    int lwork = 3*n-1;
    Vec<double> work; work.resize(lwork);
    int info;
    int itype = 1;
    dsygv_( &itype, &jobz, &uplo, &n, A_data.ptr(), &n, B_data.ptr(), &n, eig_val.ptr(), work.ptr(), &lwork, &info );
    if ( info )
        throw "pb with eig with lapack";
    //
    eig_vec.resize( n, n );
    for(int i=0,c=0;i<n;++i)
        for(int j=0;j<n;++j,++c)
            eig_vec( i, j ) = A_data[ c ];
}

}
#endif // EIGEN_PROBLEM_USING_LAPACK_H
