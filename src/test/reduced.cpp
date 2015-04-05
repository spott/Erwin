#include <iostream>
#include <petsc.h>

struct context {
    Vec H0;
    Mat D;
    std::function<double(double)> E;
};

PetscErrorCode
RHSJacobian_function( TS ts, double t_, Vec u, Mat A, Mat B, void* G_u );
PetscErrorCode
Monitor_function( TS ts_, int step, double t, Vec u, void* monitor );

PetscErrorCode
RHSJacobian_function( TS, double t_ , Vec, Mat A, Mat, void* G_u )
{
    context* c = static_cast<context*>( G_u );
    PetscScalar e = c->E( t_ );
    MatCopy( c->D, A, SAME_NONZERO_PATTERN );
    MatScale( A, e );
    MatDiagonalSet( A, c->H0, INSERT_VALUES);
    MatScale( A, std::complex<double>( 0, -1 ) );

    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    return 0;
}

PetscErrorCode
Monitor_function( TS, int step, double t, Vec u, void* )
{
    //std::function<double(double)>* E = static_cast<std::function<double(double)>*> (ef);
    PetscReal n;
    VecNorm( u, NORM_2, &n );
    std::cout << "t: " << t << " step: " << step << " norm-1: " << n - 1 << std::endl; //" ef " << (*E)(t) << std::endl;
    return 0;
}

int main( int argc, char** argv )
{
    PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
    std::cout.precision(16);
    Vec H0;
    Mat D;
    PetscViewer view;

    VecCreate(PETSC_COMM_WORLD, &H0 );
    VecSetType( H0, VECSTANDARD);
    MatCreate( PETSC_COMM_WORLD, &D );
    MatSetType( D, MATAIJ );

    PetscViewerBinaryOpen( PETSC_COMM_WORLD, "hamiltonian/energy_eigenvalues_vector.dat", FILE_MODE_READ, &view );
    VecLoad( H0, view );
    PetscViewerBinaryOpen( PETSC_COMM_WORLD, "hamiltonian/dipole_matrix.dat", FILE_MODE_READ, &view );
    MatLoad( D, view );

    //10 cycle, 800nm pulse
    std::function<double(double)> E = [](double t) {
        using namespace std;
        //return 0 * t;
        return .1 * pow( sin( .06 * 2. * 3.14159265 * t / 3), 2) * sin(.06 * t);
    };

    context c{H0, D, E};
    Mat A;
    MatDuplicate( D, MAT_SHARE_NONZERO_PATTERN, &A );
    MatDiagonalSet( A, c.H0, INSERT_VALUES);

    Vec psi0;
    MatGetVecs( A, &psi0, PETSC_NULL );
    VecSetValue(psi0, 0, 1., INSERT_VALUES);
    VecAssemblyBegin(psi0);
    VecAssemblyEnd(psi0);


    TS tss;
    TSCreate( PETSC_COMM_WORLD, &tss );
    TSSetProblemType( tss, TS_LINEAR );
    TSSetType( tss, TSCN );
    TSSetInitialTimeStep( tss, 0, 0.01 );
    TSSetDuration( tss, 5, 0.03 );
    TSSetFromOptions( tss );
    TSSetRHSFunction( tss, NULL, TSComputeRHSFunctionLinear, NULL );
    TSSetRHSJacobian( tss, A, A, &RHSJacobian_function, &c );
    TSMonitorSet(tss, &Monitor_function, &E, PETSC_NULL);
    TSSolve( tss, psi0);

    MatDestroy(&A);


    // my own solve:
    Vec tmp, psi1;
    KSP ksp;
    PetscScalar e;
    PetscReal t_,norm;
    KSPCreate( PETSC_COMM_WORLD, &ksp );
    KSPSetTolerances( ksp, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
    KSPSetFromOptions( ksp );

    MatDuplicate( D, MAT_SHARE_NONZERO_PATTERN, &A );
    MatGetVecs( A, &psi1, PETSC_NULL );
    MatGetVecs( A, &tmp, PETSC_NULL );
    VecSet(psi0, 0);
    VecSetValue(psi0, 0, 1., INSERT_VALUES);
    VecAssemblyBegin(psi0);
    VecAssemblyEnd(psi0);

    //first step:

    //A(t) = -i(H0)
    t_ = 0;
    e = c.E( t_ );
    MatCopy( c.D, A, SAME_NONZERO_PATTERN );
    MatScale( A, e );
    MatDiagonalSet( A, c.H0, INSERT_VALUES);
    MatScale( A, std::complex<double>( 0, -1 ) );
    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // 1 + dt A(t) / 2.
    MatScale(A, .01 / 2.);
    MatShift(A, 1.);

    MatMult(A, psi0, tmp);


    //Other side:
    t_ = 0.01;
    e = c.E( t_ );
    MatCopy( c.D, A, SAME_NONZERO_PATTERN );
    MatScale( A, e );
    MatDiagonalSet( A, c.H0, INSERT_VALUES);
    MatScale( A, std::complex<double>( 0, 1 ) );
    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // 1 - dt A(t) / 2.
    MatScale(A, .01 / 2.);
    MatShift(A, 1.);


    KSPSetOperators(ksp, A, A);
    KSPSolve(ksp, tmp, psi1);

    VecNorm(psi1, NORM_2, &norm);
    std::cout << "t: " << t_ << " norm-1: " << norm - 1 << " ef " << E(t_) << std::endl;

    //step 2
    t_ = 0.01;
    e = c.E( t_ );
    MatCopy( c.D, A, SAME_NONZERO_PATTERN );
    MatScale( A, e );
    MatDiagonalSet( A, c.H0, INSERT_VALUES);
    MatScale( A, std::complex<double>( 0, -1 ) );
    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // 1 + dt A(t) / 2.
    MatScale(A, .01 / 2.);
    MatShift(A, 1.);

    MatMult(A, psi1, tmp);


    //Other side:
    t_ = 0.02;
    e = c.E( t_ );
    MatCopy( c.D, A, SAME_NONZERO_PATTERN );
    MatScale( A, e );
    MatDiagonalSet( A, c.H0, INSERT_VALUES);
    MatScale( A, std::complex<double>( 0, 1 ) );
    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // 1 - dt A(t) / 2.
    MatScale(A, .01 / 2.);
    MatShift(A, 1.);

    KSPSetOperators(ksp, A, A);
    KSPSolve(ksp, tmp, psi0);

    VecNorm(psi0, NORM_2, &norm);
    std::cout << "t: " << t_ << " norm-1: " << norm - 1 << " ef " << E(t_) << std::endl;

    //step 3
    t_ = 0.02;
    e = c.E( t_ );
    MatCopy( c.D, A, SAME_NONZERO_PATTERN );
    MatScale( A, e );
    MatDiagonalSet( A, c.H0, INSERT_VALUES);
    MatScale( A, std::complex<double>( 0, -1 ) );
    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // 1 + dt A(t) / 2.
    MatScale(A, .01 / 2.);
    MatShift(A, 1.);

    MatMult(A, psi0, tmp);


    //Other side:
    t_ = 0.03;
    e = c.E( t_ );
    MatCopy( c.D, A, SAME_NONZERO_PATTERN );
    MatScale( A, e );
    MatDiagonalSet( A, c.H0, INSERT_VALUES);
    MatScale( A, std::complex<double>( 0, 1 ) );
    //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // 1 - dt A(t) / 2.
    MatScale(A, .01 / 2.);
    MatShift(A, 1.);

    KSPSetOperators(ksp, A, A);
    KSPSolve(ksp, tmp, psi1);

    VecNorm(psi1, NORM_2, &norm);
    std::cout << "t: " << t_ << " norm-1: " << norm - 1 << " ef " << E(t_) << std::endl;
}
