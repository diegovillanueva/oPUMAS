! 🚫🧠 BLAS stub for the column driver: only the neural-net emulator (warm_rain = 'emulated',
! out of scope for the Julia twin) calls dgemm; any call aborts instead of linking BLAS.
subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  implicit none
  character, intent(in) :: transa, transb
  integer, intent(in) :: m, n, k, lda, ldb, ldc
  double precision, intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
  double precision, intent(inout) :: c(ldc,*)
  error stop 'pumas_column_driver: dgemm stub called (emulated warm rain is not supported)'
end subroutine dgemm
