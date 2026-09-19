// ==> BEGIN GENERATED CODE [xckernel gauxcwriter: rows call site] <==
// Arguments are bound BY NAME from each kernel's own signature;
// regenerate rather than edit: python -m xckernel.emitters.gauxcwriter --emit-dir <dir>
if( is_mgga ) {
  xckernel::xck_gauxc_hess_rows_mgga_tau( U0, U1, U2, U3, dchi, ddchi[0], ddchi[1], ddchi[2], dden_x[ip], dden_y[ip], dden_z[ip], F_rho, F_sigma, F_tau, Gx, Gy, Gz );
  ROW(a,d,0,ip) += F_rho;
  ROW(a,d,1,ip) += F_sigma;
  ROW(a,d,2,ip) += F_tau;
  ROW(a,d,3,ip) += Gx;
  ROW(a,d,4,ip) += Gy;
  ROW(a,d,5,ip) += Gz;
} else if( is_gga ) {
  xckernel::xck_gauxc_hess_rows_gga( U0, U1, U2, U3, dchi, ddchi[0], ddchi[1], ddchi[2], dden_x[ip], dden_y[ip], dden_z[ip], F_rho, F_sigma, Gx, Gy, Gz );
  ROW(a,d,0,ip) += F_rho;
  ROW(a,d,1,ip) += F_sigma;
  ROW(a,d,2,ip) += Gx;
  ROW(a,d,3,ip) += Gy;
  ROW(a,d,4,ip) += Gz;
} else {
  xckernel::xck_gauxc_hess_rows_lda( U0, U1, U2, U3, dchi, ddchi[0], ddchi[1], ddchi[2], F_rho, Gx, Gy, Gz );
  ROW(a,d,0,ip) += F_rho;
  ROW(a,d,1,ip) += Gx;
  ROW(a,d,2,ip) += Gy;
  ROW(a,d,3,ip) += Gz;
}
// ==> END GENERATED CODE <==
