// ==> BEGIN GENERATED CODE [xckernel gauxcwriter: pulayW call site] <==
// Arguments are bound BY NAME from each kernel's own signature;
// regenerate rather than edit: python -m xckernel.emitters.gauxcwriter --emit-dir <dir>
if( is_mgga ) {
  xckernel::xck_gauxc_hess_pulayW_mgga_tau( dden_x[ip], dden_y[ip], dden_z[ip], vrho[ip], vgamma[ip], vtau[ip], weights[ip], W[0], W[1], W[2], W[3], W[4], W[5], W[8], W[10], W[12], W[15] );
} else if( is_gga ) {
  xckernel::xck_gauxc_hess_pulayW_gga( dden_x[ip], dden_y[ip], dden_z[ip], vrho[ip], vgamma[ip], weights[ip], W[0], W[1], W[2], W[3], W[4], W[8], W[12] );
} else {
  xckernel::xck_gauxc_hess_pulayW_lda( vrho[ip], weights[ip], W[0] );
}
// ==> END GENERATED CODE <==
