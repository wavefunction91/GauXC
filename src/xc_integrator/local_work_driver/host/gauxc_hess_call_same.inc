// ==> BEGIN GENERATED CODE [xckernel gauxcwriter: same call site] <==
// Arguments are bound BY NAME from each kernel's own signature;
// regenerate rather than edit: python -m xckernel.emitters.gauxcwriter --emit-dir <dir>
if( is_mgga ) {
  xckernel::xck_gauxc_hess_same_mgga_tau( xmat[k], xmat_x[k], xmat_y[k], xmat_z[k], d2c, d3c[0], d3c[1], d3c[2], dden_x[ip], dden_y[ip], dden_z[ip], vrho[ip], vgamma[ip], vtau[ip], weights[ip], sv );
} else if( is_gga ) {
  xckernel::xck_gauxc_hess_same_gga( xmat[k], xmat_x[k], xmat_y[k], xmat_z[k], d2c, d3c[0], d3c[1], d3c[2], dden_x[ip], dden_y[ip], dden_z[ip], vrho[ip], vgamma[ip], weights[ip], sv );
} else {
  xckernel::xck_gauxc_hess_same_lda( xmat[k], d2c, vrho[ip], weights[ip], sv );
}
// ==> END GENERATED CODE <==
