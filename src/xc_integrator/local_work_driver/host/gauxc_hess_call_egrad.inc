// ==> BEGIN GENERATED CODE [xckernel gauxcwriter: egrad call site] <==
// Arguments are bound BY NAME from each kernel's own signature;
// regenerate rather than edit: python -m xckernel.emitters.gauxcwriter --emit-dir <dir>
if( is_mgga ) {
  xckernel::xck_gauxc_hess_egrad_mgga_tau( ROW(a,d,0,ip), ROW(a,d,1,ip), ROW(a,d,2,ip), vrho[ip], vgamma[ip], vtau[ip], de );
} else if( is_gga ) {
  xckernel::xck_gauxc_hess_egrad_gga( ROW(a,d,0,ip), ROW(a,d,1,ip), vrho[ip], vgamma[ip], de );
} else {
  xckernel::xck_gauxc_hess_egrad_lda( ROW(a,d,0,ip), vrho[ip], de );
}
// ==> END GENERATED CODE <==
