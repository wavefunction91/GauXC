// ==> BEGIN GENERATED CODE [xckernel gauxcwriter: pair call site] <==
// Arguments are bound BY NAME from each kernel's own signature;
// regenerate rather than edit: python -m xckernel.emitters.gauxcwriter --emit-dir <dir>
if( is_mgga ) {
  xckernel::xck_gauxc_hess_pair_mgga_tau( ROW(a,dx,0,ip), ROW(b,dy,0,ip), ROW(a,dx,1,ip), ROW(b,dy,1,ip), ROW(a,dx,2,ip), ROW(b,dy,2,ip), ROW(a,dx,3,ip), ROW(a,dx,4,ip), ROW(a,dx,5,ip), ROW(b,dy,3,ip), ROW(b,dy,4,ip), ROW(b,dy,5,ip), v2rho2[ip], v2rhogamma[ip], v2rhotau[ip], v2gamma2[ip], v2gammatau[ip], v2tau2[ip], vgamma[ip], h );
} else if( is_gga ) {
  xckernel::xck_gauxc_hess_pair_gga( ROW(a,dx,0,ip), ROW(b,dy,0,ip), ROW(a,dx,1,ip), ROW(b,dy,1,ip), ROW(a,dx,2,ip), ROW(a,dx,3,ip), ROW(a,dx,4,ip), ROW(b,dy,2,ip), ROW(b,dy,3,ip), ROW(b,dy,4,ip), v2rho2[ip], v2rhogamma[ip], v2gamma2[ip], vgamma[ip], h );
} else {
  xckernel::xck_gauxc_hess_pair_lda( ROW(a,dx,0,ip), ROW(b,dy,0,ip), v2rho2[ip], h );
}
// ==> END GENERATED CODE <==
