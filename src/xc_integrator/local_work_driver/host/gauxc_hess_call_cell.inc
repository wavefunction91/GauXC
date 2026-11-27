// ==> BEGIN GENERATED CODE [xckernel gauxcwriter: cell-function call site] <==
// Arguments are bound BY NAME from each kernel's own signature;
// regenerate rather than edit: python -m xckernel.emitters.gauxcwriter --emit-dir <dir>
if( is_becke ) {
  xckernel::xck_gauxc_weight_cell_becke( mu, s, t, u );
} else {
  xckernel::xck_gauxc_weight_cell_ssf( mu, s, t, u );
}
// ==> END GENERATED CODE <==
