import os,sys,re,math
import pyexpander.lib as expander
from collocation_angular import generate_spherical_angular, generate_cartesian_angular, generate_cartesian_ls, generate_eval_lines

from io import StringIO

L_max = 6
if len(sys.argv) > 1: L_max = int(sys.argv[1])

#sphr_bf_body = []
#sphr_bf_d1_body = []

sphr_unnorm_bf_body = []
sphr_unnorm_bf_d1_body = []

cart_bf_body = []
cart_bf_d1_body = []


for L in range( L_max + 1 ):

  print( "Processing L = {} ...".format(L) )
  #sphr_ang        = generate_spherical_angular( L, False )
  sphr_unnorm_ang = generate_spherical_angular( L, True  )
  cart_ang        = generate_cartesian_angular( generate_cartesian_ls( L ) )

  #sa, sa_x, sa_y, sa_z     = generate_eval_lines( L, sphr_ang )
  sna, sna_x, sna_y, sna_z = generate_eval_lines( L, sphr_unnorm_ang )
  ca, ca_x, ca_y, ca_z     = generate_eval_lines( L, cart_ang )

  #sphr_bf_body.append( "\n  ".join(sa) )
  sphr_unnorm_bf_body.append( "\n  ".join(sna) )
  cart_bf_body.append( "\n  ".join(ca) )

  #s_d1  = "\n\n  ".join(["\n  ".join( sa_x ),  "\n  ".join(sa_y),  "\n  ".join(sa_z)])
  sn_d1 = "\n\n  ".join(["\n  ".join( sna_x ), "\n  ".join(sna_y), "\n  ".join(sna_z)])
  c_d1  = "\n\n  ".join(["\n  ".join( ca_x ),  "\n  ".join(ca_y),  "\n  ".join(ca_z)])

  #sphr_bf_d1_body.append( s_d1 )
  sphr_unnorm_bf_d1_body.append( sn_d1 )
  cart_bf_d1_body.append( c_d1 )



template_fname = 'templates/collocation_angular_template.hpp'

#sphr_var_dict = { 'L_max' : L_max, 'body' : sphr_bf_body, 'body_d1' : sphr_bf_d1_body, 'name' : 'spherical' }
sphr_unnorm_var_dict = { 'L_max' : L_max, 'body' : sphr_unnorm_bf_body, 'body_d1' : sphr_unnorm_bf_d1_body, 'name' : 'spherical_unnorm' }
cart_var_dict = { 'L_max' : L_max, 'body' : cart_bf_body, 'body_d1' : cart_bf_d1_body, 'name' : 'cartesian' }


old_sys_out = sys.stdout

sys.stdout = cart_expand = StringIO()
expander.expandFile( template_fname, external_definitions=cart_var_dict, auto_indent=True )
#sys.stdout = sphr_expand = StringIO()
#expander.expandFile( template_fname, external_definitions=sphr_var_dict, auto_indent=True )
sys.stdout = sphr_unnorm_expand = StringIO()
expander.expandFile( template_fname, external_definitions=sphr_unnorm_var_dict, auto_indent=True )

sys.stdout = old_sys_out

cart_expand = cart_expand.getvalue()
#sphr_expand = sphr_expand.getvalue()
sphr_unnorm_expand = sphr_unnorm_expand.getvalue()



# Handle Constants
constant_lines = []

# Sqrts
sqrt_regex = "sqrt\([0-9]+\)"
#sqrt_finds = re.findall( sqrt_regex, "\n".join([cart_expand,sphr_expand,sphr_unnorm_expand]) )
sqrt_finds = re.findall( sqrt_regex, "\n".join([cart_expand,sphr_unnorm_expand]) )

sqrt_finds = list(set(sqrt_finds))

for x in sqrt_finds:
  arg = x.strip('sqrt(').strip(')')
  new_str = 'sqrt_' + arg
  
  cart_expand = cart_expand.replace( x, new_str )
  #sphr_expand = sphr_expand.replace( x, new_str )
  sphr_unnorm_expand = sphr_unnorm_expand.replace( x, new_str )

  new_str = "constexpr double " + new_str + " = " + str( math.sqrt(int(arg)) ) + ";"
  constant_lines.append( new_str )

old_sys_out = sys.stdout

sys.stdout = constant_expand = StringIO()
expander.expandFile( 'templates/collocation_device_constants_template.hpp',
  external_definitions={ "const_lines" : constant_lines } )

sys.stdout = old_sys_out

constant_expand = constant_expand.getvalue()


cart_header_fname = "collocation_angular_cartesian.hpp"
#sphr_header_fname = "collocation_angular_spherical.hpp"
sphr_unnorm_header_fname = "collocation_angular_spherical_unnorm.hpp"
cons_header_fname = "collocation_device_constants.hpp"

cart_header_file = open( cart_header_fname, 'w' )
#sphr_header_file = open( sphr_header_fname, 'w' )
sphr_unnorm_header_file = open( sphr_unnorm_header_fname, 'w' )
cons_header_file = open( cons_header_fname, 'w' )

cart_header_file.write( cart_expand )
#sphr_header_file.write( sphr_expand )
sphr_unnorm_header_file.write( sphr_unnorm_expand )
cons_header_file.write( constant_expand )









# Generate Kernel Driver

#old_sys_out = sys.stdout

#sys.stdout = collocation_cartesian_kernel_expand = StringIO()
#expander.expandFile( 'collocation_kernels_template.cu', external_definitions={ 'ang_name' : 'cartesian' } )
#
#sys.stdout = collocation_spherical_kernel_expand = StringIO()
#expander.expandFile( 'collocation_kernels_template.cu', external_definitions={ 'ang_name' : 'spherical' } )
#
#sys.stdout = collocation_spherical_unnorm_kernel_expand = StringIO()
#expander.expandFile( 'collocation_kernels_template.cu', external_definitions={ 'ang_name' : 'spherical_unnorm' } )
#
#sys.stdout = old_sys_out
#
#collocation_cartesian_kernel_expand = collocation_cartesian_kernel_expand.getvalue()
#collocation_spherical_kernel_expand = collocation_spherical_kernel_expand.getvalue()
#collocation_spherical_unnorm_kernel_expand = collocation_spherical_unnorm_kernel_expand.getvalue()
#
#with open( 'collocation_kernels_cartesian.cu', 'w' ) as f:
#  f.write( collocation_cartesian_kernel_expand )
#with open( 'collocation_kernels_spherical.cu', 'w' ) as f:
#  f.write( collocation_spherical_kernel_expand )
#with open( 'collocation_kernels_spherical_unnorm.cu', 'w' ) as f:
#  f.write( collocation_spherical_unnorm_kernel_expand )
