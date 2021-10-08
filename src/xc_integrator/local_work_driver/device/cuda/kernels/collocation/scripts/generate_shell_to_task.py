import os,sys,re,math
import pyexpander.lib as expander
from collocation_angular import generate_spherical_angular, generate_cartesian_angular, generate_cartesian_ls, generate_eval_lines
import sympy
import itertools

from io import StringIO

L_max = 6
if len(sys.argv) > 1: L_max = int(sys.argv[1])

def generate_shell_to_task_lines( ang, do_grad = False, do_hess = False ):
  [x,y,z,r] = sympy.symbols('x y z r', real=True)
  [bf,bf_x,bf_y,bf_z,bf_xx,bf_xy,bf_xz,bf_yy,bf_yz,bf_zz] = sympy.symbols('radial_eval radial_eval_x radial_eval_y radial_eval_z radial_eval_xx radial_eval_xy radial_eval_xz radial_eval_yy radial_eval_yz radial_eval_zz',real=True)

  bf_eval_strs = []
  bf_x_eval_strs = []
  bf_y_eval_strs = []
  bf_z_eval_strs = []
  bf_xx_eval_strs = []
  bf_xy_eval_strs = []
  bf_xz_eval_strs = []
  bf_yy_eval_strs = []
  bf_yz_eval_strs = []
  bf_zz_eval_strs = []
  for j in range(len(ang)):
    a       = ang[j]
    a_x = sympy.diff( a, x )
    a_y = sympy.diff( a, y )
    a_z = sympy.diff( a, z )
    a_xx = sympy.diff( a, x, x )
    a_xy = sympy.diff( a, x, y )
    a_xz = sympy.diff( a, x, z )
    a_yy = sympy.diff( a, y, y )
    a_yz = sympy.diff( a, y, z )
    a_zz = sympy.diff( a, z, z )

    bf_eval = sympy.simplify( a * bf )
    bf_x_eval = sympy.simplify( a_x * bf + a * bf_x )
    bf_y_eval = sympy.simplify( a_y * bf + a * bf_y )
    bf_z_eval = sympy.simplify( a_z * bf + a * bf_z )

    bf_xx_eval = sympy.simplify( a_xx * bf + 2 * a_x * bf_x + a * bf_xx )
    bf_yy_eval = sympy.simplify( a_yy * bf + 2 * a_y * bf_y + a * bf_yy )
    bf_zz_eval = sympy.simplify( a_zz * bf + 2 * a_z * bf_z + a * bf_zz )

    bf_xy_eval = sympy.simplify( a_xy * bf + a_x * bf_y + a_y * bf_x + a * bf_xy )
    bf_xz_eval = sympy.simplify( a_xz * bf + a_x * bf_z + a_z * bf_x + a * bf_xz )
    bf_yz_eval = sympy.simplify( a_yz * bf + a_y * bf_z + a_z * bf_y + a * bf_yz )

    #bf_eval_str = 'ang_eval = {};'.format(bf_eval)
    #bf_x_eval_str = 'dang_eval_x = {};'.format(bf_x_eval)
    #bf_y_eval_str = 'dang_eval_y = {};'.format(bf_y_eval)
    #bf_z_eval_str = 'dang_eval_z = {};'.format(bf_z_eval)
    bf_eval_str   = '{}'.format(bf_eval  )
    bf_x_eval_str = '{}'.format(bf_x_eval)
    bf_y_eval_str = '{}'.format(bf_y_eval)
    bf_z_eval_str = '{}'.format(bf_z_eval)

    bf_xx_eval_str = '{}'.format(bf_xx_eval)
    bf_xy_eval_str = '{}'.format(bf_xy_eval)
    bf_xz_eval_str = '{}'.format(bf_xz_eval)
    bf_yy_eval_str = '{}'.format(bf_yy_eval)
    bf_yz_eval_str = '{}'.format(bf_yz_eval)
    bf_zz_eval_str = '{}'.format(bf_zz_eval)

    if L >= 2:
      for k in range(2,L+1):
        for X in ('x','y','z'):
          pow_str = X + '**' + str(k)
          repl_str = ''
          for K in range(k-1): repl_str = repl_str + X + '*'
          repl_str = repl_str + X
    
          bf_eval_str = bf_eval_str.replace(pow_str,repl_str)
          bf_x_eval_str = bf_x_eval_str.replace(pow_str,repl_str)
          bf_y_eval_str = bf_y_eval_str.replace(pow_str,repl_str)
          bf_z_eval_str = bf_z_eval_str.replace(pow_str,repl_str)

          bf_xx_eval_str = bf_xx_eval_str.replace(pow_str,repl_str)
          bf_xy_eval_str = bf_xy_eval_str.replace(pow_str,repl_str)
          bf_xz_eval_str = bf_xz_eval_str.replace(pow_str,repl_str)
          bf_yy_eval_str = bf_yy_eval_str.replace(pow_str,repl_str)
          bf_yz_eval_str = bf_yz_eval_str.replace(pow_str,repl_str)
          bf_zz_eval_str = bf_zz_eval_str.replace(pow_str,repl_str)

    bf_eval_strs.append( bf_eval_str )
    bf_x_eval_strs.append( bf_x_eval_str )
    bf_y_eval_strs.append( bf_y_eval_str )
    bf_z_eval_strs.append( bf_z_eval_str )
    bf_xx_eval_strs.append( bf_xx_eval_str )
    bf_xy_eval_strs.append( bf_xy_eval_str )
    bf_xz_eval_strs.append( bf_xz_eval_str )
    bf_yy_eval_strs.append( bf_yy_eval_str )
    bf_yz_eval_strs.append( bf_yz_eval_str )
    bf_zz_eval_strs.append( bf_zz_eval_str )

  ret_vals = [ bf_eval_strs ]
  if do_grad:
    ret_vals.append( bf_x_eval_strs )
    ret_vals.append( bf_y_eval_strs )
    ret_vals.append( bf_z_eval_strs )

  if do_hess:
    ret_vals.append( bf_xx_eval_strs )
    ret_vals.append( bf_xy_eval_strs )
    ret_vals.append( bf_xz_eval_strs )
    ret_vals.append( bf_yy_eval_strs )
    ret_vals.append( bf_yz_eval_strs )
    ret_vals.append( bf_zz_eval_strs )

  return ret_vals



def get_constant_lines( lines ):
  constant_lines = []

  # Sqrts
  sqrt_regex = "sqrt\([0-9]+\)"
  sqrt_finds = list(set(re.findall( sqrt_regex, "\n".join(lines) )))

  # Replace locally
  for x in sqrt_finds:
    arg = x.strip('sqrt(').strip(')')
    new_str = 'sqrt_' + arg
    new_str = "constexpr double " + new_str + " = " + str( math.sqrt(int(arg)) )
    new_str = new_str + ";"
    constant_lines.append( new_str )

  return constant_lines

def sanitize_constants( lines ):
  # Sqrts
  sqrt_regex = "sqrt\([0-9]+\)"
  sqrt_finds = list(set(re.findall( sqrt_regex, "\n".join(lines) )))

  for x in sqrt_finds:
    arg = x.strip('sqrt(').strip(')')
    new_str = 'sqrt_' + arg
    lines = [line.replace(x, new_str) for line in lines]

  return lines


# Generate the evaluation lines
cart_bf_lines = []
sph_bf_lines  = []
cart_bfx_lines = []
cart_bfy_lines = []
cart_bfz_lines = []
sph_bfx_lines  = []
sph_bfy_lines  = []
sph_bfz_lines  = []
for L in range( L_max + 1 ):
  cart_ls  = generate_cartesian_ls(L)
  cart_ang = generate_cartesian_angular(cart_ls)
  sph_ang  = generate_spherical_angular( L, True )

  cart_bf_lines.append(generate_shell_to_task_lines(cart_ang))
  sph_bf_lines.append(generate_shell_to_task_lines(sph_ang))
  [bfx,bfy,bfz] = generate_shell_to_task_lines(cart_ang,True)
  cart_bfx_lines.append(bfx)
  cart_bfy_lines.append(bfy)
  cart_bfz_lines.append(bfz)
  [bfx,bfy,bfz] = generate_shell_to_task_lines(sph_ang,True)
  sph_bfx_lines.append(bfx)
  sph_bfy_lines.append(bfy)
  sph_bfz_lines.append(bfz)


constant_lines = []
for lines in itertools.chain( cart_bf_lines, sph_bf_lines ):
  _tmp = get_constant_lines( lines )
  for line in _tmp:
    constant_lines.append(line)

  
# Sanitize wrt constants
for i,lines in enumerate(cart_bf_lines):
  cart_bf_lines[i] = sanitize_constants( lines )
for i,lines in enumerate(sph_bf_lines):
  sph_bf_lines[i] = sanitize_constants( lines )
for i,lines in enumerate(cart_bfx_lines):
  cart_bfx_lines[i] = sanitize_constants( lines )
for i,lines in enumerate(cart_bfy_lines):
  cart_bfy_lines[i] = sanitize_constants( lines )
for i,lines in enumerate(cart_bfz_lines):
  cart_bfz_lines[i] = sanitize_constants( lines )
for i,lines in enumerate(sph_bfx_lines):
  sph_bfx_lines[i] = sanitize_constants( lines )
for i,lines in enumerate(sph_bfy_lines):
  sph_bfy_lines[i] = sanitize_constants( lines )
for i,lines in enumerate(sph_bfz_lines):
  sph_bfz_lines[i] = sanitize_constants( lines )


def generate_code( eval_lines, L, eval_type, template_fname, output_fname ):
  old_sysout = sys.stdout
  var_dict = { 'eval_lines' : eval_lines, 'L' : L, 'type' : eval_type }
  sys.stdout = expand = StringIO()
  expander.expandFile( template_fname, external_definitions = var_dict, auto_indent = True )
  expand = expand.getvalue()
  sys.stdout = old_sysout

  output_file = open(output_fname, 'w')
  output_file.write(expand)

def generate_code_gradient( eval_lines, eval_lines_dx, eval_lines_dy, eval_lines_dz, L, eval_type, template_fname, output_fname ):
  old_sysout = sys.stdout
  var_dict = { 'eval_lines' : eval_lines,
               'eval_lines_dx' : eval_lines_dx, 
               'eval_lines_dy' : eval_lines_dy, 
               'eval_lines_dz' : eval_lines_dz, 
               'L' : L, 'type' : eval_type }
  sys.stdout = expand = StringIO()
  expander.expandFile( template_fname, external_definitions = var_dict, auto_indent = True )
  expand = expand.getvalue()
  sys.stdout = old_sysout

  output_file = open(output_fname, 'w')
  output_file.write(expand)



# Generate kernels
for L in range( L_max + 1 ):
  template_fname = 'templates/collocation_shell_to_task_kernels.hpp'
  cart_header_fname = "collocation_shell_to_task_kernels_cartesian_l" + str(L) + ".hpp"
  sph_header_fname  = "collocation_shell_to_task_kernels_spherical_l" + str(L) + ".hpp"
  
  generate_code( cart_bf_lines[L], L, 'cartesian', template_fname, cart_header_fname )
  generate_code( sph_bf_lines[L], L, 'spherical', template_fname, sph_header_fname )

  cart_header_fname = "collocation_shell_to_task_kernels_cartesian_l" + str(L) + "_gradient.hpp"
  sph_header_fname  = "collocation_shell_to_task_kernels_spherical_l" + str(L) + "_gradient.hpp"
  generate_code_gradient( cart_bf_lines[L], cart_bfx_lines[L], cart_bfy_lines[L], cart_bfz_lines[L],
    L, 'cartesian_gradient', template_fname, cart_header_fname )
  generate_code_gradient( sph_bf_lines[L], sph_bfx_lines[L], sph_bfy_lines[L], sph_bfz_lines[L],
    L, 'spherical_gradient', template_fname, sph_header_fname )


  #template_fname = 'templates/collocation_shell_to_task_combined_kernels.hpp'
  #cart_header_fname = "collocation_shell_to_task_combined_kernels_cartesian_l" + str(L) + ".hpp"
  #sph_header_fname  = "collocation_shell_to_task_combined_kernels_spherical_l" + str(L) + ".hpp"
  #generate_code( cart_bf_lines[L], L, 'cartesian', template_fname, cart_header_fname )
  #generate_code( sph_bf_lines[L], L, 'spherical', template_fname, sph_header_fname )

  #cart_header_fname = "collocation_shell_to_task_combined_kernels_cartesian_l" + str(L) + "_gradient.hpp"
  #sph_header_fname  = "collocation_shell_to_task_combined_kernels_spherical_l" + str(L) + "_gradient.hpp"
  #generate_code_gradient( cart_bf_lines[L], cart_bfx_lines[L], cart_bfy_lines[L], cart_bfz_lines[L],
  #  L, 'cartesian_gradient', template_fname, cart_header_fname )
  #generate_code_gradient( sph_bf_lines[L], sph_bfx_lines[L], sph_bfy_lines[L], sph_bfz_lines[L],
  #  L, 'spherical_gradient', template_fname, sph_header_fname )
