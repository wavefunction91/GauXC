import os,sys
import re
import sympy
import math
import cmath
from math import factorial as fact
from sympy import factorial as symb_fact
from sympy import factorial2 as symb_fact2
from scipy.special import binom as binomial
from sympy import exp as symb_exp
from sympy import I as symb_I

def generate_cartesian_ls( L ):
  
  l   = []
  for i in range(L+1):
    lx = L - i
    for j in range(i+1):
      ly = i - j
      lz = L - lx - ly
  
      l.append( [0, 0, 0] )
  
      for k in range( lx - 1 ):
        l[-1][0] = l[-1][0] + 1
      for k in range( ly - 1 ):
        l[-1][1] = l[-1][1] + 1
      for k in range( lz - 1 ):
        l[-1][2] = l[-1][2] + 1
  
      if lx > 0:
        l[-1][0] = l[-1][0] + 1
      if ly > 0:
        l[-1][1] = l[-1][1] + 1
      if lz > 0:
        l[-1][2] = l[-1][2] + 1
        
  return l


def generate_spherical_coeff( l, m, lx, ly, lz ):
  j = (lx + ly - abs(m))
  if j%2 == 0:
    j = int(j / 2)
  else:
    return 0.0

  prefactor = fact(2.*lx) * fact(2.*ly) * fact(2.*lz) * fact(l) 
  prefactor = prefactor * fact( l - abs(m) )
  prefactor = prefactor / (fact(2.*l) * fact(lx) * fact(ly) * fact(lz))
  prefactor = prefactor / fact( l + abs(m) )
  prefactor = math.sqrt( prefactor )


  term1 = 0.0
  for i in range( int((l - abs(m))/2)+1 ):
    term1 = term1 + binomial(l,i) * binomial(i,j) * \
                    math.pow(-1,i) * fact( 2*l - 2*i ) / \
                    fact( l - abs(m) - 2*i )
              
  term1 = term1 / math.pow(2,l) / fact(l)

  m_fact = 1.
  if m < 0:
    m_fact = -1.
                             
  term2 = 0.0 + 0.0j
  for k in range( j+1 ):
    z = cmath.exp( m_fact * math.pi / 2. * (abs(m) - lx + 2*k) * 1.j )
    term2 = term2 + binomial(j,k) * binomial( abs(m), lx - 2*k ) * z

  val = prefactor * term1 * term2

  if abs(val.real) < 1e-10:
    val = 0.0 + val.imag*1j
  if abs(val.imag) < 1e-10:
    val = val.real

  return val


def generate_spherical_coeff_symb( l, m, lx, ly, lz, unnorm = False ):
  j = (lx + ly - abs(m))
  if j%2 == 0:
    j = int(j / 2)
  else:
    return sympy.Integer(0)


  j_symb  = sympy.Integer(j)
  l_symb  = sympy.Integer(l)
  m_symb  = sympy.Integer( abs(m) )
  lx_symb = sympy.Integer(lx)
  ly_symb = sympy.Integer(ly)
  lz_symb = sympy.Integer(lz)



  prefactor = symb_fact(2*lx_symb) * symb_fact(2*ly_symb) * symb_fact(2*lz_symb) * symb_fact(l_symb) 
  prefactor = prefactor * symb_fact( l_symb - m_symb )
  prefactor = prefactor / (symb_fact(2*l_symb) * symb_fact(lx_symb) * symb_fact(ly_symb) * symb_fact(lz_symb))
  prefactor = prefactor / symb_fact( l_symb + m_symb )


  # Ed's stupid normalization convention...
  if unnorm:
    prefactor = prefactor * symb_fact2( 2*l - 1 ) / symb_fact2( 2*lx - 1 ) / symb_fact2(2*ly - 1) / symb_fact2( 2*lz - 1 )

  prefactor = sympy.sqrt( prefactor )


  term1 = sympy.Integer(0)
  for i in range( int((l - abs(m))/2)+1 ):
    term1 = term1 + sympy.Integer(binomial(l,i)) * sympy.Integer(binomial(i,j)) * \
                    sympy.Integer(math.pow(-1,i)) * symb_fact( 2*l_symb - sympy.Integer(2*i) ) / \
                    symb_fact( l_symb - m_symb - sympy.Integer(2*i) )
              
  term1 = term1 / (2**l_symb) / symb_fact(l)


  m_fact_symb = sympy.Integer(1)
  if m < 0:
    m_fact_symb = - m_fact_symb
                             
  term2 = sympy.Integer(0)
  for k in range( j+1 ):
    z = sympy.exp( m_fact_symb * sympy.pi / 2 * (m_symb - lx_symb + sympy.Integer(2*k)) * symb_I )
    term2 = term2 + sympy.Integer(binomial(j,k)) * sympy.Integer(binomial( abs(m), lx - 2*k )) * z


  return prefactor * term1 * term2



def generate_cartesian_angular( ls ):
  [x,y,z,r] = sympy.symbols('x y z r', real=True)
  
  ang = []

  for l in ls:
    ang.append(r)
    for i in range( l[0] ): ang[-1] = ang[-1] * x
    for i in range( l[1] ): ang[-1] = ang[-1] * y
    for i in range( l[2] ): ang[-1] = ang[-1] * z

    ang[-1] = ang[-1] / r

  return ang

def generate_spherical_angular( L, unnorm = False ):
  ls  = generate_cartesian_ls( L )
  angs = generate_cartesian_angular( ls )

  #r = sympy.symbols( 'r' )
  sph_angs = []
  for m in range( L + 1 ):
    tmp_p = 0
    tmp_m = 0
    for i in range(len(ls)):
      l   = ls[i]
      ang = angs[i]
      
      #c = generate_spherical_coeff( L, m, l[0],l[1],l[2] )
      c = generate_spherical_coeff_symb( L, m, l[0],l[1],l[2], unnorm )

      if m == 0:

        tmp_p = tmp_p + c * ang

      else:

        c_p = ( c + sympy.conjugate(c) ) / sympy.sqrt(2)
        c_m = ( c - sympy.conjugate(c) ) / sympy.sqrt(2) / symb_I

        tmp_p = tmp_p + c_p * ang
        tmp_m = tmp_m + c_m * ang

    sph_angs.append( (m, tmp_p) )
    if m > 0:
      sph_angs.append( (-m, tmp_m) )

  sph_angs = sorted( sph_angs, key=lambda x: x[0] ) 

  sph_angs_bare = []
  for a in sph_angs:
    sph_angs_bare.append( sympy.simplify(a[1]) )

  return sph_angs_bare


def generate_eval_lines( L, ang ):
  [x,y,z,r] = sympy.symbols('x y z r', real=True)
  [bf,bf_x,bf_y,bf_z] = sympy.symbols('bf bf_x bf_y bf_z',real=True)
  
  bf_eval_strs = []
  bf_x_eval_strs = []
  bf_y_eval_strs = []
  bf_z_eval_strs = []

  for j in range(len(ang)):
  
    a   = ang[j]
    a_x = sympy.diff( a, x )
    a_y = sympy.diff( a, y )
    a_z = sympy.diff( a, z )
  
  
  
    bf_eval = sympy.simplify( a * bf )
    bf_x_eval = sympy.simplify( a_x * bf + a * bf_x )
    bf_y_eval = sympy.simplify( a_y * bf + a * bf_y )
    bf_z_eval = sympy.simplify( a_z * bf + a * bf_z )
  
  
    bf_eval_str = 'eval[{}] = {};'.format(j,bf_eval)
    bf_x_eval_str = 'eval_x[{}] = {};'.format(j,bf_x_eval)
    bf_y_eval_str = 'eval_y[{}] = {};'.format(j,bf_y_eval)
    bf_z_eval_str = 'eval_z[{}] = {};'.format(j,bf_z_eval)
  
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
  
  
    bf_eval_strs  .append(bf_eval_str  )
    bf_x_eval_strs.append(bf_x_eval_str)
    bf_y_eval_strs.append(bf_y_eval_str)
    bf_z_eval_strs.append(bf_z_eval_str)
  
  return (bf_eval_strs, bf_x_eval_strs, bf_y_eval_strs, bf_z_eval_strs)




cart_header_fname = "gaueval_angular_cartesian.hpp"
sphr_header_fname = "gaueval_angular_spherical.hpp"
cons_header_fname = "gaueval_device_constants.hpp"

cart_header_file = open( cart_header_fname, 'w' )
sphr_header_file = open( sphr_header_fname, 'w' )
cons_header_file = open( cons_header_fname, 'w' )

L_max = 4
do_libint_norm = False
#do_libint_norm = True

preamble = """
#pragma once
#include "gaueval_device_constants.hpp"

#define GPGAUEVAL_INLINE __inline__

namespace gpgaueval {
"""



cart_header_file.write( preamble )
sphr_header_file.write( preamble )

cartesian_bf_template = """
GPGAUEVAL_INLINE __device__ void generate_cartesian_angular{}(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {{
"""

cartesian_bf_deriv1_template = """
GPGAUEVAL_INLINE __device__ void generate_cartesian_angular{}_deriv1(
  const double bf,
  const double bf_x,
  const double bf_y,
  const double bf_z,
  const double x,
  const double y,
  const double z,
  double* eval_x,
  double* eval_y,
  double* eval_z
) {{
"""

spherical_bf_template        = cartesian_bf_template.replace('cartesian','spherical')
spherical_bf_deriv1_template = cartesian_bf_deriv1_template.replace('cartesian','spherical')



constant_lines = []
for L in range( L_max + 1 ):
  
  sph_ang = generate_spherical_angular(L, do_libint_norm)
  car_ang = generate_cartesian_angular( generate_cartesian_ls(L) )
  
  sph_bf_eval_strs, sph_bf_x_eval_strs, sph_bf_y_eval_strs, sph_bf_z_eval_strs = generate_eval_lines( L, sph_ang )
  car_bf_eval_strs, car_bf_x_eval_strs, car_bf_y_eval_strs, car_bf_z_eval_strs = generate_eval_lines( L, car_ang )
  

  cartesian_bf_prototype = cartesian_bf_template.format( "_" + str(L) )
  spherical_bf_prototype = spherical_bf_template.format( "_" + str(L) )
  cartesian_bf_deriv1_prototype = cartesian_bf_deriv1_template.format( "_" + str(L) )
  spherical_bf_deriv1_prototype = spherical_bf_deriv1_template.format( "_" + str(L) )


  spherical_bf_func = spherical_bf_prototype + "\n"
  for s in sph_bf_eval_strs: spherical_bf_func = spherical_bf_func + "  " + s + "\n"
  spherical_bf_func = spherical_bf_func + "\n}\n"

  spherical_bf_deriv1_func = spherical_bf_deriv1_prototype + "\n"
  for s in sph_bf_x_eval_strs: spherical_bf_deriv1_func = spherical_bf_deriv1_func + "  " + s + "\n"
  spherical_bf_deriv1_func = spherical_bf_deriv1_func + "\n"
  for s in sph_bf_y_eval_strs: spherical_bf_deriv1_func = spherical_bf_deriv1_func + "  " + s + "\n"
  spherical_bf_deriv1_func = spherical_bf_deriv1_func + "\n"
  for s in sph_bf_z_eval_strs: spherical_bf_deriv1_func = spherical_bf_deriv1_func + "  " + s + "\n"
  spherical_bf_deriv1_func = spherical_bf_deriv1_func + "\n}\n"

  cartesian_bf_func = cartesian_bf_prototype + "\n"
  for s in car_bf_eval_strs: cartesian_bf_func = cartesian_bf_func + "  " + s + "\n"
  cartesian_bf_func = cartesian_bf_func + "\n}\n"

  cartesian_bf_deriv1_func = cartesian_bf_deriv1_prototype + "\n"
  for s in car_bf_x_eval_strs: cartesian_bf_deriv1_func = cartesian_bf_deriv1_func + "  " + s + "\n"
  cartesian_bf_deriv1_func = cartesian_bf_deriv1_func + "\n"
  for s in car_bf_y_eval_strs: cartesian_bf_deriv1_func = cartesian_bf_deriv1_func + "  " + s + "\n"
  cartesian_bf_deriv1_func = cartesian_bf_deriv1_func + "\n"
  for s in car_bf_z_eval_strs: cartesian_bf_deriv1_func = cartesian_bf_deriv1_func + "  " + s + "\n"
  cartesian_bf_deriv1_func = cartesian_bf_deriv1_func + "\n}\n"


  sqrt_regex = "sqrt\([0-9]+\)"

  sqrt_finds = re.findall( sqrt_regex, spherical_bf_func )
  sqrt_finds = sqrt_finds + ( re.findall( sqrt_regex, spherical_bf_deriv1_func ) )
  sqrt_finds = sqrt_finds + ( re.findall( sqrt_regex, cartesian_bf_func ) )
  sqrt_finds = sqrt_finds + ( re.findall( sqrt_regex, cartesian_bf_deriv1_func ) )

  sqrt_finds = list(set(sqrt_finds))

  for x in sqrt_finds:
    arg = x.strip('sqrt(').strip(')')
    new_str = 'sqrt_' + arg
    spherical_bf_func = spherical_bf_func.replace( x, new_str )
    spherical_bf_deriv1_func = spherical_bf_deriv1_func.replace( x, new_str )
    cartesian_bf_func = cartesian_bf_func.replace( x, new_str )
    cartesian_bf_deriv1_func = cartesian_bf_deriv1_func.replace( x, new_str )

    
    new_str = "constexpr double " + new_str + " = " + str( math.sqrt(int(arg)) ) + ";"
    constant_lines.append( new_str )
    

  cart_header_file.write( cartesian_bf_func )
  cart_header_file.write( cartesian_bf_deriv1_func )
  sphr_header_file.write( spherical_bf_func )
  sphr_header_file.write( spherical_bf_deriv1_func )




# Generate calling routines
cartesian_bf_calling_func = cartesian_bf_template.format('')
spherical_bf_calling_func = spherical_bf_template.format('')
cartesian_bf_deriv1_calling_func = cartesian_bf_deriv1_template.format('')
spherical_bf_deriv1_calling_func = spherical_bf_deriv1_template.format('')

am_dispatch_template = "switch( shell.l ) {{\n"
am_dispatch_template_deriv1 = "switch( shell.l ) {{\n"
for L in range( L_max + 1 ):
  bf_template = """
  case {0}:
    gaueval_{{0}}_angular_{0}(tmp, xc, yc, zc, bf_eval);
    break;
""".format(L)

  deriv1_template = """
  case {0}:
    gaueval_{{0}}_angular_{0}(tmp, xc, yc, zc, bf_eval);
    gaueval_{{0}}_angular_{0}_deriv1(tmp, tmp_x, tmp_y, tmp_z, xc, yc, zc, bf_eval, bf_x_eval, bf_y_eval, bf_z_eval);
    break;
""".format(L)


  am_dispatch_template = am_dispatch_template + bf_template
  am_dispatch_template_deriv1 = am_dispatch_template_deriv1 + deriv1_template
 


am_dispatch_template = am_dispatch_template + "}}\n"
am_dispatch_template_deriv1 = am_dispatch_template_deriv1 + "}}\n"

print(am_dispatch_template_deriv1.format('cartesian'))
print(am_dispatch_template_deriv1.format('spherical'))




footer = "} // namespace gpgaueval"
cart_header_file.write( footer )
sphr_header_file.write( footer )

constant_lines = list(set(constant_lines))
preamble = """
#pragma once

namespace gpgaueval {
"""

cons_header_file.write( preamble )
for s in constant_lines: cons_header_file.write( "  " + s + "\n" )
cons_header_file.write(footer)
