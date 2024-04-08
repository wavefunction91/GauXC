import os,sys
import pyexpander.lib as expander

from io import StringIO

L_max = 6

if len(sys.argv) > 1: L_max = int(sys.argv[1])

L_var_dict = { 'L_max' : L_max }

old_sys_out = sys.stdout

sys.stdout = col_device_expand = StringIO()
expander.expandFile( 'templates/collocation_device_template.cu',
  external_definitions=L_var_dict )
sys.stdout = old_sys_out

sys.stdout = col_shell_to_task_kernels_expand = StringIO()
expander.expandFile( 'templates/collocation_shell_to_task_kernels_template.hpp',
  external_definitions=L_var_dict )
sys.stdout = old_sys_out

col_device_expand = col_device_expand.getvalue()
col_shell_to_task_kernels_expand = col_shell_to_task_kernels_expand.getvalue()

col_device_fname = "../collocation_device.cu"
col_shell_to_task_kernels_fname = "../collocation_shell_to_task_kernels.hpp"

col_device_file = open( col_device_fname, 'w')
col_device_file.write( col_device_expand )

col_shell_to_task_kernels_file = open( col_shell_to_task_kernels_fname, 'w')
col_shell_to_task_kernels_file.write( col_shell_to_task_kernels_expand )



