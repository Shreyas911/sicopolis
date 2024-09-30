import os
import sys
import subprocess
import argparse
import numpy as np
import json

def modify_file(filename, ref_string, new_string,
		replace_or_append_or_prepend, 
		instance_number = 0, 
		instance_all = False,
		skip_line = False):

	'''
	Purpose - Modifies file using given options

	Variables - 
	filename - Name of file
	ref_string - Reference string to search within file
	new_string - String to modify ref_string with
	instance_number - Which instance of the ref_string to modify, starting from 0
	instance_all - Modify all instances if True
	skip_line - Skip one line after match with ref_string, only use for append option
	'''
	
	try :
		new_file_lines = []

		with open(filename) as f:

			file_lines = f.readlines()
			current_instance = -1
			skippy = False # A cute name for a temporary boolean
			for line in file_lines:
				
				if(skip_line is False and ref_string in line):
					current_instance = current_instance + 1
					
					if((current_instance == instance_number or instance_all is True) and replace_or_append_or_prepend == 'replace'):
						line = new_string
					elif((current_instance == instance_number or instance_all is True) and replace_or_append_or_prepend == 'append'):
						line = line + new_string
					elif((current_instance == instance_number or instance_all is True) and replace_or_append_or_prepend == 'prepend'):
						line = new_string + line
					else:
						pass

				if(skippy is True):
					current_instance = current_instance + 1
					if(current_instance == instance_number or instance_all is True):
						line = line + new_string
					skippy = False
	
				if(skip_line is True and ref_string in line and replace_or_append_or_prepend == 'append'):
					skippy = True

				new_file_lines.append(line)

		with open(filename, "w") as f:
			f.write(''.join(new_file_lines))

	except FileNotFoundError : 
		print(f'{filename} does not exist')
		sys.exit(1)

	except Exception as err :
		print(f'Some error in modifying file {filename}.')
		print(err)
		sys.exit(1)

	return None

def compile_code(mode, header, domain, 
	clean = True, dep_var=None, ind_vars = None,
	f90c = 'gfortran', cc = 'gcc', tap_adj_prof = 1):


	'''
	Purpose - Compiles code using MakefileTapenade and the given options

	Variables:
	mode - Can be grdchk or adjoint or forward or normal
	clean - If True, runs the terminal command - make -f MakefileTapenade clean
	dep_var - Specify dependent variable (generally fc)
	ind_vars - Specify independent variables (H, vx_c, vy_c, q_geo, etc.)
	cc - C compiler
	f90c - F90 compiler
	tap_adj_prof - If 1, activate profiling for adjoint
	'''

	if (mode == 'adjoint' and (dep_var is not None or ind_vars is not None)):
		pass
	elif (mode == 'forward' and (dep_var is not None or ind_vars is not None)):
		pass
	elif (mode == 'grdchk' or mode == 'normal'):
		pass
	else:
		raise TypeError('Incorrect options specified for compile_code.')

	try :
		
		if(clean is True):
			process = subprocess.run (
				['make', '-f', 'MakefileTapenade', 'clean'])

		process = subprocess.run(
			f'make -f MakefileTapenade '
			f'driver{mode} '
			f'HEADER={header} '
			f'DOMAIN_SHORT={domain} '
			f'DEP_VAR={dep_var} '
			f'IND_VARS={ind_vars} '
			f'CC={cc} '
			f'F90C={f90c} '
			f'TAP_ADJ_PROF={tap_adj_prof} ',
			shell = True, 
			check = True)

	except subprocess.CalledProcessError as error : 

		print(f'Options:')
		print(f'MODE-{mode}')
		print(f'HEADER-{header}')
		print(f'DOMAIN-{domain}')
		print(f'DEP_VAR-{dep_var}')
		print(f'IND_VARS-{ind_vars}')
		print(f'CC-{cc}')
		print(f'F90C-{f90c}')
		print(f'Error with compilation:')
		print(error)
		sys.exit(1)

	return None

def run_executable(mode):
	
	'''
	Purpose - Runs the executable after compilation

	Variables:
	mode - Can be grdchk or adjoint or forward
	'''

	try :

		process = subprocess.run(
			f'./driver{mode}',
			shell = True, 
			check = True)

	except subprocess.CalledProcessError as error : 

		print(f'Options:')
		print(f'MODE-{mode}')
		print(f'Error with running executable:')
		print(error)
		sys.exit(1)

	return None

def get_imax_jmax_kcmax_ktmax(specs_file='sico_specs.h'):
	
	'''
	Purpose - Get IMAX, JMAX, KCMAX, KTMAX from header file

	Variables:
	specs_file - Location of header file
	'''
	try :
	
		with open(specs_file) as f:
			file_lines = f.readlines()
			for line in file_lines:
				if ('#define IMAX' in line):
					IMAX = [int(s) for s in line.split() if s.isdigit()]
				elif ('#define JMAX' in line):
					JMAX = [int(s) for s in line.split() if s.isdigit()]
				elif ('#define KCMAX' in line):
					KCMAX = [int(s) for s in line.split() if s.isdigit()]
				elif ('#define KTMAX' in line):
					KTMAX = [int(s) for s in line.split() if s.isdigit()]
				else:
					pass

		return IMAX[0], JMAX[0], KCMAX[0], KTMAX[0]

	except FileNotFoundError :
		print(f'{specs_file} does not exist')
		sys.exit(1)

	except Exception as err :
		print('Some error in getting IMAX, JMAX, KCMAX, KTMAX.')
		print(err)
		sys.exit(1)

	return None

def copy_file(original_file, destination_file):

	'''
	Purpose - Copy file from source to destination

	Variables -
	original_file - Location of original file
	destination_file - Location where to copy the file
	'''
	try :
		process = subprocess.run (
			['cp', original_file, destination_file])

	except subprocess.CalledProcessError as error :
		print(f"Some issue with copying file {original_file} to {destination_file}.")
		print(error)			
		sys.exit(1)

	return None
	
def setup_grdchk(ind_var, header, domain, 
	dimension = 2,
	z_co_ord = None,
	perturbation = 1.e-3,
	limited_or_block_or_full_or_scalar = 'limited',
	block_imin = None, block_imax = None, block_jmin = None, block_jmax = None,
	grdchk_m_file = 'subroutines/tapenade/grdchk/grdchk_m.F90',
	bool_ad_after_sico_init = True,
	unit = '9999'):

	'''
	Purpose - Sets up everything correctly for compilation of grdchk run

	Variables:
	ind_var - Independent variable
	header - Name of header, for example repo_grl20_ss25ka
	domain - grl or ant
	dimension - dimension of variable being perturbed
	z_co_ord - If 3D variable give z co-ordinate
	perturbation - Amount of perturbation, generally 0.001 or 0.05
	limited_or_block_or_full_or_scalar - If limited, checks only 5 points, 
				   block runs for a block on the grid, 
				   full runs on entire grid,
				   scalar runs only for a single scalar.
	block_imin, block_imax, block_jmin, block_jmax - Specify limits of block if needed
	grdchk_m_file - Location of grdchk_m file
	bool_ad_after_sico_init - Bool that decides if perturbation is before or after sico_init
	unit - The unit to be used in FORTRAN-90 code while opening and closing the output file
	'''

	copy_file(original_file    = '../test_ad/tapenade_F90_templates/grdchk_m.F90',
		  	  destination_file = grdchk_m_file)

	IMAX, JMAX, KCMAX, KTMAX = get_imax_jmax_kcmax_ktmax()

	if(limited_or_block_or_full_or_scalar == 'full'):

		ref_string = '!@ python_automated_grdchk limited_or_block_or_full @'
		new_string = f"""   
		   do i = 0, {IMAX}
		   do j = 0, {JMAX}
		"""
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		ref_string = 'i = ipoints(p)'
		new_string = ''
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	
	
		ref_string = 'j = jpoints(p)'
		new_string = ''
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		ref_string = 'close loop over points'
		new_string = """
		   end do
		"""
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'append',
			   instance_number = 0)	

		ref_string = '!@ python_automated_grdchk scalar_or_array @'
		new_string = 'call grdchk_array'
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)
	
	elif(limited_or_block_or_full_or_scalar == 'block'
	   and float(block_imin).is_integer() 
	   and float(block_imax).is_integer()
	   and float(block_jmin).is_integer() 
	   and float(block_jmax).is_integer()):

		ref_string = '!@ python_automated_grdchk limited_or_block_or_full @'
		new_string = f"""   
		   do i = {block_imin}, {block_imax}
		   do j = {block_jmin}, {block_jmax}
		"""
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		ref_string = 'i = ipoints(p)'
		new_string = ''
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	
	
		ref_string = 'j = jpoints(p)'
		new_string = ''
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		ref_string = 'close loop over points'
		new_string = """
		   end do
		"""
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'append',
			   instance_number = 0)	

		ref_string = '!@ python_automated_grdchk scalar_or_array @'
		new_string = 'call grdchk_array'
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)
	
	elif(limited_or_block_or_full_or_scalar == 'limited'):
		ref_string = '!@ python_automated_grdchk scalar_or_array @'
		new_string = 'call grdchk_array'
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)

	elif(limited_or_block_or_full_or_scalar == 'scalar'):
		ref_string = '!@ python_automated_grdchk scalar_or_array @'
		new_string = 'call grdchk_scalar'
		modify_file(grdchk_m_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)

	else:
		print("Wrong options may have been used in grdchk setup.")
		sys.exit(1)

	if bool_ad_after_sico_init is True:

		ref_string = '!@ python_automated_grdchk_after @'
		if (dimension == 2) :
			new_string = f"""
				orig_val = {ind_var}(j,i)
				if (orig_val .ne. 0) then
					{ind_var}(j,i) = orig_val * perturbation
				else
					{ind_var}(j,i) = perturbation-1
				end if
			"""
		elif(dimension == 3 and z_co_ord is not None and z_co_ord < KCMAX):
			new_string = f"""
				orig_val = {ind_var}({z_co_ord},j,i)
				if (orig_val .ne. 0) then
					{ind_var}({z_co_ord},j,i) = orig_val * perturbation
				else
					{ind_var}({z_co_ord},j,i) = perturbation-1
				end if
			"""
		else:
			raise ValueError ("Something wrong with dimension in grdchk")
			sys.exit(1)

	else:

		ref_string = '!@ python_automated_grdchk_before @'
		if (dimension == 2) :
			new_string = f"""
				orig_val = {ind_var}(j,i)
				if (orig_val .ne. 0) then
					{ind_var}(j,i) = orig_val * perturbation
				else
					{ind_var}(j,i) = perturbation-1
				end if
			"""
		elif(dimension == 3 and z_co_ord is not None and z_co_ord < KCMAX):
			new_string = f"""
				orig_val = {ind_var}({z_co_ord},j,i)
				if (orig_val .ne. 0) then
					{ind_var}({z_co_ord},j,i) = orig_val * perturbation
				else
					{ind_var}({z_co_ord},j,i) = perturbation-1
				end if
			"""
		else:
			raise ValueError ("Something wrong with dimension in grdchk")
			sys.exit(1)

	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
	           instance_number = 0)

	if bool_ad_after_sico_init is True:

		ref_string = '!@ python_automated_grdchk_scalar_after @'
		if (dimension == 2) :
			new_string = f"""
				orig_val = {ind_var}(1,1)
				if (orig_val .ne. 0) then
					{ind_var} = orig_val * perturbation
				else
					{ind_var} = perturbation-1
				end if
			"""
		elif(dimension == 3 and z_co_ord is not None and z_co_ord < KCMAX):
			new_string = f"""
				orig_val = {ind_var}(1,1,1)
				if (orig_val .ne. 0) then
					{ind_var} = orig_val * perturbation
				else
					{ind_var} = perturbation-1
				end if
			"""
		else:
			raise ValueError ("Something wrong with dimension in grdchk")
			sys.exit(1)

	else:

		ref_string = '!@ python_automated_grdchk_scalar_before @'
		if (dimension == 2) :
			new_string = f"""
				orig_val = {ind_var}(1,1)
				if (orig_val .ne. 0) then
					{ind_var} = orig_val * perturbation
				else
					{ind_var} = perturbation-1
				end if
			"""
		elif(dimension == 3 and z_co_ord is not None and z_co_ord < KCMAX):
			new_string = f"""
				orig_val = {ind_var}(1,1,1)
				if (orig_val .ne. 0) then
					{ind_var} = orig_val * perturbation
				else
					{ind_var} = perturbation-1
				end if
			"""
		else:
			raise ValueError ("Something wrong with dimension in grdchk")
			sys.exit(1)


	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
	           instance_number = 0)

	str_ad_after_sico_init = 'after_sico_init' if bool_ad_after_sico_init else 'before_sico_init'
	
	ref_string = '!@ python_automated_grdchk IO begin @'
	new_string = f'''
	   open({unit},&
	   file=\'GradientVals_{ind_var}_{perturbation:.2E}_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat\',&
	   form="FORMATTED", status="REPLACE")
	'''
	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
	           instance_number = 0)	

	ref_string = '!@ python_automated_grdchk_scalar IO begin @'
	new_string = f'''
	   open({unit},&
	   file=\'GradientVals_{ind_var}_{perturbation:.2E}_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat\',&
	   form="FORMATTED", status="REPLACE")
	'''
	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
	           instance_number = 0)	
	
	ref_string = '!@ python_automated_grdchk IO write @'
	new_string = f'          write({unit}, fmt=\'(f40.20)\') gfd\n'
	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
	           instance_number = 0)	

	ref_string = '!@ python_automated_grdchk_scalar IO write @'
	new_string = f'          write({unit}, fmt=\'(f40.20)\') gfd\n'
	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
	           instance_number = 0)	
	
	ref_string = '!@ python_automated_grdchk IO end @'
	new_string = f'   close(unit={unit})\n'
	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
	           instance_number = 0)	

	ref_string = '!@ python_automated_grdchk_scalar IO end @'
	new_string = f'   close(unit={unit})\n'
	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
	           instance_number = 0)	
	
	ref_string = 'real(dp)                          :: orig_val, perturb_val' 
	new_string = f'   real(dp)                          :: orig_val, perturb_val = {perturbation}\n'
	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'replace',
	           instance_number = 0)

	ref_string = 'real(dp)                          :: orig_val, perturb_val' 
	new_string = f'   real(dp)                          :: orig_val, perturb_val = {perturbation}\n'
	modify_file(grdchk_m_file, ref_string, new_string,
		   replace_or_append_or_prepend = 'replace',
	           instance_number = 1)

	return None

def setup_binomial_checkpointing(status = False, number_of_steps = 20, 
	loop_file = 'subroutines/general/sico_main_loop_m.F90'):
 
	'''
	Purpose - Set up sico_main_loop_m.F90 file for binomial checkpointing
	
	Variables: 
	status - If True, setup binomial checkpointing
	number_of_steps - Number of checkpointing steps
	loop_file - Location of sico_main_loop_m.F90 file
	'''

	if status is True:
		modify_file(loop_file, 
		   'BINOMIAL-CKP', 
		   f'  !$AD BINOMIAL-CKP itercount_max+1 {number_of_steps} 1\n',
		   replace_or_append_or_prepend = 'replace',
	           instance_number = 0)	

	elif status is False: 
		modify_file(loop_file, 
		   'BINOMIAL-CKP', 
		   f'  !$NO AD BINOMIAL-CKP itercount_max+1 {number_of_steps} 1\n',
		   replace_or_append_or_prepend = 'replace',
	           instance_number = 0)	

	else: raise ValueError('Incorrect status for checkpointing.')

	return None

def setup_adjoint(ind_vars, header, domain, ckp_status,
	sicopolis_tapenade_cpp_b_file = 'sicopolis_tapenade_cpp_b.f90',
	sico_main_loop_m_cpp_b_file = 'sico_main_loop_m_cpp_b.f90',
	dimensions = [2], 
	z_co_ords = [None], limited_or_block_or_full_or_scalar = 'limited',
	bool_ad_after_sico_init = True,
	output_vars = [], output_iters = [], output_dims = [],
	output_adj_vars = [], output_adj_iters = [], output_adj_dims = []):

	'''
	Purpose - Sets up everything correctly for compilation of adjoint run

	Variables:
	ind_var - Independent variable
	header - Name of header, for example repo_grl20_ss25ka
	domain - grl or ant
	ckp_status - True if binomial checkpointing is used, else False
	sicopolis_tapenade_cpp_b_file - Location of sicopolis_tapenade_cpp_b.f90 file
	sico_main_loop_m_cpp_b_file - Location of sico_main_loop_m_cpp_b.f90 file
	dimensions - List of dimensions of independent vars
	z_co_ords - List of z co-ordinate for 3D vars
	limited_or_block_or_full_or_scalar - If limited, checks only 5 points,
				   block runs for a block on the grid,
				   full runs on entire grid,
				   scalar runs only for a single scalar.
	bool_ad_after_sico_init - Bool that decides if perturbation is before or after sico_init
	output_vars - List of normal variables to output
	output_iters - List of iter number to output the normal variables
	output_dims - List of z co-ordinates for normal variables that we output
	output_adj_vars - List of adjoint variables to output
	output_adj_iters - List of iter number to output the adjoint variables
	output_adj_dims - List of z co-ordinates for adjoint variables that we output
	'''


	if(domain == 'grl' or domain == 'ant'):
		pass
	else:
		raise ValueError('Wrong Domain')

	IMAX, JMAX, KCMAX, KTMAX = get_imax_jmax_kcmax_ktmax()

	ref_string = 'REAL(dp) :: dxi, deta, dzeta_c, dzeta_t, dzeta_r'

	if limited_or_block_or_full_or_scalar != 'scalar':

		if(domain == 'grl'):
			new_string = f'''
			INTEGER(i4b) :: i, j, p
			INTEGER(i4b), parameter :: points = 5
			INTEGER(i4b), DIMENSION(points) :: ipoints, jpoints
			DO p = 1, points
			ipoints(p) = int(real({IMAX}/2))
			jpoints(p) = int(real({JMAX}/5)) + (p-1) * points
			END DO
			'''
		elif(domain == 'ant'):
			new_string = f'''
			INTEGER(i4b) :: i, j, p
			INTEGER(i4b), parameter :: points = 5
			INTEGER(i4b), DIMENSION(points) :: ipoints, jpoints
			DO p = 1, points
			ipoints(p) = int(real({IMAX}/3)) + int(real((.85-.33)*{IMAX}/points)) * (p -1)
			jpoints(p) = int(real({JMAX}/2))
			END DO
			'''

		modify_file(sicopolis_tapenade_cpp_b_file, ref_string, new_string,
			replace_or_append_or_prepend = 'append',
				instance_all = True)
	
		ref_string = 'CALL SICO_INIT_B'
		new_string = ''

		str_ad_after_sico_init = 'after_sico_init' if bool_ad_after_sico_init else 'before_sico_init'
		
		for var_index, (ind_var, dimension, z_co_ord) in enumerate(zip(ind_vars, dimensions, z_co_ords), start = 1):

			unit = [f'{var_index}000', f'{var_index}001']


			if(dimension == 2):
				new_string = new_string + f'''
				open({unit[0]}, file=\'AdjointVals_{ind_var}b_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat\',&
					form="FORMATTED", status="REPLACE")
				open({unit[1]}, file=\'AdjointVals_{ind_var}b_{header}_{str_ad_after_sico_init}.dat\',&
					form="FORMATTED", status="REPLACE")
				do p = 1, points
				i = ipoints(p)
				j = jpoints(p)
				write ({unit[0]}, *) {ind_var}b(j,i)
				end do
				close(unit={unit[0]})
				do i = 0, {IMAX}
				do j = 0, {JMAX}
				write ({unit[1]}, *) {ind_var}b(j,i)
				end do
				end do
				close(unit={unit[1]})
				'''
			elif(dimension == 3 and z_co_ord is not None):
				new_string = new_string + f'''
				open({unit[0]}, file=\'AdjointVals_{ind_var}b_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat\',&
					form="FORMATTED", status="REPLACE")
				open({unit[1]}, file=\'AdjointVals_{ind_var}b_{header}_{str_ad_after_sico_init}.dat\',&
					form="FORMATTED", status="REPLACE")
				do p = 1, points
				i = ipoints(p)
				j = jpoints(p)
				write ({unit[0]}, *) {ind_var}b({z_co_ord},j,i)
				end do
				close(unit={unit[0]})
				do i = 0, {IMAX}
				do j = 0, {JMAX}
				write ({unit[1]}, *) {ind_var}b({z_co_ord},j,i)
				end do
				end do
				close(unit={unit[1]})
				'''
			else :	
				raise ValueError("Wrong dimensions or z coord for adjoint")
				sys.exit(1)

		if bool_ad_after_sico_init is True:
			modify_file(sicopolis_tapenade_cpp_b_file, ref_string, new_string,
				replace_or_append_or_prepend = 'prepend', instance_all = True)
		else:
			modify_file(sicopolis_tapenade_cpp_b_file, ref_string, new_string,
				replace_or_append_or_prepend = 'append', instance_all = True)
	
	else:
	
		ref_string = 'CALL SICO_INIT_B'
		new_string = ''

		str_ad_after_sico_init = 'after_sico_init' if bool_ad_after_sico_init else 'before_sico_init'
		
		for var_index, (ind_var, dimension, z_co_ord) in enumerate(zip(ind_vars, dimensions, z_co_ords), start = 1):

			unit = [f'{var_index}000', f'{var_index}001']


			if(dimension == 2):
				new_string = new_string + f'''
				open({unit[0]}, file=\'AdjointVals_{ind_var}b_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat\',&
					form="FORMATTED", status="REPLACE")
				open({unit[1]}, file=\'AdjointVals_{ind_var}b_{header}_{str_ad_after_sico_init}.dat\',&
					form="FORMATTED", status="REPLACE")
				write ({unit[0]}, *) SUM({ind_var}b)
				close(unit={unit[0]})
				write ({unit[1]}, *) SUM({ind_var}b)
				close(unit={unit[1]})
				'''
			elif(dimension == 3 and z_co_ord is not None):
				new_string = new_string + f'''
				open({unit[0]}, file=\'AdjointVals_{ind_var}b_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat\',&
					form="FORMATTED", status="REPLACE")
				open({unit[1]}, file=\'AdjointVals_{ind_var}b_{header}_{str_ad_after_sico_init}.dat\',&
					form="FORMATTED", status="REPLACE")
				write ({unit[0]}, *) SUM({ind_var}b)
				close(unit={unit[0]})
				write ({unit[1]}, *) SUM({ind_var}b)
				close(unit={unit[1]})
				'''
			else :	
				raise ValueError("Wrong dimensions or z coord for adjoint")
				sys.exit(1)

		if bool_ad_after_sico_init is True:
			modify_file(sicopolis_tapenade_cpp_b_file, ref_string, new_string,
				replace_or_append_or_prepend = 'prepend', instance_all = True)
		else:
			modify_file(sicopolis_tapenade_cpp_b_file, ref_string, new_string,
				replace_or_append_or_prepend = 'append', instance_all = True)

	if (ckp_status is True):
		ref_string = 'itercount = itercount + 1'
	elif (ckp_status is False):
		ref_string = 'END DO main_loop'
	else:
		raise ValueError("Wrong ckp_status for adjoint")
		sys.exit(1)
	
	if (output_vars and output_dims and output_iters):	
		new_string = ''
		for var_index, (output_var, dimension, iteration) in enumerate(zip(output_vars, output_dims, output_iters), start = 1):
			if iteration == -1: iteration = 'itercount_max'
			unit = [f'{var_index}9000']
			
			# 2D array diagnostics
			if dimension >= 0:
				new_string = new_string + f'''
				   if (itercount .EQ. {iteration}) THEN
				   open({unit[0]}, file=\'AdjointVals_{output_var}_{dimension}_iter_{iteration}_{header}_{str_ad_after_sico_init}.dat\',&
				       form="FORMATTED", status="REPLACE")
				   do i = 0, {IMAX}
				   do j = 0, {JMAX}
				   write ({unit[0]}, *) {output_var}({dimension},j,i)
				   end do
				   end do
				   close(unit={unit[0]})
				   end if
				'''

			# 3D array diagnostics
			elif dimension == -1:
				new_string = new_string + f'''
				   if (itercount .EQ. {iteration}) THEN
				      open({unit[0]}, file=\'AdjointVals_{output_var}_{dimension}_iter_{iteration}_{header}_{str_ad_after_sico_init}.dat\',&
				       form="FORMATTED", status="REPLACE")
				   do i = 0, {IMAX}
				   do j = 0, {JMAX}
				   write ({unit[0]}, *) {output_var}(j,i)
				   end do
				   end do
				   close(unit={unit[0]})
				   end if
				'''

			# Scalar diagnostics
			elif dimension == -2:
				new_string = new_string + f'''
				   if (itercount .EQ. {iteration}) THEN
				      open({unit[0]}, file=\'AdjointVals_{output_var}_{dimension}_iter_{iteration}_{header}_{str_ad_after_sico_init}.dat\',&
				       form="FORMATTED", status="REPLACE")
				   write ({unit[0]}, *) SUM({output_var}(j,i))
				   close(unit={unit[0]})
				   end if
				'''
		
			else:
				raise ValueError("Wrong dimensions or z coord for output_var")
				sys.exit(1)

		modify_file(sico_main_loop_m_cpp_b_file, ref_string, new_string, 
		   replace_or_append_or_prepend = 'prepend',
	           instance_number = 0)	

	else: pass

	ref_string = 'BOUNDARY_B'

	if (output_adj_vars and output_adj_dims and output_adj_iters):
		new_string = ''
		for var_index, (output_var, dimension, iteration) in enumerate(zip(output_adj_vars, output_adj_dims, output_adj_iters), start = 1):
			if iteration == -1: iteration = 'itercount_max'
			unit = [f'{var_index}8000']

			# 2D arrayb diagnostics
			if dimension >= 0:
				new_string = new_string + f'''
				   if (itercount .EQ. {iteration}) THEN
				   open({unit[0]}, file=\'AdjointVals_{output_var}b_{dimension}_iter_{iteration}_{header}_{str_ad_after_sico_init}.dat\',&
				       form="FORMATTED", status="REPLACE")
				   do i = 0, {IMAX}
				   do j = 0, {JMAX}
				   write ({unit[0]}, *) {output_var}b({dimension},j,i)
				   end do
				   end do
				   close(unit={unit[0]})
				   end if	
				'''

			# 3D arrayb diagnostics
			elif dimension == -1:
				new_string = new_string + f'''
				   if (itercount .EQ. {iteration}) THEN
				   open({unit[0]}, file=\'AdjointVals_{output_var}b_{dimension}_iter_{iteration}_{header}_{str_ad_after_sico_init}.dat\',&
				       form="FORMATTED", status="REPLACE")
				   do i = 0, {IMAX}
				   do j = 0, {JMAX}
				   write ({unit[0]}, *) {output_var}b(j,i)
				   end do
				   end do
				   close(unit={unit[0]})
				   end if	
				'''

			# scalarb diagnostics
			elif dimension == -2:
				new_string = new_string + f'''
				   if (itercount .EQ. {iteration}) THEN
				   open({unit[0]}, file=\'AdjointVals_{output_var}b_{dimension}_iter_{iteration}_{header}_{str_ad_after_sico_init}.dat\',&
				       form="FORMATTED", status="REPLACE")
				   write ({unit[0]}, *) SUM({output_var}b)
				   close(unit={unit[0]})
				   end if	
				'''
	
			else:
				raise ValueError("Wrong dimensions or z coord for output_varb")
				sys.exit(1)

		modify_file(sico_main_loop_m_cpp_b_file, ref_string, new_string,
	           replace_or_append_or_prepend = 'append',
		   instance_number = 0, skip_line = True)

	return None

def setup_forward(ind_var, header, domain,
	dimension = 2, 
	z_co_ord = None, limited_or_block_or_full_or_scalar = 'limited',
	block_imin = None, block_imax = None, block_jmin = None, block_jmax = None,
	tapenade_main_file = 'tapenade_main.F90',
	bool_ad_after_sico_init = True,
	unit = '99999'):
	

	'''
	Purpose - Setup tlm run for compilation
	
	Variables:
	ind_var - Independent variable
	header - Name of header, for example repo_grl20_ss25ka
	domain - grl or ant
	dimension - dimension of variable being perturbed
	z_co_ord - If 3D variable give z co-ordinate
	limited_or_block_or_full_or_scalar - If limited, checks only 5 points, 
				   block runs for a block on the grid, 
				   full runs on entire grid,
				   scalar runs only for a single scalar.
	block_imin, block_imax, block_jmin, block_jmax - Specify limits of block if needed
	tapenade_main_file - Location of tapenade_main file
	bool_ad_after_sico_init - Bool that decides if perturbation is before or after sico_init
	unit - The unit to be used in FORTRAN-90 code while opening and closing the output file

	'''

	if(domain == 'grl' or domain == 'ant'):
		pass
	else:
		raise ValueError('Wrong Domain')

	copy_file(original_file = '../test_ad/tapenade_F90_templates/tapenade_main.F90',
		  	  destination_file = tapenade_main_file)

	IMAX, JMAX, KCMAX, KTMAX = get_imax_jmax_kcmax_ktmax()

	if limited_or_block_or_full_or_scalar == 'full':
		ref_string = '!@ python_automated_tlm limited_or_block_or_full_or_scalar @'
		new_string = f'''
		   do i = 0, {IMAX}
		   do j = 0, {JMAX}
		'''
		modify_file(tapenade_main_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		modify_file(tapenade_main_file, 
			   'i = ipoints(p)', '',
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	
		
		modify_file(tapenade_main_file, 
			   'j = jpoints(p)', '',
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		modify_file(tapenade_main_file, 
			   'close loop over points', '   end do\n',
			   replace_or_append_or_prepend = 'append',
			   instance_number = 0)	

	elif (limited_or_block_or_full_or_scalar == 'block'
	     and float(block_imin).is_integer() 
	     and float(block_imax).is_integer() 
	     and float(block_jmin).is_integer() 
	     and float(block_jmax).is_integer()):
		ref_string = '!@ python_automated_tlm limited_or_block_or_full_or_scalar @'
		new_string = f'''
		   do i = {block_imin}, {block_imax}
		   do j = {block_jmin}, {block_jmax}
		'''

		modify_file(tapenade_main_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		modify_file(tapenade_main_file, 
			   'i = ipoints(p)', '',
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	
		
		modify_file(tapenade_main_file, 
			   'j = jpoints(p)', '',
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		modify_file(tapenade_main_file, 
			   'close loop over points', '   end do\n',
			   replace_or_append_or_prepend = 'append',
			   instance_number = 0)	

	elif(limited_or_block_or_full_or_scalar == 'limited'):
		pass

	elif(limited_or_block_or_full_or_scalar == 'scalar'):
		ref_string = '!@ python_automated_tlm limited_or_block_or_full_or_scalar @'
		new_string = f'''
		'''
		modify_file(tapenade_main_file, ref_string, new_string,
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		modify_file(tapenade_main_file, 
			   'i = ipoints(p)', '',
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	
		
		modify_file(tapenade_main_file, 
			   'j = jpoints(p)', '',
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

		modify_file(tapenade_main_file, 
			   'close loop over points', '',
			   replace_or_append_or_prepend = 'replace',
			   instance_number = 0)	

	else:
		print("Wrong options may have been used in TLM setup.")
		sys.exit(1)

	ref_string = '!@ python_automated_tlm dep_vard set 0 @'

	if (dimension == 2) :
		new_string = f'''
		            {ind_var}d = 0.0
		'''
	elif (dimension == 3 and z_co_ord is not None and z_co_ord < KCMAX):
		new_string = f'''
		            {ind_var}d = 0.0
		'''
	else:
		raise ValueError ("Something wrong with dimension in TLM")
		sys.exit(1)
	
	modify_file(tapenade_main_file, 
		   ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
		   instance_number = 0)	

	if bool_ad_after_sico_init is True:

		ref_string = '!@ python_automated_tlm_after dep_vard set 1 @'

		if limited_or_block_or_full_or_scalar != 'scalar':
			if (dimension == 2) :
				new_string = f'''
							{ind_var}d(j,i) = 1.0
				'''
			elif (dimension == 3 and z_co_ord is not None and z_co_ord < KCMAX):
				new_string = f'''
							{ind_var}d({z_co_ord},j,i) = 1.0
				'''
			else:
				raise ValueError ("Something wrong with dimension in TLM")
				sys.exit(1)
		else:
				new_string = f'''
							{ind_var}d = 1.0
				'''

	else:

		ref_string = '!@ python_automated_tlm_before dep_vard set 1 @'

		if limited_or_block_or_full_or_scalar != 'scalar':
			if (dimension == 2) :
				new_string = f'''
							{ind_var}d(j,i) = 1.0
				'''
			elif (dimension == 3 and z_co_ord is not None and z_co_ord < KCMAX):
				new_string = f'''
							{ind_var}d({z_co_ord},j,i) = 1.0
				'''
			else:
				raise ValueError ("Something wrong with dimension in TLM")
				sys.exit(1)
		else:
				new_string = f'''
							{ind_var}d = 1.0
				'''
	
	modify_file(tapenade_main_file, 
		   ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
		   instance_number = 0)	

	str_ad_after_sico_init = 'after_sico_init' if bool_ad_after_sico_init else 'before_sico_init'

	ref_string = '!@ python_automated_tlm IO begin @'
	new_string = f'''
	   open({unit}, file=\'ForwardVals_{ind_var}_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat\',&
	       form="FORMATTED", status="REPLACE")
	'''
	modify_file(tapenade_main_file, 
		   ref_string, new_string,
		   replace_or_append_or_prepend = 'append',
		   instance_number = 0)	

	modify_file(tapenade_main_file, 
		   '!@ python_automated_tlm IO write @',
		   f'          write({unit}, fmt=\'(f40.20)\') fcd\n',
		   replace_or_append_or_prepend = 'append',
		   instance_number = 0)	

	modify_file(tapenade_main_file, 
		   '!@ python_automated_tlm IO end @',
		   f'   close(unit={unit})\n',
		   replace_or_append_or_prepend = 'append',
		   instance_number = 0)

	return None

def validate_FD_AD(grdchk_file, ad_file, tolerance = 0.1, low = 0.0):

	'''
	Purpose - Validate adjoint/tlm output with grdchk output

	Variables:
	grdchk_file - Location of grdchk output data
	ad_file - Location of tlm/adjoint output data
	tolerance - Acceptable relative tolerance for discrepancies
	low - abs(values) below maxabs_value*lower_threshold are not checked
	'''

	grdchk_data = np.loadtxt(grdchk_file, dtype = float)
	ad_data = np.loadtxt(ad_file, dtype = float)

	maximum = np.max(np.abs(grdchk_data))
	max_error = np.max(np.abs(ad_data/grdchk_data-1)*(np.abs(grdchk_data)>low*maximum))

	if(max_error >= tolerance):
		raise Exception("Validation failed.")
		sys.exit(1)
	else:
		print("Validated successfully.")

def simulation(mode, header, domain, 
	      ind_var, dep_var,
	      limited_or_block_or_full_or_scalar = 'limited',
	      ind_var_dim = 2, ind_var_z_co_ord = None,
	      perturbation = 1.e-3,
		  grdchk_m_file = 'subroutines/tapenade/grdchk/grdchk_m.F90',
	      tapenade_main_file = 'tapenade_main.F90',
	      run_executable_auto = 'False',
	      bool_ad_after_sico_init = True,
	      unit = '9999',
	      block_imin = None, block_imax = None, block_jmin = None, block_jmax = None,
	      output_vars = [], output_iters = [], output_dims = [],
	      output_adj_vars = [], output_adj_iters = [], output_adj_dims = [],
	      ckp_status = False, ckp_num = 10, validation = True,
	      sicopolis_tapenade_cpp_b_file = 'sicopolis_tapenade_cpp_b.f90',
	      sico_main_loop_m_cpp_b_file = 'sico_main_loop_m_cpp_b.f90',
	      f90c = 'gfortran', cc = 'gcc', tap_adj_prof = 1, low = 0.0, tol = 0.1):

	'''
	Purpose - Sets up everything correctly for compilation of adjoint run

	Variables:
	ind_var - Independent variable
	header - Name of header, for example repo_grl20_ss25ka
	domain - grl or ant
	ind_var_z_co_ord = z_co_ord of independent var / vars (plural for adjoints only in the future)
	ind_var_dim - dimension of independent var / vars (plural for adjoints only in the future)
	perturbation - perturbation for finite difference check
	run_executable_auto - Automatically runs compiled code if True
	unit - File unit number inside FORTRA-90 code
	bool_ad_after_sico_init - Bool that decides if perturbation is before or after sico_init
	output_vars - List of normal variables to output
	output_iters - List of iter number to output the normal variables
	output_dims - List of z co-ordinates for normal variables that we output
	output_adj_vars - List of adjoint variables to output
	output_adj_iters - List of iter number to output the adjoint variables
	output_adj_dims - List of z co-ordinates for adjoint variables that we output
	ckp_status - True if binomial checkpointing is used, else False
	ckp_num - Number of checkpointing steps
	tapenade_main_file - Location of tapenade_main.F90
	sicopolis_tapenade_cpp_b_file - Location of sicopolis_tapenade_cpp_b.f90 file
	sico_main_loop_m_cpp_b_file - Location of sico_main_loop_m_cpp_b.f90 file
	cc - C compiler
	f90c - F90 compiler
	tap_adj_prof - If 1, activate profiling for adjoint
	low - abs(values) below maxabs_value*lower_threshold are not checked
    tol - How much relative error do we tolerate?
	'''

	try:
	
		os.chdir("../src/")
		print("Directory changed to ", os.getcwd())
	
	except OSError:
		print("Can't change the Current Working Directory")
		sys.exit(1)

	copy_file(f'../headers/sico_specs_{header}.h', 'sico_specs.h')

	if ind_var_dim == 3 and ind_var_z_co_ord is None:
		raise Exception ("Wrong input arguments for dimension")
		sys.exit(1)


	if mode == 'grdchk':
	
		setup_grdchk(ind_var = ind_var, header = header, domain = domain,
	        dimension = ind_var_dim,
	        z_co_ord = ind_var_z_co_ord,
	        perturbation = perturbation,
	        limited_or_block_or_full_or_scalar = limited_or_block_or_full_or_scalar,
	        block_imin = block_imin, block_imax = block_imax, block_jmin = block_jmin, block_jmax = block_jmax,
	        grdchk_m_file = grdchk_m_file,
	        bool_ad_after_sico_init = bool_ad_after_sico_init,
	        unit = unit)
	
		compile_code(mode = mode, header = header, domain = domain,
			    clean = True, dep_var=dep_var, ind_vars = ind_var,
			    f90c = f90c, cc = cc, tap_adj_prof = tap_adj_prof)
	
		print(f'grdchk compilation complete for {header}.')
	
		if run_executable_auto is True:
			run_executable('grdchk')
			print(f'grdchk execution complete for {header}.')
	
	elif mode == 'adjoint':

		if ckp_status is True:	
			setup_binomial_checkpointing(status = True, number_of_steps = ckp_num)
		else:
			setup_binomial_checkpointing(status = False)

		copy_file(original_file = '../test_ad/tapenade_F90_templates/tapenade_main.F90',
		  		  destination_file = 'tapenade_main.F90')

		compile_code(mode = mode, header = header, domain = domain,
			    clean = True, dep_var=dep_var, ind_vars = ind_var,
			    f90c = f90c, cc = cc, tap_adj_prof = tap_adj_prof)
	
		if output_dims is not None and output_iters is not None:
			output_dims = [int(x) for x in output_dims]	
			output_iters = [int(x) for x in output_iters]
		if output_adj_dims is not None and output_adj_iters is not None:
			output_adj_dims = [int(x) for x in output_adj_dims]	
			output_adj_iters = [int(x) for x in output_adj_iters]

		setup_adjoint(ind_vars = [ind_var], header = header, domain = domain, ckp_status = ckp_status,
			     sicopolis_tapenade_cpp_b_file = sicopolis_tapenade_cpp_b_file,
			     sico_main_loop_m_cpp_b_file = sico_main_loop_m_cpp_b_file,
			     dimensions = [ind_var_dim],
			     z_co_ords = [ind_var_z_co_ord], limited_or_block_or_full_or_scalar = limited_or_block_or_full_or_scalar,
			     bool_ad_after_sico_init = bool_ad_after_sico_init,
			     output_vars = output_vars, output_iters = output_iters, output_dims = output_dims,
			     output_adj_vars = output_adj_vars, output_adj_iters = output_adj_iters, output_adj_dims = output_adj_dims)

		compile_code(mode = mode, header = header, domain = domain,
			    clean = False, dep_var=dep_var, ind_vars = ind_var,
			    f90c = f90c, cc = cc, tap_adj_prof = tap_adj_prof)
		
		print(f'adjoint compilation complete for {header}.')
		
		if run_executable_auto is True:
			run_executable('adjoint')
			print(f'adjoint execution complete for {header}.')

		if validation is True:
			pert = f'{perturbation:.2e}'.upper()
			str_ad_after_sico_init = 'after_sico_init' if bool_ad_after_sico_init else 'before_sico_init'
			grdchk_file = f'GradientVals_{ind_var}_{pert}_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat'
			adjoint_file = f'AdjointVals_{ind_var}b_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat'

			if os.path.exists(grdchk_file) is False:
				raise FileNotFoundError (f'{grdchk_file} not found for validation')
			if os.path.exists(adjoint_file) is False:
				raise FileNotFoundError (f'{adjoint_file} not found for validation')

			validate_FD_AD(grdchk_file, adjoint_file, tolerance = tol, low = low)

	elif mode == 'forward':

		setup_forward(ind_var = ind_var, header = header, domain = domain,
			     dimension = ind_var_dim,
			     z_co_ord = ind_var_z_co_ord, limited_or_block_or_full_or_scalar = limited_or_block_or_full_or_scalar,
			     block_imin = block_imin, block_imax = block_imax, block_jmin = block_jmin, block_jmax = block_jmax,
			     tapenade_main_file = tapenade_main_file,
			     bool_ad_after_sico_init = bool_ad_after_sico_init,
			     unit = unit)

		compile_code(mode = mode, header = header, domain = domain,
			     clean = True, dep_var=dep_var, ind_vars = ind_var,
			     f90c = f90c, cc = cc, tap_adj_prof = tap_adj_prof)

		print(f'TLM compilation complete for {header}.')

		if run_executable_auto is True:
			run_executable('forward')
			print(f'TLM execution complete for {header}.')

		if validation is True:
			pert = f'{perturbation:.2e}'.upper()
			str_ad_after_sico_init = 'after_sico_init' if bool_ad_after_sico_init else 'before_sico_init'
			grdchk_file = f'GradientVals_{ind_var}_{pert}_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat'
			tlm_file = f'ForwardVals_{ind_var}_{header}_{limited_or_block_or_full_or_scalar}_{str_ad_after_sico_init}.dat'

			if os.path.exists(grdchk_file) is False:
				raise FileNotFoundError (f'{grdchk_file} not found for validation')
			if os.path.exists(tlm_file) is False:
				raise FileNotFoundError (f'{tlm_file} not found for validation')

			validate_FD_AD(grdchk_file, tlm_file, tolerance = tol, low = low)
	
	elif mode == 'normal':

		### The name here should match the output directory in MakefileTapenade

		process = subprocess.run (
			     ['rm', '-r', f'../sico_out/N_{header}'])
		compile_code(mode = mode, header = header, domain = domain,
			     clean = True, dep_var=dep_var, ind_vars = ind_var,
			     f90c = f90c, cc = cc, tap_adj_prof = tap_adj_prof)
		
		print(f'normal compilation complete for {header}.')
		
		if run_executable_auto is True:
			run_executable('normal')
			print(f'normal execution complete for {header}.')

	else:
		raise ValueError("Incorrect simulation mode")
		sys.exit(1)

	return None

if __name__ == "__main__":

	### Command line arguments overwrite json arguments

	# argparse.SUPPRESS ensures that empty arguments are not included in dict, and thus do not overwrite json file
	parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)

	parser.add_argument("-jsf", "--json", help="name of json data file", type=str)
	parser.add_argument("-head", "--header", help="name of header file", type=str)
	parser.add_argument("-dom", "--domain", help="short name of domain, either grl or ant", type = str)
	parser.add_argument("-dv", "--dep_var", help="name of dependent variable", type=str)
	parser.add_argument("-iv", "--ind_var", help="name of independent variable", type=str)
	parser.add_argument("-delta", "--perturbation", help="value of perturbation for grdchk", type=float)
	parser.add_argument("-ckp", "--checkpoint", help="number of steps in checkpointing", type=int)
	parser.add_argument("-dim", "--dimension", help="2D or 3D independent variable, default 2D", type=int)
	parser.add_argument("-z", "--z_co_ord", help="z co-ordinate if 3D variable", type=int)
	parser.add_argument('-ov','--output_vars', nargs='+', help='List the fields you want to output')
	parser.add_argument('-od', '--output_dims', nargs='+', help='List the z-coord of output vars, -1 if 2D')
	parser.add_argument('-oi', '--output_iters', nargs='+', help='List the iter num of output vars, -1 if itercount_max')
	parser.add_argument('-oav','--output_adj_vars', nargs='+', help='List the adjoint fields you want to output')
	parser.add_argument('-oad', '--output_adj_dims', nargs='+', help='List the z-coord of adjoint output vars, -1 if 2D')
	parser.add_argument('-oai', '--output_adj_iters', nargs='+', help='List the iter num of adjoint output vars, -1 if itercount_max')
	parser.add_argument("-run", help="Run the executables?", action="store_true")
	parser.add_argument("-cc", '--c_compiler', help="C compiler", type=str)
	parser.add_argument("-f90c", '--f90_compiler', help="F90 compiler", type=str)
	parser.add_argument("-prof", '--adj_prof', help="If 1, activate adjoint profiling", type=str)
	parser.add_argument("-lbfs", '--limited_or_block_or_full_or_scalar', help="limited or block or full or scalar?", type=str)
	parser.add_argument("-low", '--low', help="abs(values) below maxabs_value*lower_threshold are not checked", type=str)
	parser.add_argument("-asi", '--after_sico_init', help="If 1, activation of AD after sico_init", type=str)
    parser.add_argument("-tol", '--tolerance', help="maxabs relative error acceptable", type=str)

	args = parser.parse_args()

	### hasattr only correctly works for Python 3
	if hasattr(args, 'json'):
		try:
			with open(args.json, 'r') as j:
				json_dict = json.loads(j.read())
			
				# Convert Namespace to dict
				args_dict = vars(args)
				
				# Update json dict with command line arguments
				json_dict.update(args_dict)
			
				# Convert updated json dict back to args Namespace
				args = argparse.Namespace(**json_dict)
		except:
			raise Exception('Some issue with reading json file')
			sys.exit(1)

	ckp_status = False
	if hasattr(args, 'checkpoint'):
		ckp_status = True

	list_attrs = ['json', 'header', 'domain', 'ind_var', 'dep_var',
	              'dimension', 'z_co_ord', 'perturbation', 'output_vars',
	              'output_iters', 'output_dims', 'output_adj_vars',
	              'output_adj_iters', 'output_adj_dims', 'checkpoint', 'run',
	              'c_compiler', 'f90_compiler', 'adj_prof',
	              'limited_or_block_or_full_or_scalar',
	              'low', 'after_sico_init', 'tolerance']

	for attr in list_attrs:
		if not hasattr(args, attr):
			setattr(args, attr, None)

	# Default is 'limited' i.e. checks on 5 points
	if args.limited_or_block_or_full_or_scalar is None:
		args.limited_or_block_or_full_or_scalar = 'limited'

	if int(args.after_sico_init) == 1:
		args.bool_ad_after_sico_init = True
	else:
		args.bool_ad_after_sico_init = False

	for mode in ['normal', 'grdchk', 'adjoint', 'forward']:

		simulation(mode = mode, header = args.header, domain = args.domain, 
				ind_var = args.ind_var, dep_var = args.dep_var,
				limited_or_block_or_full_or_scalar = args.limited_or_block_or_full_or_scalar,
				ind_var_dim = args.dimension, ind_var_z_co_ord = args.z_co_ord,
				perturbation = args.perturbation,
				run_executable_auto = args.run, bool_ad_after_sico_init = args.bool_ad_after_sico_init,
				output_vars = args.output_vars, output_iters = args.output_iters, output_dims = args.output_dims,
				output_adj_vars = args.output_adj_vars, output_adj_iters = args.output_adj_iters, output_adj_dims = args.output_adj_dims,
				ckp_status = ckp_status, ckp_num = args.checkpoint,
				cc = args.c_compiler, f90c = args.f90_compiler, tap_adj_prof = args.adj_prof, low = float(args.low))
