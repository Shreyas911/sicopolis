import os
import sys
import subprocess
import argparse
import numpy as np

def compile_code(mode, header, domain, 
	clean = True,
	travis_ci ='', dep_var=None, ind_vars = None):

	if (mode == 'adjoint' and (dep_var is None or ind_vars is None)):
		raise TypeError('Adjoint mode but dep_var or ind_vars not specified.')

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
			f'{travis_ci}', 
			shell = True, 
			check = True)

	except subprocess.CalledProcessError as error : 

		print(f'Options:')
		print(f'MODE-{mode}')
		print(f'HEADER-{header}')
		print(f'DOMAIN-{domain}')
		print(f'DEP_VAR-{dep_var}')
		print(f'IND_VARS-{ind_vars}')
		print(f'Error with compilation:')
		print(error)
		sys.exit(1)

def run_executable(mode):

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

def get_imax_jmax(specs_file='sico_specs.h'):
	
	try :
	
		with open(specs_file) as f:
			file_lines = f.readlines()
			for line in file_lines:
				if ('#define IMAX' in line):
					IMAX = [int(s) for s in line.split() if s.isdigit()]
				if ('#define JMAX' in line):
					JMAX = [int(s) for s in line.split() if s.isdigit()]
		return IMAX[0], JMAX[0]

	except FileNotFoundError :
		print(f'{specs_file} does not exist')
		sys.exit(1)

	except Exception as err :
		print('Some error in getting IMAX and JMAX.')
		print(err)
		sys.exit(1)

def get_imax_jmax_kcmax_ktmax(specs_file='sico_specs.h'):
	
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

def copy_tapenade_m_template(template_file = '../test_ad/tapenade_m_adjoint_template.F90', 
			destination_file = 'subroutines/tapenade/tapenade_m.F90'):

	try :
		process = subprocess.run (
			['cp', template_file, destination_file])

	except subprocess.CalledProcessError as error :
		print("Some issue with copying template tapenade_m file")
		print(error)			
		sys.exit(1)
	
def setup_grdchk(ind_var, header, domain, 
	dimension = 2,
	z_co_ord = None,
	perturbation = 1.e-3,
	limited_or_full = 'limited',
	tapenade_m_file = 'subroutines/tapenade/tapenade_m.F90',
	unit = '9999'):

	copy_tapenade_m_template()

	IMAX, JMAX, KCMAX, KTMAX = get_imax_jmax_kcmax_ktmax()

	try :

		new_file_lines = []
		with open(tapenade_m_file) as f:
			file_lines = f.readlines()
	
			for line in file_lines:

				if(limited_or_full == 'full'):
					if('!@ python_automated_grdchk limited_or_full @' in line) :
						line = f'   do i = 0, {IMAX}\n' \
							+ f'   do j = 0, {JMAX}\n'

					if('i = ipoints(p)' in line) :
						line = ''
					
					if('j = jpoints(p)' in line) :
						line = ''	
			
					if('close loop over points' in line) : 
						line = line + '   end do\n'
					
				if ('!@ python_automated_grdchk @' in line):
					
					if (dimension == 2) :
						line = line \
							+ f'            orig_val = {ind_var}(j,i)\n' \
							+ f'            {ind_var}(j,i) = orig_val * perturbation\n'
					elif(dimension == 3 and z_co_ord is not None and z_co_ord < KCMAX):
						line = line \
							+ f'            orig_val = {ind_var}({z_co_ord},j,i)\n' \
							+ f'            {ind_var}({z_co_ord},j,i) = orig_val * perturbation\n'

					else:
						raise ValueError ("Something wrong with dimension in grdchk")
						sys.exit(1)
						
				if ('!@ python_automated_grdchk IO begin @' in line):
					line = line \
						+ f'   open({unit}, ' \
						+ f'file=\'GradientVals_{ind_var}_{perturbation:.2E}_\'//trim(RUNNAME)//\'_{limited_or_full}.dat\',&' \
						+ f'\n       form="FORMATTED", status="REPLACE")' 

				if ('!@ python_automated_grdchk IO write @' in line):
					line = line + f'          write({unit}, fmt=\'(f40.20)\') gfd\n'

				if ('!@ python_automated_grdchk IO end @' in line):
					line = line + f'   close(unit={unit})\n'

				if ('real(dp)                          :: orig_val, perturb_val' in line):
					line = f'   real(dp)                          :: orig_val, perturb_val = {perturbation}\n'

				new_file_lines.append(line)
	
		with open(tapenade_m_file, "w") as f:
			f.write(''.join(new_file_lines))

	except FileNotFoundError :
		print(f'{tapenade_m_file} not found.')
		sys.exit(1)

	except Exception as err:
		print("Some problem with grdchk setup.")	
		print(err)
		sys.exit(1)

def setup_binomial_checkpointing(status = False, number_of_steps = 20, 
	loop_file = 'subroutines/general/sico_main_loop_m.F90'):
 
	try:

		new_file_lines = []

		with open(loop_file, "r") as f:	

			file_lines = f.readlines()

			for line in file_lines:
				if status is True:
					if('BINOMIAL-CKP' in line):
						line = f'  !$AD BINOMIAL-CKP itercount_max+1 {number_of_steps} 1\n'
				elif status is False:
					if('BINOMIAL-CKP' in line):
						line = f'  !$NO AD BINOMIAL-CKP itercount_max+1 {number_of_steps} 1\n'
					
				else:
					raise ValueError('Incorrect status for checkpointing.')
		
				new_file_lines.append(line)
			
		with open(loop_file, "w") as f:
			f.write(''.join(new_file_lines))

	except FileNotFoundError :
		print(f'{loop_file} not found.')
		sys.exit(1)
	
	except Exception as err:
		print("Some problem with binomial checkpointing setup.")	
		print(err)
		sys.exit(1)

def setup_adjoint(ind_vars, header, domain, 
	numCore_cpp_b_file = 'numCore_cpp_b.f90',
	dimensions = [2], 
	z_co_ords = [None]):
	
	if(domain == 'grl' or domain == 'ant'):
		pass
	else:
		raise ValueError('Wrong Domain')

	IMAX, JMAX, KCMAX, KTMAX = get_imax_jmax_kcmax_ktmax()
	try :
	
		with open('numCore_cpp_b.f90') as f:
			
			file_lines = f.readlines()
			new_file_lines = []
			
			for line in file_lines:
				if ('CHARACTER(len=100) :: runname' in line):
					line = line \
						+ f'   INTEGER(i4b) :: i, j, p\n' \
						+ f'   INTEGER(i4b), parameter :: points = 5\n' \
						+ f'   INTEGER(i4b), DIMENSION(points) :: ipoints, jpoints\n' \
						+ f'   DO p = 1, points\n'
							
					if(domain == 'grl'):
						line = line \
							+ f'   ipoints(p) = int(real({IMAX}/2))\n' \
							+ f'   jpoints(p) = int(real({JMAX}/5)) + (p-1) * points\n'
		
					elif(domain == 'ant'):
						line = line \
							+ f'   ipoints(p) = int(real({IMAX}/3)) + int(real((.85-.33)*{IMAX}/points)) * (p -1)\n' \
							+ f'   jpoints(p) = int(real({JMAX}/2))\n'
					
					line = line + f'   END DO\n'
												
				if ('CALL SICO_INIT_B' in line):
					for var_index, (ind_var, dimension, z_co_ord) in enumerate(zip(ind_vars, dimensions, z_co_ords), start = 1):

						unit = [f'{var_index}000', f'{var_index}001']

						line_add = f'   open({unit[0]}, file=\'AdjointVals_{ind_var}_\'//trim(RUNNAME)//\'_limited.dat\',&\n' \
							+ f'       form="FORMATTED", status="REPLACE")\n' \
							+ f'   open({unit[1]}, file=\'AdjointVals_{ind_var}_\'//trim(RUNNAME)//\'.dat\',&\n' \
							+ f'       form="FORMATTED", status="REPLACE")\n' \
							+ f'   do p = 1, points\n' \
							+ f'   i = ipoints(p)\n' \
							+ f'   j = jpoints(p)\n'
						if (dimension == 2) :
							line_add = line_add + f'   write ({unit[0]}, *) {ind_var}b(j,i)\n'
						elif (dimension == 3 and z_co_ord is not None) :
							line_add = line_add + f'   write ({unit[0]}, *) {ind_var}b({z_co_ord},j,i)\n'
						else :
							raise ValueError("Wrong dimensions or z coord for adjoint")
							sys.exit(1)

						line_add = line_add \
							+ f'   end do\n' \
							+ f'   close(unit={unit[0]})\n' \
							+ f'   do i = 0, {IMAX}\n' \
							+ f'   do j = 0, {JMAX}\n'
						if (dimension == 2) :
							line_add = line_add + f'   write ({unit[1]}, *) {ind_var}b(j,i)\n'
						elif (dimension == 3 and z_co_ord is not None) :
							line_add = line_add + f'   write ({unit[1]}, *) {ind_var}b({z_co_ord},j,i)\n'
						else:
							raise ValueError("Wrong dimensions or z coord for adjoint")
							sys.exit(1)

						line_add = line_add \
							+ f'   end do\n' \
							+ f'   end do\n' \
							+ f'   close(unit={unit[1]})\n'

						line = line_add + line

				new_file_lines.append(line)
		
		with open('numCore_cpp_b.f90', "w") as f:
			f.write(''.join(new_file_lines))

	except FileNotFoundError :
		print(f'{numCore_cpp_b_file} not found.')
		sys.exit(1)
	except Exception as err:
		print("Some problem with adjoint setup.")
		print(err)
		sys.exit(1)

def validate_FD_AD(grdchk_file, ad_file, tolerance = 0.1):
	grdchk_data = np.loadtxt(grdchk_file, dtype = float)
	ad_data = np.loadtxt(ad_file, dtype = float)
	
	if(np.max(np.abs(ad_data/grdchk_data-1)) >= tolerance):
		raise Exception("Validation failed.")
		sys.exit(1)
	else:
		print("Validated successfully.")

if __name__ == "__main__":

	try:
	
		# Change the current working Directory
		os.chdir("../src/")
		print("Directory changed to ", os.getcwd())
	
	except OSError:
		print("Can't change the Current Working Directory")
		sys.exit(1)
	
	parser = argparse.ArgumentParser()

	parser.add_argument("-head", "--header", help="name of header file", type=str, required=True)
	parser.add_argument("-dom", "--domain", help="short name of domain, either grl or ant", type = str, required=True)
	parser.add_argument("-dv", "--dep_var", help="name of dependent variable", type=str, required=True)
	parser.add_argument("-iv", "--ind_var", help="name of independent variable", type=str, required=True)
	parser.add_argument("-delta", "--perturbation", help="value of perturbation for grdchk", type=float, required=True)
	parser.add_argument("-ckp", "--checkpoint", help="number of steps in checkpointing", type=int)
	parser.add_argument("--travis", help="travis setup", action="store_true")
	parser.add_argument("-dim", "--dimension", help="2D or 3D variable, default 2D", type=int)
	parser.add_argument("-z", "--z_co_ord", help="z co-ordinate if 3D variable", type=int)

	args = parser.parse_args()
	
	if args.dimension == 3 and args.z_co_ord is None:
		raise Exception ("Wrong input arguments for dimension")
		sys.exit(1)

	kwargs = dict(ind_var=args.ind_var, header=args.header, domain=args.domain,
		     dimension=args.dimension, perturbation=args.perturbation,z_co_ord=args.z_co_ord)

	setup_grdchk(**{k: v for k, v in kwargs.items() if v is not None},
        limited_or_full = 'limited',
        tapenade_m_file = 'subroutines/tapenade/tapenade_m.F90',
        unit = '9999')

	if args.travis:
		kwargs = dict(mode='grdchk', header=args.header, domain=args.domain, travis_ci='TRAVIS_CI=yes')
	else:
		kwargs = dict(mode='grdchk', header=args.header, domain=args.domain)

	compile_code(**{k: v for k, v in kwargs.items() if v is not None}, clean = True)	
	print(f'grdchk execution complete for {args.header}.')

	run_executable('grdchk')
	print(f'grdchk execution complete for {args.header}.')

	if args.checkpoint:	
		setup_binomial_checkpointing(status = True, number_of_steps = args.checkpoint)
	else: 
		setup_binomial_checkpointing(status = False)

	if args.travis:
		kwargs = dict(mode='adjoint', header=args.header, domain=args.domain, 
			     dep_var = args.dep_var, ind_vars = args.ind_var, travis_ci='TRAVIS_CI=yes')
	else:
		kwargs = dict(mode='adjoint', header=args.header, domain=args.domain,
			     dep_var = args.dep_var, ind_vars = args.ind_var)

	compile_code(**{k: v for k, v in kwargs.items() if v is not None}, clean = True)	

	kwargs = dict(ind_vars = [args.ind_var], header=args.header, domain=args.domain,
		     dimensions = [args.dimension], z_co_ords = [args.z_co_ord])

	setup_adjoint(**{k: v for k, v in kwargs.items() if (v is not None and v != [None])},
        	     numCore_cpp_b_file = 'numCore_cpp_b.f90')

	if args.travis:
		kwargs = dict(mode='adjoint', header=args.header, domain=args.domain, 
			     dep_var = args.dep_var, ind_vars = args.ind_var, travis_ci='TRAVIS_CI=yes')
	else:
		kwargs = dict(mode='adjoint', header=args.header, domain=args.domain,
			     dep_var = args.dep_var, ind_vars = args.ind_var)

	compile_code(**{k: v for k, v in kwargs.items() if v is not None}, clean = False)
	
	print(f'adjoint compilation complete for {args.header}.')
	
	run_executable('adjoint')
	print(f'adjoint execution complete for {args.header}.')

	validate_FD_AD(f'GradientVals_{args.ind_var}_{args.perturbation:.2E}_{args.header}_limited.dat', f'AdjointVals_{args.ind_var}_{args.header}_limited.dat')
