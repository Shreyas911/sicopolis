import os
import sys
import subprocess

def compile_code(mode, header, domain, 
	clean = True,
	lisdir = '/home/shreyas/lis-1.4.43/installation', 
	netcdf_fortran_dir = '/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4', 
	travis_ci ='', dep_var=None, ind_vars = None):

	if (mode == 'adjoint' and (dep_var is None or ind_vars is None)):
		raise TypeError('Adjoint mode but dep_var or ind_vars not specified.')

	try :
		
		if(clean is True):
			process = subprocess.run (
				['make', '-f', 'MakefileTapenade', 'clean'], 
				capture_output = True)

		process = subprocess.run(
			f'make -f MakefileTapenade '
			f'driver{mode} '
			f'HEADER={header} '
			f'DOMAIN_SHORT={domain} '
			f'LISDIR={lisdir} '
			f'NETCDF_FORTRAN_DIR={netcdf_fortran_dir} '
			f'DEP_VAR={dep_var} '
			f'IND_VARS={ind_vars} '
			f'{travis_ci}', 
			capture_output = True, 
			shell = True, 
			check = True)

	except subprocess.CalledProcessError as error : 

		print(f'Options:')
		print(f'MODE-{mode}')
		print(f'HEADER-{header}')
		print(f'DOMAIN-{domain}')
		print(f'LISDIR-{lisdir}')
		print(f'NETCDF_DIR-{netcdf_fortran_dir}')
		print(f'DEP_VAR-{dep_var}')
		print(f'IND_VARS-{ind_vars}')
		print(f'Error with compilation:')
		print(error)

def run_executable(mode):

	try :

		process = subprocess.run(
			f'./driver{mode}',
			capture_output = True, 
			shell = True, 
			check = True)

	except subprocess.CalledProcessError as error : 

		print(f'Options:')
		print(f'MODE-{mode}')
		print(f'Error with running executable:')
		print(error)

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

	except :
		print('Some error in getting IMAX and JMAX.')

def copy_tapenade_m_template(template_file = '../test_ad/tapenade_m_adjoint_template.F90', 
			destination_file = 'subroutines/tapenade/tapenade_m.F90'):

	try :
		process = subprocess.run (
			['cp', template_file, destination_file],
			capture_output = True)

	except subprocess.CalledProcessError as error :
		print("Some issue with copying template tapenade_m file")
		print(error)			
	
def setup_grdchk(ind_var, header, domain, 
	tapenade_m_file = 'subroutines/tapenade/tapenade_m.F90',
	lisdir = '/home/shreyas/lis-1.4.43/installation', 
	netcdf_fortran_dir = '/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4', 
	travis_ci ='', unit = '9999'):

	copy_tapenade_m_template()

	try :

		new_file_lines = []
		with open(tapenade_m_file) as f:
			file_lines = f.readlines()
	
			for line in file_lines:
	
				if ('!@ python_automated_grdchk @' in line):
					line = line \
						+ f'            orig_val = {ind_var}(j,i)\n' \
						+ f'            {ind_var}(j,i) = orig_val * perturbation\n'
				if ('!@ python_automated_grdchk IO begin @' in line):
					line = line \
						+ f'   open({unit}, ' \
						+ f'file=\'GradientVals_{ind_var}_\'//trim(RUNNAME)//\'.dat\',&' \
						+ f'\n       form="FORMATTED", status="REPLACE")' 
				if ('!@ python_automated_grdchk IO write @' in line):
					line = line + f'          write({unit}, fmt=\'(f40.20)\') gfd\n'
				if ('!@ python_automated_grdchk IO end @' in line):
					line = line + f'   close(unit={unit})\n'

				new_file_lines.append(line)
	
		with open(tapenade_m_file, "w") as f:
			f.write(''.join(new_file_lines))

	except FileNotFoundError :
		print(f'{tapenade_m_file} not found.')
	
	except :
		print("Some problem with grdchk setup.")	

def setup_adjoint(ind_vars, header, domain, 
	numCore_cpp_b_file = 'numCore_cpp_b.f90',
	lisdir = '/home/shreyas/lis-1.4.43/installation', 
	netcdf_fortran_dir = '/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4', 
	travis_ci =''):
	
	if(domain == 'grl' or domain == 'ant'):
		pass
	else:
		raise ValueError('Wrong Domain')

	IMAX, JMAX = get_imax_jmax()

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

					for var_index, ind_var in enumerate(ind_vars, start = 1):

						unit = [f'{var_index}000', f'{var_index}001']

						line = f'   open({unit[0]}, file=\'AdjointVals_{ind_var}_\'//trim(RUNNAME)//\'_limited.dat\',&\n' \
							+ f'       form="FORMATTED", status="REPLACE")\n' \
							+ f'   open({unit[1]}, file=\'AdjointVals_{ind_var}_\'//trim(RUNNAME)//\'.dat\',&\n' \
							+ f'       form="FORMATTED", status="REPLACE")\n' \
							+ f'   do p = 1, points\n' \
							+ f'   i = ipoints(p)\n' \
							+ f'   j = jpoints(p)\n' \
							+ f'   write ({unit[0]}, *) {ind_var}b(j,i)\n' \
							+ f'   end do\n' \
							+ f'   close(unit={unit[0]})\n' \
							+ f'   do i = 0, {IMAX}\n' \
							+ f'   do j = 0, {JMAX}\n' \
							+ f'   write ({unit[1]}, *) {ind_var}b(j,i)\n' \
							+ f'   end do\n' \
							+ f'   end do\n' \
							+ f'   close(unit={unit[1]})\n' \
							+ line

				new_file_lines.append(line)
		
		with open('numCore_cpp_b.f90', "w") as f:
			f.write(''.join(new_file_lines))

	except FileNotFoundError :
		print(f'{numCore_cpp_b_file} not found.')

	except :
		print("Some problem with adjoint setup.")

if __name__ == "__main__":

	if (len(sys.argv) != 5) :
		raise Exception('Insufficient number of command line arguments.')


	try:
	
		# Change the current working Directory
		os.chdir("../src/")
		print("Directory changed to ", os.getcwd())
	
	except OSError:

		print("Can't change the Current Working Directory")

	header_file = sys.argv[1]
	domain = sys.argv[2]
	dep_var = sys.argv[3]
	ind_var = sys.argv[4]

	setup_grdchk(ind_var, header_file, domain)
	compile_code('grdchk', header_file, domain)
	print(f'grdchk compilation complete for {header_file}.')

	run_executable('grdchk')
	print(f'grdchk execution complete for {header_file}.')

	compile_code('adjoint', header_file, domain, dep_var = dep_var, ind_vars = ind_var)
	setup_adjoint([ind_var], header_file, domain)
	compile_code('adjoint', header_file, domain, clean = False, dep_var = dep_var, ind_vars = ind_var)
	print(f'adjoint compilation complete for {header_file}.')
	
	run_executable('adjoint')
	print(f'adjoint execution complete for {header_file}.')