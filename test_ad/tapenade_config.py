import os


def grdchk(ind_vars, header, domain, lisdir = '/home/shreyas/lis-1.4.43/installation', netcdf_fortran_dir = '/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4', travis_ci ='', unit = '100'):
	for ind_var in ind_vars:
		os.system('cp ../test_ad/tapenade_m_adjoint_template.F90 subroutines/tapenade/tapenade_m.F90')
	
		new_file_lines = []
	
		with open('subroutines/tapenade/tapenade_m.F90') as f:
			file_lines = f.readlines()
	
		    ### Do stuff with the file in str form
			for line in file_lines:
	
				if ('!@ python_automated_grdchk @' in line):
					line = line + f'            orig_val = {ind_var}(j,i)\n            {ind_var}(j,i) = orig_val * perturbation\n'
					print(line)
				if ('!@ python_automated_grdchk IO begin @' in line):
					line = line + f'   open({unit}, file=\'GradientVals_{ind_var}_\'//trim(RUNNAME)//\'.dat\',&\n       form="FORMATTED", status="REPLACE")'
				if ('!@ python_automated_grdchk IO write @' in line):
					line = line + f'          write({unit}, fmt=\'(f40.20)\') gfd\n'
				if ('!@ python_automated_grdchk IO end @' in line):
					line = line + f'   close(unit={unit})\n'
				new_file_lines.append(line)

		with open('subroutines/tapenade/tapenade_m.F90', "w") as f:
			f.write(''.join(new_file_lines))
		
		os.system('make -f MakefileTapenade clean') 
		os.system(f'make -f MakefileTapenade drivergrdchk HEADER={header} DOMAIN_SHORT={domain} LISDIR={lisdir} NETCDF_FORTRAN_DIR={netcdf_fortran_dir} {travis_ci}')
		os.system('./drivergrdchk > /dev/null')

def get_imax_jmax(specs_file='sico_specs.h'):
	
	with open(specs_file) as f:
		file_lines = f.readlines()
		for line in file_lines:
			if ('#define IMAX' in line):
				IMAX = [int(s) for s in line.split() if s.isdigit()]
			if ('#define JMAX' in line):
				JMAX = [int(s) for s in line.split() if s.isdigit()]
	return IMAX[0], JMAX[0]

def adjoint(ind_vars, header, domain, run = 'separate', lisdir = '/home/shreyas/lis-1.4.43/installation', netcdf_fortran_dir = '/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4', travis_ci ='', unit = [100, 101]):
	
	IMAX, JMAX = get_imax_jmax()
	if(run == 'separate'):
		for ind_var in ind_vars:
	
			os.system('make -f MakefileTapenade clean') 
			os.system(f'make -f MakefileTapenade driveradjoint HEADER={header} DOMAIN_SHORT={domain} IND_VARS={ind_var} LISDIR={lisdir} NETCDF_FORTRAN_DIR={netcdf_fortran_dir} {travis_ci}')
			
			new_file_lines = []

			with open('numCore_cpp_b.f90') as f:
				file_lines = f.readlines()

				### Do stuff with the file in str form
				for line in file_lines:
					if ('CHARACTER(len=100) :: runname' in line):
						line = line + f'   INTEGER(i4b) :: i, j, p\n'
						line = line + f'   INTEGER(i4b), parameter :: points = 5\n'
						line = line + f'   INTEGER(i4b), DIMENSION(points) :: ipoints, jpoints\n'
						line = line + f'   DO p = 1, points\n'
						
						if(domain == 'grl'):
							line = line + f'   ipoints(p) = int(real({IMAX}/2))\n'
							line = line + f'   jpoints(p) = int(real({JMAX}/5)) + (p-1) * points\n'
						elif(domain == 'ant'):
							line = line + f'   ipoints(p) = int(real({IMAX}/3)) + int(real((.85-.33)*{IMAX}/points)) * (p -1)\n'
							line = line + f'   jpoints(p) = int(real({JMAX}/2))\n'
						else:
							raise Exception('Wrong domain')
						line = line + f'   END DO\n'
							
					if ('CALL SICO_INIT_B' in line):
						line = f'   close(unit={unit[1]})\n' + line
						line = f'   end do\n' + line
						line = f'   end do\n' + line
						line = f'   write ({unit[1]}, *) {ind_var}b(j,i)\n' + line
						line = f'   do j = 0, {JMAX}\n' + line
						line = f'   do i = 0, {IMAX}\n' + line
						line = f'   close(unit={unit[0]})\n' + line
						line = f'   end do\n' + line
						line = f'   write ({unit[0]}, *) {ind_var}b(j,i)\n' + line
						line = f'   j = jpoints(p)\n' + line
						line = f'   i = ipoints(p)\n' + line
						line = f'   do p = 1, points\n' + line
						line = f'   open({unit[1]}, file=\'AdjointVals_{ind_var}_\'//trim(RUNNAME)//\'.dat\',&\n       form="FORMATTED", status="REPLACE")\n' + line
						line = f'   open({unit[0]}, file=\'AdjointVals_{ind_var}_\'//trim(RUNNAME)//\'_limited.dat\',&\n       form="FORMATTED", status="REPLACE")\n' + line
						print(line)

					new_file_lines.append(line)
			with open('numCore_cpp_b.f90', "w") as f:
				f.write(''.join(new_file_lines))
			os.system(f'make -f MakefileTapenade driveradjoint HEADER={header} DOMAIN_SHORT={domain} IND_VARS={ind_var} LISDIR={lisdir} NETCDF_FORTRAN_DIR={netcdf_fortran_dir} {travis_ci}')
			os.system('./driveradjoint > /dev/null')
	elif(run == 'combined'):
		print("Work in progress")
	else:
		raise Exception('Choose the type of adjoint run you want.')

if __name__ == "__main__":
	try:
		# Change the current working Directory
		os.chdir("../src/")
		print("Directory changed to ", os.getcwd())
	except OSError:
		print("Can't change the Current Working Directory")
	adjoint(ind_vars = ['H', 'q_geo'], header = 'v5_grl20_ss25ka', domain = 'grl', run = 'separate')
	#grdchk(ind_vars = ['H', 'q_geo'], header = 'v5_grl20_ss25ka', domain = 'grl')
