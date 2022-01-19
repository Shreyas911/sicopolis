import numpy as np
import sys
import argparse
import os

def check_GRL_DISC(sico_specs_file):

	try:
		with open(sico_specs_file) as f:
			file_lines = f.readlines()
			for line in file_lines:
				if ('#define GRL' in line):
					GRL_flag = True
					break
				else:
					GRL_flag = False
			for line in file_lines:
				if ('#define DISC' in line):
					disc = [int(s) for s in line.split() if s.isdigit()]
					if disc[0]>0:
						DISC_flag = True
						break
					else:
						DISC_flag = False
				else:
					DISC_FLAG = False
                
		return GRL_flag and DISC_flag

	except FileNotFoundError :
		print(f'{sico_specs_file} does not exist')
		sys.exit(1)

	except Exception as err :
		print('Some error in checking the GRL and DISC case.')
		print(err)
		sys.exit(1)
	
  
  
  

if __name__ == '__main__':



	#	try:
	#		#Change the current working Directory
	#		os.chdir("../../")
	#
	#	except OSError:
	#		print("Can't change the Current Working Directory")
	#		sys.exit(1)


	
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--specsFile", help="name of specs file", type=str, required=True)
	args = parser.parse_args()

	
	result = check_GRL_DISC(args.specsFile)

	
	print(int(result))
