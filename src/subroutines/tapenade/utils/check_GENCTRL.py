#!/usr/bin/env python  
#
# Copyright 2023 Shreyas Sunil Gaikwad

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:

# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# /**********************************************************************/
# /*! \file  
#
#  \brief Check for MakefileTapenade the value of the DISC flag in case
#         the domain is GRL
# */
# /**********************************************************************/

import numpy as np
import sys
import argparse
import os

def check_GENCTRL(ad_file):

	'''
	ad_file: Location of the AD header file
	GENCTRL_flag: True if ALLOW_GENCTRL is defined, else False
	'''

	try:
		with open(ad_file) as f:
			file_lines = f.readlines()
			for line in file_lines:
				if ('#define ALLOW_GENCTRL' in line):
					return 1
		return 0

	except FileNotFoundError :
		print(f'{ad_file} does not exist')
		sys.exit(1)

	except Exception as err :
		print('Some error in checking the GENCTRL case.')
		print(err)
		sys.exit(1)

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--ADspecsFile", help="name of AD specs file", type=str, required=True)
	args = parser.parse_args()

	result = check_GENCTRL(args.ADspecsFile)

	print(int(result))
