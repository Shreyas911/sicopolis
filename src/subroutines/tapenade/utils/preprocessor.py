#!/usr/bin/env python  
#
# Copyright 2015 Kshitij Kulshreshtha
# Copyright 2015 Sri Hari Krishna Narayanan

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
#  \brief Preprocessor for generating ADOL-C external functions. 
#
# */
# /**********************************************************************/

import re
import argparse
import os.path
import sys

dep_str="""
call cost_final()
call sico_end()
"""

def readFile(name):
  alllines = ''
  inputfile= open(name, 'r')
  for line in inputfile:
    alllines+=line
  inputfile.close()
  return alllines

def writeOutput(outstr, outfile):
  outputfile= open(outfile, 'w')
  outputfile.write(outstr)
  outputfile.close()


def appendFile(outfile, outstr):
  outputfile = open(outfile,'a')
  outputfile.write(outstr)
  outputfile.write('\n')
  outputfile.close()


def extract_begin(filename):
  datalines = readFile(filename)
  datalines = re.sub(r'\r',r'',datalines)
  s = r'!@\s*begin\s+tapenade_extract\s*@(.*?)!@\s*end\s+tapenade_extract\s*@'
  p = re.compile(s,re.M|re.S)
  match = p.findall(datalines)
  extlines = ""
  for m in match:
    extlines += m
  retline = ""
  for line in extlines.splitlines():
    line = re.sub(r"sico_main_loop", "sicopolis_tapenade", line)
    retline += line + "\n"
  retline += "use cost_m \n"
  return retline


def extract_head_c(filename):
  datalines = readFile(filename)
  datalines = re.sub(r'\r',r'',datalines)
  s = r'!@\s*begin\s+tapenade_extract\s*@(.*?)!@\s*end\s+tapenade_extract\s*@'
  p = re.compile(s,re.M|re.S)
  match = p.findall(datalines)
  extlines = ""
  for m in match:
    extlines += m
  begin_str = extract_begin("subroutines/general/sico_main_loop_m.F90")
  for line in extlines.splitlines():
    line = re.sub(r"subroutines/general/sico_maths_m.F90", "subroutines/tapenade/stubs/sico_maths_m_stub.F90",line)
    line = re.sub(r"!tapenade begin subroutine sicopolis_tapenade", begin_str, line)
    line = re.sub(r"!tapenade end subroutine sicopolis_tapenade", dep_str+"end subroutine sicopolis_tapenade", line)
    appendFile('numCore.F90',line)

def preprocess(args):
  writeOutput('','numCore.F90')
  for filename in args.filenames:
    extract_head_c(filename)

sys.path = [ os.getcwd() ] + sys.path
parser = argparse.ArgumentParser(description='Preprocess input file.')
parser.add_argument('filenames', metavar='f', type=str, nargs='+',
                   help='list of files')
args = parser.parse_args()

preprocess(args)
