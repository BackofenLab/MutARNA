#!/usr/bin/env python
import sys
import argparse
from shutil import copyfile
from os.path import isfile
from subprocess import call
from subprocess import check_output

DESCRIPTION = """
dot2circ - convert ViennaRNA dotplots to circular base-pairing diagrams

Given an RNA sequence, calculate base-pairing probabilities with RNAplfold and
plot them using circos.
"""

EPILOG = """
Status: rickety prototype.

fixed relative paths, has to be called from program directory.
no checking for valid parameters.
no checking if required software is available.
folding temperature set to 25C.
folding parameters set to W150 L100.
color annotation set to my current test sequence.
"""

RNAPLFOLD_PARAMS = "-W 250 -L 200 "
RNAFOLD_PARAMS = ""

PERLBIN = ""
CIRCOSBIN = "circos"

# parse command line arguments
parser = argparse.ArgumentParser(description=DESCRIPTION, epilog=EPILOG)
parser.add_argument(
    "--sequence",
    help="input RNA sequence")
parser.add_argument(
    "--title",
    help="Set this title for the circos figure")
parser.add_argument(
    "--dp-file",
    help="Do not run fold program , use precomputed dot plot")

parser.add_argument(
    "--prefix",
    default="dot2circ",
    help="prefix of output files")
parser.add_argument("--local-fold",
                    type=bool,
                    default = False,
                    help="compute local dotplot using RNAplFold")
parser.add_argument("--temperature",
                    type=int,
                    default = 37,
                    help="folding temperature")
parser.add_argument(
    "--outputdir",
    default="./",
    help="output dir")



args = parser.parse_args()

if args.sequence is not None and args.dp_file is not None:
  print ("ERROR: Either provide a sequence or a prefolded dotplot but not both!")
  sys.exit(1)

if args.sequence is not None:
  # fold
  if args.local_fold is True:
      foldcmd = 'echo {0} | RNAplfold {1} '.format(args.sequence, RNAPLFOLD_PARAMS)
  else:
      foldcmd = 'echo {0} | RNAfold -p '.format(args.sequence)

  foldcmd += '-T {}'.format(args.temperature)

      
  foldout = check_output(foldcmd, shell=True)
else: # dotplot is provided
  if not isfile(args.dp_file):
    print ("ERROR: dotplot file {} does not exist!".format(args.dp_file))
    sys.exit(1)
  if args.local_fold is True:
    copyfile(args.dp_file, './plfold_dp.ps')
  else:
    copyfile(args.dp_file, './dot.ps')  

# create circos data
if args.local_fold is True:
    call(['./parse_plfold.sh'], shell=True)
else:
    call(['./parse_rnafold.sh'], shell=True)



# run circos
circoscmd = '{0} {1} -outputdir {2} -param image/file**="{3}.png" -param image/radius*=1000p '.format(
    PERLBIN,
    CIRCOSBIN,
    args.outputdir,
    args.prefix)

# Overide chr seq name
if (args.title):
  with open('./data/circos.karyotype.txt') as in_txt:
    line = in_txt.readline()
  assert ('chr - seq seq ' in line)
  with open('./data/circos.karyotype.txt', 'w') as out_txt:
    out_txt.write(line.replace('chr - seq seq ', 'chr - seq {} '.format(args.title)))
  circoscmd += "-param ideogram/show_label=yes "


#circoscmd += "" 
#"""\
#-param links/link/rules/annot1_1st_stem_start=139 \
#  -param links/link/rules/annot1_1st_stem_end=157 \
#-param links/link/rules/annot1_2nd_stem_start=181 \
#  -param links/link/rules/annot1_2nd_stem_end=206 \
#-param links/link/rules/annot2_1st_stem_start=179 \
#  -param links/link/rules/annot2_1st_stem_end=198 \
#-param links/link/rules/annot2_2nd_stem_start=227 \
#  -param links/link/rules/annot2_2nd_stem_end=246"""

    
    
print ('calling: "{0}"'.format(circoscmd))
print ('======circos output log========')
circosout = check_output(circoscmd, shell=True)
print (circosout.decode('ascii'))
print ('===============================')
