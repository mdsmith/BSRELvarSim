#! /usr/bin/env python

import os
import re
import glob
import argparse
import subprocess

def get_files(input, suffix):
    if os.path.isdir(input):
        input += os.sep
    file_list = glob.glob(input + "*" + suffix)
    file_list = [os.path.abspath(a) for a in file_list]
    return file_list

def batch(file_list, new_suffix):
    for file in file_list:
        bf = open(file + ".bls.bf", 'w')
        bf.write('inputRedirect = {};\n\n')
        bf.write('inputRedirect["00"] = "' + file + '";\n')
        bf.write('inputRedirect["01"] = "Y";\n')
        bf.write('inputRedirect["02"] = "BSREL";\n')
        bf.write('ExecuteAFile("'
                    + os.path.dirname(os.path.abspath(__file__))
                    + os.sep
                    + 'BLcalc.bf", inputRedirect);')
        bf.close()
        call_list = [   'HYPHYMP',
                        file + ".bls.bf"]
        output_file = open(file + ".bls.csv", 'w')
        subprocess.Popen(call_list, stdout=output_file)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input directory and general prefix")
    parser.add_argument("--suffix",
                        help="an optional suffix to identify the correct \
                        files, default is .fit",
                        default=".fit")
    parser.add_argument("--out_suffix",
                        help="an optional output suffix.  Default is \
                        .bls.csv",
                        default=".bls.csv")
    args = parser.parse_args()

    file_list = get_files(args.input, args.suffix)
    batch(file_list, args.out_suffix)
    exit(0)

