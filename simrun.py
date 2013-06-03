#! /usr/bin/env python
# XXX make fully independent

import os
import subprocess
import time

node_list = range(13,31)
local_processes = {}

def simulate(   set_file_name,
                nodeI,
                node_processes = local_processes,
                nodes = node_list):
    batchfile = open(set_file_name + ".bf", 'w')
    batchfile.write('inputRedirect = {};\n\n')
    batchfile.write('inputRedirect["00"]="'
                    + os.path.dirname(os.path.abspath(__file__))
                    + "/"
                    + set_file_name
                    + '";\n')
    batchfile.write('inputRedirect["01"]="'
                    + os.path.dirname(os.path.abspath(__file__))
                    + "/"
                    + set_file_name
                    + '.sim";\n')
    batchfile.write('ExecuteAFile' \
                    '("/data/veg/HSV/Simulation/01Initial/' \
                    'GenericSimulator.bf", inputRedirect);')
    batchfile.close()
    call_list = [   'bpsh',
                    str(nodes[nodeI]),
                    'HYPHYMP',
                    os.path.dirname(os.path.abspath(__file__))
                        + os.sep
                        + set_file_name
                        + '.bf']
    output_file = open( os.path.dirname(os.path.abspath(__file__))
                        + os.sep
                        + set_file_name
                        + '.sim'
                        + '.txt', 'w')
    node_processes[str(nodeI)] = subprocess.Popen( call_list,
                                                    stdout=output_file)
    time.sleep(1)

def run_simulation():
    return 0

if __name__ == "__main__":
    runSimulation()
