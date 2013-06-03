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
    node_processes[str(nodeI)] = subprocess.Popen(  call_list,
                                                    stdout=output_file)
    time.sleep(1)

# XXX change this to actually use the in_file argument (so you can use input
# files with different, arbitrary, names)
def run_simulation( num_taxa,
                    num_dist,
                    in_file,
                    out_file,
                    node_processes = local_processes,
                    nodes = node_list):
    dist_done = 0
    while dist_done < num_dist:
        for this_dist in range(min(num_dist - dist_done, len(nodes))):
            dist_num = this_dist + dist_done
            simulate(   out_file + "." + str(dist_num),
                        this_dist,
                        node_processes,
                        nodes)
        for sub_p in node_processes.values():
            sub_p.wait()
        dist_done += min(num_dist-dist_done, len(nodes))
    return 0

if __name__ == "__main__":
    runSimulation()
