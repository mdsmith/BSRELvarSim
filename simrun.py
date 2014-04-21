#! /usr/bin/env python

import os, sys
import subprocess
import time

node_list = range(13,31)
local_processes = {}

def simulate(   set_file_name,
                nodeI,
                alpha_sim,
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
    batchfile.write('ExecuteAFile' 
                    #'("/data/veg/HSV/Simulation/01Initial/' \
                    #'GenericSimulator.bf", inputRedirect);')
                    + '("'
                    + os.path.dirname(os.path.abspath(__file__))
                    + os.sep)
    if alpha_sim:
        batchfile.write('BSRELsimFromAlphaSettings.bf", inputRedirect);')
    else:
        batchfile.write('BSRELsimFromSettingsMod.bf", inputRedirect);')
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
    time.sleep(0.25)

# XXX change this to actually use the in_file argument (so you can use input
# files with different, arbitrary, names)
def run_simulation( num_dist,
                    in_file,
                    out_file,
                    alpha_sim,
                    node_processes = local_processes,
                    nodes = node_list):
    dist_done = 0
    while dist_done < num_dist:
        for this_dist in range(min(num_dist - dist_done, len(nodes))):
            dist_num = this_dist + dist_done
            simulate(   out_file + "." + str(dist_num),
                        this_dist,
                        alpha_sim,
                        node_processes,
                        nodes)
        for sub_p in node_processes.values():
            sub_p.wait()
        dist_done += min(num_dist-dist_done, len(nodes))
    return 0

if __name__ == "__main__":
    if len(sys.argv) == 5:
        num_dist, in_file, out_file, alpha_sim = sys.argv[1:5]
        run_simulation( int(num_dist),
                        in_file,
                        out_file,
                        alpha_sim)
    else:
        #print(sys.argv, file=sys.stderr)
        print(  "Valid usage:\n" \
                "\t- simrun <num_dist> <in_file> <out_file>\n" \
                "Where:\n" \
                "\t- <num_dist>: number of distributions\n" \
                "\t- <in_file>: input filename (may include path)\n" \
                "\t- <out_file>: output filename (may include path)\n" \
                "\t- <alpha_sim>: simulate using alpha variation\n",
                file=sys.stderr)
