#! /usr/bin/env python

import os, sys, subprocess, time
node_list = range(13,31)
local_processes = {}

def run_BSREL(  set_file_name,
                nodeI,
                rep,
                node_processes = local_processes,
                nodes = node_list):
    batchfile = open(set_file_name + ".sim." + str(rep) + ".recover.bf", 'w')
    batchfile.write('inputRedirect = {};\n\n')
    batchfile.write('inputRedirect["00"]= "Universal";\n')
    batchfile.write('inputRedirect["01"]= "No";\n')
    batchfile.write('inputRedirect["02"]="'
                    + os.path.dirname(os.path.abspath(__file__))
                    + "/"
                    + set_file_name
                    + '.sim.'
                    + str(rep)
                    + '";\n')
    batchfile.write('inputRedirect["03"]= "Y";\n')
    batchfile.write('inputRedirect["04"]= "All";\n')
    batchfile.write('inputRedirect["05"]= "";\n')
    batchfile.write('inputRedirect["06"]="'
                    + os.path.dirname(os.path.abspath(__file__))
                    + "/"
                    + set_file_name
                    + '.sim.'
                    + str(rep)
                    + '.recovered";\n')
    batchfile.write('ExecuteAFile' \
                    #'("/usr/local/lib/hyphy/TemplateBatchFiles/BranchSiteREL.bf"'
                    '("/home/martin/Software/multimodelBSREL/multiBSREL.bf"'
                    ', inputRedirect);')
    batchfile.close()
    call_list = [   'bpsh',
                    str(nodes[nodeI]),
                    'HYPHYMP',
                    os.path.dirname(os.path.abspath(__file__))
                        + os.sep
                        + set_file_name
                        + '.sim.'
                        + str(rep)
                        + '.recover.bf']
    output_file = open( os.path.dirname(os.path.abspath(__file__))
                        + os.sep
                        + set_file_name
                        + '.sim.'
                        + str(rep)
                        + '.recover.txt', 'w')
    node_processes[str(nodeI)] = subprocess.Popen( call_list,
                                                    stdout=output_file)
    time.sleep(1)

def run_all_BSREL(  num_dist,
                    out_file,
                    num_reps,
                    node_processes = local_processes,
                    nodes = node_list):
    dist_done = 0
    while dist_done < num_dist:
        for this_dist in range(min(num_dist-dist_done, len(nodes))):
            dist_num = this_dist + dist_done
            for rep in range(int(num_reps)):
                run_BSREL(  out_file + "." + str(dist_num),
                            this_dist,
                            rep,
                            node_processes,
                            nodes)
        for sub_p in node_processes.values():
            sub_p.wait()
        dist_done += min(num_dist-dist_done, len(nodes))

# XXX allow this to take an infile
if __name__ == "__main__":
    if len(sys.argv) == 4:
        num_dist, out_file, num_reps = sys.argv[1:4]
        run_all_BSREL(  int(num_dist),
                        out_file,
                        int(num_reps))
    else:
        print(  "Valid usage:\n" \
                "\t- bsrelrun <num_dist> <out_file> <num_reps>\n" \
                "Where:\n" \
                "\t- <num_dist>: number of distributions\n" \
                "\t- <out_file>: output filename (may include path)\n" \
                "\t- <num_reps>: the number of reps you simulated\n",
                file=sys.stderr)
