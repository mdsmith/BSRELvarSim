#! /usr/bin/env python
# XXX make fully independent

import os, subprocess, time
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
                    '("/usr/local/lib/hyphy/TemplateBatchFiles/BranchSiteREL.bf"'
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
