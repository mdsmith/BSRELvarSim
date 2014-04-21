#! /usr/bin/env python

import os, sys, subprocess, time
import glob, re
from subprocess import call
from threading import Thread
from queue import Queue, Empty
node_list = range(13,31)
local_processes = {}

jobs = Queue()

#HyPhy = 'HYPHYMP'
HyPhy = '/home/martin/Software/sergeiHyPhy/hyphy/HYPHYMP'

def pad(num, len):
    return str(num).zfill(len)

def run_BSREL(node, file_name, omega, alpha):
    if omega:
        response_list = []
        response_list.append("Universal")
        response_list.append("No")
        response_list.append(os.path.abspath(file_name))
        response_list.append("Y")
        response_list.append("All")
        response_list.append("")
        response_list.append(os.path.abspath(file_name) + ".recoveredOmega")
        batchfile = open(file_name + ".recoverOmega.bf", 'w')
        batchfile.write('inputRedirect = {};\n\n')
        for i, response in enumerate(response_list):
            batchfile.write('inputRedirect["'
                            + pad(i,2)
                            + '"]= "'
                            + response
                            + '";\n')
        batchfile.write('ExecuteAFile'
                        + '("'
                        + os.path.dirname(os.path.abspath(__file__))
                        + os.sep
                        + 'multiBSREL.bf"'
                        + ', inputRedirect);')
        batchfile.close()
        call_list = [   'bpsh',
                        str(node),
                        #'HYPHYMP',
                        HyPhy,
                        file_name + '.recoverOmega.bf']
        output_file = open( file_name + '.recoverOmega.txt', 'w')
        print(call_list)
        call(call_list, stdout=output_file)
        time.sleep(0.25)
    if alpha:
        response_list = []
        response_list.append("Universal")
        response_list.append("Yes")
        response_list.append(os.path.abspath(file_name))
        response_list.append("Y")
        response_list.append(os.path.abspath(file_name) + ".recoveredAlpha")
        batchfile = open(file_name + ".recoverAlpha.bf", 'w')
        batchfile.write('inputRedirect = {};\n\n')
        for i, response in enumerate(response_list):
            batchfile.write('inputRedirect["'
                            + pad(i,2)
                            + '"]= "'
                            + response
                            + '";\n')
        batchfile.write('ExecuteAFile'
                        + '("'
                        + os.path.dirname(os.path.abspath(__file__))
                        + os.sep
                        + 'multiBSREL.bf"'
                        + ', inputRedirect);')
        batchfile.close()
        call_list = [   'bpsh',
                        str(node),
                        #'HYPHYMP',
                        HyPhy,
                        file_name + '.recoverAlpha.bf']
        output_file = open( file_name + '.recoverAlpha.txt', 'w')
        print(call_list)
        call(call_list, stdout=output_file)
        time.sleep(0.25)




    '''
    #selectionTest = False
    selectionTest = True
    batchfile = open(file_name + ".recover.bf", 'w')
    batchfile.write('inputRedirect = {};\n\n')
    batchfile.write('inputRedirect["00"]= "Universal";\n')
    if selectionTest:
        batchfile.write('inputRedirect["01"]= "No";\n')
    else:
        batchfile.write('inputRedirect["01"]= "Yes";\n')
    batchfile.write('inputRedirect["02"]="'
                    + file_name
                    + '";\n')
    batchfile.write('inputRedirect["03"]= "Y";\n')
    if selectionTest:
        batchfile.write('inputRedirect["04"]= "All";\n')
        batchfile.write('inputRedirect["05"]= "";\n')
        batchfile.write('inputRedirect["06"]="'
                        + file_name
                        + '.recovered";\n')
    else:
        batchfile.write('inputRedirect["04"]="'
                        + file_name
                        + '.recovered";\n')

    batchfile.write('ExecuteAFile'
                    + '("'
                    + os.path.dirname(os.path.abspath(__file__))
                    + os.sep
                    + 'multiBSREL.bf"'
                    + ', inputRedirect);')
    batchfile.close()
    call_list = [   'bpsh',
                    str(node),
                    'HYPHYMP',
                    file_name + '.recover.bf']
    output_file = open( file_name + '.recover.txt', 'w')
    print(call_list)
    call(call_list, stdout=output_file)
    time.sleep(0.25)
    '''
    run_GTRgamma(node, file_name)

def run_GTRgamma(node, file_name):
    batchfile = open(file_name + ".recover.GTRg.bf", 'w')
    batchfile.write('inputRedirect = {};\n\n')
    batchfile.write('inputRedirect["00"]="'
                    + file_name
                    + '";\n')
    batchfile.write('inputRedirect["01"]= "CUSTOM";\n')
    batchfile.write('inputRedirect["02"]= "012345";\n')
    batchfile.write('inputRedirect["03"]= "Global w/variation";\n')
    batchfile.write('inputRedirect["04"]= "Gamma";\n')
    batchfile.write('inputRedirect["05"]= "4";\n')
    batchfile.write('inputRedirect["06"]= "Estimated";\n')
    batchfile.write('inputRedirect["07"]= "Observed";\n')
    batchfile.write('inputRedirect["08"]= "Y";\n')
    batchfile.write('inputRedirect["09"]= "Don\'t Display";\n')
    batchfile.write('ExecuteAFile'
                    + '("'
                    + '/usr/local/lib/hyphy/TemplateBatchFiles/'
                    + 'AnalyzeNucProtData.bf"'
                    + ', inputRedirect);')
    batchfile.close()
    call_list = [   'bpsh',
                    str(node),
                    #'HYPHYMP',
                    HyPhy,
                    file_name + '.recover.GTRg.bf']
    output_file = open( file_name + '.recovered.GTRg.txt', 'w')
    print(call_list)
    call(call_list, stdout=output_file)
    time.sleep(1)

def run_all_BSREL(file_list, omega=True, alpha=False):
    abspath = os.path.dirname(os.path.abspath(__file__)) + os.sep
    for file_name in file_list:
        if file_name[0] != os.sep:
            file_name = abspath + file_name
        jobs.put((file_name, omega, alpha))

def run_job(node):
    while True:
        try:
            bsrel_args = jobs.get(block=False)
            run_BSREL(node, *bsrel_args)
            jobs.task_done()
        except Empty:
            break

def get_files(indir):
    if os.path.isdir(indir):
        file_list = glob.glob(indir + os.sep + "*.sim.*")
    else:
        indir_dir, indir_file = os.path.split(os.path.abspath(indir))
        file_list = glob.glob(indir_dir + os.sep + "*"+indir_file + "*.sim.*")
    #print(file_list)
    file_list = [a for a in file_list
                if re.search("/\w+\.\d+\.sim\.\d+$", a) != None]
    return file_list

def nodes(num):
    import shlex
    from subprocess import Popen, PIPE
    cmd = shlex.split('''beomap --all-nodes --no-local --exclude
    0:1:2:3:4:5:6''')
    proc = Popen(cmd, stdout=PIPE)
    stdout, _ = proc.communicate()
    node_list = [int(node) for node in
                stdout.decode('utf-8').strip().split(':')]
    node_list.reverse()
    # XXX testing multiple jobs on one node, increase requested node number
    # above 24 to do it. Untested.
    if len(node_list) > num:
        node_list = node_list[:num]
    else:
        while len(node_list) < num:
            node_list = node_list + node_list
        node_list = node_list[:num]
    return node_list

def bsrel_main(indir, omega=True, alpha=False):
    file_list = get_files(indir)
    print(file_list)
    run_all_BSREL(file_list, omega=omega, alpha=alpha)
    for node in nodes(96):
        t = Thread(target=run_job, args=(node,))
        t.daemon = True
        t.start()
    jobs.join()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        indir = sys.argv[1]
        bsrel_main(indir)
    else:
        print(  "Valid usage:\n" \
                "\t- bsrelrun <indir>\n" \
                "Where:\n" \
                "\t- <indir>: directory containing files to run BSREL on\n",
                file=sys.stderr)
        exit(1)
