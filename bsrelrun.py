#! /usr/bin/env python

import os, sys, subprocess, time
import glob, re
from subprocess import call
from threading import Thread
from queue import Queue, Empty
node_list = range(13,31)
local_processes = {}

jobs = Queue()

def run_BSREL(node, file_name):
    selectionTest = False
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
    time.sleep(1)

def run_all_BSREL(file_list):
    abspath = os.path.dirname(os.path.abspath(__file__)) + os.sep
    for file_name in file_list:
        if file_name[0] != os.sep:
            file_name = abspath + file_name
        jobs.put(file_name)

def run_job(node):
    while True:
        try:
            bsrel_args = jobs.get(block=False)
            run_BSREL(node, bsrel_args)
            jobs.task_done()
        except Empty:
            break

def get_files(indir):
    if os.path.isdir(indir):
        file_list = glob.glob(indir + os.sep + "*.sim.*")
    else:
        indir_dir = os.path.dirname(os.path.abspath(indir))
        indir_file = os.path.basename(os.path.abspath(indir))
        file_list = glob.glob(indir_dir + os.sep + "*"+indir_file + "*.sim.*")
    print(file_list)
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
    node_list = node_list[:num]
    return node_list

def bsrel_main(indir):
    file_list = get_files(indir)
    print(file_list)
    run_all_BSREL(file_list)
    for node in nodes(24):
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
