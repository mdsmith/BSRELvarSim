#! /usr/bin/env python
# XXX have node range passed in as argument?
# XXX make fully independent

import math
import numpy as NP
from treegen import get_unrooted_tree
import os
import subprocess
import time

node = range(13,31)
local_processes = {}

def generate_settings(num_taxa, out_name, numReps, lenSeqs, rate_classes_per_branch = -1):
    this_set = {}
    if rate_classes_per_branch >= 1:
        rate_class_assignments = [  rate_classes_per_branch
                                    for _ in range((2*num_taxa) - 2)]
    else:
        rate_class_assignments = [  math.ceil(NP.random.exponential(.75))
                                    for _ in range((2*num_taxa) - 2)]
    print(rate_class_assignments)
    #rate_class_assignments = [ 1 for _ in range((2*num_taxa) - 2)]
    out = open(out_name, 'w')
    buffer = []

    #treeDef = "Tree bsrel_tree = ((1,2)Interior,3,4);\n"
    treeDef = "Tree bsrel_tree = "  + get_unrooted_tree(num_taxa) + ";\n"
    codoDef = "nuc3x4 = {4,3} [\"0.25\"];\n"
    nucbDef =   "nucleotide_bias_settings = { // relative to AG\n" \
                "\"AC\": 0.25,\n" \
                "\"AT\": 0.25,\n" \
                "\"CG\": 0.4,\n" \
                "\"CT\": 2.0,\n" \
                "\"GT\": 0.3\n" \
                "};\n"
    buffer.append(treeDef)
    buffer.append(codoDef)
    buffer.append(nucbDef)

    buffer.append("bsrel_settings = { \n")

    #for tax_name in range(1, num_taxa):
    for tax_name in range(1, (2*num_taxa) - 2):
        this_set[str(tax_name)] = generate_taxa( tax_name,
                                            rate_class_assignments[tax_name])
        buffer.append(params_to_string(this_set[str(tax_name)]))

    #this_set["Interior"] = generate_taxa( "Interior", 1)
    #buffer.append(params_to_string(this_set["Interior"]))

    buffer.append("};\n")

    repsDef = "replicates = " + str(numReps) + ";\n"
    seqlDef = "codons = " + str(lenSeqs) + ";\n"

    buffer.append(repsDef)
    buffer.append(seqlDef)

    out.writelines(buffer)
    out.close()
    return this_set

def params_to_string(param_list):
    entry = "\"" + param_list["name"] + "\" : { \"length\" : "
    entry += str(param_list["length"])
    entry += ",\n"
    entry += "\t\"omegas\" : {"
    for o,p in zip(param_list["omegas"], param_list["props"]):
        entry += "\t{ " + str(o) + ", " + str(p) + "}\n"
    entry += "}\n},\n"
    return entry

# redo these distributions properly, remove the prop vs omega sort.
def generate_taxa(tax_name, num_rate_classes):
    entry = {}
    entry["name"] = str(tax_name)
    length = NP.random.exponential(.25)
    #length = 0.5
    entry["length"] = length
    omegas = sorted([NP.random.exponential((.5 + (10*i))) for i in range(num_rate_classes)])

    # check omegas over one:
    under_one = False
    num_over_one = 0
    for omega in omegas:
        if omega < 1.0:
            under_one == True
        if omega > 1.0:
            num_over_one += 1
    # XXX this makes sure that there is at least one omega less than one
    #if under_one == False:
        #omegas[0] = omegas[0] % 1
    if num_over_one > 1:
        for omegaI in range(len(omegas) - 1):
            omegas[omegaI] = omegas[omegaI] % 1

    #omegas = [.5*i for i in range(1,num_rate_classes + 1)]
    entry["omegas"] = omegas

    props = []
    if len(omegas) > 1:
        starting_props = sorted(NP.random.rand(num_rate_classes - 1))
        last = -1
        props.append(starting_props[0])
        for p in starting_props:
            if last != -1:
                props.append(p - last)
            last = p
        props.append(1 - starting_props[-1])
    else:
        props.append(1)

    #remainingProp = 1
    #for i in range(num_rate_classes - 1):
        #props.append(max(0.1, (1 - NP.random.exponential(.5, 1)[0])) * remainingProp)
        ##props.append(0.9 * remainingProp)
        #remainingProp -= props[i]
    #props.append(1 - sum(props))
    #props = sorted(props)[::-1]
    entry["props"] = props
    return entry

def simulate(set_file_name, nodeI, node_processes = local_processes):
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
                    str(node[nodeI]),
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

def runSimulation():
    return 0

if __name__ == "__main__":
    runSimulation()