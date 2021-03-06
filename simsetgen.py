#! /usr/bin/env python
# XXX have node range passed in as argument?
# XXX make fully independent

import sys
import math
import numpy as NP
from treegen import get_unrooted_tree

def generate_settings(  num_taxa,
                        out_name,
                        numReps,
                        lenSeqs,
                        alpha_sim,
                        rate_classes_per_branch = -1):
    this_set = {}
    if rate_classes_per_branch >= 1:
        rate_class_assignments = [  rate_classes_per_branch
                                    for _ in range((2*num_taxa) - 2)]
    else:
        rate_class_assignments = NP.random.poisson( .25, (2*num_taxa) - 2)
        rate_class_assignments = [a + 1 for a in rate_class_assignments]
        #rate_class_assignments = [  math.ceil(NP.random.exponential(.75))
                                    #for _ in range((2*num_taxa) - 2)]
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
        this_set[str(tax_name)] = generate_taxa(tax_name,
                                                alpha_sim,
                                                rate_class_assignments[tax_name])
        if alpha_sim:
            buffer.append(params_to_string_alpha(this_set[str(tax_name)]))
        else:
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

def params_to_string_alpha(param_list):
    entry = "\"" + param_list["name"] + "\" : {"
    entry += "\t\"omegas\" : {"
    for n,s,p in zip(param_list["nonsyns"], param_list["syns"], param_list["props"]):
        entry += "\t{ " + str(n) + ", " + str(s) + ", "  + str(p) + "}\n"
    entry += "}\n},\n"
    return entry

# redo these distributions properly, remove the prop vs omega sort.
def generate_taxa(  tax_name,
                    alpha_sim,
                    num_rate_classes):
    entry = {}
    entry["name"] = str(tax_name)
    length = NP.random.exponential(.25)
    omegas = []
    nonsyns = []
    syns = []
    #length = 0.5
    if not alpha_sim:
        entry["length"] = length
    #omegas = sorted([NP.random.exponential((.5 + (10*i))) for i in range(num_rate_classes)])
    if alpha_sim:
        nonsyns = sorted([NP.random.exponential((.15 + (10*i))) for i in range(num_rate_classes)])
        syns = sorted([NP.random.exponential((.15 + (10*i))) for i in range(num_rate_classes)])
        # check omegas over one:
        omegas = [nonsyn/syn for nonsyn, syn in zip(nonsyns, syns)]
    else:
        omegas = sorted([NP.random.exponential((.15 + (10*i))) for i in range(num_rate_classes)])
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
    if alpha_sim:
        if num_over_one > 1:
            for nonsynI in range(len(nonsyns) - 1):
                nonsyns[nonsynI] = nonsyns[nonsynI] % 1
        entry["nonsyns"] = nonsyns
        entry["syns"] = syns
    else:
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

def generate_all_settings(  num_taxa,
                            num_dist,
                            num_reps,
                            len_seqs,
                            out_file,
                            alpha_sim,
                            rate_classes_per_branch = -1):
    inputParameterSets = {}
    for this_dist in range(num_dist):
        this_set = generate_settings(   num_taxa,
                                        out_file +  "." + str(this_dist),
                                        num_reps,
                                        len_seqs,
                                        alpha_sim,
                                        rate_classes_per_branch)
        inputParameterSets[str(this_dist)] = this_set
    # XXX this is a temporary holdover
    return inputParameterSets

if __name__ == "__main__":
    if len(sys.argv) == 6 or len(sys.argv) == 7:
        num_taxa, num_dist, num_reps, len_seqs, out_file = sys.argv[1:6]
        rate_classes_per_branch = -1
        if len(sys.argv) == 7:
            rate_classes_per_branch = sys.argv[-1]
        generate_all_settings(  int(num_taxa),
                                int(num_dist),
                                int(num_reps),
                                int(len_seqs),
                                out_file,
                                int(rate_classes_per_branch))
    else:
        print(sys.argv, file=sys.stderr)
        print(  "Valid usage:\n" \
                "\t- simsetgen <num_taxa> <num_dist> <num_reps> " \
                "<len_seqs> <out_file> [<rate_classes_per_branch>]\n" \
                "Where:\n" \
                "\t- <num_taxa>: number of taxa (power of two)\n" \
                "\t- <num_dist>: number of distributions\n" \
                "\t- <num_reps>: number of repetitions to simulate\n" \
                "\t- <len_seqs>: length of sequences to simulate\n" \
                "\t- <out_file>: output filename (may include path)\n" \
                "\t- [<rate_classes_per_branch>]: number of rate classes " \
                "per branch\n",
                file=sys.stderr)
