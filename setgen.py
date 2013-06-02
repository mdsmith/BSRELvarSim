#! /usr/bin/env python
# XXX make an output directory an optional input
# make simulation and processing separate programs/libraries

import sys
import math
import numpy as NP
import random
import os, time
import subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#node = [21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
node = range(13,31)
node_processes = {}

num_taxa = 4
#rate_class_limit = True
rate_class_limit = False

def get_tree(num_taxa):
    return  ("("
            + gen_newick([a for a in range(1, num_taxa + 1)], (2*num_taxa) - 1)
            + ")")

def gen_newick(part, number):
    if len(part) == 2:
        return  ("("
                + str(part[0])
                + ","
                + str(part[1])
                + ")"
                + str(number))
    else:
        return  ("("
                + gen_newick(   part[:int(len(part)/2)],
                                number - int(len(part)/2))
                + ","
                + gen_newick(part[int(len(part)/2):], number - 1)
                + ")"
                + str(number))

def get_unrooted_tree(num_taxa):
    return  ("("
            + gen_unrooted_newick(  [a for a in range(1, num_taxa + 1)],
                                    (2*num_taxa)-1,
                                    (2*num_taxa)-3)
            + ")")

def gen_unrooted_newick(part, number, max):
    if len(part) == 2:
        if number > max:
            return  str(part[0]) + "," + str(part[1])
        else:
            return  "(" + str(part[0]) + "," + str(part[1]) + ")" + \
                    str(number)
    else:
        if number > max:
            return  (gen_unrooted_newick(part[:int(len(part)/2)], number -
                    int(len(part)/2), max)
                    + ","
                    + gen_unrooted_newick(  part[int(len(part)/2):],
                                            number - 1,
                                            max))
        else:
            return  ("("
                    + gen_unrooted_newick(  part[:int(len(part)/2)],
                                            number - int(len(part)/2),
                                            max)
                    + ","
                    + gen_unrooted_newick(part[int(len(part)/2):], number - 1, max)
                    + ")"
                    + str(number))

def generate_settings(num_taxa, out_name):
    this_set = {}
    if rate_class_limit == True:
        rate_class_assignments = [ 1 for _ in range((2*num_taxa) - 2)]
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

def simulate(set_file_name, nodeI):
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

def run_BSREL(set_file_name, nodeI, rep):
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
                    str(node[nodeI]),
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

# return a dict of the different taxa (which are themselves dicts
def recover_fit(num_taxa, rec_file_name, dist, rep):
    results = {}
    fullfilename = rec_file_name + ".sim." + str(rep) + ".recovered.fit"
    print(fullfilename)
    resultsfile = open(fullfilename, 'r')
    lines = resultsfile.readlines()
    lines = [line for line in lines if line[:11] == "mixtureTree"]
    for taxaI in range(1, (2*num_taxa) - 2):
        results[str(taxaI)] = {}
    for line in lines:
        tokens = line.split('=')
        # check for and dispense with constraints
        if len(tokens) == 2:
            tokens[0] = tokens[0][12:]
            name, parameter = tokens[0].split('.')
            value = tokens[1].split(';')[0]
        results[name.upper()][parameter] = value
    for taxa in results.keys():
        try:
            results[taxa] = format_results(taxa, results[taxa])
        except KeyError:
            print("dist: ", dist)
            print("taxa: ", taxa)
            print(results[taxa])
            exit(1)
    return results

def recover_csv(num_taxa, rec_file_name, dist, rep):
    results = {}
    results_file = open(rec_file_name + ".sim." + str(rep) + ".recovered", 'r')
    lines = results_file.readlines()
    lines = lines[1:]
    for line in lines:
        line = line.split(',')
        results[line[0]] = {}
        results[line[0]]["name"] = line[0]
        results[line[0]]["length"] = line[-1]
        results[line[0]]["omegas"] = []
        results[line[0]]["props"] = []
        over_one = 0
        results[line[0]]["pval"] = line[7]
        if (line[3]) != "0":
            if line[3] == "inf":
                results[line[0]]["omegas"].append(10000)
            else:
                results[line[0]]["omegas"].append(line[3])
            results[line[0]]["props"].append(line[4])
            over_one = 1
        for oI in range(1, (int(line[2]) + 1 - int(over_one))):
            results[line[0]]["omegas"].append(str(line[1]))
            results[line[0]]["props"].append(str(   (1 - float(line[4]))
                                                    /float(line[2])))
    return results

def recover_simulated(num_taxa, rec_file_name, dist, rep):
    input = {}
    sim_file = open(rec_file_name + ".sim.txt", 'r')
    lines = sim_file.readlines()
    names = [line for line in lines if line[:4] == "Node"]
    names = [name.split('.')[1].rstrip() for name in names]
    lengths = [line.strip() for line in lines if line.strip()[:6] == "Length"]
    lengths = [length.split('=')[1] for length in lengths]
    for name, length in zip(names, lengths):
        input[name] = {}
        input[name]["name"] = name
        input[name]["length"] = length
    # ad hoc FSM:
    cur_name = ""
    cur_omega_list = []
    cur_prop_list = []
    for line in lines:
        if line[:4] == "Node":
            if cur_name != "":
                input[cur_name]["omegas"] = cur_omega_list
                input[cur_name]["props"] = cur_prop_list
                cur_omega_list = []
                cur_prop_list = []
            cur_name = line.split('.')[1].strip()
        if line.lstrip()[:5] == "omega":
            cur_omega_list.append(line.split('=')[1].strip())
        if line.lstrip()[:6] == "weight":
            cur_prop_list.append(line.split('=')[1].strip())
    input[cur_name]["omegas"] = cur_omega_list
    input[cur_name]["props"] = cur_prop_list
    return input

def format_results(taxa, taxa_dict):
    entry = {}
    entry["name"] = taxa
    entry["omegas"] = []
    entry["props"] = []
    omega = False
    for key in taxa_dict.keys():
        if key == 't':
            omega = True
    if omega == True:
        entry["length"] = taxa_dict["t"]
        for key in taxa_dict.keys():
            if key[:5] == "omega":
                entry["omegas"].append(taxa_dict[key])
            if key[:4] == "Paux":
                entry["props"].append(taxa_dict[key])
        entry["props"] = convolve_props(entry["props"])
    else:
        entry["length"] = taxa_dict['syn']
        if taxa_dict['syn'] == '0':
            entry["omegas"].append('20000')
        else:
            entry["omegas"].append(str(float(taxa_dict["nonsyn"])/float(taxa_dict["syn"])))
        entry["props"] = 1
    return entry

def convolve_props(prop_list):
    # go from n-1 pauxs to n props
    tbr = []
    remaining_weight = 1
    for prop in prop_list:
        this_weight = str(float(prop)*remaining_weight)
        tbr.append(this_weight)
        remaining_weight -= float(this_weight)
    tbr.append(remaining_weight)
    return tbr

# too specific, can be genericized. DEPRICATED
#def params_to_meandnds(omegas, props):
    #return sum([float(omegas[i])*float(props[i]) for i in range(len(omegas))])

def weighted_average(values, weights):
    return sum([float(values[i])*float(weights[i]) for i in range(len(values))])


def coord_to_bin(coord, bins):
    tbr = coord
    for i,bin in enumerate(bins):
        if bin > coord:
            break
        else:
            tbr = i
    #print("coord: ", coord)
    #print("binval: ", bins[tbr])
    #print("bin: ", tbr)
    return tbr

def gen_lin_bins(min, max, num):
    bin_size = (max - min)/num
    #print("binsize: ", bin_size)
    if bin_size == 0:
        return [min]
    return [min+(i*bin_size) for i in range(num)]

def gen_log_bins(min, max, num):
    log_min = math.log10(min)
    log_max = math.log10(max)
    bin_size = (log_max - log_min)/num
    if bin_size == 0:
        return [min]
    return [10**(log_min + (i * bin_size)) for i in range(num)]

# XXX incomplete
def gen_semi_log_bins(min, max, num):
    tbr = []
    first_bin = -8
    while (10**first_bin) < min:
        first_bin += 1
    first_bin -= 1
    #tbr.append(10**first_bin)
    print("start: ", 10**first_bin)
    last_bin = first_bin
    while (10**last_bin) < max:
        last_bin += 1
    #tbr = [10**i for i in range(first_bin, last_bin)]
    for i in range(first_bin, last_bin):
        if len(tbr) <= num:
            tbr.append(10**i)
    if len(tbr) < num:
        tbr.append(max)
    # XXX now go through and generate the numbers in between
    while len(tbr) < num:
        cur_index = 1
        while cur_index < len(tbr):
            if len(tbr) < num:
                tbr.insert(cur_index, ) # XXX what to insert?
                cur_index += 1
            cur_index += 1
    #if len(tbr) < num:
        #tbr.append(10**last_bin)
    print("end: ", 10**(last_bin-1))
    #while len(tbr) < num:
        #current_index = 1
        #while current_index < len(tbr) and len(tbr) < num:
            #tbr.insert(current_index, tbr[current_index]/2)
            #current_index += 2
    return tbr

# this wasn't really taking the mean... depricated
#def tuple_to_heatarray( data,
#                        x_num_bins,
#                        y_num_bins,
#                        xmax=0,
#                        xmin=1,
#                        ymax=0,
#                        ymin=1,
#                        xtype="linear",
#                        ytype="linear"):
#    xmin = min([x for x,y,p in data])
#    xmax = max([x for x,y,p in data])
#    ymin = min([y for x,y,p in data])
#    ymax = max([y for x,y,p in data])
#    if xtype == "linear":
#        xbins = gen_lin_bins(xmin, xmax, x_num_bins)
#    else:
#        xbins = gen_log_bins(xmin, xmax, x_num_bins)
#    if ytype == "linear":
#        ybins = gen_lin_bins(ymin, ymax, y_num_bins)
#    else:
#        ybins = gen_log_bins(ymin, ymax, y_num_bins)
#    print("xbins: \n", xbins)
#    heatarray = NP.ones((len(xbins), len(ybins)))
#    heatarray[:] = 99
#    heatcount = NP.zeros((len(xbins), len(ybins)))
#    #print(ybins)
#    for x,y,p in data:
#        x_bin = coord_to_bin(x, xbins)
#        y_bin = coord_to_bin(y, ybins)
#        current_val = heatcount[x_bin][y_bin] * heatarray[x_bin][y_bin]
#        heatcount[x_bin][y_bin] += 1
#        new_val = (current_val + p)/heatcount[x_bin][y_bin]
#        heatarray[x_bin][y_bin] = new_val
#        #heatarray[coord_to_bin(x, xbins)][coord_to_bin(y, ybins)] = p
#
#    for xI in range(len(heatarray[:,1])):
#        for yI in range(len(heatarray[1,:])):
#            if heatarray[xI][yI] == 99:
#                heatarray[xI][yI] = NP.NAN
#    return heatarray, xbins, ybins

# false positive to the ratio of false positives/input negative aka false
# positive/false positive + true negative
def tuple_to_false_positive_ratio_heatarray(true_negative_results,
                                            false_positive_results,
                                            x_num_bins,
                                            y_num_bins,
                                            xmax=1,
                                            xmin=0,
                                            ymax=1,
                                            ymin=0,
                                            xtype="linear",
                                            ytype="linear"):
    xmin = min([x for x,y,p in false_positive_results])
    xmax = max([x for x,y,p in false_positive_results])
    ymin = min([y for x,y,p in false_positive_results])
    ymax = max([y for x,y,p in false_positive_results])
    if xtype == "linear":
        xbins = gen_lin_bins(xmin, xmax, x_num_bins)
    else:
        xbins = gen_log_bins(xmin, xmax, x_num_bins)
    if ytype == "linear":
        ybins = gen_lin_bins(ymin, ymax, y_num_bins)
    else:
        ybins = gen_log_bins(ymin, ymax, y_num_bins)
    fp_count = NP.zeros((len(xbins), len(ybins)))
    n_count = NP.zeros((len(xbins), len(ybins)))
    for x,y,p in true_negative_results:
        x_bin = coord_to_bin(x, xbins)
        y_bin = coord_to_bin(y, ybins)
        n_count[x_bin, y_bin] += 1
    for x,y,p in false_positive_results:
        x_bin = coord_to_bin(x, xbins)
        y_bin = coord_to_bin(y, ybins)
        fp_count[x_bin, y_bin] += 1
    for xI in range(len(n_count[:,0])):
        for yI in range(len(n_count[0,:])):
            if n_count[xI,yI] == 0:
                n_count[xI,yI] = NP.NAN
            else:
                n_count[xI,yI] = fp_count[xI,yI]/n_count[xI,yI]
    return n_count, xbins, ybins

# return heatarray of true positive/(true positive + false negative) aka true
# postive/input positive
def tuple_to_true_positive_ratio_heatarray( all_positive_inputs,
                                            x_num_bins,
                                            y_num_bins,
                                            xmax=1,
                                            xmin=0,
                                            ymax=1,
                                            ymin=0,
                                            xtype="linear",
                                            ytype="linear"):
    xmin = min([x for x,y,p in all_positive_inputs])
    xmax = max([x for x,y,p in all_positive_inputs])
    ymin = min([y for x,y,p in all_positive_inputs])
    ymax = max([y for x,y,p in all_positive_inputs])
    if xtype == "linear":
        xbins = gen_lin_bins(xmin, xmax, x_num_bins)
    else:
        xbins = gen_log_bins(xmin, xmax, x_num_bins)
    if ytype == "linear":
        ybins = gen_lin_bins(ymin, ymax, y_num_bins)
    else:
        ybins = gen_log_bins(ymin, ymax, y_num_bins)
    tp_count = NP.zeros((len(xbins), len(ybins)))
    p_count = NP.zeros((len(xbins), len(ybins)))
    for x,y,p in all_positive_inputs:
        x_bin = coord_to_bin(x, xbins)
        y_bin = coord_to_bin(y, ybins)
        p_count[x_bin][y_bin] += 1
        if p < .05:
            tp_count[x_bin][y_bin] += 1
    for xI in range(len(p_count[:,0])):
        for yI in range(len(p_count[0,:])):
            if p_count[xI,yI] == 0:
                p_count[xI,yI] = NP.NAN
            else:
                p_count[xI,yI] = tp_count[xI,yI]/p_count[xI,yI]
    return p_count, xbins, ybins

def tuple_to_averaged_heatarray(tuples,
                                x_num_bins,
                                y_num_bins,
                                xmax=1,
                                xmin=1,
                                ymax=1,
                                ymin=1,
                                xtype="linear",
                                ytype="linear"):
    xmin = min([x for x,y,p in tuples])
    xmax = max([x for x,y,p in tuples])
    ymin = min([y for x,y,p in tuples])
    ymax = max([y for x,y,p in tuples])
    if xtype == "linear":
        xbins = gen_lin_bins(xmin, xmax, x_num_bins)
    else:
        xbins = gen_log_bins(xmin, xmax, x_num_bins)
    if ytype == "linear":
        ybins = gen_lin_bins(ymin, ymax, y_num_bins)
    else:
        ybins = gen_log_bins(ymin, ymax, y_num_bins)
    sum_bins = NP.zeros((len(xbins), len(ybins)))
    count_bins = NP.zeros((len(xbins), len(ybins)))
    for x,y,p in tuples:
        x_bin = coord_to_bin(x, xbins)
        y_bin = coord_to_bin(y, ybins)
        sum_bins[x_bin][y_bin] += p
        count_bins[x_bin][y_bin] += 1
    for xI in range(len(xbins)):
        for yI in range(len(ybins)):
            if count_bins[xI,yI] == 0:
                sum_bins[xI,yI] = NP.NAN
            else:
                sum_bins[xI,yI] /= count_bins[xI,yI]
    return sum_bins, xbins, ybins

# XXX generalize this to def get_input_omegas etc
def length_error(input, results):
    tbr = []
    for dist_key in input.keys():
        for repI in range(len(results[dist_key])):
            for taxa_key in input[dist_key].keys():
                input_omegas = input[dist_key][taxa_key]["omegas"]
                input_length = input[dist_key][taxa_key]["length"]
                result_length = results[dist_key][repI][taxa_key]["length"]
                difference = float(result_length) - float(input_length)
                #difference = float(input_length) - float(result_length)
                average = (float(input_length) + float(result_length))/2
                # This would be the relative difference between the two
                #relative_difference = (difference/average) * 100.0
                # This is the difference as a percent of expected
                relative_difference = (difference/float(input_length)) * 100.0
                tbr.append((float(input_omegas[0]),
                            float(input_length),
                            relative_difference))
    return tbr

def heatarray_to_heatmap(   heat_array,
                            xbins,
                            ybins,
                            xlabel,
                            ylabel,
                            figure_title,
                            name):
    masked_array = NP.ma.array(heat_array, mask=NP.isnan(heat_array))
    masked_array = NP.transpose(masked_array)
    print(masked_array)
    cmap = mpl.cm.jet
    cmap.set_bad('w',1.)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.imshow( masked_array,
                origin='lower',
                interpolation='none',
                cmap=cmap)
    plt.colorbar()
    #print("xbins:")
    #print(xbins)
    #xbin_strings = [string(x)%0.1f for x in xbins)]

    xbin_strings = ["%0.3f" % x for x in xbins]
    ybin_strings = ["%0.3f" % y for y in ybins]
    plt.xticks(range(len(xbin_strings)), xbin_strings, rotation=90)
    plt.yticks(range(len(ybin_strings)), ybin_strings)

    #plt.gcf().set_size_inches(11.5,10.5)
    plt.gcf().set_size_inches(5,5)
    plt.gcf().subplots_adjust(bottom=.20)
    plt.gcf().subplots_adjust(left=.20)
    plt.gcf().suptitle(figure_title)
    # XXX replace with cur_dir
    plt.savefig(os.path.dirname(os.path.abspath(__file__))
                        + os.sep + name, dpi=300)
    plt.clf()
    plt.cla()

def cur_dir(name):
    return (os.path.dirname(os.path.abspath(__file__)) + os.sep + name)

# Each dist is an array (for the reps) of lists where keys are taxa names and
# values are lists where keys are parameter names and values are lists or
# numbers representing parameter values
def process_csv_results(input, results):
    # pull out meandnds error
    over_one_pvals = []
    false_positive_pvals = []
    true_negative_pvals = []
    dist_errors = []
    for dist_key in input.keys():
        for repI in range(len(results[dist_key])):
            #print(input[dist_key][repI]["1"])
            # results are going to have another list layer for the reps,
            # which the input wont have
            for taxa_key in input[dist_key].keys():
                input_omegas = input[dist_key][taxa_key]["omegas"]
                results_omegas = results[dist_key][repI][taxa_key]["omegas"]
                input_props = input[dist_key][taxa_key]["props"]
                results_props = results[dist_key][repI][taxa_key]["props"]
                input_meandnds = weighted_average(  input_omegas,
                                                    input_props)
                results_meandnds = weighted_average(results_omegas,
                                                    results_props)
                dist_errors.append(math.fabs(   input_meandnds -
                                                results_meandnds))
                if max(input_omegas) > 1:
                    over_one_pvals.append(  (max(input_omegas),
                                            input_props[
                                                input_omegas.index(
                                                    max(input_omegas))],
                                            float(results [dist_key]
                                                    [repI]
                                                    [taxa_key]
                                                    ["pval"])))
                elif float(results[dist_key][repI][taxa_key]["pval"]) <= .05:
                    false_positive_pvals.append((   max(input_omegas),
                                                    input_props[
                                                        input_omegas.index(
                                                            max(input_omegas))],
                                                    float(results   [dist_key]
                                                                    [repI]
                                                                    [taxa_key]
                                                                    ["pval"])))
                else:
                    true_negative_pvals.append((   max(input_omegas),
                                                    input_props[
                                                        input_omegas.index(
                                                            max(input_omegas))],
                                                    float(results   [dist_key]
                                                                    [repI]
                                                                    [taxa_key]
                                                                    ["pval"])))
                #dist_errors.append((input_meandnds, results_meandnds))
    if len(over_one_pvals) != 0:
        print("Over one data (omega, prop, pval):\n", over_one_pvals)
        heat_array, xbins, ybins  = tuple_to_averaged_heatarray(over_one_pvals,
                                                                12,
                                                                12,
                                                                xtype="log")
                                                                #xtype="linear")
        print(heat_array)
        masked_array = NP.ma.array(heat_array, mask=NP.isnan(heat_array))
        masked_array = NP.transpose(masked_array)
        print(masked_array)
        cmap = mpl.cm.jet
        cmap.set_bad('w',1.)
        #plt.figure(figsize=(8,8), dpi=300)
        plt.xlabel('Omega value')
        plt.ylabel('Omega over one proportion')
        #plt.imshow( heat_array,
        plt.imshow( masked_array,
                    origin='lower',
                    interpolation='none',
                    cmap=cmap)
                    #extent=[2.5,4.5,2.5,4.5])
        plt.colorbar()

        xbin_strings = ["%0.3f" % x for x in xbins]
        ybin_strings = ["%0.3f" % y for y in ybins]
        plt.xticks(range(len(xbin_strings)), xbin_strings, rotation=90)
        plt.yticks(range(len(ybin_strings)), ybin_strings)

        #plt.xticks(range(len(xbins)), xbins, rotation=90)
        #plt.yticks(range(len(ybins)), ybins)
        plt.gcf().set_size_inches(5,5)
        plt.gcf().subplots_adjust(bottom=.20)
        plt.gcf().subplots_adjust(left=.20)
        plt.gcf().suptitle( "Mean p-vlaue for Estimates of Positive " \
                            "Selection")
        #ax = plt.gca()
        #for label in ax.xaxis.get_ticklabels():
            #label.set_rotation(45)
        # XXX replace with cur_dir
        plt.savefig('/home/martin/Software/Simulation/meanPval.png', dpi=300)
        plt.clf()

        # True Positives:
        heat_array, xbins, ybins  = tuple_to_true_positive_ratio_heatarray( over_one_pvals,
                                                        12,
                                                        12,
                                                        xtype="log")
                                                        #xtype="linear")
        print("true positive ratio: \n", heat_array)
        masked_array = NP.ma.array(heat_array, mask=NP.isnan(heat_array))
        masked_array = NP.transpose(masked_array)
        print(masked_array)
        cmap = mpl.cm.jet
        cmap.set_bad('w',1.)
        plt.xlabel('Omega value')
        plt.ylabel('Omega over one proportion')
        plt.imshow( masked_array,
                    origin='lower',
                    interpolation='none',
                    cmap=cmap)
                    #extent=[2.5,4.5,2.5,4.5])
        plt.colorbar()

        xbin_strings = ["%0.3f" % x for x in xbins]
        ybin_strings = ["%0.3f" % y for y in ybins]
        plt.xticks(range(len(xbin_strings)), xbin_strings, rotation=90)
        plt.yticks(range(len(ybin_strings)), ybin_strings)

        #plt.xticks(range(len(xbins)), xbins, rotation=90)
        #plt.yticks(range(len(ybins)), ybins)
        plt.gcf().suptitle("True Positive Rate")
        plt.gcf().set_size_inches(5,5)
        plt.gcf().subplots_adjust(bottom=.20)
        plt.gcf().subplots_adjust(left=.20)
        # XXX replace with cur_dir
        plt.savefig('/home/martin/Software/Simulation/truePosRatio.png',
        dpi=300)
        plt.clf()

        # False Positives:
        if len(false_positive_pvals) > 0:
            heat_array, xbins, ybins  = tuple_to_false_positive_ratio_heatarray(
                                            true_negative_pvals,
                                            false_positive_pvals,
                                            12,
                                            12,
                                            xtype="log")
                                            #xtype="linear")
            print("false positives: \n", false_positive_pvals)
            print("false positive ratio: \n", heat_array)
            masked_array = NP.ma.array(heat_array, mask=NP.isnan(heat_array))
            masked_array = NP.transpose(masked_array)
            print(masked_array)
            cmap = mpl.cm.jet
            cmap.set_bad('w',1.)
            plt.xlabel('Omega value')
            plt.ylabel('Omega proportion')
            plt.imshow( masked_array,
                        origin='lower',
                        interpolation='none',
                        cmap=cmap)
                        #extent=[2.5,4.5,2.5,4.5])
            plt.colorbar()

            xbin_strings = ["%0.3f" % x for x in xbins]
            ybin_strings = ["%0.3f" % y for y in ybins]
            plt.xticks(range(len(xbin_strings)), xbin_strings, rotation=90)
            plt.yticks(range(len(ybin_strings)), ybin_strings)

            #plt.xticks(range(len(xbins)), xbins, rotation=90)
            #plt.yticks(range(len(ybins)), ybins)
            plt.gcf().set_size_inches(5,5)
            plt.gcf().subplots_adjust(bottom=.20)
            plt.gcf().subplots_adjust(left=.20)
            # XXX replace with cur_dir
            plt.gcf().suptitle("False Positive Rate")
            plt.savefig('/home/martin/Software/Simulation/falsePosRatio.png',
            dpi=300)
            plt.clf()
    print("Input: ", input)
    print("Results: ", results)
    print("Dist Errors: ", dist_errors)
    dist_mean_error = NP.mean(dist_errors)
    dist_std_error = NP.std(dist_errors)
    print("mean: ", dist_mean_error)
    print("stddev: ", dist_std_error)
    plt.clf()
    plt.cla()
    #print(results)


    # Generate some test data
    #x = NP.random.randn(8873)
    #y = NP.random.randn(8873)

    #heatmap, xedges, yedges = NP.histogram2d(x, y, bins=50)
    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    #plt.clf()
    #plt.imshow(heatmap, extent=extent)
    #plt.show()
    #plt.savefig('/home/martin/Software/Simulation/testFig.png')

    standard_error_meandnds = 0
    tbr = standard_error_meandnds
    return tbr

def process_results(input, results):
    interlaced = {}
    for key in input.keys():
        newlist = []
        newlist.append(input[key])
        newlist.append(results[key])
        interlaced[key] = newlist
    print(interlaced)
    return 0

if len(sys.argv) != 5:
    print(  "Valid usage:\n\t- setgen [number of distributions] " \
            "[number of replicates] [length of sequences] " \
            "[output file name]\n" \
            "\t\t- [number of distributions]: The number of times " \
            "generate_settings is called\n" \
            "\t\t- [number of replicates]: The number of replicates " \
            "HyPhy is told to run\n" \
            "\t\t- [length of sequences]: In codons \n" \
            "\t\t- [output file name]: Will have dist numbers appended\n",
            file=sys.stderr)
    exit(1)
else:
    numDist = int(sys.argv[1]) # the number of times generate_settings is called
    numReps = sys.argv[2] # the number of replicates HyPhy is told to run
    lenSeqs = sys.argv[3] # the number of codons generated
    outFile = sys.argv[4] # the name of the outfile

inputParameterSets = {}
results = {}
results_with_lengths = {}
inputs_with_lengths = {}

sim_time = 0
bsrel_time = 0
python_time = 0

#for dist in range(int(numDist)):
dist_done = 0
while dist_done < numDist:
    start = time.time()
    for this_dist in range(min(numDist-dist_done, len(node))):
        dist_num = this_dist + dist_done
        this_set = generate_settings(num_taxa, outFile +  "." + str(dist_num))
        inputParameterSets[str(dist_num)] = this_set
        simulate(outFile + "." + str(dist_num), this_dist)
    for sub_p in node_processes.values():
        sub_p.wait()
    end = time.time()
    sim_time += end - start
    start = time.time()
    for this_dist in range(min(numDist-dist_done, len(node))):
        dist_num = this_dist + dist_done
        for rep in range(int(numReps)):
            run_BSREL(outFile + "." + str(dist_num), this_dist, rep)
    for sub_p in node_processes.values():
        sub_p.wait()
    end = time.time()
    bsrel_time += end - start
    start = time.time()
    for this_dist in range(min(numDist-dist_done, len(node))):
        dist_num = this_dist + dist_done
        this_sets_results = []
        this_sets_len_results = []
        for rep in range(int(numReps)):
            this_sets_len_results.append(recover_fit(   num_taxa,
                                                        outFile
                                                            + "."
                                                            + str(dist_num),
                                                        this_dist,
                                                        rep))
            this_sets_results.append(recover_csv(   num_taxa,
                                                    outFile + "."
                                                            + str(dist_num),
                                                    this_dist,
                                                    rep))
        results[str(dist_num)] = this_sets_results
        this_sets_len_inputs = recover_simulated(   num_taxa,
                                                    outFile
                                                        + "."
                                                        + str(dist_num),
                                                    this_dist,
                                                    rep)
        inputs_with_lengths[str(dist_num)] = this_sets_len_inputs
        results_with_lengths[str(dist_num)] = this_sets_len_results
    dist_done += min(numDist-dist_done, len(node))
    end = time.time()
    python_time += end - start
#process_results(inputParameterSets, results)

# recover the simulated length params
print("csv results: \n", results)
print("input with lengths: \n", inputs_with_lengths)
# recover the fit file: results_with_lengths
print("results with lengths: \n", results_with_lengths)
# XXX process the lengths
#length_errors = length_error(inputs_with_lengths, results_with_lengths)
length_errors = length_error(inputParameterSets, results)
print("length results:\n", length_errors)
length_heatarray, len_xbins, len_ybins = tuple_to_averaged_heatarray(length_errors, 12, 12)
print("\nLength heat array:\n", length_heatarray)
heatarray_to_heatmap(   length_heatarray,
                        len_xbins,
                        len_ybins,
                        "Omega value",
                        "Simulated Branch Length",
                        "Branch Length Estimation Error Relative to Real " \
                        "Length",
                        "lengthError.png")
start = time.time()
process_csv_results(inputParameterSets, results)
end = time.time()
python_time += end - start

print("Simulation time: ", sim_time)
print("BSREL time: ", bsrel_time)
print("Python time: ", python_time)
