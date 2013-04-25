#! /usr/bin/env python

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
    entry += str(param_list["len"])
    entry += ",\n"
    entry += "\t\"omegas\" : {"
    for o,p in zip(param_list["omegas"], param_list["props"]):
        entry += "\t{ " + str(o) + ", " + str(p) + "}\n"
    entry += "}\n},\n"
    return entry

def generate_taxa(tax_name, num_rate_classes):
    entry = {}
    entry["name"] = str(tax_name)
    # XXX re-enable
    #len = NP.random.exponential(.25, 1)[0]
    length = 1
    entry["len"] = length
    omegas = sorted([NP.random.exponential((.5 + (10*i))) for i in range(num_rate_classes)])

    # check omegas over one:
    under_one = False
    num_over_one = 0
    for omega in omegas:
        if omega < 1.0:
            under_one == True
        if omega > 1.0:
            num_over_one += 1
    if under_one == False:
        omegas[0] = omegas[0] % 1
    if num_over_one > 1:
        for omegaI in range(len(omegas) - 1):
            omegas[omegaI] = omegas[omegaI] % 1

    #omegas = [.5*i for i in range(1,num_rate_classes + 1)]
    entry["omegas"] = omegas
    props = []
    remainingProp = 1
    for i in range(num_rate_classes - 1):
        props.append(max(0.1, (1 - NP.random.exponential(.5, 1)[0])) * remainingProp)
        #props.append(0.9 * remainingProp)
        remainingProp -= props[i]
    props.append(1 - sum(props))
    props = sorted(props)[::-1]
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
                        + '.sim.'
                        + str(nodeI)
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
                    '("/usr/local/lib/hyphy/TemplateBatchFiles/BranchSiteREL.bf"' \
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
    resultsfile = open(rec_file_name + ".sim." + str(rep) + ".recovered.fit", 'r')
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
        results[taxa] = format_results(taxa, results[taxa])
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
        results[line[0]]["len"] = line[-1]
        results[line[0]]["omegas"] = []
        results[line[0]]["props"] = []
        over_one = 0
        results[line[0]]["pval"] = line[7]
        if (line[3]) != "0":
            results[line[0]]["omegas"].append(line[3])
            results[line[0]]["props"].append(line[4])
            over_one = 1
        for oI in range(1, (int(line[2]) + 1 - int(over_one))):
            results[line[0]]["omegas"].append(str(line[1]))
            results[line[0]]["props"].append(str(   (1 - float(line[4]))
                                                    /float(line[2])))
    return results

def format_results(taxa, taxa_dict):
    entry = {}
    entry["name"] = taxa
    entry["omegas"] = []
    entry["props"] = []
    omega = False
    for key in taxa_dict.keys():
        if key == 't':
            omega = True
    if omega:
        entry["len"] = taxa_dict["t"]
        for key in taxa_dict.keys():
            if key[:5] == "omega":
                entry["omegas"].append(taxa_dict[key])
            if key[:4] == "Paux":
                entry["props"].append(taxa_dict[key])
        entry["props"] = convolve_props(entry["props"])
    else:
        entry["len"] = taxa_dict['syn']
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

def params_to_meandnds(omegas, props):
    return sum([float(omegas[i])*float(props[i]) for i in range(len(omegas))])

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
    return [min+(i*bin_size) for i in range(num)]

#def gen_log_bins(min, max, num):
#    next_bin = .00000000001
#    while next_bin < min:
#        next_bin *= 10
#    next_bin /= 10
#    while next_bin < max:
#        yield next_bin
#        next_bin *= 10

def gen_log_bins(min, max, num):
    log_min = math.log10(min)
    log_max = math.log10(max)
    bin_size = (log_max - log_min)/num
    return [10**(log_min + (i * bin_size)) for i in range(num)]

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

def tuple_to_heatarray( data,
                        x_num_bins,
                        y_num_bins,
                        xmax=0,
                        xmin=1,
                        ymax=0,
                        ymin=1,
                        xtype="linear",
                        ytype="linear"):
    xmin = min([x for x,y,p in data])
    xmax = max([x for x,y,p in data])
    ymin = min([y for x,y,p in data])
    ymax = max([y for x,y,p in data])
    if xtype == "linear":
        xbins = gen_lin_bins(xmin, xmax, x_num_bins)
    else:
        xbins = gen_log_bins(xmin, xmax, x_num_bins)
    if ytype == "linear":
        ybins = gen_lin_bins(ymin, ymax, y_num_bins)
    else:
        ybins = gen_log_bins(ymin, ymax, y_num_bins)
    print("xbins: \n", xbins)
    heatarray = NP.ones((len(xbins), len(ybins)))
    heatarray[:] = 99
    heatcount = NP.zeros((len(xbins), len(ybins)))
    #print(ybins)
    for x,y,p in data:
        x_bin = coord_to_bin(x, xbins)
        y_bin = coord_to_bin(y, ybins)
        current_val = heatcount[x_bin][y_bin] * heatarray[x_bin][y_bin]
        heatcount[x_bin][y_bin] += 1
        new_val = (current_val + p)/heatcount[x_bin][y_bin]
        heatarray[x_bin][y_bin] = new_val
        #heatarray[coord_to_bin(x, xbins)][coord_to_bin(y, ybins)] = p

    for xI in range(len(heatarray[:,1])):
        for yI in range(len(heatarray[1,:])):
            if heatarray[xI][yI] == 99:
                heatarray[xI][yI] = NP.NAN
    return heatarray, xbins, ybins

# XXX ie correlate the results and degree of error with characteristics of
# the input
# Each dist is an array (for the reps) of lists where keys are taxa names and
# values are lists where keys are parameter names and values are lists or
# numbers representing parameter values
def process_cvs_results(input, results):
    # pull out meandnds error
    over_one_pvals = []
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
                input_meandnds = params_to_meandnds(    input_omegas,
                                                        input_props)
                results_meandnds = params_to_meandnds(  results_omegas,
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
                #dist_errors.append((input_meandnds, results_meandnds))
    if len(over_one_pvals) != 0:
        print("Over one data (omega, prop, pval):\n", over_one_pvals)
        heat_array, xbins, ybins  = tuple_to_heatarray( over_one_pvals,
                                                        20,
                                                        20,
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
        plt.xticks(range(len(xbins)), xbins, rotation=90)
        plt.yticks(range(len(ybins)), ybins)
        plt.gcf().set_size_inches(11.5,10.5)
        plt.gcf().subplots_adjust(bottom=.20)
        plt.gcf().subplots_adjust(left=.20)
        #ax = plt.gca()
        #for label in ax.xaxis.get_ticklabels():
            #label.set_rotation(45)
        plt.savefig('/home/martin/Software/Simulation/testFig.png', dpi=150)
    print("Dist Errors: ", dist_errors)
    dist_mean_error = NP.mean(dist_errors)
    dist_std_error = NP.std(dist_errors)
    print("mean: ", dist_mean_error)
    print("stddev: ", dist_std_error)
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
        for rep in range(int(numReps)):
            #this_sets_results.append(recover_fit(num_taxa, outFile + "." + str(dist), dist, rep))
            this_sets_results.append(recover_csv(   num_taxa,
                                                    outFile + "."
                                                            + str(dist_num),
                                                    this_dist,
                                                    rep))
        results[str(dist_num)] = this_sets_results
    dist_done += min(numDist-dist_done, len(node))
    end = time.time()
    python_time += end - start
#process_results(inputParameterSets, results)
start = time.time()
process_cvs_results(inputParameterSets, results)
end = time.time()
python_time += end - start

print("Simulation time: ", sim_time)
print("BSREL time: ", bsrel_time)
print("Python time: ", python_time)
