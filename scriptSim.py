#! /usr/bin/env python
# XXX make an output directory an optional input
# XXX make simulation and processing separate programs/libraries

from treegen import get_tree, get_unrooted_tree
from simsetgen import generate_all_settings
from simrun import run_simulation
from bsrelrun import run_all_BSREL
from bsrelSimParsers import (   recover_csv, recover_fit, recover_simulated,
                                recover_settings)
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
# -1 indicates no limit, draw at random
rate_classes_per_branch = -1

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

# XXX break out and make independent the graphing stuff...
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

if len(sys.argv) != 7:
    print(  "Valid usage:\n\t- scriptSim [settings/simulate/bsrel/process] \n" \
            "\t[settings/simulate/bsrel/process] [number of distributions] \n" \
            "\t[number of replicates] [length of sequences] [output file " \
            "name]\n" \
            "\t\t- [settings/simulation/bsrel/process] x2: the first and \n" \
            "\t\t  last steps to perform. Previous data must exist in the \n" \
            "\t\t  current directory. \n " \
            "\t\t- [number of distributions]: The number of times " \
            "generate_settings is called\n" \
            "\t\t- [number of replicates]: The number of replicates " \
            "HyPhy is told to run\n" \
            "\t\t- [length of sequences]: In codons \n" \
            "\t\t- [output file name]: Will have dist numbers appended\n",
            file=sys.stderr)
    exit(1)
else:
    start_step = sys.argv[1]
    end_step = sys.argv[2]
    numDist = int(sys.argv[3]) # the number of times generate_settings is called
    numReps = sys.argv[4] # the number of replicates HyPhy is told to run
    lenSeqs = sys.argv[5] # the number of codons generated
    outFile = sys.argv[6] # the name of the outfile

stages = {  "SETTINGS" : 1,
            "SIMULATE" : 2,
            "BSREL" : 3,
            "PROCESS" : 4 }

start_step = stages [start_step.upper()]
end_step = stages [end_step.upper()]
print("start: ", start_step)
print("end: ", end_step)

results = {}
results_with_lengths = {}
inputs_with_lengths = {}

sim_time = 0
bsrel_time = 0
python_time = 0

if start_step <= 1 and end_step >=1:
    print("writing settings")
    # Write settings
    start = time.time()
    generate_all_settings(  num_taxa,
                            numDist,
                            numReps,
                            lenSeqs,
                            outFile,
                            rate_classes_per_branch)
    end = time.time()
    sim_time += end - start

if start_step <= 2 and end_step >=2:
    print("running simulation")
    # Run simulator
    start = time.time()
    run_simulation(numDist, "", outFile, node_processes, node)
    end = time.time()
    sim_time += end - start

if start_step <= 3 and end_step >=3:
    print("running bsrel")
    # Run BSREL
    # XXX maybe include an option to specify these as input?
    # e.g., have 100 simulations, run BSREL on 3...
    start = time.time()
    run_all_BSREL(  numDist,
                    outFile,
                    numReps,
                    node_processes,
                    node)
    end = time.time()
    bsrel_time += end - start

if start_step <= 4 and end_step >=4:
    print("processing results")
    # Run BSREL
    # Run BSREL
    # XXX Process results
    # XXX pull out to library

    # XXX 3. separate processing and graphing
    # XXX 4. Break graphing out to library
    # XXX 5. allow specification of steps to run as argument

    # XXX if prev steps aren't done infer necessary inputs
    # XXX make a bunch of parsers in a parsing library

    parsedInputs = {}
    for dist in range(numDist):
        parsedInputs[str(dist)] = recover_settings(outFile + "." + str(dist))

    start = time.time()
    dist_done = 0
    for dist_num in range(numDist):
        this_sets_results = []
        this_sets_len_results = []
        for rep in range(int(numReps)):
            this_sets_len_results.append(recover_fit(   num_taxa,
                                                        outFile
                                                            + "."
                                                            + str(dist_num),
                                                        rep))
            this_sets_results.append(recover_csv(   outFile + "."
                                                            + str(dist_num),
                                                    rep))
        results[str(dist_num)] = this_sets_results
        this_sets_len_inputs = recover_simulated(   num_taxa,
                                                    outFile
                                                        + "."
                                                        + str(dist_num),
                                                    rep)
        inputs_with_lengths[str(dist_num)] = this_sets_len_inputs
        results_with_lengths[str(dist_num)] = this_sets_len_results
    end = time.time()
    python_time += end - start

    # recover the simulated length params
    #print("csv results: \n", results)
    #print("input with lengths: \n", inputs_with_lengths)
    # recover the fit file: results_with_lengths
    #print("results with lengths: \n", results_with_lengths)

    length_errors = length_error(parsedInputs, results)
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
    process_csv_results(parsedInputs, results)
    end = time.time()
    python_time += end - start

print("Simulation time: ", sim_time)
print("BSREL time: ", bsrel_time)
print("Python time: ", python_time)
