#! /usr/bin/env python

# Take a CSV file, fit file, etc for one taxa and convolve them as a single
# CSV file.

import csv
import os
import sys
from bsrelSimParsers import (   recover_csv, recover_fit, recover_simulated,
                                recover_settings, recover_csv_mg94)

#def append_fit(buffer, filename):
#def append_simulated(buffer, filename):

def meandnds(omegas, props):
    mean = 0
    for omega, prop in zip(omegas, props):
        mean += float(omega) * float(prop)
    return mean

def append_BSREL3(buffer, filename):
    contents = recover_csv(filename)
    branch_order = [line[0] for line in buffer[1:]]
    try:
        len_column = rep_to_column(contents, "length", branch_order)
    except KeyError:
        print("Broken file: ", filename)
        return buffer

    len_column[0] = "BSREL3_length"
    buffer = append_column(buffer, len_column)

    #print(contents)
    omegas_column = rep_to_column(contents, "omegas", branch_order)
    props_column = rep_to_column(contents, "props", branch_order)
    omegas_column[0] = "BSREL3_meandnds"
    props_column[0] = "BSREL3_propOverOne"

    mean_omegas_column = [  meandnds(omegas, props)
                            for omegas,props in
                            zip(omegas_column[1:], props_column[1:])]
    mean_omegas_column.insert(0, "BSREL3_meandns")
    #print(mean_omegas_column)

    omega_over_one_column = [   max([float(o) for o in omegas])
                                if max([float(o) for o in omegas]) > 1
                                else 0
                                for omegas in omegas_column[1:]]
    omega_over_one_column.insert(0, "BSREL3_OmegaOver1")
    #print(omega_over_one_column)

    prop_over_one_column = [props[len(props)-1]
                            if omegas[len(omegas)-1] > 1
                            else 0
                            for omegas, props in zip(omegas_column[1:],
                            props_column[1:])]
    prop_over_one_column.insert(0, "BSREL3_propOverOne")
    #print(prop_over_one_column)

    max_omega_column = [omegas[-1]
                        for omegas in omegas_column[1:]];
    max_omega_column.insert(0, "BSREL3_MaxOmega")
    #print(max_omega_column)

    max_prop_column = [props[-1]
                        for props in props_column[1:]];
    max_prop_column.insert(0, "BSREL3_MaxOmegaProp")
    #print(max_prop_column)

    buffer = append_column(buffer, omega_over_one_column)
    buffer = append_column(buffer, prop_over_one_column)
    buffer = append_column(buffer, mean_omegas_column)
    buffer = append_column(buffer, max_omega_column)
    buffer = append_column(buffer, max_prop_column)
    #print(buffer)
    return buffer

def append_MG94(buffer, filename):
# XXX may need to be configured to work with incomplete csv
    try:
        contents = recover_csv_mg94(filename)
    except IndexError:
        print("Broken file: ", filename)
        return buffer
    branch_order = [line[0] for line in buffer[1:]]
    len_column = rep_to_column(contents, "length", branch_order)
    len_column[0] = "MG94_length"
    buffer = append_column(buffer, len_column)

    #print(contents)
    omegas_column = rep_to_column(contents, "omegas", branch_order)
    omegas_column = [omegas[0][0] for omegas in omegas_column]
    omegas_column[0] = "MG94_meandnds"

    #mean_omegas_column = [  meandnds(omegas, props)
                            #for omegas,props in
                            #zip(omegas_column[1:], props_column[1:])]
    #mean_omegas_column.insert(0, "MG94_meandns")
    #print(mean_omegas_column)

    buffer = append_column(buffer, omegas_column)
    #print(buffer)
    return buffer

def append_settings(buffer, filename):
    contents = recover_settings(filename)
    num_rows = len(buffer) - 1
    branch_order = [line[0] for line in buffer[1:]]
    len_column = rep_to_column(contents, "length", branch_order)
    len_column[0] = "Settings_length"
    buffer = append_column(buffer, len_column)

    #print(contents)
    omegas_column = rep_to_column(contents, "omegas", branch_order)
    props_column = rep_to_column(contents, "props", branch_order)
    omegas_column[0] = "Settings_meandnds"
    props_column[0] = "Settings_propOverOne"

    mean_omegas_column = [  meandnds(omegas, props)
                            for omegas,props in
                            zip(omegas_column[1:], props_column[1:])]
    mean_omegas_column.insert(0, "Settings_meandns")
    #print(mean_omegas_column)

    omega_over_one_column = [   max([float(o) for o in omegas])
                                if max([float(o) for o in omegas]) > 1
                                else 0
                                for omegas in omegas_column[1:]]
    omega_over_one_column.insert(0, "Settings_OmegaOver1")
    #print(omega_over_one_column)

    prop_over_one_column = [props[len(props)-1]
                            if omegas[len(omegas)-1] > 1
                            else 0
                            for omegas, props in zip(omegas_column[1:],
                            props_column[1:])]
    prop_over_one_column.insert(0, "Settings_propOverOne")
    #print(prop_over_one_column)

    max_omega_column = [omegas[-1]
                        for omegas in omegas_column[1:]];
    max_omega_column.insert(0, "Settings_MaxOmega")
    #print(max_omega_column)

    max_prop_column = [props[-1]
                        for props in props_column[1:]];
    max_prop_column.insert(0, "Settings_MaxOmegaProp")
    #print(max_prop_column)

    buffer = append_column(buffer, omega_over_one_column)
    buffer = append_column(buffer, prop_over_one_column)
    buffer = append_column(buffer, mean_omegas_column)
    buffer = append_column(buffer, max_omega_column)
    buffer = append_column(buffer, max_prop_column)
    #print(buffer)
    return buffer

def append_csv(buffer, filename):
    file = open(filename, 'r')
    contents = file.readlines()
    buffer += contents
    return buffer

def append_column(buffer, column):
    # strip newlines, add comma, replace newlines
    for i, line in enumerate(buffer):
        buffer[i] = line.strip("\n") + "," + str(column[i]) + "\n"
    return buffer

def concat_buffers(buffer1, buffer2):
    if len(buffer1) != 0:
        buffer2 = buffer2[1:]
    return buffer1 + buffer2

def rep_to_column(rep, key, order):
    column = [key]
    values = [rep[branch][key] for branch in order]
    column += values
    return column

def rep_to_csv(rep):
    csv = [[]]
    header = []
    header.append("Branch")
    for key, value in rep.items():
        row = []
        row.append(key)
        for branch_key, branch_value in value.items():
            if header.count(branch_key) == 0:
                header.append(branch_key)
                row.append(branch_value)
            else:
                row.insert(header.index(branch_key), branch_value)
        csv.append(row)
    csv.insert(0, header)
    return csv

def write_buffer(buffer, filename):
    file = open(filename, 'w')
    file.writelines(buffer)
    return 0

# prefixes are a list of filenames that start sets "longPy.279" etc.
def run_batch(buffer, prefixes):
    # for fileset in prefixes
    for prefix in prefixes:
        buffer2 = []
        csv_filename = prefix + ".sim.0.recovered"
        settings_filename = prefix
        append_csv(buffer2, csv_filename)
        append_settings(buffer2, settings_filename)
        append_BSREL3(buffer2, settings_filename + ".sim.0.recovered.BSREL")
        append_MG94(buffer2, settings_filename + ".sim.0.recovered.mglocal.csv")
        #append_csv(buffer2, settings_filename + ".mglocal.csv")
        buffer = concat_buffers(buffer, buffer2)
    return buffer

def get_prefixes(sim_dir):
    import glob
    import re
    file_list = glob.glob(sim_dir + os.sep + "*")
    file_list = [a for a in file_list if re.search("^\w+\/\w+\.\d+$", a) != None]
    return file_list

# XXX Batch this so it can run all at once and concatenate them...
if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "-d":
        # Do dir stuff
        sim_dir = sys.argv[2]
        output_filename = sys.argv[3]
        buffer = []
        #prefixes = [os.path.join(sim_dir, "longPy." + str(number)) for number in
        #range(0,10000)]
        prefixes = get_prefixes(sim_dir)
        buffer = run_batch(buffer, prefixes)
        write_buffer(buffer, output_filename)
    elif len(sys.argv) == 4:
        csv_filename = sys.argv[1]
        settings_filename = sys.argv[2]
        output_filename = sys.argv[3]

        buffer = []
        append_csv(buffer, csv_filename)
        append_settings(buffer, settings_filename)
        # append_settings(buffer, settings_filename)
        # append_simulated(buffer, settings_filename)
        # append_fit(buffer, settings_filename)
        write_buffer(buffer, output_filename)
    else:
        print(  "Valid usage:\n\t- bsrelSimCSVconvolve.py csv_filename " \
                "settings_filename output_filename\n" \
                "\t- Or bsrelSimCSVconvolve.py -d in_dir out_file\n",
                file = sys.stderr)
        exit(1)
