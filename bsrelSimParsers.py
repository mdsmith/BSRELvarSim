
# XXX make this a reality (recover simulated doesn't recover the input
# lengths)
def recover_settings():
    return 0

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
