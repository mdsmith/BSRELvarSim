#! /usr/bin/env python

import sys

def get_tree(num_taxa):
    result = gen_newick([a for a in range(1, int(num_taxa) + 1)],
                        (2*int(num_taxa)) - 1)
    return  ("(" + str(result) + ")")

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
    result = gen_unrooted_newick(   [a for a in range(1, int(num_taxa) + 1)],
                                    (2*int(num_taxa))-1,
                                    (2*int(num_taxa))-3)
    return("(" + result + ")")

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


if __name__ == "__main__":
    if len(sys.argv) == 2:
        numBranches = sys.argv[1]
        type = "unrooted"
    elif len(sys.argv) == 3:
        numBranches = sys.argv[1]
        type = sys.argv[2]
    else:
        print(  "Valid usage:\n\t- treegen <number of taxa> [<type>]" \
                "\n\t\t-<type>: rooted or unrooted",
                file=sys.stderr)
        exit(1)
    if type == "rooted":
        print(get_tree(numBranches))
    elif type == "unrooted":
        print(get_unrooted_tree(numBranches))
