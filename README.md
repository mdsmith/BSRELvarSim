BSRELvarSim is a package for simulating and analyzing sequences as part of
the BSREL validation process.

Requirements:
    - Python 3

There are a number of ways to use it and its various components, but all of
them can be addressed through the scriptSim program. This will pull elements
from the various modules included in BSRELvarSim to do as many or as few of
the simulating and processing steps as you desire.

scriptSim.py works as follows:

    Valid usage:
        - scriptSim [settings/simulate/bsrel/process]
        [settings/simulate/bsrel/process] [number of distributions]
        [number of replicates] [length of sequences] [output file name]
            - [settings/simulation/bsrel/process] x2: the first and
              last steps to perform. Previous data must exist in the
              current directory.
            - [number of distributions]: The number of times generate_settings is called
            - [number of replicates]: The number of replicates HyPhy is told to run
            - [length of sequences]: In codons
            - [output file name]: Will have dist numbers appended

The various modules can be used independently, to run bsrel on directories of
data, simulate sequences from settings files, or generate the settings files
themselves.
