#! usr/bin/env python

import os, sys


def get_lines(file_name):
    fh = open(file_name, 'r')
    lines = fh.readlines()
    fh.close()
    if len(lines) > 1:
        lines[-1] = lines[-1] + "\n"
        return lines[0], lines[1:]
    else:
        return "", []

def write_lines(header, combined_lines):
    fh = open("combined.csv", 'w')
    fh.write(header)
    fh.writelines(combined_lines)
    fh.close()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        dir_name = sys.argv[1]
        header = ""
        combined_lines = []
        file_list = []
        for path, dirs, files in os.walk(dir_name):
            temp_files = [os.path.join(dir_name,f) for f in files]
            file_list = file_list + temp_files
        for fn in file_list:
            new_header, new_lines = get_lines(fn)
            if new_header != "":
                header = new_header
            combined_lines = combined_lines + new_lines
        write_lines(header, combined_lines)
    else:
        print("Please include an input directory name", file=sys.stderr)
        exit(1)

