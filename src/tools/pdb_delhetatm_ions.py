#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This script is an adaptation of the `pdb-tools` suite's pdb_delhetatm.py

"""
Removes all HETATM Ion records in the PDB file.

Usage:
    python pdb_delhetatm_ions.py <pdb file>

Example:
    python pdb_delhetatm_ions.py 1CTF.pdb

"""

import os
import sys
import re


def check_input(args):
    """Checks whether to read from stdin/file and validates user input/options.
    """

    # Defaults
    fh = sys.stdin  # file handle

    if not len(args):
        # Reading from pipe with default option
        if sys.stdin.isatty():
            sys.stderr.write(__doc__)
            sys.exit(1)

    elif len(args) == 1:
        if not os.path.isfile(args[0]):
            emsg = 'ERROR!! File not found or not readable: \'{}\'\n'
            sys.stderr.write(emsg.format(args[0]))
            sys.stderr.write(__doc__)
            sys.exit(1)

        fh = open(args[0], 'r')

    else:  # Whatever ...
        emsg = 'ERROR!! Script takes 1 argument, not \'{}\'\n'
        sys.stderr.write(emsg.format(len(args)))
        sys.stderr.write(__doc__)
        sys.exit(1)

    return fh


def run(fhandle):
    """
    Remove all HETATM Ions and associated records from the PDB file.

    This function is a generator.

    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.

    Yields
    ------
    str (line-by-line)
        The modified (or not) PDB line.
    """

    # CONECT 1179  746 1184 1195 1203
    char_ranges = (slice(6, 11), slice(11, 16),
                   slice(16, 21), slice(21, 26), slice(26, 31))

    het_serials = set()
    for line in fhandle:
        # In case of ions, as most of them have two-letter names
        # the ion atom name starts at position 12 all other atom names start at position 13
        # Exceptions are potassium (K), fluoride (F), vanadium (V), and tungsten (W) ions.
        # To avoid deleting two-letter atoms from ligands, check that atom name equals (part of) residue name.
        if line.startswith('HETATM') and (((line[12].isalpha()) and ((line[12:14] == line[18:20]) or (line[12:14] == line[17:19]))) or
                                          (line[13:20] == "K     K") or
                                          (line[13:20] == "F     F") or
                                          (line[13:20] == "V     V") or
                                          (line[13:20] == "W     W")):
            het_serials.add(line[6:11])
            continue
        elif line.startswith('ANISOU'):
            if line[6:11] in het_serials:
                continue
        elif line.startswith('CONECT'):
            if any(line[cr] in het_serials for cr in char_ranges):
                continue

        yield line


remove_hetatm = run


def main():
    # Check Input
    pdbfh = check_input(sys.argv[1:])

    # Do the job
    new_pdb = run(pdbfh)

    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                sys.stdout.write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        sys.stdout.write(''.join(_buffer))
        sys.stdout.flush()
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

    # last line of the script
    # We can close it even if it is sys.stdin
    pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
