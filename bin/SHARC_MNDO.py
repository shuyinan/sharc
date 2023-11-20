#!/usr/bin/env python3

# ******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2023 University of Vienna
#    Copyright (c) 2023 University of Minnesota
#
#    This file is part of SHARC.
#
#    SHARC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SHARC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
#
# ******************************************

#    ====================================================================
# ||                                                                       ||
# ||                General Remarks                                        ||
# ||                                                                       ||
#    ====================================================================
#
# This script uses several different specification for the electronic states under consideration.
# Generally, the input specs are like "3 Singlets, 0 Doublets and 3 Triplets"
# Based on this information, the states may be denoted by different systems.
#
# The most comprehensive denotation is: (mult, state, MS).
# In this case, states of higher multiplicity show up several times and the total number of states may be quite large.
# This total number of states is called nmstates in the script (equals to 12 for the above example)
#
# Since the MS components of a multiplet often share several expectation values, these need to be calculated only once.
# In this case, states can be safely denoted by (mult,state).
# The total number of these states is called nstates.
#
# In both systems, the states can be given indices. In this script and in SHARC, one first counts up the state quantum number,
# then the MS quantum number, and finally the multiplicity, like in this code snippet:
#
# i=0
# for mult in range(len(states)):
#     for MS in range(mult):
#         for state in range(states[mult]):
#             i+=1
#             print i, mult+1, state+1, MS-i/2
#
# more about this below in the docstrings of the iterator functions

# ======================================================================= #

# IMPLEMENTATION OF ADDITIONAL TASKS KEYWORDS, JOBS, ETC:
#
# A new task keyword in QMin has to be added to:
#             - readQMin (for consistency check)
#             - gettasks (planning the MNDO calculation)
#             - print QMin (optional)
#
# A new task in the Tasks array needs changes in:
#             - gettasks
#             - writeMNDOinput
#             - redotasks
#             - printtasks

# ======================================================================= #
# Modules:
# Operating system, isfile and related routines, move files, create directories
import os
import shutil
# External Calls to MNDO
import subprocess as sp
# Command line arguments
import sys
# Regular expressions
import re
# debug print for dicts and arrays
import pprint
# sqrt and other math
import math
import cmath
# runtime measurement
import datetime
# copy of arrays of arrays
from copy import deepcopy
# parallel calculations
from multiprocessing import Pool
import time
from socket import gethostname
import itertools
# write debug traces when in pool threads
import traceback
# for diabatization
import numpy as np



# ======================================================================= #

version = '3.0'
versiondate = datetime.date(2023, 10, 29)



changelogstring = '''
10.29.2023:
- creation of SHARC_MNDO.py interface, this is modified from SHARC_MOLCAS.py interface
- only semiempirical molecular orbital (SMO) followed by GUGA-CI is implemented
'''

# ======================================================================= #
# holds the system time when the script was started
starttime = datetime.datetime.now()

# global variables for printing (PRINT gives formatted output, DEBUG gives raw output)
DEBUG = True
PRINT = True

# hash table for conversion of multiplicity to the keywords used in MNDO
IToMult = {
    1: 'Singlet',
    2: 'Doublet',
    3: 'Triplet',
    4: 'Quartet',
    5: 'Quintet',
    6: 'Sextet',
    7: 'Septet',
    8: 'Octet',
    'Singlet': 1,
    'Doublet': 2,
    'Triplet': 3,
    'Quartet': 4,
    'Quintet': 5,
    'Sextet': 6,
    'Septet': 7,
    'Octet': 8
}

# hash table for conversion of polarisations to the keywords used in MNDO
IToPol = {
    0: 'X',
    1: 'Y',
    2: 'Z',
    'X': 0,
    'Y': 1,
    'Z': 2
}

# periodic table. 
periodic_table = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca',
    21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr',
    41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd',
    61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb',
    71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th',
    91: 'Pa', 92: 'U'
}

# conversion factors
au2a = 0.529177211
rcm_to_Eh = 4.556335e-6
au2eV = 27.21138386
au2kcal = 627.509608
au2debye = 2.5417469

# =============================================================================================== #
# =============================================================================================== #
# =========================================== general routines ================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #


def readfile(filename):
    try:
        f = open(filename)
        out = f.readlines()
        f.close()
    except IOError:
        print('File %s does not exist!' % (filename))
        sys.exit(12)
    return out

# ======================================================================= #


def writefile(filename, content):
    # content can be either a string or a list of strings
    try:
        f = open(filename, 'w')
        if isinstance(content, list):
            for line in content:
                f.write(line)
        elif isinstance(content, str):
            f.write(content)
        else:
            print('Content %s cannot be written to file!' % (content))
        f.close()
    except IOError:
        print('Could not write to file %s!' % (filename))
        sys.exit(13)

# ======================================================================= #


def link(PATH, NAME, crucial=True, force=True):
    # do not create broken links
    if not os.path.exists(PATH):
        print('Source %s does not exist, cannot create link!' % (PATH))
        sys.exit(14)
    if os.path.islink(NAME):
        if not os.path.exists(NAME):
            # NAME is a broken link, remove it so that a new link can be made
            os.remove(NAME)
        else:
            # NAME is a symlink pointing to a valid file
            if force:
                # remove the link if forced to
                os.remove(NAME)
            else:
                print('%s exists, cannot create a link of the same name!' % (NAME))
                if crucial:
                    sys.exit(15)
                else:
                    return
    elif os.path.exists(NAME):
        # NAME is not a link. The interface will not overwrite files/directories with links, even with force=True
        print('%s exists, cannot create a link of the same name!' % (NAME))
        if crucial:
            sys.exit(16)
        else:
            return
    os.symlink(PATH, NAME)

# ======================================================================= #


def isbinary(path):
    return (re.search(r':.* text', sp.Popen(["file", '-L', path], stdout=sp.PIPE).stdout.read()) is None)

# ======================================================================= #


def eformat(f, prec, exp_digits):
    '''Formats a float f into scientific notation with prec number of decimals and exp_digits number of exponent digits.

    String looks like:
    [ -][0-9]\\.[0-9]*E[+-][0-9]*

    Arguments:
    1 float: Number to format
    2 integer: Number of decimals
    3 integer: Number of exponent digits

    Returns:
    1 string: formatted number'''

    s = "% .*e" % (prec, f)
    mantissa, exp = s.split('e')
    return "%sE%+0*d" % (mantissa, exp_digits + 1, int(exp))

# ======================================================================= #


def measuretime():
    '''Calculates the time difference between global variable starttime and the time of the call of measuretime.

    Prints the Runtime, if PRINT or DEBUG are enabled.

    Arguments:
    none

    Returns:
    1 float: runtime in seconds'''

    endtime = datetime.datetime.now()
    runtime = endtime - starttime
    if PRINT or DEBUG:
        hours = runtime.seconds // 3600
        minutes = runtime.seconds // 60 - hours * 60
        seconds = runtime.seconds % 60
        print('==> Runtime:\n%i Days\t%i Hours\t%i Minutes\t%i Seconds\n\n' % (runtime.days, hours, minutes, seconds))
    total_seconds = runtime.days * 24 * 3600 + runtime.seconds + runtime.microseconds // 1.e6
    return total_seconds

# ======================================================================= #


def removekey(d, key):
    '''Removes an entry from a dictionary and returns the dictionary.

    Arguments:
    1 dictionary
    2 anything which can be a dictionary keyword

    Returns:
    1 dictionary'''

    if key in d:
        r = dict(d)
        del r[key]
        return r
    return d

# ======================================================================= #         OK


def containsstring(string, line):
    '''Takes a string (regular expression) and another string. Returns True if the first string is contained in the second string.

    Arguments:
    1 string: Look for this string
    2 string: within this string

    Returns:
    1 boolean'''

    a = re.search(string, line)
    if a:
        return True
    else:
        return False



# =============================================================================================== #
# =============================================================================================== #
# ============================= iterator routines  ============================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def itmult(states):

    for i in range(len(states)):
        if states[i] < 1:
            continue
        yield i + 1
    return

# ======================================================================= #


def itnmstates(states):

    for i in range(len(states)):
        if states[i] < 1:
            continue
        for k in range(i + 1):
            for j in range(states[i]):
                yield i + 1, j + 1, k - i / 2.
    return

# =============================================================================================== #
# =============================================================================================== #
# =========================================== print routines ==================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #


def printheader():
    '''Prints the formatted header of the log file. Prints version number and version date

    Takes nothing, returns nothing.'''

    print(starttime, gethostname(), os.getcwd())
    if not PRINT:
        return
    string = '\n'
    string += '  ' + '=' * 80 + '\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * 28 + 'SHARC - MNDO - Interface' + ' ' * 28 + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * 32 + 'Author: Yinan Shu' + ' ' * 31 + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '||' + ' ' * (36 - (len(version) + 1) // 2) + 'Version: %s' % (version) + ' ' * (35 - (len(version)) // 2) + '||\n'
    lens = len(versiondate.strftime("%d.%m.%y"))
    string += '||' + ' ' * (37 - lens // 2) + 'Date: %s' % (versiondate.strftime("%d.%m.%y")) + ' ' * (37 - (lens + 1) // 2) + '||\n'
    string += '||' + ' ' * 80 + '||\n'
    string += '  ' + '=' * 80 + '\n\n'
    print(string)
    if DEBUG:
        print(changelogstring)

# ======================================================================= #


def printQMin(QMin):

    if DEBUG:
        pprint.pprint(QMin)
    if not PRINT:
        return
    print('==> QMin Job description for:\n%s' % (QMin['comment']))

    string = 'Tasks:  '
    if 'h' in QMin:
        string += '\tH'
    if 'soc' in QMin:
        string += '\tSOC'
    if 'dm' in QMin:
        string += '\tDM'
    if 'grad' in QMin:
        string += '\tGrad'
    if 'nacdr' in QMin:
        string += '\tNac(ddr)'
    if 'nacdt' in QMin:
        string += '\tNac(ddt)'
    if 'overlap' in QMin:
        string += '\tOverlaps'
    if 'angular' in QMin:
        string += '\tAngular'
    if 'ion' in QMin:
        string += '\tDyson norms'
    if 'dmdr' in QMin:
        string += '\tDM-Grad'
    if 'socdr' in QMin:
        string += '\tSOC-Grad'
    if 'phases' in QMin:
        string += '\tPhases'
    print(string)

    string = 'States: '
    for i in itmult(QMin['states']):
        string += '\t%i %s' % (QMin['states'][i - 1], IToMult[i])
    print(string)

    string = 'Semiempirical Molecular Orbital (SMO) Method: %s\n' % (QMin['template']['mo_method'])
    string += 'maximum SMO SCF iterations: %i\n' % (QMin['template']['maxiter'])
    if QMin['template']['parameter']=='no':
        string += 'using default parameters in MNDO'
    elif QMin['template']['parameter']=='read':
        string += 'read parameters from fort.14 and si.dat'

    string += 'Correlation and Excited State Method: %s\n' % (QMin['template']['corr_method'])
    if QMin['template']['corr_method']=='GUGA-CI':
        string += 'number of occupied active orbital: %i\n' % (QMin['template']['nacocc'])
        string += 'number of virtual active orbital: %i\n' % (QMin['template']['nacvir'])
        string += 'excitation level: %i\n' % (QMin['template']['exlvl'])
        string += 'number of reference CSFs: %i\n' % (QMin['template']['nref'])
        if QMin['template']['refscf_method']=='auto':
            string += 'reference CSFs are determined automatically'
        elif QMin['template']['refscf_method']=='read': 
            string += 'reference CSFs are read from input:'
            for i in range(len(QMin['template']['ref_csf'])):
                for j in range(len(QMin['template']['ref_csf'][i])):
                    string += f"  {QMin['template']['ref_csf'][i][j]}"
                string += "\n"
    print(string)

    string = 'Found Geo'
    if 'veloc' in QMin:
        string += ' and Veloc! '
    else:
        string += '! '
    string += 'NAtom is %i.\n' % (QMin['natom'])
    string += 'Charge of the molecule is %i.\n' % (QMin['template']['charge'])
    print(string)

    string = '\nGeometry in Bohrs:\n'
    for i in range(QMin['natom']):
        string += '%s ' % (QMin['geo'][i][0])
        for j in range(3):
            string += '% 7.4f ' % (QMin['geo'][i][j + 1])
        string += '\n'
    print(string)

    if 'veloc' in QMin:
        string = ''
        for i in range(QMin['natom']):
            string += '%s ' % (QMin['geo'][i][0])
            for j in range(3):
                string += '% 7.4f ' % (QMin['veloc'][i][j])
            string += '\n'
        print(string)

    if 'grad' in QMin:
        string = 'Gradients:   '
        for i in range(1, QMin['nmstates'] + 1):
            if i in QMin['grad']:
                string += 'X '
            else:
                string += '. '
        string += '\n'
        print(string)

    if 'nacdr' in QMin:
        string = 'Non-adiabatic couplings:\n'
        for i in range(1, QMin['nmstates'] + 1):
            for j in range(1, QMin['nmstates'] + 1):
                if [i, j] in QMin['nacdr'] or [j, i] in QMin['nacdr']:
                    string += 'X '
                else:
                    string += '. '
            string += '\n'
        print(string)

    if 'overlap' in QMin:
        string = 'Overlaps:\n'
        for i in range(1, QMin['nmstates'] + 1):
            for j in range(1, QMin['nmstates'] + 1):
                if [i, j] in QMin['overlap'] or [j, i] in QMin['overlap']:
                    string += 'X '
                else:
                    string += '. '
            string += '\n'
        print(string)

    for i in QMin:
        if not any([i == j for j in ['h', 'dm', 'dmdr', 'geo', 'veloc', 'states', 'comment', 'LD_LIBRARY_PATH', 'grad', 'nacdr', 'ion', 'overlap', 'template']]):
            if not any([i == j for j in ['ionlist', 'ionmap']]) or DEBUG:
                string = i + ': '
                string += str(QMin[i])
                print(string)
    print('\n')
    sys.stdout.flush()

# ======================================================================= #


def printtasks(tasks):
    '''If PRINT, prints a formatted table of the tasks in the tasks list.

    Arguments:
    1 list of lists: tasks list (see gettasks for specs)'''

    # if DEBUG:
    # pprint.pprint(tasks)
    if not PRINT:
        return
    print('==> Task Queue:\n')
    for i in range(len(tasks)):
        task = tasks[i]
        if task[0] == 'smo':
            print('Semiempirical Molecular Orbtial (SMO) SCF')
        elif task[0] == 'GUGA-CI':
            print('GUGA-CI: CASCI or MR-CI, total number of states: %i' % (task[1]))
        elif task[0] == 'smo_grad':
            print('SMO gradient')
        elif task[0] == 'GUGA-CI_grad':
            print('GUGA-CI gradient of state %i' % (task[1]))
        else:
            print(task)
    print('\n')
    sys.stdout.flush()

# ======================================================================= #


def printcomplexmatrix(matrix, states):
    '''Prints a formatted matrix. Zero elements are not printed, blocks of different mult and MS are delimited by dashes. Also prints a matrix with the imaginary parts, of any one element has non-zero imaginary part.

    Arguments:
    1 list of list of complex: the matrix
    2 list of integers: states specs'''

    nmstates = 0
    for i in range(len(states)):
        nmstates += states[i] * (i + 1)
    string = 'Real Part:\n'
    string += '-' * (11 * nmstates + nmstates // 3)
    string += '\n'
    istate = 0
    for imult, i, ms in itnmstates(states):
        jstate = 0
        string += '|'
        for jmult, j, ms2 in itnmstates(states):
            if matrix[istate][jstate].real == 0.:
                string += ' ' * 11
            else:
                string += '% .3e ' % (matrix[istate][jstate].real)
            if j == states[jmult - 1]:
                string += '|'
            jstate += 1
        string += '\n'
        if i == states[imult - 1]:
            string += '-' * (11 * nmstates + nmstates // 3)
            string += '\n'
        istate += 1
    print(string)
    imag = False
    string = 'Imaginary Part:\n'
    string += '-' * (11 * nmstates + nmstates // 3)
    string += '\n'
    istate = 0
    for imult, i, ms in itnmstates(states):
        jstate = 0
        string += '|'
        for jmult, j, ms2 in itnmstates(states):
            if matrix[istate][jstate].imag == 0.:
                string += ' ' * 11
            else:
                imag = True
                string += '% .3e ' % (matrix[istate][jstate].imag)
            if j == states[jmult - 1]:
                string += '|'
            jstate += 1
        string += '\n'
        if i == states[imult - 1]:
            string += '-' * (11 * nmstates + nmstates // 3)
            string += '\n'
        istate += 1
    string += '\n'
    if imag:
        print(string)

# ======================================================================= #


def printgrad(grad, natom, geo):
    '''Prints a gradient or nac vector. Also prints the atom elements. If the gradient is identical zero, just prints one line.

    Arguments:
    1 list of list of float: gradient
    2 integer: natom
    3 list of list: geometry specs'''

    string = ''
    iszero = True
    for atom in range(natom):
        string += '%i\t%s\t' % (atom + 1, geo[atom][0])
        for xyz in range(3):
            if grad[atom][xyz] != 0:
                iszero = False
            g = grad[atom][xyz]
            if isinstance(g, float):
                string += '% .5f\t' % (g)
            elif isinstance(g, complex):
                string += '% .5f\t% .5f\t\t' % (g.real, g.imag)
        string += '\n'
    if iszero:
        print('\t\t...is identical zero...\n')
    else:
        print(string)

# ======================================================================= #


def printQMout(QMin, QMout):
    '''If PRINT, prints a summary of all requested QM output values. Matrices are formatted using printcomplexmatrix, vectors using printgrad.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout'''

    # if DEBUG:
    # pprint.pprint(QMout)
    if not PRINT:
        return
    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    print('===> Results:\n')
    # Hamiltonian matrix, real or complex
    if 'h' in QMin or 'soc' in QMin:
        eshift = math.ceil(QMout['h'][0][0].real)
        print('=> Hamiltonian Matrix:\nDiagonal Shift: %9.2f' % (eshift))
        matrix = deepcopy(QMout['h'])
        for i in range(nmstates):
            matrix[i][i] -= eshift
        printcomplexmatrix(matrix, states)
    # Dipole moment matrices
    if 'dm' in QMin:
        print('=> Dipole Moment Matrices:\n')
        for xyz in range(3):
            print('Polarisation %s:' % (IToPol[xyz]))
            matrix = QMout['dm'][xyz]
            printcomplexmatrix(matrix, states)
    # Gradients
    if 'grad' in QMin:
        print('=> Gradient Vectors:\n')
        istate = 0
        for imult, i, ms in itnmstates(states):
            print('%s\t%i\tMs= % .1f:' % (IToMult[imult], i, ms))
            printgrad(QMout['grad'][istate], natom, QMin['geo'])
            istate += 1
    # Nonadiabatic couplings
    if 'nacdr' in QMin:
        print('=> Analytical Non-adiabatic coupling vectors:\n')
        istate = 0
        for imult, i, msi in itnmstates(states):
            jstate = 0
            for jmult, j, msj in itnmstates(states):
                if imult == jmult and msi == msj:
                    print('%s\tStates %i - %i\tMs= % .1f:' % (IToMult[imult], i, j, msi))
                    printgrad(QMout['nacdr'][istate][jstate], natom, QMin['geo'])
                jstate += 1
            istate += 1
    # Overlaps
    if 'overlap' in QMin:
        print('=> Overlap matrix:\n')
        matrix = QMout['overlap']
        printcomplexmatrix(matrix, states)
        if 'phases' in QMout:
            print('=> Wavefunction Phases:\n')
            for i in range(nmstates):
                print('% 3.1f % 3.1f' % (QMout['phases'][i].real, QMout['phases'][i].imag))
            print('\n')
    # Spin-orbit coupling derivatives
    if 'socdr' in QMin:
        print('=> Spin-Orbit Gradient Vectors:\n')
        istate = 0
        for imult, i, ims in itnmstates(states):
            jstate = 0
            for jmult, j, jms in itnmstates(states):
                print('%s\t%i\tMs= % .1f -- %s\t%i\tMs= % .1f:' % (IToMult[imult], i, ims, IToMult[jmult], j, jms))
                printgrad(QMout['socdr'][istate][jstate], natom, QMin['geo'])
                jstate += 1
            istate += 1
    # Dipole moment derivatives
    if 'dmdr' in QMin:
        print('=> Dipole moment derivative vectors:\n')
        istate = 0
        for imult, i, msi in itnmstates(states):
            jstate = 0
            for jmult, j, msj in itnmstates(states):
                if imult == jmult and msi == msj:
                    for ipol in range(3):
                        print('%s\tStates %i - %i\tMs= % .1f\tPolarization %s:' % (IToMult[imult], i, j, msi, IToPol[ipol]))
                        printgrad(QMout['dmdr'][ipol][istate][jstate], natom, QMin['geo'])
                jstate += 1
            istate += 1
    # Property matrix (dyson norms)
    if 'ion' in QMin and 'prop' in QMout:
        print('=> Property matrix:\n')
        matrix = QMout['prop']
        printcomplexmatrix(matrix, states)
    sys.stdout.flush()


# =============================================================================================== #
# =============================================================================================== #
# ======================================= Matrix initialization ================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #         OK
def makecmatrix(a, b):
    '''Initialises a complex axb matrix.

    Arguments:
    1 integer: first dimension
    2 integer: second dimension

    Returns;
    1 list of list of complex'''

    mat = [[complex(0., 0.) for i in range(a)] for j in range(b)]
    return mat

# ======================================================================= #         OK


def makermatrix(a, b):
    '''Initialises a real axb matrix.

    Arguments:
    1 integer: first dimension
    2 integer: second dimension

    Returns;
    1 list of list of real'''

    mat = [[0. for i in range(a)] for j in range(b)]
    return mat


# =============================================================================================== #
# =============================================================================================== #
# =========================================== output extraction ================================= #
# =============================================================================================== #
# =============================================================================================== #
# ======================================================================= #


def getcienergy(out, mult, state, corr_method):
    '''Searches a complete MNDO output file for the MRCI energy of (mult,state).

    Arguments:
    1 list of strings: MNDO output
    2 integer: mult
    3 integer: state
    4 corr_method: correlation method

    Returns:
    1 float: total CI energy of specified state in hartree'''

    if corr_method=='none':
        energystring = '     SCF TOTAL ENERGY'
    elif corr_method=='GUGA-CI':
        energystring = ' State'
        stateindex = 1 
        enindex = 8 

    if corr_method=='none':
        for i, line in enumerate(out):
            if energystring in line:
                l = line.split()
                return float(float(l[3])/au2eV) 
        else:
            print('Energy not found!', mult, state)
            sys.exit(20)
    elif corr_method=='GUGA-CI':
        for i, line in enumerate(out):
            if energystring in line:
                l = line.split()
                l[stateindex]=l[stateindex].replace(',', '')
                if int(l[stateindex]) == state:
                    return float(float(l[enindex])/au2eV)
        else:
            print('Energy not found!', mult, state)
            sys.exit(20)
            
# ======================================================================= #


def getcidm(out, mult1, state1, state2, pol, corr_method):
    # two cases:
    # - Dipole moments are in RASSI calculation with only one JOBIPH file
    # - Dipole moments are in RASSI calculation with two JOBIPH files of same multiplicity

    if pol == 'X' or pol == 'Y' or pol == 'Z':
        pol = IToPol[pol]

    if state1>state2:
        lower_state=state2
        higher_state=state1
        diagonal=False
    elif state1<state2:
        lower_state=state1
        higher_state=state2
        diagonal=False
    elif state1==state2:
        lower_state=state1
        higher_state=state2
        diagonal=True

    if corr_method=='GUGA-CI' and diagonal:
        startstring=' State dipole moments'
        endstring=' Properties of transitions'
    elif corr_method=='GUGA-CI' and not diagonal:
        startstring=' State dipole moments'
        endstring='     TIME FOR ENERGY EVALUATION'
        if lower_state<10:
            string1=' Dipole-length electric dipole transition moments   '
            string2=' Dipole-velocity electric dipole transition moments   '
        else:
            string1=' Dipole-length electric dipole transition moments  '
            string2=' Dipole-velocity electric dipole transition moments  '
        startstring1=f"{string1}{lower_state}"
        endstring1=f"{string2}{lower_state}"



    start_line=0
    end_line=0
    start_line1=0
    end_line1=0
    if corr_method=='GUGA-CI':
        for iline, line in enumerate(out):
            if startstring in line:
                 start_line=iline+1
            if endstring in line:
                 end_line=iline+1
        if start_line==0 or end_line==0:
            print('DM element not found!', mult1, state1, state2, pol)
            sys.exit(20)
        if not diagonal:
            for jline, line in enumerate(out):
                if startstring1 in line:
                    start_line1=jline+1
                if endstring1 in line:
                    end_line1=jline



    # Now start searching at iline, looking for the requested matrix element
    state=-1
    if corr_method=='GUGA-CI' and diagonal:
        for iline, line in enumerate(out[start_line:end_line]):
            try:
                state=int(line.split()[0])
            except ValueError:
                pass
            except IndexError:
                pass
            if state==lower_state:
                if pol==0:
                    return float(float(line.split()[5])/au2debye)
                elif pol==1:
                    return float(float(line.split()[6])/au2debye)     
                elif pol==2:
                    return float(float(line.split()[7])/au2debye)
    elif corr_method=='GUGA-CI' and not diagonal:
        for iline, line in enumerate(out[start_line1:end_line1]):
            try:
                state=int(line.split()[0])
            except ValueError:
                pass
            except IndexError:
                pass
            if state==higher_state:
                if pol==0:
                    return float(float(line.split()[5])/au2debye)
                elif pol==1:
                    return float(float(line.split()[6])/au2debye)
                elif pol==2:
                    return float(float(line.split()[7])/au2debye)
    elif corr_method=='none':
        string='     SUM'
        for iline, line in enumerate(out):
            if string in line:
                if pol==0:
                    return float(float(line.split()[1])/au2debye)
                elif pol==1:
                    return float(float(line.split()[2])/au2debye)
                elif pol==2:
                    return float(float(line.split()[3])/au2debye)

# ======================================================================= #
def getgrad(out, mult, state, QMin):

 
    grad = []
    grad_state=0
    if QMin['template']['corr_method']=='GUGA-CI':
        string = 'LROOT  ='
        for iline, line in enumerate(out):
            if string in line:
                index=line.find("LROOT  =")
                if index != -1:
                    number_part = line[index + len("LROOT  ="):]
                    number_part = number_part.split()[0]
                    grad_state = int(number_part)
        if grad_state==0:
             print('Gradient not found!', mult, state) 
             sys.exit(20)

    elif QMin['template']['corr_method']=='none':
        grad_state = 1

    if grad_state==state:
        string = 'GRADIENTS (KCAL' 
        for iline, line in enumerate(out):
            if string in line:
                start_line=iline+4
                end_line=iline+5+QMin['natom']
        for atom in range(QMin['natom']):
            line=out[start_line+atom]
            grad.append([float(float(line.split()[5 + xyz])/au2kcal*au2a) for xyz in range(3)])
        return grad

# ======================================================================= #


def getQMout(out, QMin):
    '''Constructs the requested matrices and vectors using the get<quantity> routines.

    The dictionary QMout contains all the requested properties. Its content is dependent on the keywords in QMin:
    - 'h' in QMin:
                    QMout['h']: list(nmstates) of list(nmstates) of complex, the non-relaticistic hamiltonian
    - 'soc' in QMin:
                    QMout['h']: list(nmstates) of list(nmstates) of complex, the spin-orbit hamiltonian
    - 'dm' in QMin:
                    QMout['dm']: list(3) of list(nmstates) of list(nmstates) of complex, the three dipole moment matrices
    - 'grad' in QMin:
                    QMout['grad']: list(nmstates) of list(natom) of list(3) of float, the gradient vectors of every state (even if "grad all" was not requested, all nmstates gradients are contained here)
    - 'nac' in QMin and QMin['nac']==['num']:
                    QMout['nac']: list(nmstates) of list(nmstates) of complex, the non-adiabatic coupling matrix
                    QMout['mrcioverlap']: list(nmstates) of list(nmstates) of complex, the MRCI overlap matrix
                    QMout['h']: like with QMin['h']
    - 'nac' in QMin and QMin['nac']==['ana']:
                    QMout['nac']: list(nmstates) of list(nmstates) of list(natom) of list(3) of float, the matrix of coupling vectors
    - 'nac' in QMin and QMin['nac']==['smat']:
                    QMout['nac']: list(nmstates) of list(nmstates) of complex, the adiabatic-diabatic transformation matrix
                    QMout['mrcioverlap']: list(nmstates) of list(nmstates) of complex, the MRCI overlap matrix
                    QMout['h']: like with QMin['h']

    Arguments:
    1 list of strings: Concatenated MNDO output
    2 dictionary: QMin

    Returns:
    1 dictionary: QMout'''


    # get version of MNDO
    mo_method = QMin['template']['mo_method']
    corr_method = QMin['template']['corr_method']

    # Currently implemented keywords: h, grad, dm, overlap, phase.
    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    QMout = {}

    # h: get CI energies of all ci calculations and construct hamiltonian, returns a matrix(nmstates,nmstates)
    if 'h' in QMin:
        # no spin-orbit couplings, hamilton operator diagonal, only one loop
        h = makecmatrix(nmstates, nmstates)
        for istate, i in enumerate(QMin['statemap'].values()):
            mult, state, ms = tuple(i)
            h[istate][istate] = complex(getcienergy(out, mult, state, corr_method))
        QMout['h'] = h

    # DM: get vector of three dipole matrices, three nested loops, returns a list of three matrices(nmstates,nmstates)
    # Notice that we employed the dipole-length transition moment format. 
    # For relations of length and velocity format, see for example, Chem. Phys. Lett. 95, 568-572 (1983) by Roginsky, Klapisch, and Cohen
    if 'dm' in QMin:
        dm = []
        for xyz in range(3):
            dm.append(makecmatrix(nmstates, nmstates))
            for istate, i in enumerate(QMin['statemap']):
                for jstate, j in enumerate(QMin['statemap']):
                    mult1, state1, ms1 = tuple(QMin['statemap'][i])
                    mult2, state2, ms2 = tuple(QMin['statemap'][j])
                    if mult1==mult2 and ms1==ms2:
                        dm[xyz][istate][jstate] = complex(getcidm(out, mult1, state1, state2, xyz, corr_method))
                    else:
                        dm[xyz][istate][jstate] = complex(0.0)
        QMout['dm'] = dm

    # Grad:  returns a list of nmstates vectors
    if 'grad' in QMin:
        grad = []
        for i in QMin['statemap']:
            mult, state, ms = tuple(QMin['statemap'][i])
            if (mult, state) in QMin['gradmap']:
                grad.append(getgrad(out, mult, state, QMin))
            else:
                gradatom = []
                for iatom in range(natom):
                    gradatom.append([0., 0., 0.])
                grad.append(gradatom)
        QMout['grad'] = grad

    return QMout

# =============================================================================================== #
# =============================================================================================== #
# =========================================== QMout writing ===================================== #
# =============================================================================================== #
# =============================================================================================== #


# ======================================================================= #
def writeQMout(QMin, QMout, QMinfilename):
    '''Writes the requested quantities to the file which SHARC reads in. The filename is QMinfilename with everything after the first dot replaced by "out".

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout
    3 string: QMinfilename'''

    k = QMinfilename.find('.')
    if k == -1:
        outfilename = QMinfilename + '.out'
    else:
        outfilename = QMinfilename[:k] + '.out'
    if PRINT:
        print('===> Writing output to file %s in SHARC Format\n' % (outfilename))
    string = ''
    if 'h' in QMin or 'soc' in QMin:
        string += writeQMoutsoc(QMin, QMout)
    if 'dm' in QMin:
        string += writeQMoutdm(QMin, QMout)
    if 'grad' in QMin:
        string += writeQMoutgrad(QMin, QMout)
    if 'nacdr' in QMin:
        string += writeQMoutnacana(QMin, QMout)
    if 'overlap' in QMin:
        string += writeQMoutnacsmat(QMin, QMout)
    if 'socdr' in QMin:
        string += writeQMoutsocdr(QMin, QMout)
    if 'dmdr' in QMin:
        string += writeQMoutdmdr(QMin, QMout)
    if 'ion' in QMin:
        string += writeQMoutprop(QMin, QMout)
    if 'phases' in QMin:
        string += writeQmoutPhases(QMin, QMout)
    string += writeQMouttime(QMin, QMout)
    outfile = os.path.join(QMin['pwd'], outfilename)
    writefile(outfile, string)
    return

# ======================================================================= #


def writeQMoutsoc(QMin, QMout):
    '''Generates a string with the Spin-Orbit Hamiltonian in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the SOC matrix'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Hamiltonian Matrix (%ix%i, complex)\n' % (1, nmstates, nmstates)
    string += '%i %i\n' % (nmstates, nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string += '%s %s ' % (eformat(QMout['h'][i][j].real, 9, 3), eformat(QMout['h'][i][j].imag, 9, 3))
        string += '\n'
    string += '\n'
    return string

# ======================================================================= #


def writeQMoutdm(QMin, QMout):
    '''Generates a string with the Dipole moment matrices in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line. The string contains three such matrices.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the DM matrices'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Dipole Moment Matrices (3x%ix%i, complex)\n' % (2, nmstates, nmstates)
    for xyz in range(3):
        string += '%i %i\n' % (nmstates, nmstates)
        for i in range(nmstates):
            for j in range(nmstates):
                string += '%s %s ' % (eformat(QMout['dm'][xyz][i][j].real, 9, 3), eformat(QMout['dm'][xyz][i][j].imag, 9, 3))
            string += '\n'
        # string+='\n'
    return string

# ======================================================================= #


def writeQMoutgrad(QMin, QMout):
    '''Generates a string with the Gradient vectors in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates gradients are written).

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the Gradient vectors'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Gradient Vectors (%ix%ix3, real)\n' % (3, nmstates, natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        string += '%i %i ! %i %i %i\n' % (natom, 3, imult, istate, ims)
        for atom in range(natom):
            for xyz in range(3):
                string += '%s ' % (eformat(QMout['grad'][i][atom][xyz], 9, 3))
            string += '\n'
        # string+='\n'
        i += 1
    return string

# ======================================================================= #


def writeQMoutnacana(QMin, QMout):
    '''Generates a string with the NAC vectors in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. On the next line, natom and 3 are written, followed by the gradient, with one line per atom and a blank line at the end. Each MS component shows up (nmstates x nmstates vectors are written).

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the NAC vectors'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Non-adiabatic couplings (ddr) (%ix%ix%ix3, real)\n' % (5, nmstates, nmstates, natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        j = 0
        for jmult, jstate, jms in itnmstates(states):
            string += '%i %i ! %i %i %i %i %i %i\n' % (natom, 3, imult, istate, ims, jmult, jstate, jms)
            for atom in range(natom):
                for xyz in range(3):
                    string += '%s ' % (eformat(QMout['nacdr'][i][j][atom][xyz], 12, 3))
                string += '\n'
            string += ''
            j += 1
        i += 1
    return string

# ======================================================================= #


def writeQMoutnacsmat(QMin, QMout):
    '''Generates a string with the adiabatic-diabatic transformation matrix in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the transformation matrix'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Overlap matrix (%ix%i, complex)\n' % (6, nmstates, nmstates)
    string += '%i %i\n' % (nmstates, nmstates)
    for j in range(nmstates):
        for i in range(nmstates):
            string += '%s %s ' % (eformat(QMout['overlap'][j][i].real, 9, 3), eformat(QMout['overlap'][j][i].imag, 9, 3))
        string += '\n'
    string += '\n'
    return string

# ======================================================================= #


def writeQMoutdmdr(QMin, QMout):

    states = QMin['states']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Dipole moment derivatives (%ix%ix3x%ix3, real)\n' % (12, nmstates, nmstates, natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        j = 0
        for jmult, jstate, jms in itnmstates(states):
            for ipol in range(3):
                string += '%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i   pol %i\n' % (natom, 3, imult, istate, ims, jmult, jstate, jms, ipol)
                for atom in range(natom):
                    for xyz in range(3):
                        string += '%s ' % (eformat(QMout['dmdr'][ipol][i][j][atom][xyz], 12, 3))
                    string += '\n'
                string += ''
            j += 1
        i += 1
    return string

# ======================================================================= #


def writeQMoutsocdr(QMin, QMout):

    states = QMin['states']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Spin-Orbit coupling derivatives (%ix%ix3x%ix3, complex)\n' % (13, nmstates, nmstates, natom)
    i = 0
    for imult, istate, ims in itnmstates(states):
        j = 0
        for jmult, jstate, jms in itnmstates(states):
            string += '%i %i ! m1 %i s1 %i ms1 %i   m2 %i s2 %i ms2 %i\n' % (natom, 3, imult, istate, ims, jmult, jstate, jms)
            for atom in range(natom):
                for xyz in range(3):
                    string += '%s %s ' % (eformat(QMout['socdr'][i][j][atom][xyz].real, 12, 3), eformat(QMout['socdr'][i][j][atom][xyz].imag, 12, 3))
            string += '\n'
            string += ''
            j += 1
        i += 1
    return string

# ======================================================================= #


def writeQMouttime(QMin, QMout):
    '''Generates a string with the quantum mechanics total runtime in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the runtime is given

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the runtime'''

    string = '! 8 Runtime\n%s\n' % (eformat(QMout['runtime'], 9, 3))
    return string

# ======================================================================= #


def writeQMoutprop(QMin, QMout):
    '''Generates a string with the Spin-Orbit Hamiltonian in SHARC format.

    The string starts with a ! followed by a flag specifying the type of data. In the next line, the dimensions of the matrix are given, followed by nmstates blocks of nmstates elements. Blocks are separated by a blank line.

    Arguments:
    1 dictionary: QMin
    2 dictionary: QMout

    Returns:
    1 string: multiline string with the SOC matrix'''

    states = QMin['states']
    nstates = QMin['nstates']
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    string = ''
    string += '! %i Property Matrix (%ix%i, complex)\n' % (11, nmstates, nmstates)
    string += '%i %i\n' % (nmstates, nmstates)
    for i in range(nmstates):
        for j in range(nmstates):
            string += '%s %s ' % (eformat(QMout['prop'][i][j].real, 12, 3), eformat(QMout['prop'][i][j].imag, 12, 3))
        string += '\n'
    string += '\n'
    return string

# ======================================================================= #


def writeQmoutPhases(QMin, QMout):

    string = '! 7 Phases\n%i ! for all nmstates\n' % (QMin['nmstates'])
    for i in range(QMin['nmstates']):
        string += '%s %s\n' % (eformat(QMout['phases'][i].real, 9, 3), eformat(QMout['phases'][i].imag, 9, 3))
    return string


# =============================================================================================== #
# =============================================================================================== #
# =========================================== SUBROUTINES TO readQMin =========================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def checkscratch(SCRATCHDIR):
    '''Checks whether SCRATCHDIR is a file or directory. If a file, it quits with exit code 1, if its a directory, it passes. If SCRATCHDIR does not exist, tries to create it.

    Arguments:
    1 string: path to SCRATCHDIR'''

    exist = os.path.exists(SCRATCHDIR)
    if exist:
        isfile = os.path.isfile(SCRATCHDIR)
        if isfile:
            print('$SCRATCHDIR=%s exists and is a file!' % (SCRATCHDIR))
            sys.exit(28)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print('Cannot create SCRATCHDIR=%s\n' % (SCRATCHDIR))
            sys.exit(29)

# ======================================================================= #


def removequotes(string):
    if string.startswith("'") and string.endswith("'"):
        return string[1:-1]
    elif string.startswith('"') and string.endswith('"'):
        return string[1:-1]
    else:
        return string

# ======================================================================= #


def getsh2caskey(sh2cas, key):
    i = -1
    while True:
        i += 1
        try:
            line = re.sub('#.*$', '', sh2cas[i])
        except IndexError:
            break
        line = line.split(None, 1)
        if line == []:
            continue
        if key.lower() in line[0].lower():
            return line
    return ['', '']

# ======================================================================= #


def get_sh2cas_environ(sh2cas, key, environ=True, crucial=True):
    line = getsh2caskey(sh2cas, key)
    if line[0]:
        LINE = line[1]
        LINE = removequotes(LINE).strip()
    else:
        if environ:
            LINE = os.getenv(key.upper())
            if not LINE:
                if crucial:
                    print('Either set $%s or give path to %s in MNDO.resources!' % (key.upper(), key.upper()))
                    sys.exit(30)
                else:
                    return ''
        else:
            if crucial:
                print('Give path to %s in MNDO.resources!' % (key.upper()))
                sys.exit(31)
            else:
                return ''
    LINE = os.path.expandvars(LINE)
    LINE = os.path.expanduser(LINE)
    if containsstring(';', LINE):
        print("$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(), key.upper()))
        sys.exit(32)
    return LINE

# ======================================================================= #


def get_pairs(QMinlines, i):
    nacpairs = []
    while True:
        i += 1
        try:
            line = QMinlines[i].lower()
        except IndexError:
            print('"keyword select" has to be completed with an "end" on another line!')
            sys.exit(33)
        if 'end' in line:
            break
        fields = line.split()
        try:
            nacpairs.append([int(fields[0]), int(fields[1])])
        except ValueError:
            print('"nacdr select" is followed by pairs of state indices, each pair on a new line!')
            sys.exit(34)
    return nacpairs, i

# ======================================================================= #         OK


def readQMin(QMinfilename):
    '''Reads the time-step dependent information from QMinfilename. This file contains all information from the current SHARC job: geometry, velocity, number of states, requested quantities along with additional information. The routine also checks this input and obtains a number of environment variables necessary to run MNDO.

    Steps are:
    - open and read QMinfilename
    - Obtain natom, comment, geometry (, velocity)
    - parse remaining keywords from QMinfile
    - check keywords for consistency, calculate nstates, nmstates
    - obtain environment variables for path to MNDO and scratch directory, and for error handling

    Arguments:
    1 string: name of the QMin file

    Returns:
    1 dictionary: QMin'''

    # read QMinfile
    QMinlines = readfile(QMinfilename)
    QMin = {}

    # Get natom
    try:
        natom = int(QMinlines[0])
    except ValueError:
        print('first line must contain the number of atoms!')
        sys.exit(35)
    QMin['natom'] = natom
    if len(QMinlines) < natom + 4:
        print('Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task')
        sys.exit(36)

    # Save Comment line
    QMin['comment'] = QMinlines[1]

    # Get geometry and possibly velocity (for backup-analytical non-adiabatic couplings)
    QMin['geo'] = []
    QMin['veloc'] = []
    hasveloc = True
    for i in range(2, natom + 2):
        if not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
            print('Input file does not comply to xyz file format! Maybe natom is just wrong.')
            sys.exit(37)
        fields = QMinlines[i].split()
        for j in range(1, 4):
            fields[j] = float(fields[j])
        QMin['geo'].append(fields[0:4])
        if len(fields) >= 7:
            for j in range(4, 7):
                fields[j] = float(fields[j])
            QMin['veloc'].append(fields[4:7])
        else:
            hasveloc = False
    if not hasveloc:
        QMin = removekey(QMin, 'veloc')


    # Parse remaining file
    i = natom + 1
    while i + 1 < len(QMinlines):
        i += 1
        line = QMinlines[i]
        line = re.sub('#.*$', '', line)
        if len(line.split()) == 0:
            continue
        key = line.lower().split()[0]
        if 'savedir' in key:
            args = line.split()[1:]
        else:
            args = line.lower().split()[1:]
        if key in QMin:
            print('Repeated keyword %s in line %i in input file! Check your input!' % (key, i + 1))
            continue  # only first instance of key in QM.in takes effect
        if len(args) >= 1 and 'select' in args[0]:
            pairs, i = get_pairs(QMinlines, i)
            QMin[key] = pairs
        else:
            QMin[key] = args

    # Input coordinate unit for MNDO is Angstrom
    if 'unit' in QMin:
        if QMin['unit'][0] == 'angstrom':
            factor = 1. / au2a
        elif QMin['unit'][0] == 'bohr':
            factor = 1.
        else:
            print('Dont know input unit %s!' % (QMin['unit'][0]))
            sys.exit(38)
    else:
        factor = 1.

    for iatom in range(len(QMin['geo'])):
        for ixyz in range(3):
            QMin['geo'][iatom][ixyz + 1] *= factor


    if 'states' not in QMin:
        print('Keyword "states" not given!')
        sys.exit(39)
    # Calculate states, nstates, nmstates
    for i in range(len(QMin['states'])):
        QMin['states'][i] = int(QMin['states'][i])
    reduc = 0
    for i in reversed(QMin['states']):
        if i == 0:
            reduc += 1
        else:
            break
    for i in range(reduc):
        del QMin['states'][-1]
    nstates = 0
    nmstates = 0
    for i in range(len(QMin['states'])):
        nstates += QMin['states'][i]
        nmstates += QMin['states'][i] * (i + 1)
    QMin['nstates'] = nstates
    QMin['nmstates'] = nmstates


    # Various logical checks
    if 'states' not in QMin:
        print('Number of states not given in QM input file %s!' % (QMinfilename))
        sys.exit(40)

    possibletasks = ['h', 'dm', 'grad', 'dmdr', 'overlap']
    if not any([i in QMin for i in possibletasks]):
        print('No tasks found! Tasks are "h", "dm", "grad", "dmdr", "overlap".')
        sys.exit(41)

    if 'samestep' in QMin and 'init' in QMin:
        print('"Init" and "Samestep" cannot be both present in QM.in!')
        sys.exit(42)

    if 'phases' in QMin:
        QMin['overlap'] = []

    if 'overlap' in QMin and 'init' in QMin:
        print('"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"')
        sys.exit(43)

    if 'init' not in QMin and 'samestep' not in QMin:
        QMin['newstep'] = []

    if not any([i in QMin for i in ['h', 'soc', 'dm', 'grad']]) and 'overlap' in QMin:
        QMin['h'] = []

    if len(QMin['states']) > 8:
        print('Higher multiplicities than octets are not supported!')
        sys.exit(44)

    if 'soc' in QMin:
        print('Within the SHARC-MNDO interface, "soc" is not supported.')
        sys.exit(45)

    if 'nacdt' in QMin:
        print('Within the SHARC-MNDO interface, "nacdt" is not supported.')
        sys.exit(45)

    if 'molden' in QMin:
        print('Within the SHARC-MNDO interface, "molden" is not supported.')
        sys.exit(45)

    if 'ion' in QMin:
        print('Within the SHARC-MNDO interface, "ion" is not supported.')
        sys.exit(45)

    if 'overlap' in QMin:
        print('Within the SHARC-MNDO interface, "overlap" is not supported.')
        sys.exit(45)

    if 'phase' in QMin:
        print('Within the SHARC-MNDO interface, "phase" is not supported.')
        sys.exit(45)


    # Check for correct gradient list
    if 'grad' in QMin:
        if len(QMin['grad']) == 0 or QMin['grad'][0] == 'all':
            QMin['grad'] = [i + 1 for i in range(nmstates)]
            # pass
        else:
            for i in range(len(QMin['grad'])):
                try:
                    QMin['grad'][i] = int(QMin['grad'][i])
                except ValueError:
                    print('Arguments to keyword "grad" must be "all" or a list of integers!')
                    sys.exit(47)
                if QMin['grad'][i] > nmstates:
                    print('State for requested gradient does not correspond to any state in QM input file state list!')
                    sys.exit(48)

    # ===== NACdR is not available, left it here for furture use.
    # Process the overlap requests
    # identically to the nac requests
    if 'overlap' in QMin:
        if len(QMin['overlap']) >= 1:
            nacpairs = QMin['overlap']
            for i in range(len(nacpairs)):
                if nacpairs[i][0] > nmstates or nacpairs[i][1] > nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(49)
        else:
            QMin['overlap'] = [[j + 1, i + 1] for i in range(nmstates) for j in range(i + 1)]

    # type conversion has already been done
    if 'nacdr' in QMin:
        if len(QMin['nacdr']) >= 1:
            nacpairs = QMin['nacdr']
            for i in range(len(nacpairs)):
                if nacpairs[i][0] > nmstates or nacpairs[i][1] > nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(50)
        else:
            QMin['nacdr'] = [[j + 1, i + 1] for i in range(nmstates) for j in range(i)]

    # Turn off nacdr and overlap currently.
    if 'nacdr' in QMin:
        prin('Overlap or NACdR is currently not supported in MNDO')
        sys.exit(51)

    # obtain the statemap
    statemap = {}
    i = 1
    for imult, istate, ims in itnmstates(QMin['states']):
        statemap[i] = [imult, istate, ims]
        i += 1
    QMin['statemap'] = statemap

    # get the set of states for which gradients actually need to be calculated
    gradmap = set()
    if 'grad' in QMin:
        for i in QMin['grad']:
            gradmap.add(tuple(statemap[i][0:2]))
    gradmap = sorted(gradmap)
    QMin['gradmap'] = gradmap

    # get the list of statepairs for NACdr calculation
    nacmap = set()
    if 'nacdr' in QMin:
        for i in QMin['nacdr']:
            s1 = statemap[i[0]][0:2]
            s2 = statemap[i[1]][0:2]
            if s1[0] != s2[0] or s1 == s2:
                continue
            if s1[1] > s2[1]:
                continue
            nacmap.add(tuple(s1 + s2))
    nacmap = list(nacmap)
    nacmap.sort()
    QMin['nacmap'] = nacmap







    # open MNDO.resources
    filename = 'MNDO.resources'
    if os.path.isfile(filename):
        sh2cas = readfile(filename)
    else:
        print('Warning: No MNDO.resources found!')
        print('Reading resources from SH2CAS.inp')
        sh2cas = readfile('SH2CAS.inp')

    QMin['pwd'] = os.getcwd()

    QMin['mndo'] = get_sh2cas_environ(sh2cas, 'MNDO')
    os.environ['MNDO'] = QMin['mndo']

    driver = get_sh2cas_environ(sh2cas, 'driver', crucial=False)
    if driver == '':
        driver = os.path.join(QMin['mndo'], 'mndo2020')
        if not os.path.isfile(driver):
            print('No driver (mndo2020) found in $MNDO/. Please add the path to the driver via the "driver" keyword.')
            sys.exit(52)
    QMin['driver'] = driver

    # ===== ION is not currently implemented, left it here for furture use.
    if 'ion' in QMin:
        QMin['wfoverlap'] = get_sh2cas_environ(sh2cas, 'wfoverlap', crucial=False)
        if not QMin['wfoverlap']:
            ciopath = os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')), 'wfoverlap.x')
            if os.path.isfile(ciopath):
                QMin['wfoverlap'] = ciopath
            else:
                print('Give path to wfoverlap.x in MNDO.resources!')
                sys.exit(51)
        # read ncore and ndocc from resources
        line = getsh2caskey(sh2cas, 'numfrozcore')
        if line[0]:
            try:
                QMin['ncore'] = max(0, int(line[1]))
            except ValueError:
                print('numfrozcore does not evaluate to numerical value!')
                sys.exit(52)
        line = getsh2caskey(sh2cas, 'numocc')
        if line[0]:
            try:
                QMin['ndocc'] = int(line[1])
            except ValueError:
                print('numocc does not evaluate to numerical value!')
                sys.exit(53)


    # Set up scratchdir
    line = get_sh2cas_environ(sh2cas, 'scratchdir', False, False)
    if line is None:
        line = QMin['pwd'] + '/SCRATCHDIR/'
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    # checkscratch(line)
    QMin['scratchdir'] = line


    # Set up savedir
    if 'savedir' in QMin:
        # savedir may be read from QM.in file
        line = QMin['savedir'][0]
    else:
        line = get_sh2cas_environ(sh2cas, 'savedir', False, False)
        if line is None or line == '':
            line = os.path.join(QMin['pwd'], 'SAVEDIR')
    line = os.path.expandvars(line)
    line = os.path.expanduser(line)
    line = os.path.abspath(line)
    if 'init' in QMin:
        checkscratch(line)
    QMin['savedir'] = line


    line = getsh2caskey(sh2cas, 'debug')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global DEBUG
            DEBUG = True

    line = getsh2caskey(sh2cas, 'no_print')
    if line[0]:
        if len(line) <= 1 or 'true' in line[1].lower():
            global PRINT
            PRINT = False

    QMin['memory'] = 500
    line = getsh2caskey(sh2cas, 'memory')
    if line[0]:
        try:
            QMin['memory'] = int(line[1])
        except ValueError:
            print('MNDO memory does not evaluate to numerical value!')
            sys.exit(54)
    else:
        print('WARNING: Please set memory for MNDO in MNDO.resources (in MB)! Using 1000 MB default value!')
    os.environ['MNDOMEM'] = str(QMin['memory'])
    os.environ['MNDO_MEM'] = str(QMin['memory'])

    QMin['ncpu'] = 1
    line = getsh2caskey(sh2cas, 'ncpu')
    if line[0]:
        try:
            QMin['ncpu'] = int(line[1])
        except ValueError:
            print('Number of CPUs does not evaluate to numerical value!')
            sys.exit(55)

    QMin['mpi_parallel'] = False
    line = getsh2caskey(sh2cas, 'mpi_parallel')
    if line[0]:
        QMin['mpi_parallel'] = True


    QMin['schedule_scaling'] = 0.6
    line = getsh2caskey(sh2cas, 'schedule_scaling')
    if line[0]:
        try:
            x = float(line[1])
            if 0 < x <= 2.:
                QMin['schedule_scaling'] = x
        except ValueError:
            print('Schedule scaling does not evaluate to numerical value!')
            sys.exit(56)

    QMin['Project'] = 'MNDO'
    os.environ['Project'] = QMin['Project']

    QMin['delay'] = 0.0
    line = getsh2caskey(sh2cas, 'delay')
    if line[0]:
        try:
            QMin['delay'] = float(line[1])
        except ValueError:
            print('Submit delay does not evaluate to numerical value!')
            sys.exit(57)

    QMin['Project'] = 'MNDO'
    os.environ['Project'] = QMin['Project']
    os.environ['MNDO_OUTPUT'] = 'PWD'


    # open template
    template = readfile('MNDO.template')

    QMin['template'] = {}
    #nacocc: number of occupied active orbital - ici1
    #nacvir: number of virtual active orbital - ici2
    #exlvl: excitation level in GUGA-CI - levexc
    #nref: number of reference configuration state functions - nciref
    #mo_method: semiempirical molecular orbital theory used - iop
    #corr_method: correlation method used - kci, currently only support kci=0 or kci=5
    #             none - kci=0; GUGA-CI - kci=5
    #refscf_method: GUGA-CI reference configuration state function method - mciref
    #               auto - mciref=0, SCF for nciref=1; SCF+D for nciref=2; 
    #                      SCF+SD for nciref=3
    #               read - mciref=1: read an additional line for reference orbitals
    #parameter: if uses parameter file in fort.14 and si.dat - iparok
    #           no - use default parameters 
    #           yes - set iparok=3, and read fort.14 and si.dat files (this will be improved later) 
    #escf: energy convergence criteria for SCF calculations - iscf
    #dscf: density convergence criteria for SCF calculations - iplscf
    #egcut: cutoff for three-center orthogonalization correction (both energy and gradient) - icuts
    #maxiter: maximum number of SCF iterations - kitscf
    integers = ['nacocc', 'nacvir', 'exlvl', 'nref', 'escf', 'dscf', 'egcut', 'gcut', 'maxiter', 'charge']
    strings = ['mo_method', 'corr_method', 'refscf_method', 'parameter']
    floats = ['displ']
    #set defaults:
    QMin['template']['roots'] = [0 for i in range(8)]
    QMin['template']['rootpad'] = [0 for i in range(8)]
    QMin['template']['mo_method'] = 'ODM3'
    QMin['template']['corr_method'] = 'GUGA-CI'
    QMin['template']['refscf_method'] = 'auto'
    QMin['template']['parameter'] = 'no'
    QMin['template']['nref'] = 1 #single SCF reference, ended up with SMO-CASCI
    QMin['template']['escf'] = 9 #10^-9 eV for GUGA-CI, default used in MNDO
    QMin['template']['dscf'] = 9 #10^-9 for GUGA-CI, defaul used in MNDO 
    QMin['template']['egcut'] = -1 #no cutoff, defauled used in MNDO
    QMin['template']['maxiter'] = 200 #200 iterations, default used in MNDO
    QMin['template']['charge'] = 0 #neutral
    QMin['template']['displ'] = 0.005 #displacement for numerical gradient

    for line in template:
        orig = re.sub('#.*$', '', line).split(None, 1)
        line = re.sub('#.*$', '', line).lower().split()
        if len(line) == 0:
            continue
        if 'spin' in line[0]:
            QMin['template']['roots'][int(line[1]) - 1] = int(line[3])
        elif 'roots' in line[0]:
            for i, n in enumerate(line[1:]):
                QMin['template']['roots'][i] = int(n)
        elif 'rootpad' in line[0]:
            for i, n in enumerate(line[1:]):
                QMin['template']['rootpad'][i] = int(n)
        elif line[0] in integers:
            QMin['template'][line[0]] = int(line[1])
        elif line[0] in strings:
            QMin['template'][line[0]] = line[1]

    # roots must be larger or equal to states
    for i, n in enumerate(QMin['template']['roots']):
        if i == len(QMin['states']):
            break
        if not n >= QMin['states'][i]:
            print('Too few states in state-averaging in multiplicity %i! %i requested, but only %i given' % (i + 1, QMin['states'][i], n))
            sys.exit(59)

    # check rootpad
    for i, n in enumerate(QMin['template']['rootpad']):
        if i == len(QMin['states']):
            break
        if not n >= 0:
            print('Rootpad must not be negative!')
            sys.exit(60)

    # condense roots list
    for i in range(len(QMin['template']['roots']) - 1, 0, -1):
        if QMin['template']['roots'][i] == 0:
            QMin['template']['roots'].pop(i)
        else:
            break
    QMin['template']['rootpad'] = QMin['template']['rootpad'][:len(QMin['template']['roots'])]


    if QMin['template']['corr_method'] == 'GUGA-CI':
        necessary = ['nacocc', 'nacvir', 'exlvl']
        for i in necessary:
            if i not in QMin['template']:
                print('Key %s missing in template file!' % (i))
                sys.exit(61)

    if QMin['template']['corr_method']=='none':
        if QMin['nstates']>1:
            print('Multiple electronic states calculation requires GUGA-CI, should set corr_method to GUGA-CI')
            sys.exit(61)

    # logic checks:
    if QMin['template']['nref']>3 and QMin['template']['refscf_method']=='auto':
        print('Automatic assignment of reference configuration state functions can not be used when nref > 3')
        sys.exit(62)

    # find semiempirical molecular orbital method
    allowed_mo_methods = ['MNDO', 'MNDO/d', 'MNDO/H', 'MNDO/dH', 'MNDOC', 'AM1', 'PM3', 'OM1', 'OM2', 'OM3', 'ODM2', 'ODM3', 'mndo', 'mndo/d', 'mndo/h', 'mndo/dh', 'mndoc', 'am1', 'pm3', 'om1', 'om2', 'om3', 'odm2', 'odm3']
    for i, m in enumerate(allowed_mo_methods):
        if QMin['template']['mo_method'] == m:
            QMin['mo_method'] = i
            break
    else:
        print('Unknown semiempirical molecular orbital method "%s" given in MNDO.template' % (QMin['template']['mo_method']))
        sys.exit(64)

    # now read reference CSFs
    if QMin['template']['corr_method'] == 'GUGA-CI':
        nref=QMin['template']['nref']
        nacocc=QMin['template']['nacocc']
        nacvir=QMin['template']['nacvir']
        if QMin['template']['nref']>3 and QMin['template']['refscf_method']=='read':
            QMin['template']['ref_csf'] = [[0] * (nacocc + nacvir) for _ in range(nref)] #reference CSFs. 
            ref_csf_lines = []
            for line in template:
                data = line.strip().split()
                if len(data)==(nacocc+nacvir+1) and data[0] == 'ref_csf':
                    ref_csf_lines.append(data[1:])
                elif len(data)!=(nacocc+nacvir+1):
                    print('provided reference CSFs do not seem correct, should have occupation number for each active orbital')
                    sys.exit(65)
            for i, data in enumerate(ref_csf_lines):
                for j, n in enumerate(data):
                    QMin['template']['ref_csf'][i][j]=int(n)


    # decide which type of gradients to do:
    # 0 = analytical MNDO gradient in one MNDO input (only works for corr_method=0, i.e. SMO only)
    # 1 = analytical MNDO gradients in separate MNDO inputs, possibly distributed over several CPUs (DEFAULT)
    # 2 = numerical MNDO gradients (dmdr, socdr), possibly distributed over several CPUs
    if 'grad' in QMin:
        if QMin['template']['corr_method']=='none' and QMin['nstates']==1:
            QMin['gradmode'] = 0
        elif QMin['template']['corr_method']=='GUGA-CI':
            QMin['gradmode'] = 1

    if 'dmdr' in QMin:
        print('Numerical gradients due to dmdr...')
        QMin['gradmode'] = 2
        QMin['displ'] = QMin['template']['displ']         

    QMin['ncpu'] = max(1, QMin['ncpu'])
 
    # Check the save directory
    try:
        ls = os.listdir(QMin['savedir'])
        err = 0
    except OSError:
        print('Problems reading SCRADIR=%s' % (QMin['savedir']))
        sys.exit(68)

    if PRINT:
        printQMin(QMin)

    return QMin


# =============================================================================================== #
# =============================================================================================== #
# =========================================== gettasks and setup routines ======================= #
# =============================================================================================== #
# =============================================================================================== #

def gettasks(QMin):

    # Currently implemented keywords: soc, dm, grad,
    tasks = []

    list_to_do = [(i, j) for i, j in enumerate(QMin['states'])]

    # for imult,nstates in enumerate(QMin['states']):
    for imult, nstates in list_to_do:
        if nstates == 0:
            continue

        # SMO 
        if QMin['template']['corr_method']=='none':
            #now you need to be careful here, because in MNDO, the following definition is used
            #imult=0: closed-shell singlet; imult=1: open-shell singlet; imult=2: doublet; imult=3: triplet 
            if imult==0:
                tasks.append(['smo', imult, QMin['template']['roots'][imult]])
            else:
                tasks.append(['smo', imult + 1, QMin['template']['roots'][imult]])

        # SMO followed by GUGA-CI
        if QMin['template']['corr_method']=='GUGA-CI':
            if imult==0:
                tasks.append(['smo', imult, QMin['template']['roots'][imult]])
            else:
                tasks.append(['smo', imult + 1, QMin['template']['roots'][imult]])
            tasks.append(['GUGA-CI', QMin['template']['roots'][imult]])

        # Gradients
        if QMin['gradmode'] == 0:
            for i in QMin['gradmap']:
                if i[0] == imult + 1:
                    if QMin['template']['corr_method']=='none':
                        tasks.append(['smo_grad'])
                    elif QMin['template']['corr_method']=='GUGA-CI':
                        tasks.append(['GIGA-CI_grad', i[1]])

        #space and comment 
        tasks.append(['space'])
        tasks.append(['comments'])

        #Geometry 
        tasks.append(['geom'])

        tasks.append(['space'])

        #reference CSFs
        if QMin['template']['nref']>4:
            tasks.append(['ref_csfs'])


    if DEBUG:
        printtasks(tasks)

    return tasks

# ======================================================================= #


def writeMNDOinput(tasks, QMin):

    string = ''

    for task in tasks:


        if task[0] == 'smo':
            #iop
            if QMin['template']['mo_method'] == 'MNDO' or QMin['template']['mo_method'] == 'mndo':
                string += 'iop=0 '
            elif QMin['template']['mo_method'] == 'MNDO/d' or QMin['template']['mo_method'] == 'mndo/d':
                string += 'iop=-10 '
            elif QMin['template']['mo_method'] == 'MNDO/H' or QMin['template']['mo_method'] == 'mndo/h':
                string += 'iop=-3 '
            elif QMin['template']['mo_method'] == 'MNDO/dH' or QMin['template']['mo_method'] == 'mndo/dh':
                string += 'iop=-13 '
            elif QMin['template']['mo_method'] == 'MNDOC' or QMin['template']['mo_method'] == 'mndoc':
                string += 'iop=-1 '
            elif QMin['template']['mo_method'] == 'AM1' or QMin['template']['mo_method'] == 'am1':
                string += 'iop=-2 '
            elif QMin['template']['mo_method'] == 'PM3' or QMin['template']['mo_method'] == 'pm3':
                string += 'iop=-7 '
            elif QMin['template']['mo_method'] == 'OM1' or QMin['template']['mo_method'] == 'om1':
                string += 'iop=-5 '
            elif QMin['template']['mo_method'] == 'OM2' or QMin['template']['mo_method'] == 'om2':
                string += 'iop=-6 '
            elif QMin['template']['mo_method'] == 'OM3' or QMin['template']['mo_method'] == 'om3':
                string += 'iop=-8 '
            elif QMin['template']['mo_method'] == 'ODM2' or QMin['template']['mo_method'] == 'odm2':
                string += 'iop=-22 '
            elif QMin['template']['mo_method'] == 'ODM3' or QMin['template']['mo_method'] == 'odm3':
                string += 'iop=-23 '
            elif QMin['template']['mo_method'] == 'MINDO/3' or QMin['template']['mo_method'] == 'mindo/3':
                string += 'iop=1 '
            elif QMin['template']['mo_method'] == 'CNDO/2' or QMin['template']['mo_method'] == 'cndo/2':
                string += 'iop=2 '
            elif QMin['template']['mo_method'] == 'SCC-DFTB' or QMin['template']['mo_method'] == 'scc-dftb':
                string += 'iop=5 '
            elif QMin['template']['mo_method'] == 'SCC-DFTB/J' or QMin['template']['mo_method'] == 'scc-dftb/j':
                string += 'iop=6 '
            #iparok 
            if QMin['template']['parameter']=='no':
                string += 'iparok=0 '
            elif QMin['template']['parameter']=='yes':
                string += 'iparok=3 '
            string += 'kitscf=%i ' % (QMin['template']['maxiter'])
            string += 'jop=-1 '
            string += 'igeom=1 +\n'
            string += 'kharge=%i ' % (QMin['template']['charge'])
            string += 'imult=%i ' % (task[1])
            string += 'iscf=%i ' % (QMin['template']['escf'])
            string += 'iplscf=%i ' % (QMin['template']['dscf'])
            string += 'icuts=%i \n' % (QMin['template']['egcut'])

        elif task[0] == 'GUGA-CI':
            string = string.replace('icuts=%i \n' % (QMin['template']['egcut']), 'icuts=%i +\n' % (QMin['template']['egcut']))
            string += 'kci=5 '
            string += 'ioutci=2 '
            string += 'ici1=%i ' % (QMin['template']['nacocc'])
            string += 'ici2=%i ' % (QMin['template']['nacvir'])
            string += 'levexc=%i ' % (QMin['template']['exlvl'])
            string += 'nciref=%i ' % (QMin['template']['nref'])
            string += 'iroot=%i ' % (task[1])
            string += 'iuvcd=3 \n'

 
        elif task[0] == 'smo_grad':
            string = string.replace("jop=-1 ", "jop=-2 ")


        elif task[0] == 'GIGA-CI_grad':
            string = string.replace('iuvcd=3 \n', 'iuvcd=3 +\n')
            string = string.replace("jop=-1 ", "jop=-2 ")
            string += 'lroot=%i \n' % (task[1])


        elif task[0] == 'space':
            string += '    \n'


        elif task[0] == 'comments':
            string += 'input from SHARC-MNDO interface\n'


        elif task[0] == 'geom':
            for iatom, atom in enumerate(QMin['geo']):
                for number, symbol in periodic_table.items():
                    if atom[0]==symbol:
                        atomic_number = number
                string += f'{atomic_number:2d}{"        "}{atom[1] * au2a:10.5f}{"          "}{atom[2] * au2a:10.5f}{"          "}{atom[3] * au2a:10.5f}'
                string += '\n' 


        elif task[0] == 'ref_csfs':
            for i in range(len(QMin['template']['ref_csf'])):
                for j in range(len(QMin['template']['ref_csf'][i])):
                    string += f"  {a[i][j]}" 
                string += '\n'


        else:
            print('Unknown task keyword %s found in writeMNDOinput!' % task[0])
            print(task)
            sys.exit(70)

    return string

# ======================================================================= #


def setupWORKDIR(WORKDIR, tasks, QMin):
    # mkdir the WORKDIR, or clean it if it exists, then copy all necessary JobIph files from pwd and savedir
    # then put the geom.xyz and MNDO.input files

    # setup the directory
    if os.path.exists(WORKDIR):
        if os.path.isfile(WORKDIR):
            print('%s exists and is a file!' % (WORKDIR))
            sys.exit(72)
        elif os.path.isdir(WORKDIR):
            if DEBUG:
                print('Remake\t%s' % WORKDIR)
            shutil.rmtree(WORKDIR)
            os.makedirs(WORKDIR)
    else:
        try:
            if DEBUG:
                print('Make\t%s' % WORKDIR)
            os.makedirs(WORKDIR)
        except OSError:
            print('Can not create %s\n' % (WORKDIR))
            sys.exit(73)

    # write MNDO.input
    inputstring = writeMNDOinput(tasks, QMin)
    filename = os.path.join(WORKDIR, 'MNDO.input')
    writefile(filename, inputstring)
    if DEBUG:
        print(inputstring)
        print('MNDO input written to: %s' % (filename))

    # make subdirs
    if QMin['mpi_parallel']:
        for i in range(QMin['ncpu'] - 1):
            subdir = os.path.join(WORKDIR, 'tmp_%i' % (i + 1))
            os.makedirs(subdir)

    return


# ======================================================================= #
def runMNDO(WORKDIR, MNDO, driver, ncpu, strip=False):
    prevdir = os.getcwd()
    os.chdir(WORKDIR)
    os.environ['WorkDir'] = WORKDIR
    os.environ['MNDO_NPROCS'] = str(ncpu)
    path = driver
    if not os.path.isfile(path):
        print('ERROR: could not find Molcas driver ("mndo2020") in $MNDO/!')
        sys.exit(74)

    string = path + ' < MNDO.input > MNDO.out'
    stdoutfile = open(os.path.join(WORKDIR, 'MNDO.out'), 'w')
    stderrfile = open(os.path.join(WORKDIR, 'MNDO.err'), 'w')
    if PRINT or DEBUG:
        starttime = datetime.datetime.now()
        sys.stdout.write('START:\t%s\t%s\t"%s"\n' % (WORKDIR, starttime, string))
        sys.stdout.flush()
    try:
        runerror = sp.call(string, shell=True, stdout=stdoutfile, stderr=stderrfile)
        # pass
    except OSError:
        print('Call have had some serious problems:', OSError)
        sys.exit(75)
    stdoutfile.close()
    stderrfile.close()
    if PRINT or DEBUG:
        endtime = datetime.datetime.now()
        sys.stdout.write('FINISH:\t%s\t%s\tRuntime: %s\tError Code: %i\n' % (WORKDIR, endtime, endtime - starttime, runerror))
        sys.stdout.flush()
    os.chdir(prevdir)
    if strip and not DEBUG:
        stripWORKDIR(WORKDIR)
    return runerror

# ======================================================================= #


def doDisplacement(QMin, idir, displ):
    iatom, ixyz, isign = tuple(idir)
    QMin1 = deepcopy(QMin)
    QMin1['geo'][iatom][ixyz + 1] += isign * displ
    return QMin1


# ======================================================================= #

def parallel_speedup(N, scaling):
    # computes the parallel speedup from Amdahls law
    # with scaling being the fraction of parallelizable work and (1-scaling) being the serial part
    return 1. / ((1 - scaling) + scaling / N)

# ======================================================================= #


def divide_slots(ncpu, ntasks, scaling):
    # this routine figures out the optimal distribution of the tasks over the CPU cores
    #   returns the number of rounds (how many jobs each CPU core will contribute to),
    #   the number of slots which should be set in the Pool,
    #   and the number of cores for each job.
    minpar = 1
    ntasks_per_round = ncpu // minpar
    if ncpu == 1:
        ntasks_per_round = 1
    ntasks_per_round = min(ntasks_per_round, ntasks)
    optimal = {}
    for i in range(1, 1 + ntasks_per_round):
        nrounds = int(math.ceil(float(ntasks) // i))
        ncores = ncpu // i
        optimal[i] = nrounds / parallel_speedup(ncores, scaling)
    # print optimal
    best = min(optimal, key=optimal.get)
    nrounds = int(math.ceil(float(ntasks) // best))
    ncores = ncpu // best

    cpu_per_run = [0 for i in range(ntasks)]
    if nrounds == 1:
        itask = 0
        for icpu in range(ncpu):
            cpu_per_run[itask] += 1
            itask += 1
            if itask >= ntasks:
                itask = 0
        nslots = ntasks
    else:
        for itask in range(ntasks):
            cpu_per_run[itask] = ncores
        nslots = ncpu // ncores
    # print(nrounds,nslots,cpu_per_run)
    return nrounds, nslots, cpu_per_run


# ======================================================================= #
def generate_joblist(QMin):
    '''split the full job into subtasks, each with a QMin dict, a WORKDIR
    structure: joblist = [ {WORKDIR: QMin, ..}, {..}, .. ]
    each element of the joblist is a set of jobs,
    and all jobs from the first set need to be completed before the second set can be processed.'''

    joblist = []

    if QMin['gradmode'] == 0:
        # case of serial gradients on one cpu
        QMin1 = deepcopy(QMin)
        QMin1['master'] = []
        if QMin['mpi_parallel']:
            QMin1['ncpu'] = QMin['ncpu']
        else:
            QMin1['ncpu'] = 1
        QMin['nslots_pool'] = [1]
        joblist.append({'master': QMin1})

    elif QMin['gradmode'] == 1:
        # case of analytical gradients for several states on several cpus

        # we will do wavefunction and dm, overlap always first
        # afterwards we will do all gradients and nacdr asynchonously
        QMin1 = deepcopy(QMin)
        QMin1['master'] = []
        QMin1['gradmap'] = []
        #QMin1['nacmap'] = []
        if QMin['mpi_parallel']:
            QMin1['ncpu'] = QMin['ncpu']
        else:
            QMin1['ncpu'] = 1
        QMin['nslots_pool'] = [1]
        joblist.append({'master': QMin1})

        QMin2 = deepcopy(QMin)
        remove = ['h', 'soc', 'dm', 'comment', 'ncpu', 'init', 'veloc', 'overlap', 'ion']
        for r in remove:
            QMin2 = removekey(QMin2, r)
        QMin2['gradmode'] = 0
        QMin2['samestep'] = []
        ntasks = len(QMin['gradmap'])
        if QMin['mpi_parallel']:
            # nrounds,nslots,cpu_per_run=divide_slots(QMin['ncpu'],ntasks,QMin['schedule_scaling'])
            nrounds = ntasks
            nslots = 1
            cpu_per_run = [QMin['ncpu']] * ntasks
        else:
            nrounds = 1
            nslots = QMin['ncpu']
            cpu_per_run = [1] * ntasks
        joblist.append({})

        icount = 0
        for grad in QMin['gradmap']:
            QMin3 = deepcopy(QMin2)
            QMin3['gradmap'] = [grad]
            QMin3['ncpu'] = cpu_per_run[icount]
            icount += 1
            joblist[-1]['grad_%i_%i' % grad] = QMin3
        QMin['nslots_pool'].append(nslots)


    elif QMin['gradmode'] == 2:
        # case of numerical gradients for ALL states, plus optionally gradients of DM and SOC
        # if only energy gradients:
        # -> do central point first, and n-1 displacements in parallel
        # -> do all other displacements afterwards
        QMin1 = deepcopy(QMin)
        QMin1['master'] = []
        QMin1['gradmap'] = []
        QMin1['master_displacement'] = []
        if 'h' not in QMin:
            QMin1['h'] = []
        remove = ['grad', 'dmdr']
        for r in remove:
            QMin1 = removekey(QMin1, r)
        if 'dmdr' in QMin:
            QMin1['dm'] = []
        if QMin['mpi_parallel']:
            QMin1['ncpu'] = QMin['ncpu']
        else:
            QMin1['ncpu'] = 1
        QMin['nslots_pool'] = [1]
        joblist.append({'master': QMin1})

        QMin2 = deepcopy(QMin)
        remove = ['comment', 'ncpu', 'veloc', 'grad', 'h', 'dm', 'dmdr']
        for r in remove:
            QMin2 = removekey(QMin2, r)
        QMin2['newstep'] = []
        QMin2['gradmap'] = []
        ntasks = 6 * QMin['natom']
        if QMin['mpi_parallel']:
            # nrounds,nslots,cpu_per_run=divide_slots(QMin['ncpu'],ntasks,QMin['schedule_scaling'])
            nrounds = ntasks
            nslots = 1
            cpu_per_run = [QMin['ncpu']] * ntasks
        else:
            nrounds = 1
            nslots = QMin['ncpu']
            cpu_per_run = [1] * ntasks
        QMin['nslots_pool'].append(nslots)

        #if 'socdr' in QMin or 'dmdr' in QMin:
        #    idispl=QMin['ncpu']-1
        #else:
        #    idispl=0
        icount = 0
        joblist.append({})
        for iatom in range(QMin['natom']):
            for ixyz in range(3):
                for isign in [-1., 1.]:
                    # idispl+=1
                    # if idispl==QMin['ncpu']:
                    #     joblist.append({})
                    QMin3 = deepcopy(QMin2)
                    QMin3 = doDisplacement(QMin3, [iatom, ixyz, isign], QMin['displ'])

                    # if 'socdr' in QMin or 'dmdr' in QMin or idispl>QMin['ncpu']:
                    QMin3['displacement'] = []
                    remove = ['init']
                    for r in remove:
                        QMin3 = removekey(QMin3, r)
                    if 'grad' in QMin:
                        QMin3['h'] = []
                    if 'dmdr' in QMin:
                        QMin3['dm'] = []
                    # if 'socdr' in QMin or 'dmdr' in QMin:
                    #QMin3['overlap'] = [[j + 1, i + 1] for i in range(QMin['nmstates']) for j in range(i + 1)]
                    QMin3['ncpu'] = cpu_per_run[icount]
                    icount += 1

                    jobname = 'displ_%i_%i_%s' % (iatom, ixyz, {-1.: 'p', 1.: 'n'}[isign])
                    joblist[-1][jobname] = QMin3

    if DEBUG:
        pprint.pprint(joblist, depth=3)
    return QMin, joblist

# ======================================================================= #


def run_calc(WORKDIR, QMin):
    err = 96
    irun = -1
    while err == 96:
        irun += 1
        try:
            Tasks = gettasks(QMin)
            setupWORKDIR(WORKDIR, Tasks, QMin)
            strip = 'keepintegrals' not in QMin
            err = runMNDO(WORKDIR, QMin['mndo'], QMin['driver'], QMin['ncpu'], strip)
        except Exception as problem:
            print('*' * 50 + '\nException in run_calc(%s)!' % (WORKDIR))
            traceback.print_exc()
            print('*' * 50 + '\n')
            raise problem
    return err

# ======================================================================= #


def runjobs(joblist, QMin):

    print('>>>>>>>>>>>>> Starting the MNDO job execution')

    errorcodes = {}
    for ijobset, jobset in enumerate(joblist):
        if not jobset:
            continue
        pool = Pool(processes=QMin['nslots_pool'][ijobset])
        for job in jobset:
            QMin1 = jobset[job]
            WORKDIR = os.path.join(QMin['scratchdir'], job)

            errorcodes[job] = pool.apply_async(run_calc, [WORKDIR, QMin1])
            # errorcodes[job]=run_calc(WORKDIR,QMin1)
            time.sleep(QMin['delay'])
        pool.close()
        pool.join()

        print('')

    for i in errorcodes:
        errorcodes[i] = errorcodes[i].get()

    if PRINT:
        string = '  ' + '=' * 40 + '\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '||' + ' ' * 10 + 'All Tasks completed!' + ' ' * 10 + '||\n'
        string += '||' + ' ' * 40 + '||\n'
        string += '  ' + '=' * 40 + '\n'
        print(string)
        j = 0
        string = 'Error Codes:\n\n'
        for i in errorcodes:
            string += '\t%s\t%i' % (i + ' ' * (10 - len(i)), errorcodes[i])
            j += 1
            if j == 4:
                j = 0
                string += '\n'
        print(string)

    if any((i != 0 for i in errorcodes.values())):
        print('Some subprocesses did not finish successfully!')
        # sys.exit(76)

    return errorcodes

# ======================================================================= #


def collectOutputs(joblist, QMin, errorcodes):

    QMout = {}

    for jobset in joblist:
        for job in jobset:
            if errorcodes[job] == 0:
                outfile = os.path.join(QMin['scratchdir'], job, 'MNDO.out')
                print('Reading %s' % (outfile))
                out = readfile(outfile)
                QMout[job] = getQMout(out, jobset[job])
            else:
                if 'master' in job or 'grad' in job:
                    print('Job %s did not finish sucessfully!' % (job))
                    sys.exit(77)
                elif 'displ' in job:
                    QMout[job] = get_zeroQMout(jobset[job])

    # if DEBUG:
        # pprint.pprint(QMout,width=130)
    if DEBUG:
        for i in sorted(QMout):
            QMout1 = QMout[i]
            for j in joblist:
                if i in j:
                    QMin1 = j[i]
                    break
            print('==============================> %s <==============================' % (i))
            printQMout(QMin1, QMout1)

    return QMout

# ======================================================================= #


def phase_correction(matrix):
    length = len(matrix)
    phase_corrected_matrix = [[.0 for x in range(length)] for x in range(length)]

    for i in range(length):
        diag = matrix[i][i].real

        # look if diag is significant and negative & switch phase
        if diag ** 2 > 0.5 and diag < 0:
            for j in range(length):
                phase_corrected_matrix[j][i] = matrix[j][i] * -1
        # otherwise leave values as is
        else:
            for j in range(length):
                phase_corrected_matrix[j][i] = matrix[j][i]

    return phase_corrected_matrix
# ======================================================================= #


def loewdin_orthonormalization(A):
    '''
    returns loewdin orthonormalized matrix
    '''

    # S = A^T * A
    S = np.dot(A.T, A)

    # S^d = U^T * S * U
    S_diag_only, U = np.linalg.eigh(S)

    # calculate the inverse sqrt of the diagonal matrix
    S_diag_only_inverse_sqrt = [1. / (float(d) ** 0.5) for d in S_diag_only]
    S_diag_inverse_sqrt = np.diag(S_diag_only_inverse_sqrt)

    # calculate inverse sqrt of S
    S_inverse_sqrt = np.dot(np.dot(U, S_diag_inverse_sqrt), U.T)

    # calculate loewdin orthonormalized matrix
    A_lo = np.dot(A, S_inverse_sqrt)

    # normalize A_lo
    A_lo = A_lo.T
    length = len(A_lo)
    A_lon = np.zeros((length, length), dtype=complex)

    for i in range(length):
        norm_of_col = np.linalg.norm(A_lo[i])
        A_lon[i] = [e / (norm_of_col ** 0.5) for e in A_lo[i]][0]

    return A_lon.T

# ======================================================================= #


def calculate_W_dQi(H, S, e_ref):
    '''
    Return diabatized H and S
    '''

    # get diagonal of Hamiltonian
    H = np.diag([e - e_ref for e in np.diag(H)])

    # do phase correction if necessary
    if any([x for x in np.diag(S) if x < 0]):
        S = phase_correction(S)

    # do loewdin orthonorm. on overlap matrix
    U = loewdin_orthonormalization(np.matrix(S))

    return np.dot(np.dot(U.T, H), U), U



# ======================================================================= #
def overlapsign(x):
    overlapthreshold = 0.8
    if abs(x) < overlapthreshold:
        return 0.0
    else:
        return math.copysign(1, x)

# ======================================================================= #


def numdiff(enp, enn, enc, displ, o1p, o2p, o1n, o2n, iatom, idir):
    o1p = overlapsign(o1p)
    o2p = overlapsign(o2p)
    o1n = overlapsign(o1n)
    o2n = overlapsign(o2n)

    enp *= o1p * o2p
    enn *= o1n * o2n

    if (o1p == 0.0 or o2p == 0.0) and (o1n == 0.0 or o2n == 0.0):
        print('Numerical differentiation failed, both displacements have bad overlap! iatom=%i, idir=%i' % (iatom, idir))
        sys.exit(78)
    if o1p == 0.0 or o2p == 0.0:
        print('Using one-sided NumDiff for iatom=%i, idir=%i. Retaining only negative displacement.' % (iatom, idir))
        g = (enc - enn) / displ
    elif o1n == 0.0 or o2n == 0.0:
        print('Using one-sided NumDiff for iatom=%i, idir=%i. Retaining only positive displacement.' % (iatom, idir))
        g = (enp - enc) / displ
    else:
        g = (enp - enn) / 2. / displ
    return -g

# ======================================================================= #


def arrangeQMout(QMin, QMoutall):

    # sys.exit(0)

    QMout = {}
    if 'h' in QMin or 'soc' in QMin:
        QMout['h'] = QMoutall['master']['h']
    if 'dm' in QMin:
        QMout['dm'] = QMoutall['master']['dm']

    if 'grad' in QMin:
        if QMin['gradmode'] == 0:
            QMout['grad'] = QMoutall['master']['grad']

        elif QMin['gradmode'] == 1:
            zerograd = [[0.0 for xyz in range(3)] for iatom in range(QMin['natom'])]
            grad = []
            for i in sorted(QMin['statemap']):
                mult, state, ms = tuple(QMin['statemap'][i])
                if (mult, state) in QMin['gradmap']:
                    name = 'grad_%i_%i' % (mult, state)
                    grad.append(QMoutall[name]['grad'][i - 1])
                else:
                    grad.append(zerograd)
            QMout['grad'] = grad

        elif QMin['gradmode'] == 2:
            grad = [[[0.0 for xyz in range(3)] for iatom in range(QMin['natom'])] for istate in range(QMin['nmstates'])]
            for iatom in range(QMin['natom']):
                for xyz in range(3):
                    namep = 'displ_%i_%i_p' % (iatom, xyz)
                    namen = 'displ_%i_%i_n' % (iatom, xyz)
                    displ = QMin['displ']

                    # diabatization
                    if QMin['template']['diab_num_grad']:
                        Hmaster = QMoutall['master']['h']
                        Hpos = deepcopy(QMoutall[namep]['h'])
                        Spos = deepcopy(QMoutall[namep]['overlap'])
                        Hneg = deepcopy(QMoutall[namen]['h'])
                        Sneg = deepcopy(QMoutall[namen]['overlap'])
                        Hpos, Spos = calculate_W_dQi(Hpos, Spos, QMout['h'][0][0])
                        Hneg, Sneg = calculate_W_dQi(Hneg, Sneg, QMout['h'][0][0])
                        QMoutall[namep]['h_diab'] = Hpos
                        QMoutall[namep]['s_diab'] = Spos
                        QMoutall[namen]['h_diab'] = Hneg
                        QMoutall[namen]['s_diab'] = Sneg
                    else:
                        Hmaster = QMoutall['master']['h']
                        Hpos = QMoutall[namep]['h']
                        Hneg = QMoutall[namen]['h']
                        Spos = QMoutall[namep]['overlap']
                        Sneg = QMoutall[namen]['overlap']


                    for istate in range(QMin['nmstates']):

                        enc = Hmaster[istate][istate].real

                        enp = Hpos[istate][istate].real
                        ovp = Spos[istate][istate].real

                        enn = Hneg[istate][istate].real
                        ovn = Sneg[istate][istate].real

                        g = numdiff(enp, enn, enc, displ, ovp, ovp, ovn, ovn, iatom, xyz)
                        grad[istate][iatom][xyz] = g
            QMout['grad'] = grad

    if 'dmdr' in QMin:
        dmdr = [[[[[0.0 for xyz in range(3)] for iatom in range(QMin['natom'])] for istate in range(QMin['nmstates'])] for jstate in range(QMin['nmstates'])] for ipol in range(3)]
        displ = QMin['displ']
        for iatom in range(QMin['natom']):
            for xyz in range(3):
                namep = 'displ_%i_%i_p' % (iatom, xyz)
                namen = 'displ_%i_%i_n' % (iatom, xyz)
                for ipol in range(3):
                    for istate in range(QMin['nmstates']):
                        for jstate in range(QMin['nmstates']):

                            # diabatization
                            if QMin['template']['diab_num_grad']:
                                Hmaster = QMoutall['master']['dm'][ipol]

                                Hpos = deepcopy(QMoutall[namep]['dm'][ipol])
                                Hpos = np.array(Hpos) - np.diag([e - QMout['h'][0][0] for e in np.diag(Hpos)])
                                Spos = QMoutall[namep]['s_diab']
                                Hpos = np.dot(np.dot(Spos.T, Hpos), Spos)

                                Hneg = deepcopy(QMoutall[namen]['dm'][ipol])
                                Hneg = np.array(Hneg) - np.diag([e - QMout['h'][0][0] for e in np.diag(Hneg)])
                                Sneg = QMoutall[namen]['s_diab']
                                Hneg = np.dot(np.dot(Sneg.T, Hneg), Sneg)
                            else:
                                Hmaster = QMoutall['master']['dm'][ipol]
                                Hpos = QMoutall[namep]['dm'][ipol]
                                Hneg = QMoutall[namen]['dm'][ipol]
                                Spos = QMoutall[namep]['overlap']
                                Sneg = QMoutall[namen]['overlap']

                            enc = Hmaster[istate][jstate].real

                            enp = Hpos[istate][jstate].real
                            o1p = Spos[istate][istate].real
                            o2p = Spos[jstate][jstate].real

                            enn = Hneg[istate][jstate].real
                            o1n = Sneg[istate][istate].real
                            o2n = Sneg[jstate][jstate].real

                            g = numdiff(enp, enn, enc, displ, o1p, o2p, o1n, o2n, iatom, xyz)
                            dmdr[ipol][istate][jstate][iatom][xyz] = g
        QMout['dmdr'] = dmdr

    if PRINT:
        print('\n===================================================================')
        print('========================= Final Results ===========================')
        print('===================================================================')
        printQMout(QMin, QMout)

    return QMout

# ======================================================================= #
def get_zeroQMout(QMin):
    nmstates = QMin['nmstates']
    natom = QMin['natom']
    QMout = {}
    if 'h' in QMin or 'soc' in QMin:
        QMout['h'] = [[complex(0.0) for i in range(nmstates)] for j in range(nmstates)]
    if 'dm' in QMin:
        QMout['dm'] = [[[complex(0.0) for i in range(nmstates)] for j in range(nmstates)] for xyz in range(3)]
    if 'overlap' in QMin:
        QMout['overlap'] = [[complex(0.0) for i in range(nmstates)] for j in range(nmstates)]
    if 'grad' in QMin:
        QMout['grad'] = [[[0., 0., 0.] for i in range(natom)] for j in range(nmstates)]
    return QMout







# ======================================================================= #
def cleanupSCRATCH(SCRATCHDIR):
    ''''''
    if PRINT:
        print('===> Removing directory %s\n' % (SCRATCHDIR))
    try:
        if True:
            shutil.rmtree(SCRATCHDIR)
        else:
            print('not removing anything. SCRATCHDIR is %s' % SCRATCHDIR)
    except OSError:
        print('Could not remove directory %s' % (SCRATCHDIR))

# ======================================================================= #

def stripWORKDIR(WORKDIR):
    ls = os.listdir(WORKDIR)
    keep = ['MNDO.out']
    for ifile in ls:
        delete = True
        for k in keep:
            if containsstring(k, ifile):
                delete = False
        if delete:
            rmfile = os.path.join(WORKDIR, ifile)
            if not DEBUG:
                if os.path.isdir(rmfile):
                    cleanupSCRATCH(rmfile)
                else:
                    os.remove(rmfile)

# ======================================================================= #



# =============================================================================================== #
# =============================================================================================== #
# =========================================== Main routine  ===================================== #
# =============================================================================================== #
# =============================================================================================== #

# ========================== Main Code =============================== #
def main():

    # Retrieve PRINT and DEBUG
    try:
        envPRINT = os.getenv('SH2CAS_PRINT')
        if envPRINT and envPRINT.lower() == 'false':
            global PRINT
            PRINT = False
        envDEBUG = os.getenv('SH2CAS_DEBUG')
        if envDEBUG and envDEBUG.lower() == 'true':
            global DEBUG
            DEBUG = True
    except ValueError:
        print('PRINT or DEBUG environment variables do not evaluate to numerical values!')
        sys.exit(90)

    # Process Command line arguments
    if len(sys.argv) != 2:
        print('Usage:\n./SHARC_MNDO.py <QMin>\n')
        print('version:', version)
        print('date:', versiondate)
        print('changelog:\n', changelogstring)
        sys.exit(91)
    QMinfilename = sys.argv[1]

    # Print header
    printheader()

    # Read QMinfile
    QMin = readQMin(QMinfilename)

    # make list of jobs
    QMin, joblist = generate_joblist(QMin)

    # run all MNDO jobs
    errorcodes = runjobs(joblist, QMin)

    # get output
    QMoutall = collectOutputs(joblist, QMin, errorcodes)

    # format final output
    QMout = arrangeQMout(QMin, QMoutall)

    # Measure time
    runtime = measuretime()
    QMout['runtime'] = runtime

    # Write QMout
    writeQMout(QMin, QMout, QMinfilename)

    # Remove Scratchfiles from SCRATCHDIR
    if not DEBUG:
        cleanupSCRATCH(QMin['scratchdir'])
        if 'cleanup' in QMin:
            cleanupSCRATCH(QMin['savedir'])
    if PRINT or DEBUG:
        print('#================ END ================#')


if __name__ == '__main__':
    main()






# kate: indent-width 4
