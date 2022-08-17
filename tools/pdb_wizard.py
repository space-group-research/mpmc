#!/usr/bin/python3
#
# PDB Wizard v0.3.0
# copyright Adam Hogan 2021-2022
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import copy
import os
import sys

import numpy as np
from numpy import array, cos, pi, sin, sqrt

list_of_elements = [
    "Ac",
    "Ag",
    "Al",
    "Am",
    "Ar",
    "As",
    "At",
    "Au",
    "B",
    "Ba",
    "Be",
    "Bh",
    "Bi",
    "Bk",
    "Br",
    "C",
    "Ca",
    "Cd",
    "Ce",
    "Cf",
    "Cl",
    "Cm",
    "Co",
    "Cr",
    "Cs",
    "Cu",
    "Db",
    "Dy",
    "Er",
    "Es",
    "Eu",
    "F",
    "Fe",
    "Fm",
    "Fr",
    "Ga",
    "Gd",
    "Ge",
    "H",
    "He",
    "Hf",
    "Hg",
    "Ho",
    "Hs",
    "I",
    "In",
    "Ir",
    "K",
    "Kr",
    "La",
    "Li",
    "Lr",
    "Lu",
    "Md",
    "Mg",
    "Mn",
    "Mo",
    "Mt",
    "N",
    "Na",
    "Nb",
    "Nd",
    "Ne",
    "Ni",
    "No",
    "Np",
    "O",
    "Os",
    "P",
    "Pa",
    "Pb",
    "Pd",
    "Pm",
    "Po",
    "Pr",
    "Pt",
    "Pu",
    "Ra",
    "Rb",
    "Re",
    "Rf",
    "Rh",
    "Rn",
    "Ru",
    "S",
    "Sb",
    "Sc",
    "Se",
    "Sg",
    "Si",
    "Sm",
    "Sn",
    "Sr",
    "Ta",
    "Tb",
    "Tc",
    "Te",
    "Th",
    "Ti",
    "Tl",
    "Tm",
    "U",
    "V",
    "W",
    "Xe",
    "Y",
    "Yb",
    "Zn",
    "Zr",
]

element_masses = {
    "H": 1.00797,
    "He": 4.0026,
    "Li": 6.941,
    "Be": 9.01218,
    "B": 10.81,
    "C": 12.011,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.998403,
    "Ne": 20.179,
    "Na": 22.98977,
    "Mg": 24.305,
    "Al": 26.98154,
    "Si": 28.0855,
    "P": 30.97376,
    "S": 32.06,
    "Cl": 35.453,
    "K": 39.0983,
    "Ar": 39.948,
    "Ca": 40.08,
    "Sc": 44.9559,
    "Ti": 47.9,
    "V": 50.9415,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.847,
    "Ni": 58.7,
    "Co": 58.9332,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.72,
    "Ge": 72.59,
    "As": 74.9216,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.8,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.9059,
    "Zr": 91.22,
    "Nb": 92.9064,
    "Mo": 95.94,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.9055,
    "Pd": 106.4,
    "Ag": 107.868,
    "Cd": 112.41,
    "In": 114.82,
    "Sn": 118.69,
    "Sb": 121.75,
    "I": 126.9045,
    "Te": 127.6,
    "Xe": 131.3,
    "Cs": 132.9054,
    "Ba": 137.33,
    "La": 138.9055,
    "Ce": 140.12,
    "Pr": 140.9077,
    "Nd": 144.24,
    "Pm": 145,
    "Sm": 150.4,
    "Eu": 151.96,
    "Gd": 157.25,
    "Tb": 158.9254,
    "Dy": 162.5,
    "Ho": 164.9304,
    "Er": 167.26,
    "Tm": 168.9342,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.9479,
    "W": 183.85,
    "Re": 186.207,
    "Os": 190.2,
    "Ir": 192.22,
    "Pt": 195.09,
    "Au": 196.9665,
    "Hg": 200.59,
    "Tl": 204.37,
    "Pb": 207.2,
    "Bi": 208.9804,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226.0254,
    "Ac": 227.0278,
    "Pa": 231.0359,
    "Th": 232.0381,
    "Np": 237.0482,
    "U": 238.029,
    "Pu": 242,
    "Am": 243,
    "Bk": 247,
    "Cm": 247,
    "No": 250,
    "Cf": 251,
    "Es": 252,
    "Hs": 255,
    "Mt": 256,
    "Fm": 257,
    "Md": 258,
    "Lr": 260,
    "Rf": 261,
    "Bh": 262,
    "Db": 262,
    "Sg": 263,
}


class Atom:
    def __init__(self, x, y, z, name):
        self.name = name.strip()
        self.x = array([float(x), float(y), float(z)])
        element = "".join([i for i in self.name[:2] if i.isalpha()])
        element = element.lower().capitalize()

        if element not in list_of_elements:
            element = element[0]
            if element not in list_of_elements:
                print("!!! Invalid element {} !!!".format(name))

        if element == "H":
            self.bond_r = 0.8
            self.vdw = 1.2
        elif element == "O":
            self.bond_r = 1.3
            self.vdw = 1.8
        elif element == "N" or element == "C":
            self.bond_r = 1.6
            self.vdw = 2.0
        else:
            self.bond_r = 2.0
            self.vdw = 3.0

        self.element = element
        self.charge = 0.0
        self.alpha = 0.0
        self.epsilon = 0.0
        self.sigma = 0.0
        self.c6 = 0.0
        self.c8 = 0.0
        self.c10 = 0.0
        self.mass = element_masses[self.element]
        self.id = 0

    def __str__(self):
        return "Atom instance {} {}".format(self.element, self.id)


class PBC:
    def __init__(self, a, b, c, alpha, beta, gamma):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        basis00 = a
        basis01 = 0.0
        basis02 = 0.0
        basis10 = b * cos(pi / 180.0 * gamma)
        basis11 = b * sin(pi / 180.0 * gamma)
        basis12 = 0.0
        basis20 = c * cos(pi / 180.0 * beta)
        basis21 = ((b * c * cos(pi / 180.0 * alpha)) - (basis10 * basis20)) / basis11
        basis22 = sqrt(c * c - basis20 * basis20 - basis21 * basis21)

        self.basis_matrix = array(
            [
                [basis00, basis01, basis02],
                [basis10, basis11, basis12],
                [basis20, basis21, basis22],
            ]
        )

        self.volume = basis00 * (basis11 * basis22 - basis12 * basis21)
        self.volume += basis01 * (basis12 * basis20 - basis10 * basis22)
        self.volume += basis02 * (basis10 * basis21 - basis11 * basis20)

        self.inverse_volume = 1.0 / self.volume

        reciprocal_basis00 = self.inverse_volume * (
            basis11 * basis22 - basis12 * basis21
        )
        reciprocal_basis01 = self.inverse_volume * (
            basis02 * basis21 - basis01 * basis22
        )
        reciprocal_basis02 = self.inverse_volume * (
            basis01 * basis12 - basis02 * basis11
        )
        reciprocal_basis10 = self.inverse_volume * (
            basis12 * basis20 - basis10 * basis22
        )
        reciprocal_basis11 = self.inverse_volume * (
            basis00 * basis22 - basis02 * basis20
        )
        reciprocal_basis12 = self.inverse_volume * (
            basis02 * basis10 - basis00 * basis12
        )
        reciprocal_basis20 = self.inverse_volume * (
            basis10 * basis21 - basis11 * basis20
        )
        reciprocal_basis21 = self.inverse_volume * (
            basis01 * basis20 - basis00 * basis21
        )
        reciprocal_basis22 = self.inverse_volume * (
            basis00 * basis11 - basis01 * basis10
        )

        self.reciprocal_basis_matrix = array(
            [
                [reciprocal_basis00, reciprocal_basis10, reciprocal_basis20],
                [reciprocal_basis01, reciprocal_basis11, reciprocal_basis21],
                [reciprocal_basis02, reciprocal_basis12, reciprocal_basis22],
            ]
        )

    def update(self, a, b, c, alpha, beta, gamma):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        basis00 = a
        basis01 = 0.0
        basis02 = 0.0
        basis10 = b * cos(pi / 180.0 * gamma)
        basis11 = b * sin(pi / 180.0 * gamma)
        basis12 = 0.0
        basis20 = c * cos(pi / 180.0 * beta)
        basis21 = ((b * c * cos(pi / 180.0 * alpha)) - (basis10 * basis20)) / basis11
        basis22 = sqrt(c * c - basis20 * basis20 - basis21 * basis21)

        self.basis_matrix = array(
            [
                [basis00, basis01, basis02],
                [basis10, basis11, basis12],
                [basis20, basis21, basis22],
            ]
        )

        self.volume = basis00 * (basis11 * basis22 - basis12 * basis21)
        self.volume += basis01 * (basis12 * basis20 - basis10 * basis22)
        self.volume += basis02 * (basis10 * basis21 - basis11 * basis20)

        self.inverse_volume = 1.0 / self.volume

        reciprocal_basis00 = self.inverse_volume * (
            basis11 * basis22 - basis12 * basis21
        )
        reciprocal_basis01 = self.inverse_volume * (
            basis02 * basis21 - basis01 * basis22
        )
        reciprocal_basis02 = self.inverse_volume * (
            basis01 * basis12 - basis02 * basis11
        )
        reciprocal_basis10 = self.inverse_volume * (
            basis12 * basis20 - basis10 * basis22
        )
        reciprocal_basis11 = self.inverse_volume * (
            basis00 * basis22 - basis02 * basis20
        )
        reciprocal_basis12 = self.inverse_volume * (
            basis02 * basis10 - basis00 * basis12
        )
        reciprocal_basis20 = self.inverse_volume * (
            basis10 * basis21 - basis11 * basis20
        )
        reciprocal_basis21 = self.inverse_volume * (
            basis01 * basis20 - basis00 * basis21
        )
        reciprocal_basis22 = self.inverse_volume * (
            basis00 * basis11 - basis01 * basis10
        )

        self.reciprocal_basis_matrix = array(
            [
                [reciprocal_basis00, reciprocal_basis10, reciprocal_basis20],
                [reciprocal_basis01, reciprocal_basis11, reciprocal_basis21],
                [reciprocal_basis02, reciprocal_basis12, reciprocal_basis22],
            ]
        )

    def min_image(self, dx):
        img = np.matmul(dx, self.reciprocal_basis_matrix)
        img = np.round(img)
        di = np.matmul(img, self.basis_matrix)
        dx_return = dx - di
        r = np.sqrt(np.dot(dx_return, dx_return))
        return r

    def wrap(self, dx):
        img = np.matmul(dx, self.reciprocal_basis_matrix)
        img = np.round(img)
        di = np.matmul(img, self.basis_matrix)
        dx_return = dx - di
        return dx_return

    def wrap_forward(self, dx):
        img = np.matmul(dx, self.reciprocal_basis_matrix)
        img = np.floor(img)
        di = np.matmul(img, self.basis_matrix)
        dx_return = dx - di
        return dx_return


def progressbar(it, prefix="", size=60, out=sys.stdout):
    count = len(it)

    def show(j):
        x = int(size * j / count)
        print(
            "{}[{}{}] {}/{}".format(prefix, "█" * x, "." * (size - x), j, count),
            end="\r",
            file=out,
            flush=True,
        )

    show(0)

    for i, item in enumerate(it):
        yield item
        show(i + 1)

    print("", flush=True, file=out)


def get_forcefield(name):
    ffs = []
    opls_aa_uff = {
        "H": {
            "mass": 1.0079,
            "alpha": 0.41380,
            "sigma": 2.42,
            "epsilon": 15.11,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "C": {
            "mass": 12.011,
            "alpha": 1.2866,
            "sigma": 3.55,
            "epsilon": 35.25,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "N": {
            "mass": 14.0067,
            "alpha": 0.97157,
            "sigma": 3.25000,
            "epsilon": 85.60000,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "O": {
            "mass": 15.999,
            "alpha": 0.852,
            "sigma": 3.118,
            "epsilon": 30.19,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "P": {
            "mass": 30.974,
            "alpha": 3.35,
            "sigma": 3.69456,
            "epsilon": 153.48197,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Cl": {
            "mass": 35.4527,
            "alpha": 2.40028,
            "sigma": 3.516377,
            "epsilon": 114.23084,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Zn": {
            "mass": 65.38,
            "alpha": 1.98870,
            "sigma": 2.46155,
            "epsilon": 62.39923,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Pt": {
            "mass": 195.084,
            "alpha": 8.56281,
            "sigma": 2.45353,
            "epsilon": 40.25756,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Pd": {
            "mass": 106.42,
            "alpha": 5.25926,
            "sigma": 2.5827,
            "epsilon": 24.1545,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "B": {
            "mass": 10.811,
            "alpha": 0.6634,
            "sigma": 3.63754,
            "epsilon": 90.57952,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Ni": {
            "mass": 58.69340,
            "alpha": 0.38980,
            "sigma": 2.52480,
            "epsilon": 7.55330,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Cr": {
            "mass": 51.99610,
            "alpha": 3.50740,
            "sigma": 2.69320,
            "epsilon": 7.55330,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Cu": {
            "mass": 63.54630,
            "alpha": 2.19630,
            "sigma": 3.11400,
            "epsilon": 2.51600,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Co": {
            "mass": 58.93320,
            "alpha": 3.26440,
            "sigma": 2.55870,
            "epsilon": 7.04980,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "F": {
            "mass": 18.99800,
            "alpha": 0.444747,
            "sigma": 2.996983,
            "epsilon": 25.160979,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Si": {
            "mass": 28.08500,
            "alpha": 2.133000,
            "sigma": 3.826410,
            "epsilon": 202.294269,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Br": {
            "mass": 79.90400,
            "alpha": 3.493000,
            "sigma": 3.732000,
            "epsilon": 126.392600,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "Ti": {
            "mass": 47.86710,
            "alpha": 3.24280,
            "sigma": 2.82860,
            "epsilon": 8.56050,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "In": {
            "mass": 114.8180,
            "alpha": 2.00000,
            "sigma": 3.97600,
            "epsilon": 301.40000,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
        "W": {
            "mass": 183.84000,
            "alpha": 3.65453,
            "sigma": 2.73420,
            "epsilon": 33.73830,
            "c6": 0.0,
            "c8": 0.0,
            "c10": 0.0,
        },
    }
    ffs.append(opls_aa_uff)

    phahst = {
        "Cu": {
            "mass": 63.54630,
            "alpha": 0.29252,
            "sigma": 2.73851,
            "epsilon": 8.82345,
            "c6": 6.96956,
            "c8": 262.82938,
            "c10": 13951.49740,
        },
        "C": {
            "mass": 12.01100,
            "alpha": 0.71317,
            "sigma": 3.35929,
            "epsilon": 4.00147,
            "c6": 11.88969,
            "c8": 547.51694,
            "c10": 27317.97855,
        },
        "O": {
            "mass": 15.99900,
            "alpha": 1.68064,
            "sigma": 3.23867,
            "epsilon": 3.89544,
            "c6": 27.70093,
            "c8": 709.36452,
            "c10": 19820.89339,
        },
        "H": {
            "mass": 1.00790,
            "alpha": 0.02117,
            "sigma": 1.87446,
            "epsilon": 3.63874,
            "c6": 0.16278,
            "c8": 5.03239,
            "c10": 202.99322,
        },
    }
    ffs.append(phahst)
    return ffs[name]


def set_atom_ids(system):
    for ind, atom in enumerate(system):
        atom.id = ind + 1


def gcd_list(my_list):
    result = my_list[0]
    for x in my_list[1:]:
        if result < x:
            temp = result
            result = x
            x = temp
        while x != 0:
            temp = x
            x = result % x
            result = temp
    return result


def read_pdb(system, filename):
    file = open(filename, "r")
    line = file.readline()
    pbc = None
    try:
        for line in progressbar(file.readlines(), "Reading file {} ".format(filename)):
            tokens = line.split()

            # ignore these lines
            if (
                tokens[0] == "REMARK"
                or tokens[0] == "HEADER"
                or tokens[0] == "MODEL"
                or tokens[0] == "TER"
                or tokens[0] == "TITLE"
                or tokens[0] == "AUTHOR"
                or tokens[0] == "EXPDTA"
                or tokens[0] == "SEQRES"
                or tokens[0] == "ANISOU"
                or tokens[0] == "SCALE1"
                or tokens[0] == "SCALE2"
                or tokens[0] == "SCALE3"
            ):
                continue

            # stop at END or ENDMDL
            if tokens[0] == "END" or tokens[0] == "ENDMDL":
                break

            # create atoms
            if tokens[0] == "ATOM" or tokens[0] == "HETATM":
                atom = Atom(line[30:38], line[38:46], line[46:54], line[12:16])
                system.append(atom)

            # load box
            if tokens[0] == "CRYST1":
                a = float(tokens[1])
                b = float(tokens[2])
                c = float(tokens[3])
                alpha = float(tokens[4])
                beta = float(tokens[5])
                gamma = float(tokens[6])
                pbc = PBC(a, b, c, alpha, beta, gamma)

    except ValueError:
        sys.exit("Error reading line:\n{}".format(line))

    print("")
    set_atom_ids(system)
    file.close()
    return pbc


def read_xyz(system, filename):
    file = open(filename, "r")
    file.readline()
    line = file.readline()
    pbc = None

    try:
        tokens = line.split()
        if len(tokens) != 6:
            raise ValueError
        a = float(tokens[0])
        b = float(tokens[1])
        c = float(tokens[2])
        alpha = float(tokens[3])
        beta = float(tokens[4])
        gamma = float(tokens[5])
        pbc = PBC(a, b, c, alpha, beta, gamma)
    except ValueError:
        print("Couldn't locate a b c alpha beta gamma on second line of .xyz file")

    try:
        for line in progressbar(file.readlines(), "Reading file {} ".format(filename)):

            # ignore blank lines
            if line == "" or line == "\n":
                continue

            tokens = line.split()
            atom = Atom(tokens[1], tokens[2], tokens[3], tokens[0])

            # attempt to read the fifth column as the charge, pass otherwise
            try:
                charge = tokens[4]
                atom.charge = float(charge)
            except (ValueError, IndexError):
                pass

            system.append(atom)

    except ValueError:
        sys.exit("Error reading line:\n{}".format(line))

    set_atom_ids(system)
    file.close()
    return pbc


def overlap_detector(system, pbc):
    messages = []
    overlapping_atoms = True
    while overlapping_atoms:
        overlapping_atoms = False
        for atom in progressbar(system, "Detecting overlapping atoms "):
            for atom2 in system:
                if atom.id != atom2.id and not overlapping_atoms:
                    dx = atom.x - atom2.x
                    r = pbc.min_image(dx)
                    if r < 0.05:
                        overlapping_atoms = True
                        messages.append(
                            "Deleting overlapping atoms\n{:>3} {:>5} --- {:>3} {:>5}\n {}\n {}".format(
                                atom.element,
                                atom.id,
                                atom2.element,
                                atom2.id,
                                atom.x,
                                atom2.x,
                            )
                        )
                        system.remove(atom2)
    set_atom_ids(system)
    for message in messages:
        print(message)


def list_close_contacts(system, pbc):
    messages = []
    for atom in progressbar(system):
        for atom2 in system:
            if atom2.id > atom.id:
                dx = atom.x - atom2.x
                r = pbc.min_image(dx)
                vdw = 0.5 * (atom.vdw + atom2.vdw)
                bond_r = 0.5 * (atom.bond_r + atom2.bond_r)
                if vdw > r > bond_r:
                    element_string = "{}-{}".format(atom.element, atom2.element)
                    messages.append(
                        "{:<5} {:>5} {:>5}   r = {}".format(
                            element_string, atom.id, atom2.id, round(r, 6)
                        )
                    )
    print("\nClose contacts:\n")
    for message in messages:
        print(message)


def list_bonds(system, pbc):
    messages = []
    for atom in progressbar(system):
        for atom2 in system:
            if atom2.id > atom.id:
                dx = atom.x - atom2.x
                r = pbc.min_image(dx)
                bond_r = 0.5 * (atom.bond_r + atom2.bond_r)
                if r < bond_r:
                    element_string = "{}-{}".format(atom.element, atom2.element)
                    messages.append(
                        "{:<5} {:>5} {:>5}   r = {}".format(
                            element_string, atom.id, atom2.id, round(r, 6)
                        )
                    )
    print("\nBonded atoms:\n")
    for message in messages:
        print(message)


def list_angles(system, pbc):
    messages = []
    mols = find_molecules(system, pbc)
    for mol in progressbar(mols):
        for atom in mol:
            for atom2 in mol:
                for atom3 in mol:
                    if (
                        atom2.id == atom.id
                        or atom3.id == atom2.id
                        or atom3.id == atom.id
                    ):
                        continue
                    dx1 = atom2.x - atom.x
                    r1 = pbc.min_image(dx1)
                    bond_r1 = 0.5 * (atom.bond_r + atom2.bond_r)
                    dx2 = atom2.x - atom3.x
                    r2 = pbc.min_image(dx2)
                    bond_r2 = 0.5 * (atom2.bond_r + atom3.bond_r)
                    if r1 < bond_r1 and r2 < bond_r2:
                        u_dx1 = pbc.wrap(dx1)
                        u_dx2 = pbc.wrap(dx2)
                        u_dx1 = u_dx1 / np.sqrt(np.dot(u_dx1, u_dx1))
                        u_dx2 = u_dx2 / np.sqrt(np.dot(u_dx2, u_dx2))
                        bond_angle = np.degrees(
                            np.arccos(np.clip(np.dot(u_dx1, u_dx2), -1.0, 1.0))
                        )
                        element_string = "{}-{}-{}".format(
                            atom.element, atom2.element, atom3.element
                        )
                        messages.append(
                            "{:<7} {:>5} {:>5} {:>5}   angle = {:>6} r1, r2 = {:>6}, {:>6}".format(
                                element_string,
                                atom.id,
                                atom2.id,
                                atom3.id,
                                np.round(bond_angle, 2),
                                np.round(r1, 3),
                                np.round(r2, 3),
                            )
                        )
    print("\nAngles:\n")
    for message in messages:
        print(message)


def list_lone_atoms(system, pbc):
    lone_atoms = []
    for atom in progressbar(system):
        lone_atom = True
        for atom2 in system:
            if atom2.id != atom.id:
                dx = atom.x - atom2.x
                r = pbc.min_image(dx)
                vdw = 0.5 * (atom.vdw + atom2.vdw)
                if r < vdw:
                    lone_atom = False
        if lone_atom:
            lone_atoms.append(atom)
    if len(lone_atoms) == 0:
        print("\nNo lone atoms found\n")
    else:
        print("\nLone atoms:\n")
        for atom in lone_atoms:
            print("{:>3} {:>5} {}".format(atom.element, atom.id, atom.x))


def delete_lone_atoms(system, pbc):
    lone_atoms = []
    for atom in progressbar(system):
        lone_atom = True
        for atom2 in system:
            if atom2.id != atom.id:
                dx = atom.x - atom2.x
                r = pbc.min_image(dx)
                vdw = 0.5 * (atom.vdw + atom2.vdw)
                if r < vdw:
                    lone_atom = False
        if lone_atom:
            lone_atoms.append(atom)
    if len(lone_atoms) == 0:
        print("\nNo lone atoms found\n")
    else:
        print("\nDeleting lone atoms:\n")
        for atom in lone_atoms:
            print("{:>3} {:>5} {}".format(atom.element, atom.id, atom.x))
            system.remove(atom)


def write_xyz(system, pbc, filename):
    out = open(filename, "w")
    out.write(str(len(system)))
    out.write(
        "\n{} {} {} {} {} {}\n".format(
            pbc.a, pbc.b, pbc.c, pbc.alpha, pbc.beta, pbc.gamma
        )
    )
    for atom in system:
        out.write("{} {} {} {}\n".format(atom.element, atom.x[0], atom.x[1], atom.x[2]))
    out.close()
    print("Wrote {}".format(filename))


def write_standard_pdb(system, pbc, filename):
    mols = find_molecules(system, pbc)

    out = open(filename, "w")
    out.write("MODEL        1\n")
    out.write("COMPND    {:<69}\n".format(filename))
    out.write("AUTHOR    GENERATED BY PDB WIZARD\n")
    out.write(
        "CRYST1  {:>7}  {:7}  {:7} {:>6} {:>6} {:>6} P 1           1\n".format(
            round(pbc.a, 3),
            round(pbc.b, 3),
            round(pbc.c, 3),
            round(pbc.alpha, 2),
            round(pbc.beta, 2),
            round(pbc.gamma, 2),
        )
    )

    atom_id = 1
    for idx, mol in enumerate(mols):

        mol.sort(key=lambda atom: atom.element)
        mol_name = "UNK"
        mol_elements = [atom.element for atom in mol]
        if mol_elements == ["H", "H", "O"]:
            mol_name = "HOH"
        elif mol_elements == ["He"]:
            mol_name = "HE"
        elif mol_elements == ["Ne"]:
            mol_name = "NE"
        elif mol_elements == ["Ar"]:
            mol_name = "AR"
        elif mol_elements == ["Kr"]:
            mol_name = "KR"
        elif mol_elements == ["Xe"]:
            mol_name = "XE"
        elif mol_elements == ["He"]:
            mol_name = "HE"
        elif mol_elements == ["Zn"]:
            mol_name = "ZNA"
        elif mol_elements == ["Cl", "Cl", "Cl", "Cl", "Zn"]:
            mol_name = "ZNC"

        base_atom = mol[-1]
        additional_tag = 1

        for atom in mol:
            atom.id = atom_id

            other_elements = [
                other_atom.element for other_atom in mol if atom is not other_atom
            ]
            if atom.element in other_elements:
                atom.name = atom.element + str(additional_tag)

            dx = atom.x - base_atom.x

            dx = pbc.wrap(dx)

            atom.x = base_atom.x + dx

            out.write(
                "HETATM {:>4}  {:<3} {:>3} A {:>4}    {:>7} {:>7} {:>7}  1.00  0.00          {:>2}\n".format(
                    atom.id,
                    atom.name,
                    mol_name,
                    idx + 1,
                    round(atom.x[0], 3),
                    round(atom.x[1], 3),
                    round(atom.x[2], 3),
                    atom.element,
                )
            )

            atom_id += 1
            additional_tag += 1

    for idx, mol in enumerate(mols):

        mol.sort(key=lambda atom: atom.element)

        mol_elements = [atom.element for atom in mol]

        if mol_elements == ["Cl", "Cl", "Cl", "Cl", "Zn"]:
            out.write("CONECT {:>4} {:>4}\n".format(mol[0].id, mol[4].id))
            out.write("CONECT {:>4} {:>4}\n".format(mol[1].id, mol[4].id))
            out.write("CONECT {:>4} {:>4}\n".format(mol[2].id, mol[4].id))
            out.write("CONECT {:>4} {:>4}\n".format(mol[3].id, mol[4].id))

    out.write("END\n")
    out.close()
    set_atom_ids(system)
    print("Wrote {}".format(filename))


def find_molecules(system, pbc):
    set_atom_ids(system)
    merge_lists = []
    for atom1 in progressbar(system, "Finding molecules "):
        for atom2 in system:
            if atom2.id != atom1.id:
                dx = atom1.x - atom2.x
                r = pbc.min_image(dx)
                bond_r = 0.5 * (atom1.bond_r + atom2.bond_r)
                if r < bond_r:
                    found = False
                    for merge_list in merge_lists:
                        if atom1.id in merge_list:
                            if atom2.id not in merge_list:
                                merge_list.append(atom2.id)
                            found = True
                        elif atom2.id in merge_list:
                            if atom1.id not in merge_list:
                                merge_list.append(atom1.id)
                            found = True
                    if not found:
                        merge_lists.append([atom1.id, atom2.id])

    for atom in system:
        lone = True
        for merge_list in merge_lists:
            if atom.id in merge_list:
                lone = False
        if lone:
            merge_lists.append([atom.id])

    mols = [
        [system[atom_id - 1] for atom_id in merge_list] for merge_list in merge_lists
    ]

    return mols


def apply_ff_to_system(system, ff):
    for atom in system:
        try:
            atom.mass = ff[atom.element]["mass"]
            atom.alpha = ff[atom.element]["alpha"]
            atom.epsilon = ff[atom.element]["epsilon"]
            atom.sigma = ff[atom.element]["sigma"]
            atom.c6 = ff[atom.element]["c6"]
            atom.c8 = ff[atom.element]["c8"]
            atom.c10 = ff[atom.element]["c10"]
        except KeyError:
            print(
                "!!! atom {} not found in forcefield, parameters set to all zeros !!!".format(
                    atom.element
                )
            )


def write_mpmc_pdb(system, pbc, filename, write_charges=False, write_params=False):
    out = open(filename, "w")
    out.write(
        "CRYST1  {:>7}  {:7}  {:7} {:>6} {:>6} {:>6} P 1           1\n".format(
            round(pbc.a, 3),
            round(pbc.b, 3),
            round(pbc.c, 3),
            round(pbc.alpha, 2),
            round(pbc.beta, 2),
            round(pbc.gamma, 2),
        )
    )
    for atom in progressbar(system):
        out.write(
            "ATOM {:>6} {:<4} MOF F    1    {:>7} {:>7} {:>7}".format(
                atom.id,
                atom.element,
                round(atom.x[0], 3),
                round(atom.x[1], 3),
                round(atom.x[2], 3),
                atom.element,
            )
        )
        if write_params is True:
            out.write(" {:>9.6}".format(atom.mass))
        if write_charges is True or write_params is True:
            out.write(" {:>8.4}".format(atom.charge))
        if write_params is True:
            out.write(
                " {:>8.4} {:>8.4} {:>8.4} 0.0 0.0 {:>8.4} {:>10.4} {:>10.2}\n".format(
                    atom.alpha, atom.epsilon, atom.sigma, atom.c6, atom.c8, atom.c10
                )
            )
        else:
            out.write(" xxx{}xxx\n".format(atom.name.strip()))
    borders = [
        array([0, 0, 0]),
        array([1, 0, 0]),
        array([0, 1, 0]),
        array([0, 0, 1]),
        array([1, 1, 0]),
        array([1, 0, 1]),
        array([0, 1, 1]),
        array([1, 1, 1]),
    ]
    for ind, pos in enumerate(borders):
        border_pos = np.matmul(pos, pbc.basis_matrix)
        out.write(
            "ATOM {:>6} {:<4} BOX F    2    {:>7} {:>7} {:>7} 0.0 0.0 0.0 0.0 0.0"
            "\n".format(
                len(system) + ind,
                "X",
                round(border_pos[0], 3),
                round(border_pos[1], 3),
                round(border_pos[2], 3),
            )
        )
    connections = [
        [0, 1],
        [0, 2],
        [0, 3],
        [1, 4],
        [1, 5],
        [2, 4],
        [2, 6],
        [4, 7],
        [5, 7],
        [6, 7],
        [3, 6],
        [3, 5],
    ]
    for connection in connections:
        out.write(
            "CONECT {:>4} {:>4}\n".format(
                len(system) + connection[0], len(system) + connection[1]
            )
        )
    out.write(
        "REMARK BOX BASIS[0] = {:20.14f} {:20.14f} {:20.14f}\n".format(
            pbc.basis_matrix[0][0], pbc.basis_matrix[0][1], pbc.basis_matrix[0][2]
        )
    )
    out.write(
        "REMARK BOX BASIS[1] = {:20.14f} {:20.14f} {:20.14f}\n".format(
            pbc.basis_matrix[1][0], pbc.basis_matrix[1][1], pbc.basis_matrix[1][2]
        )
    )
    out.write(
        "REMARK BOX BASIS[2] = {:20.14f} {:20.14f} {:20.14f}\n".format(
            pbc.basis_matrix[2][0], pbc.basis_matrix[2][1], pbc.basis_matrix[2][2]
        )
    )
    out.write("END\n")
    out.close()
    print("Wrote {}".format(filename))


def print_info(system, pbc, filename):
    print("")
    print(
        r"   ___  ___  ___   __    __ _                  _ " + "\n"
        r"  / _ \/   \/ __\ / / /\ \ (_)______ _ _ __ __| |" + "\n"
        r" / /_)/ /\ /__\// \ \/  \/ / |_  / _` | '__/ _` |" + "\n"
        r"/ ___/ /_// \/  \  \  /\  /| |/ / (_| | | | (_| |" + "\n"
        r"\/  /___,'\_____/   \/  \/ |_/___\__,_|_|  \__,_|"
    )
    print("\nfilename: {}".format(filename))
    print(
        "\nCell:\n{:>7}  {:7}  {:7} {:>6} {:>6} {:>6}\n".format(
            round(pbc.a, 3),
            round(pbc.b, 3),
            round(pbc.c, 3),
            round(pbc.alpha, 2),
            round(pbc.beta, 2),
            round(pbc.gamma, 2),
        )
    )
    print(
        "{:20.14f} {:20.14f} {:20.14f}".format(
            pbc.basis_matrix[0][0], pbc.basis_matrix[0][1], pbc.basis_matrix[0][2]
        )
    )
    print(
        "{:20.14f} {:20.14f} {:20.14f}".format(
            pbc.basis_matrix[1][0], pbc.basis_matrix[1][1], pbc.basis_matrix[1][2]
        )
    )
    print(
        "{:20.14f} {:20.14f} {:20.14f}".format(
            pbc.basis_matrix[2][0], pbc.basis_matrix[2][1], pbc.basis_matrix[2][2]
        )
    )
    print(
        "Volume: {:10.2f} A^3 Density: {:10.4} g/cm^3".format(
            pbc.volume, sum([atom.mass for atom in system]) * 1.66054 / pbc.volume
        )
    )
    print_formula_unit(system)


def print_formula_unit(system):
    atom_dict = {}
    for atom in system:
        if atom.element not in atom_dict:
            atom_dict[atom.element] = 1
        else:
            atom_dict[atom.element] += 1
    print("\nTotal number of atoms:\n")
    for ele in atom_dict:
        print("{} {}".format(ele, atom_dict[ele]))
    atom_n = [atom_dict[i] for i in atom_dict]
    atoms_gcd = gcd_list(atom_n)
    print("\nFormula unit\n")
    for ele in atom_dict:
        print("{} {}".format(ele, int(atom_dict[ele] / atoms_gcd)))


def wrapall_forward(system, pbc):
    for atom in progressbar(system):
        atom.x = pbc.wrap_forward(atom.x)
    print("\nWrapped atoms forward of origin")


def wrapall(system, pbc):
    for atom in progressbar(system):
        atom.x = pbc.wrap(atom.x)
    print("\nWrapped atoms around origin")


def extend_axis(system, pbc):
    print(
        "\nCurrent cell:\n{:>7}  {:7}  {:7} {:>6} {:>6} {:>6}\n".format(
            round(pbc.a, 3),
            round(pbc.b, 3),
            round(pbc.c, 3),
            round(pbc.alpha, 2),
            round(pbc.beta, 2),
            round(pbc.gamma, 2),
        )
    )
    print(
        "{:20.14f} {:20.14f} {:20.14f}".format(
            pbc.basis_matrix[0][0], pbc.basis_matrix[0][1], pbc.basis_matrix[0][2]
        )
    )
    print(
        "{:20.14f} {:20.14f} {:20.14f}".format(
            pbc.basis_matrix[1][0], pbc.basis_matrix[1][1], pbc.basis_matrix[1][2]
        )
    )
    print(
        "{:20.14f} {:20.14f} {:20.14f}".format(
            pbc.basis_matrix[2][0], pbc.basis_matrix[2][1], pbc.basis_matrix[2][2]
        )
    )
    while True:
        try:
            axis = input(
                "\nWhat axis would you like to extend? 0, 1, 2 or x, y, z or q(uit)\n\n> "
            )
            if (
                axis == "q"
                or axis == "quit"
                or axis == "Q"
                or axis == "Quit"
                or axis == "QUIT"
            ):
                return
            if axis == "x" or axis == "X":
                axis = 0
            elif axis == "y" or axis == "Y":
                axis = 1
            elif axis == "z" or axis == "Z":
                axis = 2
            axis = int(axis)
            if axis < 0 or axis > 2:
                raise ValueError
            times = input("\nHow many times would you like to extend it?\n\n> ")
            times = int(times)
            if times < 1:
                raise ValueError
            break
        except ValueError:
            print("!!! Error converting input to int or x, y, z !!!")

    new_atoms = []
    for i in np.arange(times):
        for atom in system:
            new_atom = copy.deepcopy(atom)
            if axis == 0:
                new_atom.x += (i + 1) * pbc.basis_matrix[0]
            elif axis == 1:
                new_atom.x += (i + 1) * pbc.basis_matrix[1]
            elif axis == 2:
                new_atom.x += (i + 1) * pbc.basis_matrix[2]
            new_atoms.append(new_atom)

    for atom in new_atoms:
        system.append(atom)

    if axis == 0:
        pbc.update((times + 1) * pbc.a, pbc.b, pbc.c, pbc.alpha, pbc.beta, pbc.gamma)
    elif axis == 1:
        pbc.update(pbc.a, (times + 1) * pbc.b, pbc.c, pbc.alpha, pbc.beta, pbc.gamma)
    elif axis == 2:
        pbc.update(pbc.a, pbc.b, (times + 1) * pbc.c, pbc.alpha, pbc.beta, pbc.gamma)
    set_atom_ids(system)


def list_coords(system):
    for atom in system:
        print("{} {}".format(atom.element, atom.x))


def vmd_preview(system, pbc):
    write_mpmc_pdb(system, pbc, "pdb_wizard.tmp.pdb")
    os.system("vmd pdb_wizard.tmp.pdb")
    os.system("rm pdb_wizard.tmp.pdb")


def edit_h_dist(system, pbc):
    while True:
        second_element = "XX"
        distance = 0
        try:
            second_element = input(
                "\nLook for hydrogens bonded with which element? (e.g. C, O, N, etc)\n\n> "
            )
            if second_element not in list_of_elements:
                raise ValueError
            distance = input(
                "\nWhat distance (in angstroms) shall {}-H bonds be set to?\n\n> ".format(
                    second_element
                )
            )
            distance = float(distance)
            break
        except ValueError:
            print("!!! Error finding element or reading distance !!!")
    messages = []
    for atom in progressbar(system):
        for atom2 in system:
            if atom2.id > atom.id:
                if (atom.element == "H" and atom2.element == second_element) or (
                    atom2.element == "H" and atom.element == second_element
                ):
                    if atom.element == "H":
                        h_atom = atom
                        other_atom = atom2
                    else:
                        h_atom = atom2
                        other_atom = atom
                    dx = h_atom.x - other_atom.x
                    r = pbc.min_image(dx)
                    bond_r = 0.5 * (atom.bond_r + atom2.bond_r)
                    if r < bond_r:
                        element_string = "{}-{}".format(atom.element, atom2.element)
                        messages.append(
                            "{:<5} {:>5} {:>5}".format(
                                element_string, atom.id, atom2.id
                            )
                        )
                        dx = pbc.wrap(dx)
                        dx *= distance / r
                        h_atom.x = other_atom.x + dx
    for message in messages:
        print(message)


def write_mpmc_options(system, pbc):
    write_charges = 0
    write_force_field = 0
    force_field = 0

    while True:
        try:
            write_charges = input(
                "\nWould you like to read in charges?\n"
                "('yes', 'y', 1 or 'no', 'n', 0)\n\n> "
            )
            if write_charges == "yes" or write_charges == "y":
                write_charges = 1
            if write_charges == "no" or write_charges == "n":
                write_charges = 0
            write_charges = int(write_charges)
            if write_charges > 1 or write_charges < 0:
                raise ValueError
            break
        except ValueError:
            print("!!! Error reading input !!!")

    if write_charges == 1:
        while True:
            charges_filename = input("\ncharges file name > ")
            try:
                charges = []
                for line in open(charges_filename, "r").readlines():
                    charges.append(float(line))
                if len(charges) == len(system):
                    print("Applying charges ...")
                    for ind, atom in enumerate(system):
                        atom.charge = charges[ind]
                elif len(system) % len(charges) == 0:
                    print(
                        "Number of charges a multiple of the number of atoms, applying charges recursively ..."
                    )
                    for ind, atom in enumerate(system):
                        i = ind % len(charges)
                        atom.charge = charges[i]
                else:
                    raise ValueError
                break
            except TypeError:
                print("!!! Something went wrong reading charges file !!!")
            except ValueError:
                print(
                    "!!! Number of charges doesn't match (a multiple) of the number of atoms !!!\n"
                    "(or something else went wrong reading charges file)"
                )
            except FileNotFoundError:
                print("!!! File not found !!!")

    while True:
        try:
            write_force_field = input(
                "\nWould you like to automatically apply a forcefield to this MPMC .pbd file?\n"
                "('yes', 'y', 1 or 'no', 'n', 0)\n\n> "
            )
            if write_force_field == "yes" or write_force_field == "y":
                write_force_field = 1
            if write_force_field == "no" or write_force_field == "n":
                write_force_field = 0
            write_force_field = int(write_force_field)
            if write_force_field > 1 or write_force_field < 0:
                raise ValueError
            break
        except ValueError:
            print("!!! Error reading input !!!")

    if write_force_field == 1:
        while True:
            try:
                force_field = input(
                    "\nWhich force field?\n"
                    "valid answers are 'OPLSAA' (0) or 'PHAHST' (1)\n\n> "
                )
                if force_field == "OPLSAA":
                    force_field = 0
                if force_field == "PHAHST":
                    force_field = 1
                force_field = int(force_field)
                if force_field > 1 or force_field < 0:
                    raise ValueError
                apply_ff_to_system(system, get_forcefield(force_field))
                break
            except ValueError:
                print("!!! Error reading input !!!")

    filename = input("\noutput filename > ")

    if write_charges == 1:
        bool_charges = True
    else:
        bool_charges = False
    if write_force_field == 1:
        bool_ff = True
    else:
        bool_ff = False

    write_mpmc_pdb(
        system, pbc, filename, write_charges=bool_charges, write_params=bool_ff
    )
    return


def geom_analysis(system, pbc):
    while True:
        option = 0
        try:
            option = input(
                "\nWhat would you like to do?\n\n\
                1 = list bonds\n\
                2 = list close vdw contacts\n\
                3 = list angles\n\
                4 = list lone atoms\n\
                5 = delete lone atoms\n\
                6 = list coordinates\n\
                7 = edit hydrogen bond distances\n\
                8 = preview with VMD\n\
                9 = back to main menu\n\n> "
            )
            option = int(option)
        except ValueError:
            print("!!! Error converting input to int !!!")
        if option == 1:
            list_bonds(system, pbc)
        elif option == 2:
            list_close_contacts(system, pbc)
        elif option == 3:
            list_angles(system, pbc)
        elif option == 4:
            list_lone_atoms(system, pbc)
        elif option == 5:
            delete_lone_atoms(system, pbc)
        elif option == 6:
            list_coords(system)
        elif option == 7:
            edit_h_dist(system, pbc)
        elif option == 8:
            vmd_preview(system, pbc)
        elif option == 9:
            return


def main_loop():

    if len(sys.argv) != 2:
        sys.exit("Usage: python3 pdb_wizard.py <filename.[xyz|pdb]>")

    system = []

    input_filename = sys.argv[1]

    if input_filename.split(".")[-1] == "xyz":
        pbc = read_xyz(system, input_filename)
    elif (
        input_filename.split(".")[-1] == "pdb" or input_filename.split(".")[-1] == "ent"
    ):
        pbc = read_pdb(system, input_filename)
    elif input_filename.split(".")[-1] == "cif":
        sys.exit("Use Mercury to save a .cif as a .xyz or .pdb first")
    else:
        sys.exit("Unable to determine if input file is xyz or pdb (please rename)")

    if pbc is None:
        while True:
            try:
                a = input("Enter cell information\na>     ")
                a = float(a)
                b = input("b>     ")
                b = float(b)
                c = input("c>     ")
                c = float(c)
                alpha = input("alpha> ")
                alpha = float(alpha)
                beta = input("beta>  ")
                beta = float(beta)
                gamma = input("gamma> ")
                gamma = float(gamma)
                break
            except ValueError:
                print("!!! Error converting input to float !!!\n")
        pbc = PBC(a, b, c, alpha, beta, gamma)

    overlap_detector(system, pbc)

    while True:
        print_info(system, pbc, input_filename)
        option = 0
        try:
            option = input(
                "\nWhat would you like to do?\n\n\
                1 = geometry analysis\n\
                2 = extend along axis\n\
                3 = wrap atoms from (0, 0, 0) to (1, 1, 1)\n\
                4 = wrap atoms from (-1/2, -1/2, -1/2) to (1/2, 1/2, 1/2)\n\
                5 = update cell dimensions\n\
                6 = write .xyz\n\
                7 = write MPMC .pdb\n\
                8 = write standardized .pdb\n\
                9 = quit\n\n> "
            )
            option = int(option)
        except ValueError:
            print("!!! Error converting input to int !!!")
        if option == 1:
            geom_analysis(system, pbc)
        elif option == 2:
            extend_axis(system, pbc)
        elif option == 3:
            wrapall_forward(system, pbc)
        elif option == 4:
            wrapall(system, pbc)
        elif option == 5:
            while True:
                try:
                    a = input("Enter cell information\na>     ")
                    a = float(a)
                    b = input("b>     ")
                    b = float(b)
                    c = input("c>     ")
                    c = float(c)
                    alpha = input("alpha> ")
                    alpha = float(alpha)
                    beta = input("beta>  ")
                    beta = float(beta)
                    gamma = input("gamma> ")
                    gamma = float(gamma)
                    break
                except ValueError:
                    print("!!! Error converting input to float !!!\n")
            pbc.update(a, b, c, alpha, beta, gamma)
        elif option == 6:
            filename = input("\noutput filename > ")
            write_xyz(system, pbc, filename)
        elif option == 7:
            write_mpmc_options(system, pbc)
        elif option == 8:
            filename = input("\noutput filename > ")
            write_standard_pdb(system, pbc, filename)
        elif option == 9:
            break
        else:
            print("\nInvalid option!")


if __name__ == "__main__":
    main_loop()
