# -*- coding: utf-8 -*-
"""Crystallographic classes and methods for pypospack.

This module represents simulation cells are crystallographic representions and
classes and methods to support this represention.

Attributes:
    iso_chemical_symbols(list): list of chemical symbols
    atom_info(dict): dictionary of chemical symbols, and a dictionary of their values
"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017,2018,2019"
__license__ = "MIT License"
__version__ = "1.0"

import subprocess, yaml, os
from copy import copy, deepcopy
import numpy as np
import numpy.linalg as linalg


from mexm.structure.utils import get_atomic_weight
from mexm.structure.lattice import Lattice
from mexm.structure.atom import Atom
from mexm.structure.simulationcell import SimulationCell
from mexm.structure.structuredb import StructureDatabase

# functions
from mexm.structure.simulationcell import convert_coordinates

def make_super_cell(structure, sc):
    """makes a supercell from a given cell

    Args:
        structure(pypospack.crystal.SimulationCell): the base structure from
            which the supercell will be made from.
        sc (:obj:`list` of :obj:`int`): the number of repeat units in the h1, h2, and h3
            directions
    """
    assert isinstance(structure, SimulationCell)
    assert isinstance(sc, list)
    assert all([isinstance(k,int) for k in sc])

    supercell = SimulationCell()
    supercell.structure_comment = "{}x{}x{}".format(sc[0],sc[1],sc[2])

    supercell.lattice = Lattice(
        a0 = structure.a0,
        H = structure.H
    )
    supercell.h1 = structure.h1 * sc[0]
    supercell.h2 = structure.h2 * sc[1]
    supercell.h3 = structure.h3 * sc[2]


    # add supercell atoms
    for i in range(sc[0]):
        for j in range(sc[1]):
            for k in range(sc[2]):
                for atom in structure.atomic_basis:
                    sc_atom  = deepcopy(atom)
                    sc_atom.position = [
                        (i+atom.position[0])/sc[0],
                        (j+atom.position[1])/sc[1],
                        (k+atom.position[2])/sc[2]\
                    ]
                    supercell.atomic_basis.append(deepcopy(sc_atom))

    # return a copy of the supercell
    return supercell

class RadialDistributionFunction(object):

    def pairCorrelationFunction_3D(x, y, z, S, rMax, dr):
        """Compute the three-dimensioanl pair correlation function

        Compute the three-dimensional pair correlation function for a set of
        spherical particles contained in a cube with side length S.  This simple
        function finds reference particles such that a sphere of radius rMax drawn
        around the particle will fit entirely within the cube, eliminating the need
        to compensate for edge effects.  If no such particles exist, an error is
        returned.  Try a smaller rMax...or write some code to handle edge effects! ;)

        Arguments:
            x(numpy.ndarray): an array of x positions of centers of particles
            y(numpy.ndarray): an array of y positions of centers of particles
            z(numpy.ndarray): an array of z positions of centers of particles
            S(list):length of each side of the cube in space
            rMax(float): outer diameter of largest spherical shell
            dr(float): increment for increasing radius of spherical shell

        Returns:
            tuple: first element is a numpy array containing the correlation
                function g(r) radii. the second element is a numpy array
                containing the radii of the spherical shells used to compute
                g(r).  the final element are the indices of the refence particles
        """
        from numpy import zeros, sqrt, where, pi, mean, arange, histogram

        # Find particles which are close enough to the cube center that a sphere of radius
        # rMax will not cross any face of the cube
        bools1 = x > rMax
        bools2 = x < (S - rMax)
        bools3 = y > rMax
        bools4 = y < (S - rMax)
        bools5 = z > rMax
        bools6 = z < (S - rMax)

        interior_indices, = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
        num_interior_particles = len(interior_indices)

        if num_interior_particles < 1:
            raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                    will lie entirely within a cube of side length S.  Decrease rMax\
                    or increase the size of the cube.")

        edges = arange(0., rMax + 1.1 * dr, dr)
        num_increments = len(edges) - 1
        g = zeros([num_interior_particles, num_increments])
        radii = zeros(num_increments)
        numberDensity = len(x) / S**3

        # Compute pairwise correlation for each interior particle
        for p in range(num_interior_particles):
            index = interior_indices[p]
            d = sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2)
            d[index] = 2 * rMax

            (result, bins) = histogram(d, bins=edges, normed=False)
            g[p,:] = result / numberDensity

        # Average g(r) for all interior particles and compute radii
        g_average = zeros(num_increments)
        for i in range(num_increments):
            radii[i] = (edges[i] + edges[i+1]) / 2.
            rOuter = edges[i + 1]
            rInner = edges[i]
            g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))

        return (g_average, radii, interior_indices)
        # Number of particles in shell/total number of particles/volume of shell/number density
        # shell volume = 4/3*pi(r_outer**3-r_inner**3)


def get_nearest_neighbor_information(a, lattice_type):
    lattice_types = {
            'fcc':FaceCenteredCubic }
    return lattice_types[lattice_type].get_nearest_neighbor_information(a)

def get_fcc_nearest_neighbor_distance(a0,NN):
    assert isinstance(a0,float)
    assert isinstance(NN,int) or isinstance(NN,float)

    nn_distances = [
        0,
        0.707 * a0,
        1.000 * a0,
        1.225 * a0,
        1.414 * a0,
        1.581 * a0
    ]

    import math
    return nn_distances[math.ceil(NN)] \
            + (NN%1)*(nn_distances[math.floor(NN)]-nn_distances[math.ceil(NN)])
