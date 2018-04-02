# [[file:~/Workspace/Programming/ovito-scripts/carbon-trjactory.note::27877ab9-40e3-47f4-a1a9-d2dafb69fc27][27877ab9-40e3-47f4-a1a9-d2dafb69fc27]]
#! /usr/bin/env python
# -*- coding: utf-8 -*-
#====================================================================#
#   DESCRIPTION:  ovito script for counting carbon atom number
#                 using the built-in slice modifier
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#         NOTES:  ovitos test.py
#        AUTHOR:  Wenping Guo (ybyygu)
#         EMAIL:  winpng@gmail.com
#       LICENCE:  GPL version 2 or upper
#       CREATED:  2016-04-17 14:39
#       UPDATED:  2018-04-01 18:19
#====================================================================#

TARGET_ID = 2                   # <= change the target particle type
IRON_TYPE_ID = 1                # <= change Fe particle type ID
SLICE_WIDTH = 1.2               # <= change the slice modifier width

import numpy as np
import ovito

from ovito import *
from ovito.io import *
from ovito.modifiers import *

def load_lammps_trjfile(trjfile):
    node = import_file(trjfile, multiple_frames=True)
    return node

def apply_slice_modifier(node, distance, slice_width, select=False):
    if select:
        node.modifiers.append(ClearSelectionModifier())

    node.modifiers.append(SliceModifier(
        normal=(0, 0, 1),       # on the z direction
        distance=distance,
        slice_width=slice_width,
        inverse=select,
        select=select
    ))

    # apply the slice modifier
    node.compute()

def count_carbon_atoms(node):
    # got the id number of C element
    items = node.output.particle_properties.particle_type.array
    carbon_id_list = [x for x in items if x == TARGET_ID]
    n_carbon = len(carbon_id_list)
    return n_carbon

def count_using_slice(node, distance, slice_width):
    """count carbon atoms using the slice modifier"""

    node.modifiers.clear()
    apply_slice_modifier(node,
                         distance=distance,
                         slice_width=slice_width,
                         select=False)
    n_carbon = count_carbon_atoms(node)
    n_sliced = node.output.number_of_particles
    return n_carbon, n_sliced


def show_positions(node, filename, particle):
    """save particle positions (z-component only) in numpy native format, so we can
    load them in external python environment for plotting using seaborn
    """
    modifier = SelectParticleTypeModifier(property="Particle Type")
    modifier.types = {particle}
    node.modifiers.append(modifier)
    # delete other particles
    node.modifiers.append(InvertSelectionModifier())
    node.modifiers.append(DeleteSelectedParticlesModifier())
    node.compute()
    positions = node.output.particle_properties.position.array
    # save the z component
    np.save(filename, positions.T[-1])

def find_first_layer_of_iron(node):
    """scan for iron atoms from the top"""

    zmax = 16                   # <= change
    zmin = -16                  # <= change
    nlayer = 225                # <= change the expected number of iron atoms per layer
    ntarget = nlayer*0.3        # <= change the expected number of iron atoms for the first layer

    for z in range(16, -16, -1):
        node.modifiers.clear()
        apply_slice_modifier(node,
                             distance=z,
                             slice_width=SLICE_WIDTH,
                             select=True)
        selection = node.output.particle_properties.selection.array
        ptypes  = node.output.particle_properties.particle_type.array
        cond = np.all([ptypes == IRON_TYPE_ID, selection == 1], axis=0)
        selected = selection[cond]
        n = np.count_nonzero(selected)
        if n > ntarget:
            return z

def main():
    trjfile = "tests/fe100_c_600K.lammpstrj"  # <= change
    node = load_lammps_trjfile(trjfile)
    total_frames = node.source.num_frames
    print("loaded number of atoms: {}".format(node.source.number_of_particles))
    print("loaded number of frames: {}".format(total_frames))

    # forward the time frame
    print("{:^10}{:^10}{:^10}{:^10}{:^10}".format("N_carbon", "N_sliced", "frame", "nlayer", "z"))

    jump_step = int(total_frames / 10.0)  # <= change the frame step size
    for frame in range(0, total_frames, jump_step):
        # Jump to the animation frame.
        ovito.dataset.anim.current_frame = frame
        # show_positions(node, particle=1, filename="Fe/positions-{:05}".format(frame))
        z = find_first_layer_of_iron(node)
        # z = 0
        n, ntotal = count_using_slice(node, distance=z, slice_width=SLICE_WIDTH)
        print("{:^10}{:^10}{:^10}{:^10}{:^10}".format(n, ntotal, frame, 1, z))
        n, ntotal = count_using_slice(node, distance=z+SLICE_WIDTH, slice_width=SLICE_WIDTH)
        print("{:^10}{:^10}{:^10}{:^10}{:^10}".format(n, ntotal, frame, 2, z))
        n, ntotal = count_using_slice(node, distance=z+2*SLICE_WIDTH, slice_width=SLICE_WIDTH)
        print("{:^10}{:^10}{:^10}{:^10}{:^10}".format(n, ntotal, frame, 3, z))

if __name__ == '__main__':
    main()
# 27877ab9-40e3-47f4-a1a9-d2dafb69fc27 ends here
