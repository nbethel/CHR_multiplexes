import pyrosetta
import sys
from pyrosetta.rosetta.core.select import get_residues_from_subset

import argparse
import glob

def clean_file_name(input_file):
    tag= ""
    for i in range(0, len(input_file)):# removes .pdb from end of file name
        if input_file[i] == "." and len(input_file) - i <= 5:
            break
        tag += input_file[i]
        if input_file[i] == "/": #removes path before file name
           tag = ""
    return tag

def getMutRes(filename):
     pdbname = clean_file_name(filename)
     res_string = pdbname.split('_')[0]
     nums = []
     lets = []
     for n in range(len(res_string)):
         if res_string[n].isdigit():
             if res_string[n-1].isdigit():
                 nums[-1] = nums[-1] + res_string[n]
             else:
                 nums.append(res_string[n])
         else:
             lets.append(res_string[n])
     return [lets,nums]

def get_resnum_string(mutres):
    res_string = ""
    for i, res in enumerate(mutres[1]):
        if i < len(mutres[1]) - 1:
            res_string += res + ","
        else:
            res_string += res
    return res_string

def get_residue_list(pose, selector):
    residue_list = get_residues_from_subset(selector.apply(pose))
    return residue_list

def initialize_xmls():
    obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string("""
    <RESIDUE_SELECTORS>
            <Layer name="coreN2" select_core="true" select_boundary="false" select_surface="false" ball_radius="2.0" use_sidechain_neighbors="true" sc_neighbor_dist_exponent="0.7" core_cutoff="4.8"/>
            <Chain name="chainA" chains="A"/>
            <And name="cA" selectors="coreN2,chainA"/>
    </RESIDUE_SELECTORS>
    """)
    return obj

def select_ligand_neighborhood_old(pose, input_files):
    coord_res = getMutRes(input_files)
    coord_res = get_resnum_string(coord_res)
    coord_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(coord_res)
    ligand = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector("PCA")
    motif = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(coord_selector, ligand)
    iron = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand, 5.0)
    iron.set_atom_names_for_distance_measure(["FE1"])
    end = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand, 5.0)
    end.set_atom_names_for_distance_measure(["C8"])
    neighborhood = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(iron,end)
    return neighborhood

def select_ligand_neighborhood(pose, obj):
    
    neighborhood = obj.get_residue_selector("neighborhood")
    return neighborhood

def select_interface(pose, obj):
    interface = obj.get_residue_selector("sel_interface")
    return interface

def select_pore(pose, obj):
    interface = obj.get_residue_selector("pore")
    return interface

def select_not_surface_sheet(pose):
    surface = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    surface.set_layers(False, False, True)
    sheet = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector()
    sheet.set_selected_ss('E')
    surface_sheet = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(surface, sheet)
    not_surface_sheet = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(surface_sheet)
    return not_surface_sheet

def select_wide_boundary(pose):
    not_surface = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    not_surface.set_layers(False, True, False)
    not_surface.set_cutoffs(4.5,1.2)
    return not_surface

def select_wide_boundary_chainA(pose):
    boundary_selector = select_wide_boundary(pose)
    chainA = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('A')
    chainA_boundary = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(boundary_selector, chainA)
    return chainA_boundary

def select_all(pose):
    pose = pyrosetta.pose_from_file(pose)
    select_all = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    return select_all

def make_pos_file(residue_list, pdb, prefix = ""):
    pos_list = ''
    if '3mer' in pdb:
        divi=3
    elif '4mer' in pdb:
        divi=4
    elif '5mer' in pdb:
        divi=5
    llen=len(residue_list)
    spll=pdb.split('/')[-1].split('_')
    cnt=0
    for i in residue_list:
      #  if i>repL:
      #      break
        pos_list += str(i) + "A,"
        cnt+=1
    filename = 'pos/' + prefix + clean_file_name(pdb) + ".pos"
    pos_file = open(filename, "w")
    pos_file.write(pos_list)
    print("Wrote " + filename)

def pymol_list(residue_list, pdb_file):
    pymol_list = 'color red, res '
    for i, residue in enumerate(residue_list):
        if i < len(residue_list) - 1:
            pymol_list += str(residue) + '+'
        else:
            pymol_list += str(residue)
    pos_file = open("pymol_list.txt", "a")
    pos_file.write(pdb_file + "\n")
    pos_file.write(pymol_list + "\n\n")



pyrosetta.init("")

obj = initialize_xmls()


for pdb_file in glob.glob('../5_carve_bb/*pdb'):
    pose = pyrosetta.pose_from_file(pdb_file)
 
    residue_selector = obj.get_residue_selector("cA") ### MOD make a function that returns your residue selector
    residue_list = get_residue_list(pose, residue_selector)
    pymol_list(residue_list, pdb_file)

    make_pos_file(residue_list, pdb_file)
