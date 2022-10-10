#!/home/rdkibler/.conda/envs/pyroml/bin/python3

import pyrosetta
from glob import glob
pyrosetta.init("--beta_nov16 true --mute all")
simple_threading_mover = pyrosetta.rosetta.protocols.simple_moves.SimpleThreadingMover()
detect_symmetry = pyrosetta.rosetta.protocols.symmetry.DetectSymmetry()

for fs in glob("../output/alignments/*fa"):
    parent_pose = pyrosetta.pose_from_file("../../5_carve_bb/"+fs.split('/')[-1][:-3]+'.pdb')
    detect_symmetry.apply(parent_pose)
    fin=open(fs,'r')
    fin.readline()
    line=fin.readline()[:-1]
    simple_threading_mover.set_sequence(line,1)
    mutant_pose = parent_pose.clone()
    simple_threading_mover.apply(mutant_pose)
    mutant_pose.dump_pdb(fs.split('/')[-1][:-3]+'_1.pdb')
    line=fin.readline()[:-1]
    simple_threading_mover.set_sequence(line,1)
    mutant_pose = parent_pose.clone()
    simple_threading_mover.apply(mutant_pose)
    mutant_pose.dump_pdb(fs.split('/')[-1][:-3]+'_2.pdb')
    line=fin.readline()[:-1]
    simple_threading_mover.set_sequence(line,1)
    mutant_pose = parent_pose.clone()
    simple_threading_mover.apply(mutant_pose)
    mutant_pose.dump_pdb(fs.split('/')[-1][:-3]+'_3.pdb')
    line=fin.readline()[:-1]
    simple_threading_mover.set_sequence(line,1)
    mutant_pose = parent_pose.clone()
    simple_threading_mover.apply(mutant_pose)
    mutant_pose.dump_pdb(fs.split('/')[-1][:-3]+'_4.pdb')
