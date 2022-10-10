import os
import numpy as np
from os import path
import subprocess
for Nreps in range(4,7):
   # lline2='../../4_reloop_v2/two_helix_scaffolds_2mers/'+lline2
    #print spl[-3]
    pdbname='reordered.pdb'
    fin3=open(pdbname,'r')
    lline3=fin3.readline()
    origLen=0
    while len(lline3)>0:
        if len(lline3)>4 and lline3[:4]=="ATOM" and lline3[21]=='A' and lline3[13:15]=='CA':
            origLen+=1
        lline3=fin3.readline()
    fin3.close()
    origLen=origLen/2.0
    startID=1

    #if repLen!=66:
    #    print pdbname
    fout=open('propagate.tcl','w')
    fout.write('mol new reordered.pdb\n')
    fout.write('set sel1 [atomselect top "residue %d to %d and chain A and backbone"]\n'%((startID-1),(origLen-1)))
    fout.write('set sel2 [atomselect top "residue %d to %d and chain A and backbone"]\n'%((startID+origLen-1),(2*origLen-1)))
    fout.write('set transformation_matrix [measure fit $sel1 $sel2]\n')
    fout.write('set move_sel [atomselect top "residue %d to %d and chain A"]\n'%(startID-1,origLen-1))
    for i in range(Nreps*2):
        fout.write("$move_sel move $transformation_matrix\n") 
        fout.write("$move_sel writepdb out%d.pdb\n"%(i))
    fout.write("quit\n")
    fout.close()
    proc = subprocess.Popen(['/home/nbethel/Documents/vmd/run_vmd_tmp', 'dispdev', 'text', '-e', 'propagate.tcl'],
                            stdin =  subprocess.PIPE,
                            stdout = subprocess.PIPE,
                            stderr = subprocess.PIPE
                        )

    (out, err) = proc.communicate()
    chains=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S']
    cwd=os.getcwd()
    mer=float(cwd.split("_")[-2][:-3])
    print(mer)
    fout=open('%drpt.pdb'%Nreps,'w' )
    for i in range(0,int(mer)):
     sj=0
     ej=Nreps
     rot=i
     resid=0
     for j in range(sj,ej):
      fin=open("out%d.pdb"%j,'r')
      lline=fin.readline()
      xxx=np.zeros((1,3))
      while len(lline)>0:
        if len(lline)>4 and lline[:4]=="ATOM":
            if  lline[13:15]=='N ':
                resid+=1
            xxx[0,0]=float(lline[30:38])
            xxx[0,1]=float(lline[38:46])
            xxx[0,2]=float(lline[46:54])
            Rz=np.array([[np.cos(np.pi*2.0*float(rot)/mer),np.sin(np.pi*2.0*float(rot)/mer),0],[-np.sin(np.pi*2.0*float(rot)/mer),np.cos(np.pi*2.0*float(rot)/mer),0],[0.0,0.0,1.0]])
            xxx=np.transpose(np.matmul(Rz,np.transpose(xxx)))
            fout.write(lline[:21]+chains[i]+'%4d'%resid+lline[26:30]+'%8.3f%8.3f%8.3f'%(xxx[0,0],xxx[0,1],xxx[0,2])+lline[54:])
        lline=fin.readline()
      fin.close()
    fout.close()
