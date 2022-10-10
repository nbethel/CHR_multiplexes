import os
import subprocess
import numpy as np
from os import path
import sys
from glob import glob
rosetta_path='/software/rosetta/latest/'

fout=open('repack_tasks2','w')
cnt=0

for lline2 in glob('h*pdb'):
    cnt+=1
    
    spl=lline2.split('_')
    pdbname=lline2
    fin=open(pdbname,'r')
    lline=fin.readline()
    repLen=0
    while len(lline)>0:
        if len(lline)>4 and lline[:4]=="ATOM" and lline[21]=='A' and lline[13:15]=='CA':
            repLen+=1
        lline=fin.readline()
    fin.close()

    rep1=repLen/4.0
    outDir=pdbname[:-4]
    offS=int(spl[0][1:])+int(spl[1][1:])-1
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
        fin3=open(pdbname,'r')
        fout3=open(outDir+'/rechained.pdb','w')
        lline3=fin3.readline()
        xxx=np.zeros((1,3))
        chainIDs=['A','B','C','D','E','F']
        while len(lline3)>0:
            if len(lline3)>3 and lline3[:4]=='ATOM' and lline3[12:15] not in ['2H ','3H ' ]:
                if (int(lline3[22:26])>(2*rep1+offS) or int(lline3[22:26])<=offS):
                    lline3=fin3.readline()
                    continue
                if lline3[21]=='A' and   int(lline3[22:26])<=rep1+offS:
                    outChain=chainIDs[0]
                elif lline3[21]=='A' and  int(lline3[22:26])>rep1+offS:
                    outChain=chainIDs[1]
                elif lline3[21]=='B' and  int(lline3[22:26])<=rep1+offS:
                    outChain=chainIDs[2]
                elif lline3[21]=='B' and  int(lline3[22:26])>rep1+offS:
                    outChain=chainIDs[3]
                elif lline3[21]=='C' and  int(lline3[22:26])<=rep1+offS:
                    outChain=chainIDs[4]
                elif lline3[21]=='C' and  int(lline3[22:26])>rep1+offS:
                    outChain=chainIDs[5]
#                else:
#                    print('break point not found')
                if lline3[12:15]=='1H ':
                    fout3.write(lline3[:12]+' H '+lline3[15:21]+outChain+lline3[22:])
                else:
                    fout3.write(lline3[:21]+outChain+lline3[22:])
            lline3=fin3.readline()
        fin3.close()
        fout3.close()

        fout2=open(outDir+'/helix.symm', "w+")
        ferr=open(outDir+'/helix.err', "w+")
        proc = subprocess.Popen([rosetta_path+'/src/apps/public/symmetry/make_symmdef_file.pl', '-p '+outDir+'/rechained.pdb','-m HELIX','-r 20.0','-a A','-b B','-i C'],
                                stdin = subprocess.PIPE,
                                stdout = fout2,
                                stderr = ferr
                            )

        (out, err) = proc.communicate()
        fout2.close()

    pdbname0=pdbname.split('/')[-1]

    fout.write('cd "'+outDir+'" ; ')
    fout.write(rosetta_path + '/bin/rosetta_scripts -s ')
    fout.write('rechained_INPUT.pdb'+' -parser:script_vars ')
    fout.write('s1c=%d e1c=%d s2c=%d e2c=%d  '
            %((offS+1),(offS+3),(rep1+offS-2),(rep1+offS)))
    fout.write('-nstruct 1 -parser:protocol ../design.xml @../flags -beta_nov16 -beta_nov16_cart  -holes:dalphaball ../../helper_scripts/DAlphaBall.gcc; ')
    fout.write('python ../reorder.py ; ')
    fout.write('python ../rearrange_cyc.py ; ')
    fout.write('rm out*pdb\n')
fout.close()
