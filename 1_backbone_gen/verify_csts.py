from glob import glob
import os
import numpy as np
import shutil
foutW=open('../2_align_poses/cyl_tasks','w')
for dire in glob("h*"):
    for pdb in glob(dire+"/*dummy*pdb"):
        fin=open(dire+'/design.csts','r')
        lline=fin.readline()
        anyunsat=False
        while len(lline)>0:
            if len(lline)>8 and lline[:4]=='Ambi':
              satis=False
              while lline[:3]!='END':
                if len(lline)>8 and lline[:4]=='Atom':
                    ind1=lline.split()[2]
                    ind2=lline.split()[4]
                    dist=float(lline.split()[-3])
                    fin2=open(pdb,'r')
                    lline2=fin2.readline()
                    while len(lline2)>0:
                        if len(lline2)>4 and lline2[:4]=='ATOM' and lline2[13:15]=='CA':
                            if lline2.split()[5]==ind1:
                                vec1=np.array([float(lline2[30:38]),float(lline2[38:46]),float(lline2[46:54])])
                            elif lline2.split()[5]==ind2:
                                vec2=np.array([float(lline2[30:38]),float(lline2[38:46]),float(lline2[46:54])])
                                distT=np.linalg.norm(vec2-vec1)
                                if distT<=dist:
                                    satis=True
                                break
                        lline2=fin2.readline()
                    fin2.close()
                lline=fin.readline()
              if satis==False:
                  anyunsat=True
                  break
            lline=fin.readline()            
        fin.close()
        if anyunsat==False:
            shutil.copyfile(pdb,'selects/'+pdb.split('/')[-2]+'_'+pdb.split('/')[-1])
        else:
            os.remove(pdb)
