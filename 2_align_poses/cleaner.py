import os
from glob import glob
for lline in glob('rescore*pdb'):
    spl=lline.split("_")
    mer=int(spl[-2][:-3])
    if os.path.isfile(spl[0]+"_"+spl[1]+"_"+spl[2]+"_"+spl[3]+"_"+spl[4]+"_"+spl[5]+"_"+spl[6]+"_"+spl[7]+'_%dmer_0001.pdb'%(mer+2)):
        os.remove(lline)
        continue
