import os
from glob import glob
import shutil
cntB=0
cntC=0
cntBoth=0
for lline2 in glob('h*mer'):
 if not os.path.exists(lline2+'/sasa_score.sc'):
    continue
 for ins in ['4rpt']:
  fin=open(lline2+'/sasa_score.sc','r')
  lline=fin.readline()
  lline=fin.readline()
  lline=fin.readline()
  flag1=False
  flag2=False
  while len(lline)>0:
    if (len(lline)>6 or lline[:4]=='SCOR') and ins in lline:
        spl=lline.split()
        splM=lline2.split('_')
        mer=float(splM[-1][:-3])
        if len(spl)>7:
            if mer==2 and float(spl[2])<7.0 and (float(spl[4])/mer)>10:
                #fout.write(spl[-1]+'\n')
                cc=(float(spl[4]))
                flag1=True
            elif mer>2 and float(spl[2])<7.0 and (float(spl[4])/mer)>20:
                cc=(float(spl[4]))
                flag1=True
            if float(spl[2])>=7.0:
                cntB+=1
            if mer>2 and (float(spl[4])/mer)<=20:
                cntC+=1
            if float(spl[2])<7.0 and mer>2 and (float(spl[4])/mer)<=20 and (float(spl[4])/mer)>=10 :
                print(lline2)
            if mer==2 and (float(spl[4])/mer)<=10:
                cntC+=1
            if (float(spl[2])>=7.0) or (mer>2 and (float(spl[4])/mer)<=20) or (mer==2 and (float(spl[4])/mer)<=10):
                cntBoth+=1
        else:
            print(lline[:-1],mer,lline2)
    if (len(lline)>6 or lline[:4]=='SCOR') and '6rpt' in lline:
        spl=lline.split()
        score=float(spl[1])
        if flag1 and score<100.0:
            shutil.copyfile(lline2+'/'+ins+'.pdb','../4_repack_side2/'+lline2+'_'+ins+'.pdb')
            break
    lline=fin.readline()
 fin.close()
print(cntB,cntC,cntBoth)
