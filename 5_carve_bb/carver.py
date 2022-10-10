import os
import numpy as np
from os import path
import shutil
import subprocess
from glob import glob
def makeBreak(bNum,repL,baseF,lline):
    spl=baseF.split("_")
    #spl=bp.split("_")
    bbounds=[[0,4],[1,3],[0,3],[1,2]]
    break1=0
    nH=2
    endL=3
    repLen=0
    for i in range(4):
        repLen+=int(spl[i][1:])
    for i in range(0,2*(bbounds[bNum][0]),2):
        break1+=int(spl[i][1:])+int(spl[i+1][1:])
        print(break1 )
    if (bbounds[bNum][1]>=nH):
        break2=5*repLen
        bcnt=bbounds[bNum][1]-(nH)
    else:
        break2=6*repLen
        bcnt=bbounds[bNum][1]
    for i in range(0,(2*(bcnt)),2):
        if (endL-1-i)<0:
            print('uhoh!!')
        break2-=(int(spl[endL-i][1:])+int(spl[endL-1-i][1:]))
    for i in range(2*(bcnt%nH),2*(bcnt%nH)+1):
        break2-=(int(spl[endL-i][1:]))


    fin2=open(lline+'/6rpt.pdb','r')
    fout=open(baseF+'_bk%d.pdb'%bNum,'w')
    lline2=fin2.readline()
    resid=0
    while len(lline2)>0:
        if len(lline2)>3 and lline2[:4]=='ATOM':
            if ((int(lline2[22:26])<=break1) or (int(lline2[22:26])>=break2)):
                lline2=fin2.readline()
                continue
            if  lline2[13:15]=='N ':
                resid+=1
             # xxxN=fittedfunc(xxx,i)
            fout.write(lline2[:22]+'%4d'%resid+lline2[26:])
        lline2=fin2.readline()
    fin2.close()


for lline in glob('../4_repack_side2/h*rpt'):
    baseF=lline.split('/')[-1]
    print(baseF)
    if path.exists(lline+'/5rpt.pdb') and not path.exists(baseF+'_bk1.pdb'):
      print('here')
      fout2=open('test.out', "w+")
      ferr=open('test.err', "w+")
      proc = subprocess.Popen(['../helper_scripts/blue.sh', lline+'/4rpt.pdb'],
                                stdin = subprocess.PIPE,
                                stdout = fout2,
                                stderr = ferr
                            )

      (out, err) = proc.communicate()
      fout2.close()
      with open('test.out') as f:
          for line in f:
             pass
          last_line = line
      Tline=last_line.split()[-1][:-1]
      outStr=''
      tmp='L'
      cnt=0
      for i in Tline:
          if i==tmp:
               cnt+=1
          else:
              outStr+=tmp.lower()+'%d_'%cnt
              cnt=1
              tmp=i
      spl=outStr.split('_')
      if spl[0][0]!='l' or spl[1][0]!='h' or spl[2][0]!='l' or spl[3][0]!='h' or spl[4][0]!='l':
         print("ERROR")
      bp='h%d_l%d_h%d_l%d'%((int(spl[0][1:])+int(spl[1][1:])),int(spl[2][1:]),int(spl[3][1:]),int(spl[4][1:]))
      splB=baseF.split('_')
      baseF=bp
      for b in range(4,len(splB)):
          baseF+='_'+splB[b]
      fin2=open(lline+'/5rpt.pdb','r')
      lline2=fin2.readline()
      repLen=0
      while len(lline2)>0:
        if len(lline2)>4 and lline2[:4]=="ATOM" and lline2[21]=='A' and lline2[13:15]=='CA':
            repLen+=1
        lline2=fin2.readline()
      fin2.close()
      repLen/=5.0
      print(baseF)
      for i in range(0,4):
        makeBreak(i,repLen,baseF,lline)
