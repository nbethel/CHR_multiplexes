import os
import numpy as np
from os import path
import subprocess
from scipy.optimize import leastsq,minimize
from scipy.optimize import brute,basinhopping

def fitfunc(p,xyz,frag):
        Rx=np.array([[1,0,0],[0,np.cos(p[2]),np.sin(p[2])],[0,-np.sin(p[2]),np.cos(p[2])]])
        Ry=np.array([[np.cos(p[3]),0,-np.sin(p[3])],[0,1,0],[np.sin(p[3]),0,np.cos(p[3])]])
        xyzt=np.copy(xyz)
        xyzt[:,0]=xyzt[:,0]-p[0]
        xyzt[:,1]=xyzt[:,1]-p[1]
        xyzt2=np.transpose(np.matmul(Ry,np.matmul(Rx,np.transpose(xyzt))))
        zdiffs=np.zeros(len(xyzt2[:,2])-1)
        for i in range(0,len(zdiffs)):
            zdiffs[i]=xyzt2[i,2]-xyzt2[i+1,2]
        #return zdiffs
        return np.sqrt(np.square(xyzt2[:,0])+np.square(xyzt2[:,1])),zdiffs

def fittedfunc(p,xyz,flip,mer):
        Rx=np.array([[1,0,0],[0,np.cos(p[2]),np.sin(p[2])],[0,-np.sin(p[2]),np.cos(p[2])]])
        Ry=np.array([[np.cos(p[3]),0,-np.sin(p[3])],[0,1,0],[np.sin(p[3]),0,np.cos(p[3])]])
        Rz=np.array([[np.cos(3.14*float(flip)*2.0/mer),np.sin(3.14*float(flip)*2.0/mer),0],[-np.sin(3.14*float(flip)*2.0/mer),np.cos(3.14*float(flip)*2.0/mer),0],[0.0,0.0,1.0]])
        xyzt=np.copy(xyz)
        xyzt[:,0]=xyzt[:,0]-p[0]
        xyzt[:,1]=xyzt[:,1]-p[1]
        xyzt2=np.transpose(np.matmul(Rz,np.matmul(Ry,np.matmul(Rx,np.transpose(xyzt)))))
        zdiffs=np.zeros(len(xyzt2[:,2])-1)
        for i in range(0,len(zdiffs)):
            zdiffs[i]=xyzt2[i,2]-xyzt2[i+1,2]
        return xyzt2
def errfunc(p,x,y,z,frag):
        xyzt=np.zeros((len(x),3))
        xyzt[:,0]=x
        xyzt[:,1]=y
        xyzt[:,2]=z
        outT=0
        if frag:
            print(p,np.std(fitfunc(p, xyz,frag)[0])+np.std(fitfunc(p, xyz,frag)[1]))
        if frag and (p[0]>=1 or p[0]<=-1 ):
            outT+= 10*(np.abs(p[0]))
        if frag and (p[1]>=1 or p[1]<=-1):
            outT+= 10*(np.abs(p[1]))
        if frag and (p[2] <-3.14 or p[2]>3.14) :
            outT+= 10*(np.abs(p[2]))
        return np.std(fitfunc(p, xyz,frag)[0]) #error functiono
def radfunc(p,x,y,z):
        xyzt=np.zeros((len(x),3))
        xyzt[:,0]=x
        xyzt[:,1]=y
        xyzt[:,2]=z
        return np.mean(fitfunc(p, xyz,False)[0]) #error functiono

pdbname='rechained_INPUT_0001.pdb'
cwd=os.getcwd()
mer=int(cwd.split("_")[-1][:-3])
if (path.exists('score.sc')) and (path.exists('rechained_INPUT_0001.pdb')) and (path.exists('rechained_symm.pdb')):
    fin2=open('score.sc','r')
    lline2=fin2.readline()
    lline2=fin2.readline()
    lline2=fin2.readline()
    fin2.close()
    spl2=lline2.split()
    fin3=open('rechained_INPUT_0001.pdb','r')
    lline3=fin3.readline()
    repLen=0
    while len(lline3)>0:
        if len(lline3)>4 and lline3[:4]=="ATOM" and lline3[21]=='A' and lline3[13:15]=='CA':
            repLen+=1
        lline3=fin3.readline()
    fin3.close()
    fin2=open('rechained_INPUT_0001.pdb','r')
    lline2=fin2.readline()
    cntt=0
    spl=cwd.split("/")[-1].split('_')
    offS=int(spl[0][1:])+int(spl[1][1:])-1
    print(offS)
    while len(lline2)>0:
        if len(lline2)>3 and lline2[:4]=='ATOM':
          if lline2[13:15]=='N ' and int(lline2[22:26])==(offS+1):
              cntt+=1
        lline2=fin2.readline()
    fin2.close()
    fin2=open('rechained_INPUT_0001.pdb','r')
    lline2=fin2.readline()
    xyz=np.zeros((cntt,3))
    cntt=0
    while len(lline2)>0:
        if len(lline2)>3 and lline2[:4]=='ATOM':
          if lline2[13:15]=='N ' and int(lline2[22:26])==(offS+10):
              xyz[cntt,0]=float(lline2[30:38])
              xyz[cntt,1]=float(lline2[38:46])
              xyz[cntt,2]=float(lline2[46:54])
              cntt+=1
        lline2=fin2.readline()
    fin2.close()
    p = np.array([0,0,0,0])
    aa=np.mean(xyz[:,0])
    bb=np.mean(xyz[:,1])
    cc=np.mean(xyz[:,2])
    xyz[:,0]=xyz[:,0]#-np.mean(xyz[:,0])
    xyz[:,1]=xyz[:,1]#-np.mean(xyz[:,1])
    xyz[:,2]=xyz[:,2]#-np.mean(xyz[:,2])
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    est_p  = minimize(errfunc,p,args=(x,y,z,False)).x
    print(pdbname[:-4],errfunc(p,x,y,z,False),errfunc(est_p,x,y,z,False))
    fout=open('reordered.pdb','w')
    print(mer)
    for i in range(0,1):
      chainIDs=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O']
      fin2=open('rechained_symm.pdb','r')
      lline2=fin2.readline()
      h1=0
      h1c=np.array([0.0,0.0,0.0])
      htmin=100.0
      while len(lline2)>0:
        if len(lline2)>3 and lline2[:4]=='ATOM':
          if lline2[13:15]=='N ' and int(lline2[22:26])==(offS+1):
                htc=np.array([float(lline2[30:38]),float(lline2[38:46]),float(lline2[46:54])])
                if lline2[21]=='A':
                    h1c=np.copy(htc)
                elif lline2[21]==chainIDs[mer]:
                    h2c=np.copy(htc)
        lline2=fin2.readline()
      htmin=np.linalg.norm(h1c-h2c)
      fin2.close()
      fin2=open('rechained_INPUT_0001.pdb','r')
      lline2=fin2.readline()
      chainIDs=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O']
      h1=0
      h1Z=0.0
      h2=5
      h3=10
      h1c=np.array([0.0,0.0,0.0])
      h2c=np.array([0.0,0.0,0.0])
      h3c=np.array([0.0,0.0,0.0])
      while len(lline2)>0:
        if len(lline2)>3 and lline2[:4]=='ATOM':
          if lline2[13:15]=='N ' and int(lline2[22:26])==(offS+1):
     #           print(np.abs(np.linalg.norm(h1c-htc)-htmin))
                htc=np.array([float(lline2[30:38]),float(lline2[38:46]),float(lline2[46:54])])
                if h1==0:
                    h1c=np.copy(htc)
                    outChain=chainIDs[h1]
                    h1+=1
                elif (np.abs(np.linalg.norm(h2c-htc)-htmin)<(0.01) or h2c[0]==0.0):
                    h2c=np.copy(htc)
                    outChain=chainIDs[h2]
                    h2+=1
                elif (np.abs(np.linalg.norm(h3c-htc)-htmin)<(0.01) or h3c[0]==0.0):
                    h3c=np.copy(htc)
                    outChain=chainIDs[h3]
                    h3+=1
                elif np.abs(np.linalg.norm(h1c-htc)-htmin)<(0.01):
                    h1c=np.copy(htc)
                    outChain=chainIDs[h1]
                    h1+=1
          #      else:
          #          print lline
          if h1>=5 or h2>=10  or h3>=15:
              break
          xxx=np.zeros((1,3))
       # if h1>=6:
        #    print 'uhoh'
        lline2=fin2.readline()
      fin2.close()
      if h1>=3:
          pickChain=1
      elif h2>=8:
          pickChain=2
      elif h3>=13:
          pickChain=3
      else:
          print(h1,h2,h3,htmin,'here')
          #lline=fin.readline()
          continue
      fin2=open('rechained_INPUT_0001.pdb','r')
      lline2=fin2.readline()
      chainIDs=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R']
      h1=0
      h1Z=0.0
      h2=5
      h3=10
      h1c=np.array([0.0,0.0,0.0])
      h2c=np.array([0.0,0.0,0.0])
      h3c=np.array([0.0,0.0,0.0])
      resid=0
      while len(lline2)>0:
        if len(lline2)>3 and lline2[:4]=='ATOM' and not (lline2[12:16] in ['1H  ','2H  ','3H  ',' OXT']):
          
          if lline2[13:15]=='N ' and int(lline2[22:26])==(offS+1):
                htc=np.array([float(lline2[30:38]),float(lline2[38:46]),float(lline2[46:54])])
                if h1==0:
                    h1c=np.copy(htc)
                    outChain=chainIDs[h1]
                    h1+=1
                elif (np.abs(np.linalg.norm(h2c-htc)-htmin)<(0.01) or h2c[0]==0.0):
                    h2c=np.copy(htc)
                    outChain=chainIDs[h2]
                    h2+=1
                elif (np.abs(np.linalg.norm(h3c-htc)-htmin)<(0.01) or h3c[0]==0.0):
                    h3c=np.copy(htc)
                    outChain=chainIDs[h3]
                    h3+=1
                elif np.abs(np.linalg.norm(h1c-htc)-htmin)<(0.01):
                    h1c=np.copy(htc)
                    outChain=chainIDs[h1]
                    h1+=1
                else:
                    outChain='Z'
          if h1>=6 or h2>=11 or h3>=16:
              break
          xxx=np.zeros((1,3))
          
          if (outChain in ['A','B','C','D','E'] and pickChain==1) or (outChain in ['F','G','H','I','J'] and pickChain==2) or (outChain in ['K','L','M','N','O'] and pickChain==3):
              if ((int(lline2[22:26])<(repLen+1) and h1==1 and pickChain==1) or ((int(lline2[22:26])<(repLen+1) and h2==6 and pickChain==2)) or ((int(lline2[22:26])<(repLen+1) and h3==11 and pickChain==3))):
                  lline2=fin2.readline()
                  continue
              if  ((int(lline2[22:26])>(repLen) and h1>=3 and pickChain==1) or ((int(lline2[22:26])>(repLen) and h2>=8 and pickChain==2)) or ((int(lline2[22:26])>(repLen) and h3>=13 and pickChain==3)) ):
                  lline2=fin2.readline()
                  continue
              if  ((h1>3 and pickChain==1) or ((h2>8 and pickChain==2)) or ((h3>13 and pickChain==3)) ):
                  lline2=fin2.readline()
                  continue
            #  xxx[0,0]=float(lline2[30:38])-np.mean(xs)
            #  xxx[0,1]=float(lline2[38:46])-np.mean(ys)
            #  xxx[0,2]=float(lline2[46:54])
              xxx[0,0]=float(lline2[30:38])#-aa
              xxx[0,1]=float(lline2[38:46])#-bb
              xxx[0,2]=float(lline2[46:54])#-cc
              if pickChain==1:
                  htt=i
              elif pickChain==2:
                  htt=i
              else:
                  htt=i
              if i==2:
                  ddd=-1
              else:
                  ddd=i
              if  lline2[13:15]=='N ':
                resid+=1
              xxxN=fittedfunc(est_p,xxx,i,float(mer))
             # xxxN=fittedfunc(xxx,i)
              fout.write(lline2[:21]+chainIDs[htt]+'%4d'%resid+lline2[26:30]+'%8.3f%8.3f%8.3f'%(xxxN[0,0],xxxN[0,1],xxxN[0,2])+lline2[54:])
        #if h1>=6:
        #    break
        lline2=fin2.readline()
      fin2.close()
#    print spl[-9]+'_'+spl[-8]+'_'+spl[-7]+'_'+spl[-6]+'_'+spl[-5]+'_'+spl[-3]+'_'+spl[-2],h1,h2
    fout.close()
