import numpy as np
from scipy.optimize import leastsq
from scipy.optimize import brute,basinhopping
from os import path
from glob import glob


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
    return np.sqrt(np.square(xyzt2[:,0])+np.square(xyzt2[:,1])),zdiffs

def fittedfunc(p,xyz,flip,tot):
    Rx=np.array([[1,0,0],[0,np.cos(p[2]),np.sin(p[2])],[0,-np.sin(p[2]),np.cos(p[2])]])
    Ry=np.array([[np.cos(p[3]),0,-np.sin(p[3])],[0,1,0],[np.sin(p[3]),0,np.cos(p[3])]])
    Rz=np.array([[np.cos(np.pi*2.0*flip/tot),np.sin(np.pi*2.0*flip/tot),0],[-np.sin(np.pi*2.0*flip/tot),np.cos(np.pi*2.0*flip/tot),0],[0.0,0.0,1.0]])
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
    if frag and (p[0]>=1 or p[0]<=-1 ):
        outT+= 10*(np.abs(p[0]))
    if frag and (p[1]>=1 or p[1]<=-1):
        outT+= 10*(np.abs(p[1]))
    if frag and (p[2] <-3.14 or p[2]>3.14) :
        outT+= 10*(np.abs(p[2]))
    return outT+np.std(fitfunc(p, xyz,frag)[0])+np.std(fitfunc(p, xyz,frag)[1]) #error functiono

def radfunc(p,x,y,z):
    xyzt=np.zeros((len(x),3))
    xyzt[:,0]=x
    xyzt[:,1]=y
    xyzt[:,2]=z
    return np.mean(fitfunc(p, xyz,False)[0]) #error functiono

def findaxis(xyz,p,th,pdbname,aa,bb,cc,rl):

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    est_p  = brute(errfunc, (slice(p[0]-5.0, p[0]+50, 1.0)                    , slice(p[1]-1.0, p[1]+1, 1.0),
                             slice(p[2]-3.5, p[2]+3.5, 0.2)                , slice(p[3]-3.5,p[3]+ 3.6, 0.2)                      ), args=(x,y,z,False))
    print(pdbname[:-4],errfunc(est_p,x,y,z,False),radfunc(est_p,x,y,z))
    if errfunc(est_p,x,y,z,False)>0.005:
        return est_p
    splPDB=pdbname.split('/')
    chains=['A','B','C','D','E','F','G','H','I','J','K','L','M','O','P','Q','R','S']
    for mm in range(2,9):
      tot=float(mm)
      pdbOut=splPDB[-1][:-4]+'_%dmer.pdb'%mm
      fout=open(pdbOut,'w')
      fout.write('REMARK '+ pdbname[:-4]+ ' %6.4f %8.4f\n'%(errfunc(est_p,x,y,z,False),radfunc(est_p,x,y,z)))
      for ll in range(0,int(tot)):
        fin=open(pdbname,'r')
        lline=fin.readline()
        xxx=np.zeros((1,3))
        while len(lline)>0:
          if len(lline)>3 and lline[:4]=='ATOM':
              if int(lline[22:26])>rl*4:
                  lline=fin.readline()
                  continue
              xxx[0,0]=float(lline[30:38])-aa
              xxx[0,1]=float(lline[38:46])-bb
              xxx[0,2]=float(lline[46:54])-cc
              xxxN=fittedfunc(est_p,xxx,ll,tot)
              fout.write(lline[:21]+chains[ll]+lline[22:30]+'%8.3f%8.3f%8.3f'%(xxxN[0,0],xxxN[0,1],xxxN[0,2])+lline[54:])
          lline=fin.readline()
        fin.close()
      fout.close()
    xyz[:,0]=xyz[:,0]-np.mean(xyz[:,0])
    xyz[:,1]=xyz[:,1]-np.mean(xyz[:,1])
    xyz[:,2]=xyz[:,2]-np.mean(xyz[:,2])

for pdbname in glob("../1_backbone_gen/selects/*pdb"):

   # pdbname=pdbnames[int(sys.argv[1])]
    xyz=np.zeros((12,3))
    splPDB=pdbname.split('/')

    if not path.exists(pdbname):
        continue

    fin=open(pdbname,'r')
    spl=pdbname[:-1].split('_')
    rl=0
    lline=fin.readline()

    while len(lline)>0:
        if len(lline)>3 and lline[:4]=='ATOM':
            if lline[13:15]=='CA' and lline[21]=='A':
                rl+=1
        lline=fin.readline()
    fin.close()
    rl=rl/12

    fin=open(pdbname,'r')
    i=0
    lline=fin.readline()

    while len(lline)>0:
      if len(lline)>3 and lline[:4]=='ATOM':
          if lline[13:15]=='CA' and lline[21]=='A':
            xyz[i,0]+=float(lline[30:38])/rl
            xyz[i,1]+=float(lline[38:46])/rl
            xyz[i,2]+=float(lline[46:54])/rl
            if int(lline[22:26])%rl==0:
                i+=1
                if i==12:
                    break
      lline=fin.readline()
    fin.close()

    aa=np.mean(xyz[:,0])
    bb=np.mean(xyz[:,1])
    cc=np.mean(xyz[:,2])
    xyz[:,0]=xyz[:,0]-aa
    xyz[:,1]=xyz[:,1]-bb
    xyz[:,2]=xyz[:,2]-cc

    np.set_printoptions(suppress=True)    
    p = np.array([0,0,0,0])

    findaxis(xyz,p,0.00001,pdbname,aa,bb,cc,rl)

