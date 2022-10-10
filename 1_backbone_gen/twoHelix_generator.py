from shutil import copyfile
import os
foutT=open("two_helix_tasks",'w')
for i in range(12,23):
    for j in range(2,5):
        for k in range(12,23):
            for l in range(2,5):
                        os.mkdir("h%d_l%d_h%d_l%d"%(i,j,k,l)) 
                        os.chdir("h%d_l%d_h%d_l%d"%(i,j,k,l))
                        fout=open('design.csts','w')
                        fout.write('AtomPair CA %d CA %d SCALARWEIGHTEDFUNC 50 BOUNDED 9.0 18.0 1 tag\n'%((i+j+k+1),(i+j+k+1+i+j+k+l)))
                        fout.write('AtomPair CA %d CA %d SCALARWEIGHTEDFUNC 50 BOUNDED 9.0 18.0 1 tag\n'%((i+1),(i+1+i+j+k+l)))
                        fout.close()
                        fout=open("design.blueprint",'w')
                        fout.write("1 A HA\n")
                        for ii in range(2,i+1):
                            fout.write("0 x HA\n")
                        for ii in range(1,j+1):
                            fout.write("0 x LD\n")
                        for ii in range(1,k+1):
                            fout.write("0 x HA\n")
                        for ii in range(1,l):
                            fout.write("0 x LD\n")
                        fout.write("0 x LD")
                        fout.close()
                        os.chdir("..")
                        for ii in range(5):
                            foutT.write("cd h%d_l%d_h%d_l%d\n"%(i,j,k,l))
foutT.close()
