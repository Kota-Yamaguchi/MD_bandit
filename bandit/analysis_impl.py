import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.lib.util import get_weights, deprecate
import matplotlib.pyplot as plt
#plt.switch_backend("tkagg")
import os
from MDAnalysis.analysis.pca import PCA
from argparse import ArgumentParser
import sys
from .analysis import Analysis

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument("-o","--option",action='store_true')
    argparser.add_argument("-s","--setting", help="gro file.")
    argparser.add_argument("-f","--force", nargs="*", help="trajectory.")
    argparser.add_argument("-r","--reference", help="referense .")
    argparser.add_argument("-fi","--fitting", help="fitting trajectory.")
    argparser.add_argument("--rmsd",help="what atoms you want to calculate ??")
    argparser.add_argument("--pca",help="what atoms you want to calculate ??")
    argparser.add_argument("--path", default=".")
    args = argparser.parse_args()
    return args

#cudir = os.getcwd()

class Analysis_impl2(Analysis):
    def __init__(self): 
      self.option = get_option()
      path = self.option.path
      self.cudir = path
      if self.option.option == True:
        self.top = self.option.setting
        obje = self.option.reference

        if len(self.option.force)==1:
            traj = self.option.force
        elif len(self.option.force) != 1:
            traj = self.option.force

        print(np.array(traj).shape)
        self.ref = Universe(self.top, obje)
        #traj = traj[0]
        print(traj)
        trj = Universe(self.top, traj)
        print("fitting")
        align.AlignTraj(trj, self.ref, select="{}".format(self.option.fitting),filename="alin.xtc").run()
        print("fited")
        self.trj = Universe(self.top, "alin.xtc")
        print("align")
        print(self.trj.trajectory)     
      else:
        self.top = None
        self.trj = None
        self.ref = None  
      self.pca_eigen_vec = None
       # a=map(lambda i : np.array(i), self.trj.select_atoms("name CA").trajectory)
       # trj = np.array(list(a))
       # b=map(lambda i : trj[i].reshape(trj.shape[1]*3),range(trj.shape[0]))
       # self.traj_row = np.array(list(b))
       # print(self.traj_row.shape)
    
    def root_mean_square_deviation(self):
        rmsd1=np.array([])
        print("calculating rmsd1")
        for ts in self.trj.trajectory:
            a =rmsd(self.trj.select_atoms("{}".format(self.option.rmsd)).positions, self.ref.select_atoms("{}".format(self.option.rmsd)).positions)
            rmsd1 = np.append(rmsd1, a)
        print(rmsd1)
        return rmsd1


    def PCA_(self, n_components=30, select_="name CA"):
      #if self.option.option != True: 
        #pc_space = np.array([])
        if self.option.option == True:
          if os.path.exists(self.cudir+"/eigvec.npy") ==True:
              pca = PCA(self.trj ,select="{}".format(select_)).run()
              eig_vec = np.load(self.cudir+"/eigvec.npy")
              #eig_vec = pre_eigen_vec
              PC=pca.transform(self.trj.select_atoms("{}".format(select_)), n_components=n_components,eigen_vec = eig_vec)
              print("axis",PC.shape)
              #pc_space = np.append(pc_space, PC, axis =0)
              #PC = np.dot(self.traj_row, pre_eigen_vec)
          
          else:
              print("You don't have eigen vector, make eigvec ")
              pca = PCA(self.trj ,select="{}".format(select_)).run()
              PC=pca.transform(self.trj.select_atoms("{}".format(select_)), n_components=n_components)
              eig_vec = pca.p_components
              np.save(self.cudir+"/eigvec.npy", eig_vec)
              #pc_space = np.append(pc_space, PC, axis =1)
              
          print("tar",PC.shape)       
          PC_ref=pca.transform(self.ref.select_atoms("{}".format(select_)), n_components=n_components, eigen_vec = eig_vec)
          print("ref",PC_ref.shape)
          PC_ref=PC_ref.reshape(len(PC_ref[0]))
          print("ref reshape",PC_ref.shape)
          
          pcnorm = np.linalg.norm(PC, axis =1)
          PCnormalize = [PC[i]/pcnorm[i] for i in range(len(pcnorm))]
          PCnormalize = np.array(PCnormalize)

          pcrefnorm = np.linalg.norm(PC_ref) 
          PCnormalize_ref = PC_ref/pcrefnorm
          PC_dif = np.dot(PCnormalize,PCnormalize_ref)
          
          PC_dif = np.abs(PC_dif-1)
          print("dif",PC_dif)
          return PC, PC_dif, eig_vec 



    def ranking(self, RC, rank = 10, reverse = True):
        hist = RC
        rank_argv = []
        top_rank = sorted(hist, reverse = reverse)[:int(rank)]
        for n in range(len(Series(top_rank))):
            for i in range(len(hist)):
                if hist[i] == Series(top_rank)[n]:
                    rank_argv.append(i)
        score = np.array([rank_argv, top_rank])
        return score

    def translater(self, top, xtc, frame, number):
        traj = md.load_xtc(xtc, top=top)
        traj = traj[frame]
        traj.save_gro("candi{0}_{1}.gro".format(number,frame))

    def main(self):
        if self.option.rmsd != None:
            rmsd=self.root_mean_square_deviation()
            return rmsd
        if self.option.pca != None:
            pca ,pca_dif, eigvec = self.PCA_(select_=self.option.pca)
            return pca_dif
if __name__ == "__main__":
 a = analysis()
 RC=a.main()
 np.save("RC.npy",RC)
 #pca,eigvec = a.PCA_()
 
 #plt.plot(rmsd,color="r")
 #plt.plot(pca.T[0],color="b")
 #plt.show()
