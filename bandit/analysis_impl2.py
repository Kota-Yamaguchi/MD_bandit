import numpy as np
#from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis as mda
from MDAnalysis import Universe
from pandas import Series
#from MDAnalysis.lib.util import get_weights, deprecate
#import matplotlib.pyplot as plt
#plt.switch_backend("tkagg")
import os
from MDAnalysis.analysis.pca_mod import PCA_MOD as PCA
#from argparse import ArgumentParser
import sys
from analysis import Analysis
from scipy.spatial import distance

class Analysis_impl2(Analysis):
	def __init__(self , trj, top, ref): 
		self.cudir = os.getcwd()
		self.ref = Universe(top, ref)
		self.trj = Universe(top, trj)


	def rootMeanSquareDeviation(self, rank_argv=None ,part = "name CA"):
		rmsd1=np.array([])
		print("calculating rmsd1")
		for ts in self.trj.trajectory:
			a =rmsd(self.trj.select_atoms("{}".format(part)).positions, self.ref.select_atoms("{}".format(part)).positions)
			rmsd1 = np.append(rmsd1, a)
		print(rmsd1)

		if rank_argv != None:
			return rmsd1[rank_argv]
		
		return rmsd1
	
	def contactMap(self, rank_argv=None):
		print("calculate ContactMap ")
		d = map(lambda i : np.array(i), self.trj.trajectory)
		d = np.array(list(d))
		
		dist_M = map(lambda i : distance.cdist(d[i], d[i] ,metric="euclidean"), range(d.shape[0]))
		dist_M = np.array(list(dist_M))
		
		d_ref = map(lambda i : np.array(i), self.ref.trajectory)
		d_ref = np.array(list(d_ref))
		ref_dist_M = map(lambda i : distance.cdist(d_ref[i], d_ref[i] ,metric="euclidean"), range(d_ref.shape[0]))
		ref_dist_M = np.array(list(ref_dist_M))
		
		hist_dif = []
		for i in range(dist_M.shape[0]):
			dist_dif = np.sum((dist_M[i] -ref_dist_M)**2)
			hist_dif.append(dist_dif)
		hist_dif = np.array(hist_dif)
		hist_dif = self.min_max(hist_dif)	
		if np.all(rank_argv) != None:
			return hist_dif[rank_argv]				
	
		return hist_dif

	def pca(self, n_components=30, select_="name CA",rank_argv=None):
      #if self.option.option != True: 
        #pc_space = np.array([])
		if os.path.exists(self.cudir+"/eigvec.npy") ==True:
			pca = PCA(self.trj ,select="{}".format(select_)).run()
			eig_vec = np.load(self.cudir+"/eigvec.npy")
              #eig_vec = pre_eigen_vec
			PC=pca.transform(self.trj.select_atoms("{}".format(select_)), n_components=n_components,eigen_vec = eig_vec)
              #pc_space = np.append(pc_space, PC, axis =0)
              #PC = np.dot(self.traj_row, pre_eigen_vec)
          
		else:
			print("You don't have eigen vector, make eigvec ")
			pca = PCA(self.trj ,select="{}".format(select_)).run()
			PC=pca.transform(self.trj.select_atoms("{}".format(select_)), n_components=n_components)
			eig_vec = pca.p_components
			np.save(self.cudir+"/eigvec.npy", eig_vec)
              #pc_space = np.append(pc_space, PC, axis =1)
              
		PC_ref=pca.transform(self.ref.select_atoms("{}".format(select_)), n_components=n_components, eigen_vec = eig_vec)
		PC_ref=PC_ref.reshape(len(PC_ref[0]))
          
		pcnorm = np.linalg.norm(PC, axis =1)
		PCnormalize = [PC[i]/pcnorm[i] for i in range(len(pcnorm))]
		PCnormalize = np.array(PCnormalize)

		pcrefnorm = np.linalg.norm(PC_ref) 
		PCnormalize_ref = PC_ref/pcrefnorm
		PC_dif = np.dot(PCnormalize,PCnormalize_ref)
          
		PC_dif = np.abs(PC_dif-1)
#		if rank_argv != None:
 #                       return PC_dif[rank_argv]
		return PC_dif #eig_vec 
	
	def min_max(self, x, axis=None):
		min = x.min(axis=axis, keepdims=True)
		max = x.max(axis=axis, keepdims=True)
		result = (x-min)/(max-min)
		return result


	def ranking(self, RC, rank = 10, reverse = True):
		hist = RC
		rank_argv = []
		top_rank = sorted(hist, reverse = reverse)[:int(rank)]
		for n in range(len(Series(top_rank))):
			for i in range(len(hist)):
				if hist[i] == Series(top_rank)[n]:
					rank_argv.append(int(i))
		score = np.array([rank_argv, top_rank])
		return score
