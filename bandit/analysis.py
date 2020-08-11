import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis as mda
from MDAnalysis import Universe
#from MDAnalysis.lib.util import get_weights, deprecate
import matplotlib.pyplot as plt
#plt.switch_backend("tkagg")
import os
from MDAnalysis.analysis.pca import PCA
#from argparse import ArgumentParser
#import sys

class Analysis():
	def rootMeanSquareDeviation(self):
		return null


	def pca(self, n_components=30, select_="name CA"):
		return null, null, null 

	def contactMap(self):
		return null

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

