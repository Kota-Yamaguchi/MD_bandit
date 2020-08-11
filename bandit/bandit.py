import numpy as np
import os
from analysis_impl2 import Analysis_impl2 
import sys
class Bandit:
	def __init__(self, trj, top, ref):
		self.cudir = os.getcwd()
		self.learningRate = 0.8
		self.RC_list =np.array( ["RMSD", "ContactMap", "PCA"])
		
		if os.path.exists(self.cudir+"/prob.npy") == False:
			self.RC_prob = np.zeros(self.RC_list.shape)
			self.RC_prob += 1/self.RC_list.shape[0]
		
		else:
			self.RC_prob = np.load(self.cudir+"/prob.npy")
			
		if os.path.exists(self.cudir+"/Score.npy") == True:
			self.score = np.load(self.cudir+"/Score.npy")
			self.choice = max(self.score)
		else:
			self.score = {"RMSD": 0,
				"ContactMap":0,
				"PCA":0}
			choice = np.random.choice(self.RC_list, p = self.RC_prob) 
	
		self.RC = Analysis_impl2(trj, top, ref)		

		self.calcRC = {"RMSD" : (lambda i  : self.RC.rootMeanSquareDeviation(i)) }  
		self.calcRC["ContactMap"] = lambda i  : self.RC.contactMap(i)
		self.calcRC["PCA"] = lambda i : self.RC.pca(i)
		


	def run(self,banditScoreName="ContactMap"):
		
		choice = np.random.choice(self.RC_list, p = self.RC_prob) 
		result = self.calcRC[choice](None)
		rank = self.RC.ranking(result)
		rank_argv = rank[0]
		rank_argv = [int(n) for n in rank_argv ]
		rankerValue = self.calcRC[banditScoreName](rank_argv)
		topRanker = np.max(rankerValue)
		if os.path.exists(self.cudir+"/preBanditScore.npy") == False:
			
			np.save(self.cudir+"/preBanditScore.npy", topRanker)
			delta = topRanker
			self.score[choice] = delta
			np.save(self.cudir+"/Score.npy", self.score)
		else:
			preBanditScore = np.load(self.cudir+"/preBanditScore.npy")
			delta = topRanker - preBanditScore
			preScore = np.load(self.cudir+"/Score.npy")
			self.score[choice] = delta
			for i in preScore:
				self.score[i]+=(self.learningRate * preScore[i])	
			np.save(self.cudir+"/Score.npy", self.score)
		
		self._updateProbabilitySoftmax(self.score, tau=10)
			
	def _updateProbabilitySoftmax(self, score, tau = 10):
		for i in self.RC_list:
			self.prob = np.exp(score[i])/(np.sum(np.exp(score[i]/10)))
		np.save(self.cudir+"/prob.npy", self.prob)
		return self.prob
		

	#def _score(self, ):
		
		#current_rankerのscore計算
	#	currentRC = Analysis_inple2(pre_traj, top, ref)
	#	scoreRC = {"RMSD" : (lambda : RC.rootMeanSquareDeviation()) } 
	#	scoreRC["ContactMap"] = lambda : RC.contactMap()
         #       scoreRC["PCA"] = lambda : RC.pca()
		
		
if __name__=="__main__":
	from argparse import ArgumentParser
	def get_option():
		argparser = ArgumentParser()
		argparser.add_argument("-s","--setting", help="gro file.")
		argparser.add_argument("-t","--traj", nargs="*", help="trajectory.")
		argparser.add_argument("-r","--reference", help="referense .")
		args = argparser.parse_args()
		return args
	option = get_option()
	bandit = Bandit(option.traj, option.setting ,option.reference)
	bandit.run()