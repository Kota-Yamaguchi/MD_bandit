import numpy as np
import os
from analysis_impl2 import Analysis_impl2 
import sys
class Bandit:
	def __init__(self, trj, top, ref):
		self.cudir = os.getcwd()
		self.forgettingRate = 0.8
		self.rankNumber = 100
		self.banditScoreName = "ContactMap"
		self.RC_list =np.array( ["RMSD", "ContactMap", "PCA"])
		
		if os.path.exists(self.cudir+"/prob.npy") == False:
			self.RC_prob = np.zeros(self.RC_list.shape)
			self.RC_prob += 1/self.RC_list.shape[0]
		
		else:
			self.RC_prob = np.load(self.cudir+"/prob.npy",allow_pickle=True)
			
		if os.path.exists(self.cudir+"/Score.npy") == True:
			self.score = np.load(self.cudir+"/Score.npy",allow_pickle=True)[0]
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
		


	def run(self,banditScoreName=self.banditScoreName):
		print("Probability  of one reaction coordinate {}".format(self.RC_prob))		
		choice = np.random.choice(self.RC_list, p = self.RC_prob) 
		result = self.calcRC[choice](None)
		rank = self.RC.ranking(result,rank=self.rankNumber)

		with open(self.cudir+"/rank.txt", "w") as f:
				f.write(str(rank.T))	



		rank_argv = rank[0]
		rank_argv = [int(n) for n in rank_argv ]
		rankerValue = self.calcRC[banditScoreName](rank_argv)
		topRanker = np.max(rankerValue)

		self._score(topRanker, choice)
		#if os.path.exists(self.cudir+"/preBanditScore.npy") == False:
		#	
		#	np.save(self.cudir+"/preBanditScore.npy", topRanker)
		#	delta = topRanker
		#	self.score[choice] = delta
		#	np.save(self.cudir+"/Score.npy", [self.score])
		#else:
		#	preBanditScore = np.load(self.cudir+"/preBanditScore.npy",allow_pickle=True)
		#	delta = topRanker - preBanditScore
		#	preScore = np.load(self.cudir+"/Score.npy",allow_pickle=True)[0]
		#	self.score[choice] = delta
		#	for i in preScore:
		#		self.score[i]+=(self.learningRate * preScore[i])	
		#	np.save(self.cudir+"/Score.npy", [self.score])
		
		self._updateProbabilitySoftmax(self.score, tau=10)
			
	def _updateProbabilitySoftmax(self, score, tau = 10):
		print("Update Probability ")
		n = 0
		sigma = 0
		for i in self.RC_list:
			sigma += np.exp(score[i]/tau)
		for i in self.RC_list:
			print("{} : {}".format(i ,np.exp(score[i]/tau)/sigma))
			self.RC_prob[n] = np.exp(score[i]/tau)/sigma
			n+=1
		np.save(self.cudir+"/prob.npy", self.RC_prob)
		return self.RC_prob
		

	def _score(self, Ranker, choiceRC):
		if os.path.exists(self.cudir+"/preBanditScore.npy") == False:
			
			np.save(self.cudir+"/preBanditScore.npy", Ranker)
			delta = Ranker
			self.score[choiceRC] = delta
			np.save(self.cudir+"/Score.npy", [self.score])
		else:
			preBanditScore = np.load(self.cudir+"/preBanditScore.npy",allow_pickle=True)
			delta = Ranker - preBanditScore
			preScore = np.load(self.cudir+"/Score.npy",allow_pickle=True)[0]
			self.score[choiceRC] = delta
			for i in preScore:
				self.score[i]+=(self.forgettingRate * preScore[i])	
			np.save(self.cudir+"/Score.npy", [self.score])

		
		
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
