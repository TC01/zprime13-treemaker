#!/usr/bin/env python

import os

import multiprocessing

from Zprime_Inclusive_Treemaker import Zprime_Inclusive_Treemaker

# Define a timeout.
mpTimeout = 43200

triggers = ["triggers/MuonID_Z_RunCD_Reco74X_Dec1.root", "triggers/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root"]

def run(source, output, isData):
	print "*** Running over " + source

	# Initialize pileup:
	PUF = TFile("pileup/PUnom.root")
	PUn = PUF.Get("pileup")
	PUU = TFile("pileup/PUup.root")
	PUu = PUU.Get("pileup")
	PUD = TFile("pileup/PUdn.root")
	PUd = PUD.Get("pileup")
	pileup = [PUn, PUu, PUd]


	treemaker = Zprime_Inclusive_Treemaker(output, source, isData, triggers)
	treemaker.Fill('B2GTTreeMaker/B2GTree', pileup)
	print "*** Cleaning up..."

def multiprocess(directory):
	# I *was* going to specify everything we needed but I got lazy...

	# Set up multiprocessing.
	pool = multiprocessing.Pool()
	results = []

	for path, dirs, files in os.walk(directory):
		if path == directory:
			for dir in dirs:
				outputName = dir
				source = os.path.join(directory, dir) + "/"

				# Decide which things are data.
				isData = False
				if "SingleElectron" in outputName or "SingleMu" in outputName:
					isData = True
				
				# Decide which things to exclude-- we only want the PromptReco data, I think?
				if isData and '17Jul2015' in outputName:
					continue

				#run(source, outputName, isData)
				result = pool.apply_async(run, (source, outputName, isData, ))
				results.append(result)

	# Close and join the pool.
	pool.close()
	pool.join()
	for result in results:
		result.get(timeout=mpTimeout)


if __name__ == '__main__':
	# Test code.
	#run("/uscms_data/d3/jkarancs/B2GTTreeNtuple/Aug13/SingleElectron_Run2015B-PromptReco/", "test", False)

	# Things to actually run.

	multiprocess("/eos/uscms/store/user/bjr/b2g/zprime-trees-round2/")
#	multiprocess("/eos/uscms/store/user/osherson/b2g/QCD_SAFE/")
#	multiprocess("/eos/uscms/store/user/bjr/b2g/backgrounds/")
#	multiprocess("/eos/uscms/store/user/anovak/b2g/")

#	run("/eos/uscms/store/user/bjr/b2g-particles/trees_data/SingleElectron_Run2015D-PromptReco_v4-decosa", "SingleElectron_Run2015D-PromptReco_v4-decosa", True)
