#!/usr/bin/env python

import os

import multiprocessing

from Zprime_Inclusive_Treemaker import Zprime_Inclusive_Treemaker

# Define a timeout.
mpTimeout = 43200

def run(source, output, isData):
	print "*** Running over " + source
	treemaker = Zprime_Inclusive_Treemaker(output, source, isData)
	treemaker.Fill('B2GTTreeMaker/B2GTree')
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
				if "SingleElectron" or "SingleMu" in outputName:
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
	multiprocess("/uscms_data/d3/jkarancs/B2GTTreeNtuple/Aug13/")
	multiprocess("/eos/uscms/store/user/bjr/b2g/trees/Sep22")
