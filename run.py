#!/usr/bin/env python

import os

from Zprime_Inclusive_Treemaker import Zprime_Inclusive_Treemaker

def run(source, output, isData):
	print "*** Running over " + source
	treemaker = Zprime_Inclusive_Treemaker(output, source, isData)
	treemaker.Fill('B2GTTreeMaker/B2GTree')
	print "*** Cleaning up..."

def main(directory):
	# I *was* going to specify everything we needed but I got lazy...
	for path, dirs, files in os.walk(directory):
		if path == directory:
			for dir in dirs:
				outputName = dir
				source = os.path.join(directory, dir)

				# Decide which things are data.
				isData = False
				if "SingleElectron" or "SingleMu" in outputName:
					isData = True
				
				# Decide which things to exclude-- we only want the PromptReco data, I think?
				if isData and '17Jul2015' in outputName:
					continue

				run(source + "/", outputName, isData)


if __name__ == '__main__':
	# Test code.
	#run("/uscms_data/d3/jkarancs/B2GTTreeNtuple/Aug13/SingleElectron_Run2015B-PromptReco/", "test", False)
	
	main("/uscms_data/d3/jkarancs/B2GTTreeNtuple/Aug13/")
