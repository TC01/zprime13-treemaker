#!/usr/bin/env python

from Zprime_Inclusive_Treemaker import Zprime_Inclusive_Treemaker

def run(source, output, isData):
	print "*** Running over " + source
	treemaker = Zprime_Inclusive_Treemaker(output, source, isData)
	treemaker.Fill('B2GTTreeMaker/B2GTree')
	print "*** Cleaning up..."

if __name__ == '__main__':
	run("/uscms_data/d3/jkarancs/B2GTTreeNtuple/Aug13/SingleElectron_Run2015B-PromptReco/", "test", False)
