# Treemaker for Zprime --> tT analysis (but could in theory also be used for tt, thus "inclusive")
# Runs on TTrees from Janos' Treemaker (syntax based on those vars).
#### Marc Osherson, oshersonmarc@gmail.com Summer 2015

# include standard libraries:
import os
import sys
import math
from array import array
import ROOT
from ROOT import *

# The leppt.
lepPtCut = 25.

# Class Definition:
class Zprime_Inclusive_Treemaker:
	# Because the normalizations are now stored as event_weights directly in the TTrees, all we need is the folder containint the trees itself and a flag for looking at MC quantities.
	def __init__(self, OutFileName, TreeFolder, isData, lepFiles): 
		self.outname = OutFileName
		self.files = []
		self.TF = TreeFolder
		print "Files in " + self.TF + " :"
		for path, dir, files in os.walk(TreeFolder):
			for file in files:
				if file.endswith(".root"):
					self.files.append(os.path.join(path, file))
					print(file)
		print "IS THIS DATA? --- " + str(isData)
		self.lepFiles = lepFiles
		self.isData = isData
		self.__book__()

	def getUncertainty(self, rootfile, name, outputName):
		dataEff = rootfile.Get(name + "/efficienciesMC/abseta_pt_MC")
		mcEff = rootfile.Get(name + "/efficienciesDATA/abseta_pt_DATA")
		dataEff.Divide(mcEff)
		return dataEff.Clone(outputName)

	def GetLepUncertainties(self,lF):
		# lepID:
		muFile_ID = TFile(lF[0])
		self.muFile_ID = muFile_ID
		self.MC_2D_muID_eff = self.muFile_ID.Get("NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/efficienciesMC/abseta_pt_MC")
		self.DATA_2D_muID_eff = self.muFile_ID.Get("NUM_LooseID_DEN_genTracks_PAR_pt_spliteta_bin1/efficienciesDATA/abseta_pt_DATA")
		# now divide them and keep the result for later:
		self.DATA_2D_muID_eff.Divide(self.MC_2D_muID_eff)
		self.muID_eff = self.DATA_2D_muID_eff.Clone("ID_eff")

		# Triggers. This could be made nicer but I'm feeling lazy...
		self.triggerFile = ROOT.TFile(lF[1])
		self.runC_rereco_mu45_eff = self.getUncertainty(self.triggerFile, "runCreRECO_Mu45_eta2p1_PtEtaBins", "RunC_Mu45_Eff")
		self.runD_mu45_eff = self.getUncertainty(self.triggerFile, "runD_Mu45_eta2p1_PtEtaBins", "RunD_Mu45_Eff")

	# Create Branches for new trees (Error values will generally be -99.9... you should never see -99.9 except for the third-jet, which isn't required in all events.
	def __book__(self):
		print "booking..."
		## Writing new Tree:
		self.f = TFile( self.outname + ".root", "RECREATE" )
		self.tree = TTree("tree", "tree")
		# Branches of that tree:
		self.weight = array('f', [-99.9])
		self.addBranch('weight', self.weight)		
		self.isEl = array('i', [-9])
		self.addBranch('isEl', self.isEl)	
		self.isMu = array('i', [-9])
		self.addBranch('isMu', self.isMu)
		self.metPt = array('f', [-99.9])
		self.addBranch('metPt', self.metPt)
		self.metPhi = array('f', [-99.9])
		self.addBranch('metPhi', self.metPhi)
		self.lepPt = array('f', [-99.9])
		self.addBranch('lepPt', self.lepPt)
		self.lepPhi = array('f', [-99.9])
		self.addBranch('lepPhi', self.lepPhi)
		self.lepEta = array('f', [-99.9])
		self.addBranch('lepEta', self.lepEta)
		self.lepIsLoose = array('f', [-9])
		self.addBranch('lepIsLoose', self.lepIsLoose)
		self.lepIsTight = array('f', [-9])
		self.addBranch('lepIsTight', self.lepIsTight)
		# Lepton 2D and Triangle variables (based off of nearest "nice" jet)
		self.lep2Ddr = array('f', [-99.9])
		self.addBranch('lep2Ddr', self.lep2Ddr)
		self.lep2Drel = array('f', [-99.9])
		self.addBranch('lep2Drel', self.lep2Drel)
		
		# Leptonic isolation.
		self.lepIso = array('f', [-1.0])
		self.addBranch('lepIso', self.lepIso)
		self.lepMiniIso = array('f', [-1.0])
		self.addBranch('lepMiniIso', self.lepMiniIso)
		
		# Also lep2D vars based off of tag jet (see below).
		self.tagLep2Ddr = array('f', [-99.9])
		self.addBranch('tagLep2Ddr', self.tagLep2Ddr)
		self.tagLep2Drel = array('f', [-99.9])
		self.addBranch('tagLep2Drel', self.tagLep2Drel)

		# Also keep Janos's leptonic 2D cuts, but only one of them depending on the lepton.
		self.ak4jetV1dr = array('f', [-99.9])
		self.ak4jetV1rel = array('f', [-99.9])
		self.ak4jetV2dr = array('f', [-99.9])
		self.ak4jetV2rel = array('f', [-99.9])
		self.ak4jetV3dr = array('f', [-99.9])
		self.ak4jetV3rel = array('f', [-99.9])
		
		self.addBranch('ak4jetV1dr', self.ak4jetV1dr)
		self.addBranch('ak4jetV1rel', self.ak4jetV1rel)
		self.addBranch('ak4jetV2dr', self.ak4jetV2dr)
		self.addBranch('ak4jetV2rel', self.ak4jetV2rel)
		self.addBranch('ak4jetV3dr', self.ak4jetV3dr)
		self.addBranch('ak4jetV3rel', self.ak4jetV3rel)

		# We also want the actual 4 vector quantities for tag + lep.
		self.tagLepPt = array('f', [-99.9])
		self.addBranch('tagLepPt', self.tagLepPt)
		self.tagLepPhi = array('f', [-99.9])
		self.addBranch('tagLepPhi', self.tagLepPhi)
		self.tagLepEta = array('f', [-99.9])
		self.addBranch('tagLepEta', self.tagLepEta)
		self.tagLepMass = array('f', [-99.9])
		self.addBranch('tagLepMass', self.tagLepMass)

		# Monte Carlo truth info.
		self.MCtoppt = array('f', [-1.0])
		self.MCantitoppt = array('f', [-1.0])
		self.addBranch('MCtoppt', self.MCtoppt)
		self.addBranch('MCantitoppt', self.MCantitoppt)

		self.isLeptonic = array('f', [0.0])
		self.isTauEvent = array('f', [0.0])
		self.addBranch('isLeptonic', self.isLeptonic)
		self.addBranch('isTauEvent', self.isTauEvent)

		# Keep systematics info and stuff.
		self.XSec = array('f', [-1.0])
		self.NEventCorr = array('f', [-1.0])
		self.LumiWeight = array('f', [-1.0])
		self.GenWeight = array('f', [-1.0])
		self.muIDWeight = array('f', [-1.0])
		self.muIDWeightErr = array('f', [-1.0])
		self.runCreReco_mu45_Weight = array('f', [-1.0])
		self.runCreReco_mu45_WeightErr = array('f', [-1.0])
		self.runD_mu45_Weight = array('f', [-1.0])
		self.runD_mu45_WeightErr = array('f', [-1.0])

		self.addBranch('XSec', self.XSec)
		self.addBranch('NEventCorr', self.NEventCorr)
		self.addBranch('LumiWeight', self.LumiWeight)
		self.addBranch('GenWeight', self.GenWeight)
		self.addBranch('muIDWeight', self.muIDWeight)
		self.addBranch('muIDWeightErr', self.muIDWeightErr)
		self.addBranch('runCreReco_mu45_Weight', self.runCreReco_mu45_Weight)
		self.addBranch('runCreReco_mu45_WeightErr', self.runCreReco_mu45_WeightErr)
		self.addBranch('runD_mu45_Weight', self.runD_mu45_Weight)
		self.addBranch('runD_mu45_WeightErr', self.runD_mu45_WeightErr)

		# TAG JET
		self.tagJetPt = array('f', [-99.9])
		self.addBranch('tagJetPt', self.tagJetPt)
		self.tagJetPhi = array('f', [-99.9])
		self.addBranch('tagJetPhi', self.tagJetPhi)
		self.tagJetEta = array('f', [-99.9])
		self.addBranch('tagJetEta', self.tagJetEta)
		self.tagJetPrMass = array('f', [-99.9])
		self.addBranch('tagJetPrMass', self.tagJetPrMass)
		self.tagJetSDMass = array('f', [-99.9])
		self.addBranch('tagJetSDMass', self.tagJetSDMass)
		self.tagJetTau1 = array('f', [-99.9])
		self.addBranch('tagJetTau1', self.tagJetTau1)
		self.tagJetTau2 = array('f', [-99.9])
		self.addBranch('tagJetTau2', self.tagJetTau2)
		self.tagJetTau3 = array('f', [-99.9])
		self.addBranch('tagJetTau3', self.tagJetTau3)

		# W RECONSTRUCTION
		self.wPt = array('f', [-99.9])
		self.addBranch('wPt', self.wPt)
		self.wPhi = array('f', [-99.9])
		self.addBranch('wPhi', self.wPhi)
		self.wEta = array('f', [-99.9])
		self.addBranch('wEta', self.wEta)
		self.wAltPt = array('f', [-99.9])
		self.addBranch('wAltPt', self.wAltPt)
		self.wAltPhi = array('f', [-99.9])
		self.addBranch('wAltPhi', self.wAltPhi)
		self.wAltEta = array('f', [-99.9])
		self.addBranch('wAltEta', self.wAltEta)

		# TWO LIGHT JETS AND RELATED QUANTITIES:
		# - jets themselves:	
		self.numLightJets = array('i', [-9])
		self.addBranch('numLightJets', self.numLightJets)	
		self.leadJetPt = array('f', [-99.9])
		self.addBranch('leadJetPt', self.leadJetPt)
		self.leadJetPhi = array('f', [-99.9])
		self.addBranch('leadJetPhi', self.leadJetPhi)
		self.leadJetEta = array('f', [-99.9])
		self.addBranch('leadJetEta', self.leadJetEta)
		self.leadJetMass = array('f', [-99.9])
		self.addBranch('leadJetMass', self.leadJetMass)
		self.offJetPt = array('f', [-99.9])
		self.addBranch('offJetPt', self.offJetPt)
		self.offJetPhi = array('f', [-99.9])
		self.addBranch('offJetPhi', self.offJetPhi)
		self.offJetEta = array('f', [-99.9])
		self.addBranch('offJetEta', self.offJetEta)
		self.offJetMass = array('f', [-99.9])
		self.addBranch('offJetMass', self.offJetMass)	
		self.leadJetCSV = array('f', [-99.9])
		self.addBranch('leadJetCSV', self.leadJetCSV)	
		self.offJetCSV = array('f', [-99.9])
		self.addBranch('offJetCSV', self.offJetCSV)
		# - with W:
		self.lepTopPt = array('f', [-99.9])
		self.addBranch('lepTopPt', self.lepTopPt)
		self.lepTopPt2 = array('f', [-99.9])
		self.addBranch('lepTopPt2', self.lepTopPt2)
		self.lepTopMass = array('f', [-99.9])
		self.addBranch('lepTopMass', self.lepTopMass)
		self.lepTopMass2 = array('f', [-99.9])
		self.addBranch('lepTopMass2', self.lepTopMass2)
		self.lepTopCSV2 = array('f', [-99.9])
		self.addBranch('lepTopCSV2', self.lepTopCSV2)
		# - with TAGJET:
		self.hadTopPt = array('f', [-99.9])
		self.addBranch('hadTopPt', self.hadTopPt)
		self.hadTopPt2 = array('f', [-99.9])
		self.addBranch('hadTopPt2', self.hadTopPt2)
		self.hadTopMass = array('f', [-99.9])
		self.addBranch('hadTopMass', self.hadTopMass)
		self.hadTopMass2 = array('f', [-99.9])
		self.addBranch('hadTopMass2', self.hadTopMass2)
		# TOTAL EVENT QUANTITIES:
		self.eventMass = array('f', [-99.9])
		self.addBranch('eventMass', self.eventMass)
		self.eventMass2 = array('f', [-99.9])
		self.addBranch('eventMass2', self.eventMass2)

		# Triggers that we care about.
		self.HLT_Mu24_eta2p1 = array('f', [-1.0])
		self.HLT_Mu45_eta2p1 = array('f', [-1.0])
		self.HLT_Mu40_eta2p1_PFJet200_PFJet50 = array('f', [-1.0])
		self.addBranch('HLT_Mu24_eta2p1', self.HLT_Mu24_eta2p1)
		self.addBranch('HLT_Mu45_eta2p1', self.HLT_Mu45_eta2p1)
		self.addBranch('HLT_Mu40_eta2p1_PFJet200_PFJet50', self.HLT_Mu40_eta2p1_PFJet200_PFJet50)
		self.HLT_Ele33_CaloIdM_TrackIdM_PFJet30 = array('f', [-1.0])
		self.HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50 = array('f', [-1.0])
		self.addBranch('HLT_Ele33_CaloIdM_TrackIdM_PFJet30', self.HLT_Ele33_CaloIdM_TrackIdM_PFJet30)
		self.addBranch('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50', self.HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50 )

	def LoadBranch(self, Tree, var):
		Tree.SetBranchAddress(var[0], var[1])

	def Fill(self, TreeName): # Loop through events and fill them. Actual Fill step is done at the end, allowing us to make a few quality control cuts.

		if not self.isData:
			self.GetLepUncertainties(self.lepFiles)

		total = 0
		print "filling..."
		for i in self.files:
			print "Reading from " + i
			File = TFile(i)

			# Rather than halting print an error message...
			try:
				Tree = File.Get(TreeName)
				n = Tree.GetEntries()
			except:
				print "Error: unable to read a tree from " + i

			total += n-1
			for j in range(0, n): # Here is where we loop over all events.
				if j % 5000 == 0:
	      				percentDone = float(j) / float(n) * 100.0
	       				print 'Processing {0:10.0f}/{1:10.0f} : {2:5.2f} %'.format(j, n, percentDone )
				Tree.GetEntry(j) # populates all our self. arrays with the tree's values

				# Get the weight!
				# This is apparently evt_lumi_weight now...
				if not self.isData:
					self.weight[0] = Tree.evt_Lumi_Weight

			############# FAT JET PART ################
				tagJetIndex = -1
				nTagJets = 0
				# Find if we have (exactly) one heavy jet.
				for i in range(min(Tree.jetAK8_size,4)):		
					if Tree.jetAK8_prunedMass[i] > 50 and Tree.jetAK8_Pt[i] > 200 and math.fabs(Tree.jetAK8_Eta[0]) < 2.1:
						tagJetIndex = i
						nTagJets += 1
				if nTagJets != 1:
					continue
				self.tagJetPt[0] = Tree.jetAK8_Pt[tagJetIndex]
				self.tagJetEta[0] = Tree.jetAK8_Eta[tagJetIndex]
				self.tagJetPhi[0] = Tree.jetAK8_Phi[tagJetIndex]
				self.tagJetPrMass[0] = Tree.jetAK8_prunedMass[tagJetIndex]
				self.tagJetSDMass[0] = Tree.jetAK8_softDropMass[tagJetIndex]
				self.tagJetTau1[0] = Tree.jetAK8_tau1[tagJetIndex]
				self.tagJetTau2[0] = Tree.jetAK8_tau2[tagJetIndex]
				self.tagJetTau3[0] = Tree.jetAK8_tau3[tagJetIndex]
				TAGJET = ROOT.TLorentzVector()
				TAGJET.SetPtEtaPhiE(self.tagJetPt[0],self.tagJetEta[0],self.tagJetPhi[0],Tree.jetAK8_E[tagJetIndex])

			############# MET PART ################
				if Tree.met_Pt[0] > 25:
					self.metPt[0] = Tree.met_Pt[0]
					self.metPhi[0] = Tree.met_Phi[0]
				else:
					continue

			############# LEPTON PART ################

				lepE = -1
				self.isMu[0] = 0
				self.isEl[0] = 0	
				if len(Tree.el_Pt) == 0 and len(Tree.mu_Pt) > 0:
					if Tree.mu_Pt[0] > lepPtCut and math.fabs(Tree.mu_Eta[0]) < 2.1:
						self.isMu[0] = 1
						self.isEl[0] = 0
				elif len(Tree.el_Pt) > 0 and len(Tree.mu_Pt) == 0:
					if Tree.el_Pt[0] > lepPtCut and math.fabs(Tree.el_Eta[0]) < 2.1:
						self.isEl[0] = 1
						self.isMu[0] = 0
				elif len(Tree.el_Pt) > 0 and len(Tree.mu_Pt) > 0:
					if Tree.mu_Pt[0] > lepPtCut and math.fabs(Tree.mu_Eta[0]) < 2.1 and 2 * Tree.el_Pt[0] < Tree.mu_Pt[0]:
						self.isMu[0] = 1
						self.isEl[0] = 0
					elif Tree.el_Pt[0] > lepPtCut and math.fabs(Tree.el_Eta[0]) < 2.1 and 2 * Tree.mu_Pt[0] < Tree.el_Pt[0]:
						self.isEl[0] = 1
						self.isMu[0] = 0

				# Code reuse = bad! Only do this once for significantly cleaner code...
				if self.isMu[0] == 1:
					self.lepPt[0] = Tree.mu_Pt[0]
					self.lepPhi[0] = Tree.mu_Phi[0]
					self.lepEta[0] = Tree.mu_Eta[0]
					self.lepIsLoose[0] = Tree.mu_IsLooseMuon[0]
					self.lepIsTight[0] = Tree.mu_IsTightMuon[0]
					self.lepMiniIso[0] = Tree.mu_MiniIso[0]
					self.lepIso[0] = Tree.mu_Iso04[0]
					self.ak4jetV1dr = Tree.mu_AK4JetV1DR[0]
					self.ak4jetV2dr = Tree.mu_AK4JetV2DR[0]
					self.ak4jetV3dr = Tree.mu_AK4JetV3DR[0]
					self.ak4jetV1rel = Tree.mu_AK4JetV1PtRel[0]
					self.ak4jetV2rel = Tree.mu_AK4JetV2PtRel[0]
					self.ak4jetV3rel = Tree.mu_AK4JetV3PtRel[0]
					lepE = Tree.mu_E[0]
				elif self.isEl[0] == 1:
					self.lepPt[0] = Tree.el_Pt[0]
					self.lepPhi[0] = Tree.el_Phi[0]
					self.lepEta[0] = Tree.el_Eta[0]
					self.lepIsLoose[0] = Tree.el_isLoose[0]
					self.lepIsTight[0] = Tree.el_isTight[0]
					self.lepMiniIso[0] = Tree.el_MiniIso[0]
					self.lepIso[0] = Tree.el_Iso03[0]
					self.ak4jetV1dr = Tree.el_AK4JetV1DR[0]
					self.ak4jetV2dr = Tree.el_AK4JetV2DR[0]
					self.ak4jetV3dr = Tree.el_AK4JetV3DR[0]
					self.ak4jetV1rel = Tree.el_AK4JetV1PtRel[0]
					self.ak4jetV2rel = Tree.el_AK4JetV2PtRel[0]
					self.ak4jetV3rel = Tree.el_AK4JetV3PtRel[0]
					lepE = Tree.el_E[0]
				if lepE == -1:
					continue

				# Keep whichever of Janos's 2D cuts is correct for this lepton.
				if self.isEl[0] == 1:
					self.ak4jetV1dr = Tree.el_AK4JetV1DR[0]
					self.ak4jetV2dr = Tree.el_AK4JetV2DR[0]
					self.ak4jetV3dr = Tree.el_AK4JetV3DR[0]
					self.ak4jetV1rel = Tree.el_AK4JetV1PtRel[0]
					self.ak4jetV2rel = Tree.el_AK4JetV2PtRel[0]
					self.ak4jetV3rel = Tree.el_AK4JetV3PtRel[0]
				elif self.isMu[0] == 1:
					self.ak4jetV1dr = Tree.mu_AK4JetV1DR[0]
					self.ak4jetV2dr = Tree.mu_AK4JetV2DR[0]
					self.ak4jetV3dr = Tree.mu_AK4JetV3DR[0]
					self.ak4jetV1rel = Tree.mu_AK4JetV1PtRel[0]
					self.ak4jetV2rel = Tree.mu_AK4JetV2PtRel[0]
					self.ak4jetV3rel = Tree.mu_AK4JetV3PtRel[0]

				lep = ROOT.TLorentzVector()
				lep.SetPtEtaPhiE(self.lepPt[0],self.lepEta[0],self.lepPhi[0],lepE)
				met = ROOT.TLorentzVector()
				met.SetPtEtaPhiM(self.metPt[0],0.,self.metPhi[0],0.)
				Ws = make_W(met,lep)
				self.wPt[0] = Ws[0].Pt()
				self.wPhi[0] = Ws[0].Phi()
				self.wEta[0] = Ws[0].Eta()
				self.wAltPt[0] = Ws[1].Pt()
				self.wAltPhi[0] = Ws[1].Phi()
				self.wAltEta[0] = Ws[1].Eta()

			############# LIGHT JET PART ################
				lightJetList = []
				lightJetIndex = []
				d2dcutIndex = -1
				d2dcutDR = 9999.9

				# Also compute the DR between the tagged jet and the lepton.
				self.tagLep2Ddr[0] = lep.DeltaR(TAGJET)
				self.tagLep2Drel[0] = lep.Perp(TAGJET.Vect())

				tagLep = TAGJET + lep
				self.tagLepMass[0] = tagLep.M()
				self.tagLepPt[0] = tagLep.Pt()
				self.tagLepEta[0] = tagLep.Eta()
				self.tagLepPhi[0] = tagLep.Phi()

				# Find the light jet (if a good candidate exists).
				for i in range(min(Tree.jetAK4_size,6)):
					iJet = 	ROOT.TLorentzVector()
					iJet.SetPtEtaPhiE(Tree.jetAK4_Pt[i],Tree.jetAK4_Eta[i],Tree.jetAK4_Phi[i],Tree.jetAK4_E[i])
					if iJet.Pt() > 100. and iJet.DeltaR(TAGJET) > 0.6 and math.fabs(iJet.Eta()) < 2.1:
						lightJetList.append(iJet)
						lightJetIndex.append(i)
						if iJet.DeltaR(lep) < d2dcutDR:
							d2dcutDR = iJet.DeltaR(lep)
							#self.lep2Ddr[0] = iJet.DeltaR(lep)
							#self.lep2Drel[0] = iJet.Perp(lep.Vect())
							self.lep2Ddr[0] = lep.DeltaR(iJet)
							self.lep2Drel[0] = lep.Perp(iJet.Vect())

							# Store the pt

				if len(lightJetList) < 1:
					continue
				self.numLightJets[0] = len(lightJetList)
				# Always fill with lead jet
				self.leadJetPt[0] = lightJetList[0].Pt()
				self.leadJetEta[0] = lightJetList[0].Eta()
				self.leadJetPhi[0] = lightJetList[0].Phi()
				self.leadJetMass[0] = lightJetList[0].M()
				self.leadJetCSV[0] = Tree.jetAK4_CSV[lightJetIndex[0]]
				if len(lightJetList) > 1: # Fill offjet vars 
					self.offJetPt[0] = lightJetList[1].Pt()
					self.offJetEta[0] = lightJetList[1].Eta()
					self.offJetPhi[0] = lightJetList[1].Phi()
					self.offJetMass[0] = lightJetList[1].M()
					self.offJetCSV[0] = Tree.jetAK4_CSV[lightJetIndex[1]]

			############# EVENT QUANTITIES ################
				lepTop = Ws[0] + lightJetList[0]
				self.lepTopPt[0] = lepTop.Pt()
				self.lepTopMass[0] =  lepTop.M()
				self.eventMass[0] = (lightJetList[0] + TAGJET + Ws[0]).M()


				# Hadronic top assignments.
				hadTop = TAGJET + lightJetList[0]
				self.hadTopPt[0] = hadTop.Pt()
				self.hadTopMass[0] =  hadTop.M()
				if len(lightJetList) > 1:
					if lep.DeltaR(lightJetList[0]) < lep.DeltaR(lightJetList[1]):
						hadTop2 = TAGJET + lightJetList[1]
						lepTop2 = Ws[0] + lightJetList[0]
						self.lepTopCSV2[0] = Tree.jetAK4_CSV[lightJetIndex[0]]
					elif lep.DeltaR(lightJetList[1]) < lep.DeltaR(lightJetList[0]):
						hadTop2 = TAGJET + lightJetList[0]
						lepTop2 = Ws[0] + lightJetList[1]
						self.lepTopCSV2[0] = Tree.jetAK4_CSV[lightJetIndex[1]]

					self.lepTopPt2[0] = lepTop2.Pt()
					self.lepTopMass2[0] =  lepTop2.M()
					self.hadTopPt2[0] = hadTop2.Pt()
					self.hadTopMass2[0] =  hadTop2.M()
					self.eventMass2[0] = (lightJetList[0] + lightJetList[1] + TAGJET + Ws[0]).M()
				
				# MC Filling:
				######## Monte Carlo stuff.
				if not self.isData:
					# LEP SF
					self.muIDWeight[0] = get_lep_SF_and_Unc(self.muID_eff, min(lep.Pt(),119.9), math.fabs(lep.Eta()))[0]
					self.muIDWeightErr[0] = get_lep_SF_and_Unc(self.muID_eff, min(lep.Pt(),119.9), math.fabs(lep.Eta()))[1]
					self.runCreReco_mu45_Weight[0] = get_lep_SF_and_Unc(self.runC_rereco_mu45_eff, min(lep.Pt(),119.9), math.fabs(lep.Eta()))[0]
					self.runCreReco_mu45_WeightErr[0] = get_lep_SF_and_Unc(self.runC_rereco_mu45_eff, min(lep.Pt(),119.9), math.fabs(lep.Eta()))[1]
					self.runD_mu45_Weight[0] = get_lep_SF_and_Unc(self.runD_mu45_eff, min(lep.Pt(),119.9), math.fabs(lep.Eta()))[0]
					self.runD_mu45_WeightErr[0] = get_lep_SF_and_Unc(self.runD_mu45_eff, min(lep.Pt(),119.9), math.fabs(lep.Eta()))[1]
					# JET SF
					# TRUTH
					count = 0
					foundT = False
					foundTbar = False
					for id in Tree.gen_ID:
						if id == 6:
							foundT = True
							self.MCtoppt[0] = Tree.gen_Pt[count]
						if id == -6:
							foundTbar = True
							self.MCantitoppt[0] = Tree.gen_Pt[count]
						count += 1
					if not foundT:
						self.MCtoppt[0] = -1.0
					if not foundTbar:
						self.MCantitoppt[0] = -1.0

					# More Monte Carlo Stuff.
					self.isLeptonic[0] = 0.0
					self.isTauEvent[0] = 0.0
					leptonIDs = [11, 13, 15]
					for mcIndex in xrange(len(Tree.gen_ID)):
						id = Tree.gen_ID[mcIndex]
						# Exciting. This is different...
						try:
							momID = Tree.gen_MomID[mcIndex]
						except:
							momID = Tree.gen_Mom0ID[mcIndex]
						if id in leptonIDs and momID == abs(24):
							self.isLeptonic[0] += 1.0
						if id == 15:
							self.isTauEvent[0] = 1.0
							if momID == abs(24):
								self.isTauEvent = 2.0

					# Monte Carlo - systematics stuff.
					self.XSec[0] = Tree.evt_XSec
					self.NEventCorr[0] = Tree.evt_NEvent_Corr
					self.LumiWeight[0] = Tree.evt_Lumi_Weight
					self.GenWeight[0] = Tree.evt_Gen_Weight


					# Monte Carlo - triggers.
					self.HLT_Mu24_eta2p1[0] = Tree.HLT_Mu24_eta2p1
					self.HLT_Mu45_eta2p1[0] = Tree.HLT_Mu45_eta2p1
					self.HLT_Mu40_eta2p1_PFJet200_PFJet50[0] = Tree.HLT_Mu40_eta2p1_PFJet200_PFJet50
					self.HLT_Ele33_CaloIdM_TrackIdM_PFJet30[0] = Tree.HLT_Ele33_CaloIdM_TrackIdM_PFJet30
					self.HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50[0] = Tree.HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50


#	#	#	#	# Fill things:
				self.tree.Fill()

			File.Close()
		self.f.cd()
		print str(total) + " events processed"
		self.tree.Write()
		print "Tree written."
		self.f.Write()
		print "File written."
		self.f.Close()
		print "File closed?"
	def addBranch(self, name, var): # Just to make things easier to read (and reduce typos)
		self.tree.Branch(name, var, name+'/F')
	def __del__(self):
	        self.f.Close()
def get_lep_SF_and_Unc(scan, pt, eta):
	SF = scan.GetBinContent(scan.FindBin(eta,pt))
	SF_Err = scan.GetBinError(scan.FindBin(eta,pt))
	return [SF, SF_Err]
def make_W(met, lep): #both should be TLor vectors.
	newmet = ROOT.TLorentzVector()
	newmet_m = ROOT.TLorentzVector()
	newmet_p = ROOT.TLorentzVector()
	newmet.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
	newmet_m.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
	newmet_p.SetPtEtaPhiM(met.Pt(),0,met.Phi(),0)
	phivec = [math.cos(met.Phi()), math.sin(met.Phi())]
	P_lep = math.sqrt((lep.Px()*lep.Px())+(lep.Py()*lep.Py())+(lep.Pz()*lep.Pz()))
	P_phi = (lep.Px()*phivec[0])+(lep.Py()*phivec[1])
	b = (80.4*80.4) + (P_lep*P_lep) - (lep.E()*lep.E()) + (2*met.Pt()*P_phi)
	arg = (lep.E()*lep.E()) * ((4*met.Pt()*met.Pt()*((lep.Pz()*lep.Pz())-(lep.E()*lep.E())))+(b*b))
	if arg <= 0:
		Pz_met = lep.Pz()*b/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
		newmet.SetPz(Pz_met)
		newmet.SetE(math.sqrt(newmet.Px()*newmet.Px()+newmet.Py()*newmet.Py()+newmet.Pz()*newmet.Pz()))
		return [newmet+lep, newmet+lep]
	else:
		Pz_met_p = ((lep.Pz()*b)+math.sqrt(arg))/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
		Pz_met_m = ((lep.Pz()*b)-math.sqrt(arg))/(2*((lep.E()*lep.E()) -(lep.Pz()*lep.Pz())))
		newmet_p.SetPz(Pz_met_p)
		newmet_p.SetE(math.sqrt(newmet_p.Px()*newmet_p.Px()+newmet_p.Py()*newmet_p.Py()+newmet_p.Pz()*newmet_p.Pz()))
		newmet_m.SetPz(Pz_met_m)
		newmet_m.SetE(math.sqrt(newmet_m.Px()*newmet_m.Px()+newmet_m.Py()*newmet_m.Py()+newmet_m.Pz()*newmet_m.Pz()))
		return [newmet_p+lep, newmet_m+lep]

#### TEST THE ABOVE FUNCTIONS:
if __name__ == '__main__':
	F = "/eos/uscms/store/user/bjr/b2g/zprime-trees/ZprimeToTprimeT_TprimeToWB_MZp-2000Nar_MTp-1200Nar_RH_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_ZprimeToTprimeT_TprimeToWB_MZp-2000Nar_MTp-1200Nar/160105_194555/0000/"
	test = Zprime_Inclusive_Treemaker("test", F, False, ["triggers/MuonID_Z_RunCD_Reco74X_Dec1.root", "triggers/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root"])
	test.Fill("B2GTTreeMaker/B2GTree")
	print "Cleaning up..."

