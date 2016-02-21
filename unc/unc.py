#!/usr/bin/env python

import optparse
import os,sys
import json
import pickle
import ROOT
from subprocess import Popen, PIPE

jetuncLabel = []
jetuncLabel.append("JES_down")
jetuncLabel.append("AbsoluteStat")
jetuncLabel.append("AbsoluteScale")
jetuncLabel.append("AbsoluteFlavMap")
jetuncLabel.append("AbsoluteMPFBias")
jetuncLabel.append("Fragmentation")
jetuncLabel.append("SinglePionECAL")
jetuncLabel.append("SinglePionHCAL")
jetuncLabel.append("TimeEta")
jetuncLabel.append("TimePt")
jetuncLabel.append("RelativeJEREC1")
jetuncLabel.append("RelativeJEREC2")
jetuncLabel.append("RelativeJERHF")
jetuncLabel.append("RelativePtBB")
jetuncLabel.append("RelativePtEC1")
jetuncLabel.append("RelativePtEC2")
jetuncLabel.append("RelativePtHF")
jetuncLabel.append("RelativeFSR")
jetuncLabel.append("RelativeStatEC")
jetuncLabel.append("RelativeStatHF")
jetuncLabel.append("PileUpDataMC")
jetuncLabel.append("PileUpPtRef")
jetuncLabel.append("PileUpPtBB")
jetuncLabel.append("PileUpPtEC1")
jetuncLabel.append("PileUpPtEC2")
jetuncLabel.append("PileUpPtHF")
jetuncLabel.append("FlavorPureGluon")
jetuncLabel.append("FlavorPureQuark")
jetuncLabel.append("FlavorPureCharm")
jetuncLabel.append("FlavorPureBottom")

def unc():
    inFileURL = "../porting/MC13TeV_TTJets_m169v5.root"
    fIn=ROOT.TFile.Open(inFileURL)
    tree=fIn.Get('data')
    totalEntries=tree.GetEntriesFast()
    tree.GetEntry(5)
    print len( tree.Jet_pt )
    print len( tree.Jet_uncs )     
    for x in range(tree.nJet) :
      print tree.Jet_pt[x], tree.Jet_eta[x]
      for y in range(29*x, 29*(x+1)) :
        print jetuncLabel[y%29], tree.Jet_uncs[y]
    

def main():
    unc()

if __name__ == "__main__":
    sys.exit(main())
