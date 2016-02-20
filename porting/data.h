//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb 21 02:11:53 2016 by ROOT version 6.04/14
// from TTree ntupleData/ntupleData
// found on file: MC13TeV_TTJets_m169v5.root
//////////////////////////////////////////////////////////

#ifndef ntupleData_h
#define ntupleData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ntupleData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Evt;
   Int_t           LumiBlock;
   Int_t           nPV;
   Float_t         nPUtrue;
   Float_t         PUWeights[3];
   Float_t         LepSelEffWeights[3];
   Float_t         TopPtWgt;
   Int_t           TrigWord;
   Int_t           nLepton;
   Float_t         Lepton_pt[2];   //[nLepton]
   Float_t         Lepton_eta[2];   //[nLepton]
   Float_t         Lepton_phi[2];   //[nLepton]
   Int_t           Lepton_id[2];   //[nLepton]
   Int_t           Lepton_gid[2];   //[nLepton]
   Int_t           Lepton_ch[2];   //[nLepton]
   Float_t         MET_pt;
   Float_t         MET_phi;
   Int_t           nGenWeight;
   Float_t         GenWeights[223];   //[nGenWeight]
   Int_t           nJet;
   Float_t         Jet_uncs[13][29];   //[nJet]
   Float_t         Jet_pt[13];   //[nJet]
   Float_t         Jet_genpt[13];   //[nJet]
   Float_t         Jet_eta[13];   //[nJet]
   Float_t         Jet_phi[13];   //[nJet]
   Float_t         Jet_mass[13];   //[nJet]
   Float_t         Jet_CombIVF[13];   //[nJet]
   Int_t           Jet_flavour[13];   //[nJet]

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Evt;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_nPUtrue;   //!
   TBranch        *b_PUWeights;   //!
   TBranch        *b_LepSelEffWeights;   //!
   TBranch        *b_TopPtWgt;   //!
   TBranch        *b_TrigWord;   //!
   TBranch        *b_nLepton;   //!
   TBranch        *b_Lepton_pt;   //!
   TBranch        *b_Lepton_eta;   //!
   TBranch        *b_Lepton_phi;   //!
   TBranch        *b_Lepton_id;   //!
   TBranch        *b_Lepton_gid;   //!
   TBranch        *b_Lepton_ch;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_nGenWeight;   //!
   TBranch        *b_GenWeights;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_uncs;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_genpt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_CombIVF;   //!
   TBranch        *b_Jet_flavour;   //!

   ntupleData(TTree *tree=0);
   virtual ~ntupleData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
#ifdef ntupleData_cxx
ntupleData::ntupleData(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MC13TeV_TTJets_m169v5.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MC13TeV_TTJets_m169v5.root");
      }
      f->GetObject("ntupleData",tree);

   }
   Init(tree);
}

ntupleData::~ntupleData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntupleData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntupleData::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ntupleData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Evt", &Evt, &b_Evt);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("nPUtrue", &nPUtrue, &b_nPUtrue);
   fChain->SetBranchAddress("PUWeights", PUWeights, &b_PUWeights);
   fChain->SetBranchAddress("LepSelEffWeights", LepSelEffWeights, &b_LepSelEffWeights);
   fChain->SetBranchAddress("TopPtWgt", &TopPtWgt, &b_TopPtWgt);
   fChain->SetBranchAddress("TrigWord", &TrigWord, &b_TrigWord);
   fChain->SetBranchAddress("nLepton", &nLepton, &b_nLepton);
   fChain->SetBranchAddress("Lepton_pt", Lepton_pt, &b_Lepton_pt);
   fChain->SetBranchAddress("Lepton_eta", Lepton_eta, &b_Lepton_eta);
   fChain->SetBranchAddress("Lepton_phi", Lepton_phi, &b_Lepton_phi);
   fChain->SetBranchAddress("Lepton_id", Lepton_id, &b_Lepton_id);
   fChain->SetBranchAddress("Lepton_gid", Lepton_gid, &b_Lepton_gid);
   fChain->SetBranchAddress("Lepton_ch", Lepton_ch, &b_Lepton_ch);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("nGenWeight", &nGenWeight, &b_nGenWeight);
   fChain->SetBranchAddress("GenWeights", GenWeights, &b_GenWeights);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_uncs", Jet_uncs, &b_Jet_uncs);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_genpt", Jet_genpt, &b_Jet_genpt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_CombIVF", Jet_CombIVF, &b_Jet_CombIVF);
   fChain->SetBranchAddress("Jet_flavour", Jet_flavour, &b_Jet_flavour);
   Notify();
}

Bool_t ntupleData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntupleData::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntupleData::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif
