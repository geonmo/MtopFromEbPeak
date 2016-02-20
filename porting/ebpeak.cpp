#include<iostream>
#include<TH1F.h>
#include<TFile.h>
#include <unordered_map>
#include"OptionParser.h"
#include<fstream>
#include"json/json.h"
#include"data.h"
#include<TLorentzVector.h>


//namespace po = boost::program_options;
using namespace std;
using optparse::OptionParser;

class sampleInfo{
  
private : 
  double xsec_;
  bool isMC_;
  string origin_; 
  string group_;
  int color_;
public :
  double getXsec() { return xsec_; }
  bool   isMC() { return isMC_; }
  string getOrigin() { return origin_; }
  string getGroup() { return group_; }
  int    getColor() { return color_; }

  sampleInfo() {}
  sampleInfo( const Json::Value database ) {
    xsec_ = database[0].asFloat();
    if ( database[1].asInt() ==0 ) isMC_ = true;
    else isMC_ = false;
    origin_ = database[2].asString();
    group_  = database[3].asString();
    color_ = database[4].asInt();
  }
};
           

  
class runBJetEnergyPeak{

private :
std::string inFileURL, outFileURL;
double xsec ;

std::unordered_map<std::string, TH1F*> hist;

public :
runBJetEnergyPeak( std::string inFileURL, std::string outFileURL, double xsec = -1.0) {
	this->inFileURL=inFileURL;
	this->outFileURL=outFileURL;
	this->xsec = xsec;
  // Booking the histogram.
  bookHist1D("nvtx",";Vertex multiplicity; Events",30,0,30 );
  bookHist1D("nbtags",";b-tag multiplicity; Events",5,0,5);
  bookHist1D("bjeten",";Energy [GeV]; Jets",30,0,300);
  bookHist1D("bjetenls",";log(E);  1/E dN_{b jets}/dlog(E)",20,3.,7.);

  ana();
}

void bookHist1D( std::string name, std::string title, int nbin, int xmin, int xmax) {
	hist[name] = new TH1F(name.c_str(), title.c_str(), nbin, xmin, xmax);
	hist[name]->Sumw2();
	hist[name]->SetDirectory(0);
}

void ana(){
    cout<<"...analysing "<<inFileURL<<endl;

    //open file and loop over events tree
    TFile* fIn = TFile::Open( inFileURL.c_str() );
    TTree* tree = (TTree*)fIn->Get("data");
    int totalEntries=tree->GetEntriesFast();

    ntupleData* ntuple = new ntupleData(tree) ;  
    for ( int i=0 ; i < totalEntries ; ++i){

        ntuple->GetEntry(i);
        //if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
        //#require at least two jets
        int nJets=0, nBtags=0;
        vector<TLorentzVector> taggedJetsP4;

        for ( int ij=0 ; ij < ntuple->nJet ; ++ij) {
            //#get the kinematics and select the jet
            auto jp4=TLorentzVector();
            jp4.SetPtEtaPhiM(ntuple->Jet_pt[ij],ntuple->Jet_eta[ij],ntuple->Jet_phi[ij],ntuple->Jet_mass[ij]);
            if (jp4.Pt()<30 || TMath::Abs(jp4.Eta())>2.4)  continue;
            //#count selected jet
            nJets +=1;

            //#save P4 for b-tagged jet
            if ( ntuple->Jet_CombIVF[ij]>0.605) {
                nBtags+=1;
                taggedJetsP4.push_back(jp4);
            }
        }
        if ( nJets<2 ) continue; 
        //#generator level weight only for MC
        float evWgt=1.0;
        if ( (int)xsec != 1  )          evWgt  = xsec*ntuple->LepSelEffWeights[0]*ntuple->PUWeights[0];
        if (ntuple->nGenWeight>0)  evWgt *= ntuple->GenWeights[0];

        //#ready to fill the histograms
        hist["nvtx"]->Fill(ntuple->nPV,evWgt);
        hist["nbtags"]->Fill(nBtags,evWgt);

        //#use up to two leading b-tagged jets
        for ( int ij=0 ; ij < taggedJetsP4.size() ; ++ij) {
            if (ij>1) break;
            hist["bjeten"]->Fill(taggedJetsP4[ij].E(),evWgt);
            hist["bjetenls"]->Fill(TMath::Log(taggedJetsP4[ij].E()),evWgt/taggedJetsP4[ij].E());
        }
    }
    //#all done with this file
    fIn->Close();

    //#save histograms to file
    TFile* fOut=TFile::Open(outFileURL.c_str(),"RECREATE");
    for (auto ptr = hist.begin() ; ptr != hist.end() ; ++ptr) hist[ptr->first]->Write();
    fOut->Close();
}

};

int main( int argc, char* argv[ ]) {

  OptionParser parser = OptionParser() .description("As like runBJetEnergyPeak.py");
  parser.add_option("-j", "--json").dest("json").help("json with list of files").set_default(""); //        type='string')
  parser.add_option("-i", "--inDir").dest("inDir").help("input directory with files").set_default("");//        type='string')
  parser.add_option("-o", "--outDir").dest("outDir").help("output directory").set_default("analysis"); // type='string')
  parser.add_option("-n", "--njobs").dest("njobs").help("# jobs to run in parallel").set_default(0).type("int");

  optparse::Values options = parser.parse_args(argc, argv);
  vector<string> args = parser.args();

  runBJetEnergyPeak ebpeak( options["inDir"], options["outDir"] ) ;

  //read list of samples
  ifstream infile; 
  infile.open(options["json"]); 

  Json::Value root;   // will contains the root value after parsing.
  Json::Reader reader;
  bool parsedSuccess = reader.parse(infile, root);
  if(not parsedSuccess)
  {
     // Report failures and their locations in the document.
     cout<<"Failed to parse JSON"<<endl
     <<reader.getFormatedErrorMessages()
     <<endl;
     return 1;
  }

  std::unordered_map<string, sampleInfo> sampleLists;
  auto datafiles = root.getMemberNames();
  for( string data : datafiles ) {
    std::cout<<data<<std::endl;
    const Json::Value subdata = root[data];
    sampleLists[data] = sampleInfo( subdata ) ;
  }

    /*
    //#run the analysis jobs
    if opt.njobs == 0:
        for inFileURL, outFileURL, xsec in taskList:
            runBJetEnergyPeak(inFileURL=inFileURL, outFileURL=outFileURL, xsec=xsec)
    else:
        from multiprocessing import Pool
        pool = Pool(opt.njobs)
        pool.map(runBJetEnergyPeakPacked,taskList)
    */
  for( auto sample : sampleLists ) {
        string inFileURL  = TString::Format("%s/%s.root",options["inDir"].c_str(),sample.first.c_str()).Data();  //# For NTU cluster,
        //#if not os.path.isfile(inFileURL): continue
        double xsec = sample.second.getXsec();
        string outFileURL = TString::Format("%s/%s.root",options["outDir"].c_str(),sample.first.c_str()).Data();
        //taskList.append( (inFileURL,outFileURL,xsec) )
        cout<<inFileURL<<endl;
        cout<<outFileURL<<endl; 
        cout<<xsec<<endl;
  } 

	return 0;
}

/*    
def runBJetEnergyPeak(inFileURL, outFileURL, xsec=None):



"""
Wrapper to be used when run in parallel
"""
def runBJetEnergyPeakPacked(args):
    
    try:
        return runBJetEnergyPeak(inFileURL=args[0],
                                 outFileURL=args[1],
                                 xsec=args[2])
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],args[0])
        print 50*'<'
        return False


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
        
    #create the analysis jobs
    taskList = []
    for sample, sampleInfo in samplesList: 
        #inFileURL  = 'root://cmsxrootd.fnal.gov//%s/%s.root' % (opt.inDir,sample)  ## For other sites,
        inFileURL  = '%s/%s.root' % (opt.inDir,sample)  # For NTU cluster,
        #if not os.path.isfile(inFileURL): continue
        xsec=sampleInfo[0] if sampleInfo[1]==0 else None        
        outFileURL = '%s/%s.root' % (opt.outDir,sample)
        taskList.append( (inFileURL,outFileURL,xsec) )

    #run the analysis jobs
    if opt.njobs == 0:
        for inFileURL, outFileURL, xsec in taskList:
            runBJetEnergyPeak(inFileURL=inFileURL, outFileURL=outFileURL, xsec=xsec)
    else:
        from multiprocessing import Pool
        pool = Pool(opt.njobs)
        pool.map(runBJetEnergyPeakPacked,taskList)

    #all done here
    print 'Analysis results are available in %s' % opt.outDir
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())

*/
