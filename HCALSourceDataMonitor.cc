// -*- C++ -*-
//
// Package:    HCALSourceDataMonitor
// Class:      HCALSourceDataMonitor
// 
/**\class HCALSourceDataMonitor HCALSourceDataMonitor.cc HCALSourcing/HCALSourceDataMonitor/src/HCALSourceDataMonitor.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seth Cooper,32 4-B03,+41227675652,
//         Created:  Tue Jul  2 10:47:48 CEST 2013
// $Id: HCALSourceDataMonitor.cc,v 1.7 2013/07/23 08:41:57 scooper Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/HcalRecHit/interface/HcalSourcePositionData.h"
#include "TBDataFormats/HcalTBObjects/interface/HcalTBTriggerData.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TTree.h"

//
//	Basic Constants
//
#define MAXCHPEREVENT 500
#define NUMPHIS 36
#define NUMETAS 13
#define NUMBINS 32
#define NUMDEPTHS 2
#define NUMTUBETYPES 2

//
// class declaration
//
class HCALSourceDataMonitor : public edm::EDAnalyzer {
   public:
      explicit HCALSourceDataMonitor(const edm::ParameterSet&);
      ~HCALSourceDataMonitor();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	  //
	  //	OLD!!!
	  //	Isn't being used anymore!!!
	  //
      bool isDigiAssociatedToSourceTube(const HcalDetId& detId, std::string tubeName);
	  
	  //
	  //	These are the ANALYSIS Methods. Names are self-explanatory
	  //
	  int analyze_NoSrc(const edm::Event&, const edm::EventSetup&);
	  int analyze_wSrc(const edm::Event&, const edm::EventSetup&);

	  //
	  //	Members
	  //
	  std::string _rootFileName;
	  TFile *_rootFile;
	  bool printRawHistograms_;
      bool checkIntegrals_;
      bool selectDigiBasedOnTubeName_;
      int naiveEvtNum_;
	  int _iEvent;
	  int _analysisType;
      std::set<HcalDetId> emptyHistogramChannelsSet;
	  int _verbosity;

	  //
	  //	ROOT Containers/Output Related
	  //
	  TTree *_tree;
	  Int_t _numChannels;
	  Int_t _msgCounter;
	  Int_t _mtrI;
	  Int_t _mtrV;
	  Int_t _reelPos;
	  Int_t _orbitNum;
	  Int_t _index;
	  uint32_t _driverStatus;
	  char _tubeName[100];
	  Int_t _phi[MAXCHPEREVENT];
	  Int_t _eta[MAXCHPEREVENT];
	  Int_t _depth[MAXCHPEREVENT];
	  uint16_t _rawC0[MAXCHPEREVENT][32];
	  uint16_t _rawC1[MAXCHPEREVENT][32];
	  uint16_t _rawC2[MAXCHPEREVENT][32];
	  uint16_t _rawC3[MAXCHPEREVENT][32];
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

std::string intToString(int num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}

//
// constructors and destructor
//
HCALSourceDataMonitor::HCALSourceDataMonitor(const edm::ParameterSet& iConfig) :
	_rootFileName(iConfig.getUntrackedParameter<std::string>("RootFileName")),
  printRawHistograms_ (iConfig.getUntrackedParameter<bool>("PrintRawHistograms",false)),
  checkIntegrals_ (iConfig.getUntrackedParameter<bool>("CheckHistogramIntegrals",false)),
  selectDigiBasedOnTubeName_ (iConfig.getUntrackedParameter<bool>("SelectDigiBasedOnTubeName",true)),
	_analysisType(iConfig.getUntrackedParameter<int>("AnalysisType")),
	_verbosity(iConfig.getUntrackedParameter<int>("Verbosity"))
{
	_iEvent = 0;
	_rootFile = new TFile(_rootFileName.c_str(), "RECREATE");
	_rootFile->cd();
	
	//
	//	ROOT Initializations
	//
	_tree = new TTree("Events", "Events");
	_tree->Branch("motorCurrent", &_mtrI);
	_tree->Branch("motorVoltage", &_mtrV);
	_tree->Branch("reelPos", &_reelPos);
	_tree->Branch("orbitNumber", &_orbitNum);
	_tree->Branch("index", &_index);
	_tree->Branch("driverStatus", &_driverStatus);
	_tree->Branch("tubeName", &_tubeName, "tubeName/C");
	_tree->Branch("numChannels", &_numChannels, "numChannels/I");
	_tree->Branch("phi", _phi, "phi[numChannels]/I");
	_tree->Branch("eta", _eta, "eta[numChannels]/I");
	_tree->Branch("depth", _depth, "depth[numChannels]/I");

	_tree->Branch("rawC0", _rawC0, "rawC0[numChannels][32]/s");
	_tree->Branch("rawC1", _rawC1, "rawC1[numChannels][32]/s");
	_tree->Branch("rawC2", _rawC2, "rawC2[numChannels][32]/s");
	_tree->Branch("rawC3", _rawC3, "rawC3[numChannels][32]/s");
}


HCALSourceDataMonitor::~HCALSourceDataMonitor()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
//
bool HCALSourceDataMonitor::isDigiAssociatedToSourceTube(const HcalDetId& detId, std::string tubeName)
{
  return false;
}

int HCALSourceDataMonitor::analyze_wSrc(const edm::Event& iEvent,

		const edm::EventSetup& iSetup)
{
	using namespace std;
	using namespace edm;

	//
	//	HCAL-related
	//	-> get the map
	//	-> get the source pos info
	//	-> get trigger info
	//
	vector<Handle<HcalHistogramDigiCollection> > hcalHistDigiCollHandleVec;
	iEvent.getManyByType(hcalHistDigiCollHandleVec);

	//
	//	Map
	//
	edm::ESHandle<HcalDbService> pSetup;
	iSetup.get<HcalDbRecord>().get(pSetup);

	Handle<HcalSourcePositionData> hcalSourcePositionDataHandle;
	iEvent.getByType(hcalSourcePositionDataHandle);
	const HcalSourcePositionData *spd = hcalSourcePositionDataHandle.product();

	Handle<HcalTBTriggerData> hcalTBTriggerDataHandle;
	iEvent.getByType(hcalTBTriggerDataHandle);
	const HcalTBTriggerData *triggerData = hcalTBTriggerDataHandle.product();

	string currentSourceTubeName = spd->tubeNameFromCoord();
	sprintf(_tubeName, "%s", currentSourceTubeName.c_str());
	_mtrI = spd->motorCurrent();
	_mtrV = spd->motorVoltage();
	_reelPos = spd->reelCounter();
	_orbitNum = triggerData->orbitNumber();
	_index = spd->indexCounter();
	_driverStatus = (uint32_t)spd->status();

	
	//
	//	Go over every channel
	//
	int iCh = 0;
	for (vector<Handle<HcalHistogramDigiCollection> >::iterator
			it=hcalHistDigiCollHandleVec.begin(); 
			it!=hcalHistDigiCollHandleVec.end(); ++it)
	{
		if (!it->isValid())
		{
			cout << "### Invalid Digi Collection Found;" << endl;
			continue;
		}

		const HcalHistogramDigiCollection &hcalHistDigiColl=*(*it);
		HcalHistogramDigiCollection::const_iterator idigi;
		for (idigi=hcalHistDigiColl.begin(); idigi!=hcalHistDigiColl.end();
				idigi++)
		{
			const HcalDetId detId = idigi->id();
			int eta = detId.ieta();
			int phi = detId.iphi();
			int depth = detId.depth();

			_eta[iCh] = eta;
			_phi[iCh] = phi;
			_depth[iCh] = depth;

			//
			//	For Debugging
			//
			if (_verbosity>1)
				cout << "### eta=" << eta << "  phi=" << phi << "  depth=" 
					<< depth << endl;

			for (int iBin=0; iBin<NUMBINS; iBin++)
			{
				_rawC0[iCh][iBin] = idigi->get(0, iBin);
				_rawC1[iCh][iBin] = idigi->get(1, iBin);
				_rawC2[iCh][iBin] = idigi->get(2, iBin);
				_rawC3[iCh][iBin] = idigi->get(3, iBin);

//				int allCapsSum = idigi->get(0, iBin) + idigi->get(1, iBin) + 
//				idigi->get(2, iBin) + idigi->get(3, iBin);
//				if (allCapsSum==0 && iBin<15)
//					cout << "$$$$$$$$$$$$$$$$$$$$$ WARNING: ZERO allCapsSum for Bin: "
//						<< iBin
//						<< endl;
			}
			iCh++;
		}
	}

	_iEvent++;
	_numChannels = iCh;

	_rootFile->cd();
	_tree->Fill();
	if (_verbosity>1)
		cout << "### Number of Channels Per Event=" << _numChannels << endl;

	return 0;
}


// ------------ method called for each event  ------------
void
HCALSourceDataMonitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

	analyze_wSrc(iEvent, iSetup);

}


// ------------ method called once each job just before starting event loop  ------------
void 
HCALSourceDataMonitor::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HCALSourceDataMonitor::endJob() 
{
	using namespace std;

	_rootFile->cd();
	_tree->Write();
	_rootFile->Close();

  std::cout << "The following " << emptyHistogramChannelsSet.size() << " channels had at least one empty histogram:" << std::endl;
  // print out list of empty hist channels
  for(std::set<HcalDetId>::const_iterator itr = emptyHistogramChannelsSet.begin(); itr != emptyHistogramChannelsSet.end(); ++itr)
  {
    std::cout << *itr << std::endl;
  }
  std::cout << "End of list of empty histogram channels." << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
HCALSourceDataMonitor::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HCALSourceDataMonitor::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HCALSourceDataMonitor::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HCALSourceDataMonitor::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HCALSourceDataMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HCALSourceDataMonitor);
