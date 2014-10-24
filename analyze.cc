





#include <signal.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sstream>
#include <vector>
#include <math.h>
#include <string>

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TMath.h"
#include "TGraph.h"

//
//	Constants
//
#define MAXCHPEREVENT 500
#define NUMPHIS 36
#define NUMETAS 13
#define NUMBINS 32
#define NUMDEPTHS 2
#define NUMTUBETYPES 2
#define PEDCUT_1TS 10
#define PEDCUT_2TS 18

double qBins_1TS[NUMBINS+1] = {
	0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.,
	17., 19., 21., 23., 25., 27., 29., 32., 35., 38., 41., 45., 49., 53., 58.,
	63., 68};
double qBins_2TS[NUMBINS+1] = {
	0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18., 20.,
	22., 26., 30., 34., 38., 42., 50., 58., 66., 74., 82., 96., 114., 130.,
	162., 194., 256.};

//
//	Global Vars
//
bool runAll = true;

using namespace std;
using namespace ROOT;

struct RawInput
{
	TTree *tree;
	Int_t numChannels;
	Int_t phi[MAXCHPEREVENT];
	Int_t eta[MAXCHPEREVENT];
	Int_t depth[MAXCHPEREVENT];

	char tubeName[100];
	Int_t msgCounter;
	Int_t mtrI;
	Int_t mtrV;
	Int_t reelPos;
	Int_t orbitNum;
	Int_t index;
	uint32_t driverStatus;

	uint16_t rawC0[MAXCHPEREVENT][NUMBINS];
	uint16_t rawC1[MAXCHPEREVENT][NUMBINS];
	uint16_t rawC2[MAXCHPEREVENT][NUMBINS];
	uint16_t rawC3[MAXCHPEREVENT][NUMBINS];
};

struct ServiceVariables
{
	double *qBins_toUse;
	double cutToUse;
	int verbosity;
};

struct HistoOutput
{
	TFile *rootFile;
	TH1D *hRAW1D_wSrc[NUMPHIS][NUMETAS][NUMDEPTHS][NUMTUBETYPES];
	TH1D *hRAW1D_F_wSrc[NUMPHIS][NUMETAS][NUMDEPTHS][NUMTUBETYPES];
	TH1D *hRAW1D_B_wSrc[NUMPHIS][NUMETAS][NUMDEPTHS][NUMTUBETYPES];
	TH1D *hRAW1D_NoSrc[NUMPHIS][NUMETAS][NUMDEPTHS];

	TH1D *hMeansOfPED[NUMPHIS][NUMETAS][NUMDEPTHS];
	TH1D *hRMSsOfPED[NUMPHIS][NUMETAS][NUMDEPTHS];

	TGraph *gSPESigVsOrbit[NUMPHIS][NUMETAS][NUMDEPTHS][NUMTUBETYPES];
	TGraph *gSPESigVsReel[NUMPHIS][NUMETAS][NUMDEPTHS][NUMTUBETYPES];

	TH2D *hRMSofRMS;
	TH2D *hMeansOfRMS;
	TH2D *hRMSofMeans;
	TH2D *hMeansOfMeans;
};

struct SourceTubeInfo
{
	bool isSourced;
	double tubeStart;
	double tubeEnd;
	int wedge;
	string tubeName;
	char groove;
};

struct TubesMap
{
	SourceTubeInfo tubes[NUMPHIS][NUMETAS][NUMTUBETYPES];
};

//
//	Sub Functions
//
void initOutput(HistoOutput&, ServiceVariables&, string);
vector<string>& readFiles(string, int verbosity=0);
void branch(RawInput&);
void analyze(int&, RawInput&, HistoOutput&, TubesMap&, ServiceVariables&);
int readTubesMap(string fileNameMap, TubesMap&, 
		int verbosity=0);
double addHisto(int, RawInput& raw, TH1D*, ServiceVariables&, int v=0);
void sigHandler(int);
double computeSPESignal(int, RawInput&, ServiceVariables&, int v=0);


//
//	Main Entry Point
//
int main(int argc, char** argv)
{
	//
	//	Register Signals
	//
	signal(SIGABRT, &sigHandler);
	signal(SIGTERM, &sigHandler);
	signal(SIGINT, &sigHandler);

	int verbosity = atoi(argv[1]);
	string listOfFiles = argv[2];
	string rootOutFileName = argv[3];
	string fileNameMap = argv[4];
	int iTS = atoi(argv[5]);


	int globalEvents = 0;
	RawInput raw;
	HistoOutput out;
	TubesMap tubesMap;
	ServiceVariables service;
	service.qBins_toUse = qBins_1TS;
	service.cutToUse = PEDCUT_2TS;
	service.verbosity = verbosity;
	if (iTS==1)
		service.qBins_toUse = qBins_1TS;
	else if (iTS==2)
		service.qBins_toUse = qBins_2TS;
	readTubesMap(fileNameMap, tubesMap, verbosity);
	initOutput(out, service, rootOutFileName);
	vector<string> &vInputFiles = readFiles(listOfFiles, verbosity);

	for (vector<string>::iterator it=vInputFiles.begin(); 
			(it!=vInputFiles.end()) && runAll==true; ++it)
	{
		cout << "### Processing File: " << *it << endl;
		TFile *inFile = TFile::Open(it->c_str());
		raw.tree = (TTree*)inFile->Get("Events");
		branch(raw);
		analyze(globalEvents, raw, out, tubesMap, service);
		inFile->Close();
	}

	
	//
	//	Weigh 1D histos and Create 2D Noise Maps
	//
	for (int iphi=0; iphi<NUMPHIS; iphi++)
		for (int ieta=0; ieta<NUMETAS; ieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				//
				//	Noise
				//
				double rmsOfrms, meanOfrms, rmsOfMean, meanOfMean;
				rmsOfrms = out.hRMSsOfPED[iphi][ieta][idepth]->GetRMS();
				meanOfrms = out.hRMSsOfPED[iphi][ieta][idepth]->GetMean();
				rmsOfMean = out.hMeansOfPED[iphi][ieta][idepth]->GetRMS();
				meanOfMean = out.hMeansOfPED[iphi][ieta][idepth]->GetMean();
				out.hRMSofRMS->Fill(2*iphi+1, ieta+29, rmsOfrms);
				out.hMeansOfRMS->Fill(2*iphi+1, ieta+29, meanOfrms);
				out.hRMSofMeans->Fill(2*iphi+1, ieta+29, rmsOfMean);
				out.hMeansOfMeans->Fill(2*iphi+1, ieta+29, meanOfMean);

				//
				//	Weighing
				//
				for (int iTubeType=0; iTubeType<NUMTUBETYPES; iTubeType++)
				{
					for (int iBin=0; iBin<NUMBINS; iBin++)
					{
						double newContent = 
						out.hRAW1D_wSrc[iphi][ieta][idepth][iTubeType]->GetBinContent(iBin+1)/out.hRAW1D_wSrc[iphi][ieta][idepth][iTubeType]->GetBinWidth(iBin+1);
						out.hRAW1D_wSrc[iphi][ieta][idepth][iTubeType]->SetBinContent(iBin+1, newContent);

						if (iTubeType==0)
						{
							newContent = 
							out.hRAW1D_NoSrc[iphi][ieta][idepth]->GetBinContent(
									iBin+1)/out.hRAW1D_NoSrc[iphi][ieta][idepth]->GetBinWidth(iBin+1);
							out.hRAW1D_NoSrc[iphi][ieta][idepth]->SetBinContent(
									iBin+1, newContent);
						}
					}

					//	
					//	Save Graphs
					//
					out.rootFile->cd("SRC/SPESignals/VsOrbit");
					out.gSPESigVsOrbit[iphi][ieta][idepth][iTubeType]->Write();

					out.rootFile->cd("SRC/SPESignals/VsReelPos");
					out.gSPESigVsReel[iphi][ieta][idepth][iTubeType]->Write();
				}
			}


	out.rootFile->Write();
	out.rootFile->Close();
	return 0;
}

//
//
//
void branch(RawInput& raw)
{
	raw.tree->SetBranchAddress("numChannels", &raw.numChannels);
	raw.tree->SetBranchAddress("phi", &raw.phi);
	raw.tree->SetBranchAddress("eta", &raw.eta);
	raw.tree->SetBranchAddress("depth", &raw.depth);

	raw.tree->SetBranchAddress("motorCurrent", &raw.mtrI);
	raw.tree->SetBranchAddress("motorVoltage", &raw.mtrV);
	raw.tree->SetBranchAddress("reelPos", &raw.reelPos);
	raw.tree->SetBranchAddress("orbitNumber", &raw.orbitNum);
	raw.tree->SetBranchAddress("index", &raw.index);
	raw.tree->SetBranchAddress("driverStatus", &raw.driverStatus);
	raw.tree->SetBranchAddress("tubeName", &raw.tubeName);

	raw.tree->SetBranchAddress("rawC0", raw.rawC0);
	raw.tree->SetBranchAddress("rawC1", raw.rawC1);
	raw.tree->SetBranchAddress("rawC2", raw.rawC2);
	raw.tree->SetBranchAddress("rawC3", raw.rawC3);

	return;
}

//
//
//
void analyze(int &globalEvents, RawInput& raw, HistoOutput& out,
		TubesMap &map, 	ServiceVariables &service)
{
	cout << "### Analyzing a New Tree..." << endl;
	string someName = "UNKNOWN";
	int pointD1 = 0;
	int pointD2 = 0;

	//
	//	For each event
	//
	for (int iEvent=0; (iEvent<raw.tree->GetEntries()) && runAll==true;
			iEvent++)
	{
		raw.tree->GetEntry(iEvent);

		//
		//	Print #events processed.
		//
//		if (iEvent%1000 == 0)
//			cout << "### Processing Event=" << iEvent << endl;

		//
		//	Print the TubeName
		//
		string sourceTubeName(raw.tubeName);
		if (someName != sourceTubeName)
		{
			someName = sourceTubeName;
			cout << pointD1 << "  " << pointD2 << endl;
			pointD1 = 0;
			pointD2 = 0;
			cout << "### Current Source Tube: " << someName << endl;
		}

		int srcEta, srcPhi, srciTubeType, srcWedge, srcTubeNumber;
		char srcTubeType;
		sscanf(sourceTubeName.c_str(), "HFP%d_ETA%d_PHI%d_T%d%c%*[^]",
				&srcWedge, &srcEta, &srcPhi, &srcTubeNumber,
				&srcTubeType);
		if (srcTubeType=='A' || srcTubeType=='_')
			srciTubeType = 0;
		else if (srcTubeType=='B')
			srciTubeType = 1;
		else 
			cout << "### ERROR: Wrong Tube Type!" << endl;
		int srciEta = abs(srcEta) - 29;
		int srciPhi = (srcPhi-1)/2;

		//
		//	Check the reel position validity
		//
		if (service.verbosity>2)
			cout << "### ReelPos=" << raw.reelPos << "  start="
				<< map.tubes[srciPhi][srciEta][srciTubeType].tubeStart
				<< "  end="
				<< map.tubes[srciPhi][srciEta][srciTubeType].tubeEnd
				<< endl;
		

		//
		//	Iterate through all the channels
		//
		for (int iCh=0; iCh<raw.numChannels; iCh++)
		{
			int iphi = raw.phi[iCh];	int iiphi=(iphi-1)/2;
			int ieta = raw.eta[iCh]; int iieta=abs(ieta)-29;
			int depth = raw.depth[iCh]; int idepth=depth-1;

			if (service.verbosity>0)
				cout << "### phi=" << iphi << "  eta=" << ieta << "  depth="
					<< depth << endl;

			//
			//	If this is a channel with Source
			//
			if (srciEta==iieta && srciPhi==iiphi)
			{
				double speSignal = computeSPESignal(iCh, raw, service, 
						service.verbosity);
				if (depth==1)
				{
					if (map.tubes[srciPhi][srciEta][srciTubeType].tubeStart <= 
							raw.reelPos && 
							map.tubes[srciPhi][srciEta][srciTubeType].tubeEnd >= 
							raw.reelPos)
						addHisto(iCh, raw, 
							out.hRAW1D_wSrc[iphi][ieta][idepth][srciTubeType],
							service, service.verbosity);
					out.gSPESigVsOrbit[iphi][ieta][idepth][srciTubeType]
						->SetPoint(pointD1, raw.orbitNum, speSignal);
					out.gSPESigVsReel[iphi][ieta][idepth][srciTubeType]->SetPoint(
						pointD1, raw.orbitNum, speSignal);
					pointD1++;
				}
				else
				{
					if ((map.tubes[iphi][ieta][srciTubeType].tubeEnd -
								300)>=raw.reelPos &&
							map.tubes[iphi][ieta][srciTubeType].tubeStart <= 
							raw.reelPos)
						addHisto(iCh, raw,
								out.hRAW1D_wSrc[iphi][ieta][idepth][srciTubeType],
								service, service.verbosity);
					out.gSPESigVsOrbit[iphi][ieta][idepth][srciTubeType]
						->SetPoint(pointD2, raw.orbitNum, speSignal);
					out.gSPESigVsReel[iphi][ieta][idepth][srciTubeType]->SetPoint(
						pointD2, raw.orbitNum, speSignal);
					pointD2++;
				}

			}
			else
			{
				//
				//	This is a channel without a source,
				//	apply geometric isolation and fill the NoSrc Histo
				//
				if ((srciPhi==1 && abs(srciPhi-iphi)<5) || (srciPhi==71 && 
							abs(srciPhi - iphi)<5))
					continue;
/*				if (abs(srciEta - ieta)>5 && abs(srciPhi - iphi)>5)
				{
					addHisto(iCh, raw, 
							out.hRAW1D_NoSrc[iphi][ieta][idepth], service, 0);

					//	
					//	To produce 2D maps of RMSs/Means Event by Event
					//
					TH1D *h = new TH1D("h1", "h1", 32, service.qBins_toUse);
					for (int iBin=0; iBin<NUMBINS; iBin++)
					{
						int allCapsSum = raw.rawC0[iCh][iBin] + 
							raw.rawC1[iCh][iBin] + 
							raw.rawC2[iCh][iBin] + raw.rawC3[iCh][iBin];
						h->Fill(service.qBins_toUse[iBin], allCapsSum);
					}
					h->GetXaxis()->SetRange(0, 15);
					out.hMeansOfPED[iphi][ieta][idepth]->Fill(h->GetMean());
					out.hRMSsOfPED[iphi][ieta][idepth]->Fill(h->GetRMS());
					delete h;
				}
				*/
			}
		}

		globalEvents++;
	}
}

//
//	Compute SPE signal
//
double computeSPESignal(int iCh, RawInput& raw, ServiceVariables& service, int v)
{
	double totSum = 0;
	double speSum = 0;
	for (int iBin=0; iBin<NUMBINS; iBin++)
	{		
		int allCapsSum = raw.rawC0[iCh][iBin] + raw.rawC1[iCh][iBin] + 
			raw.rawC2[iCh][iBin] + raw.rawC3[iCh][iBin];
		totSum += allCapsSum;
		double charge = (service.qBins_toUse[iBin] + 
				service.qBins_toUse[iBin+1])/2.;
		if (charge>service.cutToUse)
			speSum += allCapsSum*charge;
	}

	return speSum/totSum;
}

//
//	Signal Handler
//
void sigHandler(int sig)
{
	cout << "### Signal: " << sig << " caughter. Exiting..." << endl;
	runAll = false;
}

//
//	Add RAW histo to a Total Histogram
//	Returns SPE Signal
//
double addHisto(int iCh, RawInput& raw, TH1D *hTotal, ServiceVariables &service, 
		int v)
{
	if (v>0)
		cout << "### Adding Histos..." << endl;

	double totSum = 0;
	double speSum = 0;
	for (int iBin=0; iBin<NUMBINS; iBin++)
	{
		int allCapsSum = raw.rawC0[iCh][iBin] + raw.rawC1[iCh][iBin] + 
			raw.rawC2[iCh][iBin] + raw.rawC3[iCh][iBin];
		totSum += allCapsSum;
		if (hTotal->GetBinCenter(iBin+1)>=service.cutToUse)
		{
			double charge = hTotal->GetBinCenter(iBin+1);
			speSum += allCapsSum*charge;
		}

		if (v>1)
			cout << "### AllCapsSum=" << allCapsSum << endl;

		hTotal->Fill(service.qBins_toUse[iBin], allCapsSum);
	}


	return speSum/totSum;
}

//
//	
//
void initOutput(HistoOutput &out, ServiceVariables &service,
		string fileName)
{
	out.rootFile = new TFile(fileName.c_str(), "RECREATE");
	char histName[200];

	cout << "### Start HF Output Initialization..." << endl;
	
	//
	//	Init HF Histos
	//
	out.rootFile->mkdir("SRC");
	out.rootFile->mkdir("NOSRC");

	out.rootFile->cd("SRC");
	gDirectory->mkdir("1D");
	gDirectory->mkdir("1D_F");
	gDirectory->mkdir("1D_B");

	gDirectory->mkdir("SPESignals");
	gDirectory->cd("SPESignals");
	gDirectory->mkdir("VsOrbit");
	gDirectory->mkdir("VsReelPos");

	out.rootFile->cd("NOSRC");
	out.hRMSofRMS = new TH2D("RMSofRMS", "RMSofRMS", 72, 0, 72,
			13, 29, 42);
	out.hMeansOfRMS = new TH2D("MeansOfRMS", "MeansOfRMS", 72, 0, 72,
			13, 29, 42);
	out.hRMSofMeans = new TH2D("RMSofMeans", "RMSofMeans", 72, 0, 72,
			13, 29, 42);
	out.hMeansOfMeans = new TH2D("MeansOfMeans", "MeansOfMeans", 72, 0, 72,
			13, 29, 42);

	gDirectory->mkdir("1D");
	gDirectory->mkdir("PEDMEANS");
	gDirectory->mkdir("PEDRMSS");
//	gDirectory->mkdir()

	for (int iphi=0; iphi<NUMPHIS; iphi++)
		for (int ieta=0; ieta<NUMETAS; ieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
				for (int iTubeType=0; iTubeType<NUMTUBETYPES; 
						iTubeType++)
				{
					out.rootFile->cd("SRC/1D");
					sprintf(histName, "1D_wSrc_%d_PHI%d_ETA%d_D%d",
							iTubeType, 2*iphi+1, ieta+29, idepth+1);
					out.hRAW1D_wSrc[iphi][ieta][idepth][iTubeType] = new 
						TH1D(histName, histName, 32, service.qBins_toUse);

					out.rootFile->cd("SRC/1D_F");
					sprintf(histName, "1D_F_wSrc_%d_PHI%d_ETA%d_D%d",
							iTubeType, 2*iphi+1, ieta+29, idepth+1);
					out.hRAW1D_F_wSrc[iphi][ieta][idepth][iTubeType] = new 
						TH1D(histName, histName, 32, service.qBins_toUse);

					out.rootFile->cd("SRC/1D_B");
					sprintf(histName, "1D_B_wSrc_%d_PHI%d_ETA%d_D%d",
							iTubeType, 2*iphi+1, ieta+29, idepth+1);
					out.hRAW1D_B_wSrc[iphi][ieta][idepth][iTubeType] = new 
						TH1D(histName, histName, 32, service.qBins_toUse);

					out.rootFile->cd("SRC/SPESignals/VsOrbit");
					sprintf(histName, "SPEsigVsOrbit_%d_PHI%d_ETA%d_D%d",
							iTubeType, 2*iphi+1, ieta+29, idepth+1);
					out.gSPESigVsOrbit[iphi][ieta][idepth][iTubeType] = new 
						TGraph();
					out.gSPESigVsOrbit[iphi][ieta][idepth][iTubeType]->SetName(
							histName);

					out.rootFile->cd("SRC/SPESignals/VsReelPos");
					sprintf(histName, "SPEsigVsReelPos_%d_PHI%d_ETA%d_D%d",
							iTubeType, 2*iphi+1, ieta+29, idepth+1);
					out.gSPESigVsReel[iphi][ieta][idepth][iTubeType] = new 
						TGraph();
					out.gSPESigVsReel[iphi][ieta][idepth][iTubeType]->SetName(
							histName);
				}

	for (int iphi=0; iphi<NUMPHIS; iphi++)
		for (int ieta=0; ieta<NUMETAS; ieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				out.rootFile->cd("NOSRC/1D");
				sprintf(histName, "1D_NoSrc_PHI%d_ETA%d_D%d",
						2*iphi+1, ieta+29, idepth+1);
				out.hRAW1D_NoSrc[iphi][ieta][idepth] = new 
					TH1D(histName, histName, 32, service.qBins_toUse);

				out.rootFile->cd("NOSRC/PEDMEANS");
				sprintf(histName, "PEDMeans_PHI%d_ETA%d_D%d",
						2*iphi+1, ieta+29, idepth+1);
				out.hMeansOfPED[iphi][ieta][idepth] = new TH1D(histName, 
						histName, 1000, 0, 20);

				out.rootFile->cd("NOSRC/PEDRMSS");
				sprintf(histName, "PEDRMSs_PHI%d_ETA%d_D%d",
						2*iphi+1, ieta+29, idepth+1);
				out.hRMSsOfPED[iphi][ieta][idepth] = new TH1D(histName,
						histName, 1000, 0, 5);
			}

	cout << "### Initialization is DONE" << endl;
}

//
//
//
vector<string>& readFiles(string listFiles, int v)
{
	vector<string> *vFileNames = new vector<string>;
	ifstream in(listFiles.c_str());
	string fileN;

	while (in >> fileN)
		vFileNames->push_back(fileN);

	return *vFileNames;
}

//
//	Read in Goemetry Coefficients and Tubes Map
//
int readTubesMap(string fileNameMap, TubesMap &map,
		int verbosity)
{
	cout << "### Reading in Tubes Map..." << endl;

	ifstream mapFile(fileNameMap.c_str());
	if (!mapFile)
	{
		cout << "### ERROR: No Tubes Map File..." << endl;
		return -1;
	}

	string tubeName, wedgeName, swaste;
	int ieta, iphi, tubeStart, tubeEnd, iwaste;
	char groove;

	int counter = 0;
	while (mapFile >> tubeName)
	{
		mapFile >> wedgeName >> iwaste >> swaste 
			>> ieta >> ieta >> iphi >> iwaste >> iwaste 
			>> tubeStart >> tubeEnd >> groove;

		int iiphi = (iphi-1)/2;
		int iieta = abs(ieta) - 29;
		int checkRBX, checkPhi, checkEta, checkTubeNumber;
		char tubeType;
		sscanf(tubeName.c_str(), "HFP%d_ETA%d_PHI%d_T%d%c%*[^]",
				&checkRBX, &checkEta, &checkPhi, &checkTubeNumber,
				&tubeType);
		int iTubeType = 0;
		if (tubeType=='A' || tubeType=='_')
			iTubeType = 0;
		else if (tubeType=='B')
			iTubeType = 1;
		else
			cout << "### ERROR: Wrong Tube Type..." << endl;

		map.tubes[iiphi][iieta][iTubeType].tubeStart = tubeStart;
		map.tubes[iiphi][iieta][iTubeType].tubeEnd = tubeEnd;
		map.tubes[iiphi][iieta][iTubeType].tubeName = tubeName;
		map.tubes[iiphi][iieta][iTubeType].groove = groove;

		if (verbosity > 0)
			cout << tubeName << "  " << ieta << "  " << iphi << "  "
				<< tubeStart << "  " << tubeEnd << "  " << groove << endl;

		counter++;
	}

	cout << "### DONE Reading Tubes Map File!" << endl;
	return 0;
}

















