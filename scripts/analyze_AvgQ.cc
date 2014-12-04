


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
#include "TF1.h"

using namespace std;
using namespace ROOT;

#define NUMPHIS 36
#define NUMETAS 13
#define NUMTUBETYPES 2
#define NUMDEPTHS 2
#define NUMBINS 32
#define PEDCUT_1TS 10
#define PEDCUT_2TS 18

double GEVPER25NS[NUMDEPTHS] = { 744E-6, 706E-6};
//	
//	Define Source Activity Correction. 
//	This value shows how much source strength decreased from November 2013
//	till the month the sourcing data has been taken
//	[0] - HFP
//	[1] - HFM
//
double SOURCEACTIVITYCORRECTION[2] = {0.91, 0.876};
//#define SOURCEACTIVITYCORRECTION 0.92398
//#define SOURCEACTIVITYCORRECTION 0.88117

//
//	Adapt this for either HFP or HFM
//	#swaps in HFP is 34
//	#swaps in HFM is 32
//
//#define NUMSWAPS 32
#define NUMSWAPS 34

#include "map.cc"
#include "readSwaps.cc"

Double_t qBins_1TS[NUMBINS+1] = {
	0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 17.,
	19., 21., 23., 25., 27., 29., 32., 35., 38., 41., 45., 49., 53., 58.,
	63., 68};

Double_t qBins_2TS[NUMBINS+1] = {
	0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 14., 16., 18., 20., 22.,
	26., 30., 34., 38., 42., 50., 58., 66., 74., 82., 96., 114., 130., 
	162., 194., 256.};

Double_t *qBins_toUse;
Double_t cut_toUse;

struct Coordinates
{
	int iphi;
	int ieta;
	int idepth;
};

Coordinates swappedChs[NUMSWAPS];

struct SourceTubeInfo
{
	bool isSourced;
	double tubeStart;
	double tubeEnd;
	string tubeName;
	char groove;
	double em;
	double had;
	double uncert;
};

struct TubesMap
{
	SourceTubeInfo tubes[NUMPHIS][NUMETAS][NUMTUBETYPES];
};

struct TOut
{
	TFile *rootFile;
	ofstream txtFile;

	//
	//	The last dimension differentiates between w/ OF and w/o OF
	//	[0] -> w/ OF
	//	[1] -> w/o OF
	//
	TProfile *pSignalvsETA_UNCORR[NUMTUBETYPES][NUMDEPTHS][2];
	TProfile *pSignalvsETA_CORR[NUMTUBETYPES][NUMDEPTHS][2];
	TProfile *pSFSBvsETA[NUMTUBETYPES][NUMDEPTHS][2];
	TProfile *pSB23toSvsETA[NUMTUBETYPES][NUMDEPTHS][2];

	TH2D *hPedMean_S[NUMTUBETYPES][NUMDEPTHS];
	TH2D *hPedMean_B[NUMTUBETYPES][NUMDEPTHS];
	TH2D *hPedSigma_S[NUMTUBETYPES][NUMDEPTHS];
	TH2D *hPedSigma_B[NUMTUBETYPES][NUMDEPTHS];
	TH2D *hSwapRatio[NUMTUBETYPES][NUMDEPTHS];
//	TH2D *hADC2GeV_OV2[NUMTUBETYPES][NUMDEPTHS][2];

	TH1D *hSignal_UNCORR[NUMTUBETYPES][NUMDEPTHS][2];
	TH1D *hSignal_CORR[NUMTUBETYPES][NUMDEPTHS][2];
	TH1D *hSignal_OV2[NUMTUBETYPES][NUMDEPTHS][2];
	TH1D *hADC2GeV_OV2_UNCORR[NUMTUBETYPES][NUMDEPTHS][2];
	TH1D *hADC2GeV_OV2_CORR[NUMTUBETYPES][NUMDEPTHS][2];
	TH1D *hSFSB[NUMTUBETYPES][NUMDEPTHS][2];
	TH1D *hSB23toS[NUMTUBETYPES][NUMDEPTHS][2];
	
	TH1D *hGainsRatio;
	TH1D *hRatio_NoOF2wOF;

//	TGraph *gSignalvsGain;

//	TProfile *
};


void analyze_AvgQ(int hfSide, string inFileName, string outFileName, 
		string outTxtFileName, 
		string mapFileName, string gcfFile, string eMapFileName,
		string gainMapFileName, string swapsFileName, int numTS)
{
	analyze(hfSide, inFileName, outFileName, outTxtFileName,
			mapFileName, gcfFile, eMapFileName, gainMapFileName, 
			swapsFileName, numTS);
}

double getSwapRatio(TH1D *h, int iTS);
int initOut(TOut&, int iTS, string rootFileName, string txtFileName);

void analyze(int hfSide, string inFileName, string outFileName, string outTxtFileName, 
		string mapFileName, string gcfFile, 
		string eMapFileName, string gainMapFileName,
		string swapsFileName, int numTS)
{
	
	//	Open up the In File
	//
	TFile *in = new TFile(inFileName.c_str());

	//
	//	Get Electronics and Gain Maps + Define Output
	//
	TEMap eMap;
	getMap(hfSide, eMapFileName, eMap);
	TGainMap gainMap;
	getGainMap(hfSide, gainMapFileName, gainMap, eMap);	
	TSwaps swaps;
	readSwaps(swapsFileName, swaps, 0);
	TOut vOut;
	initOut(vOut, numTS, outFileName, outTxtFileName);

	bool ooo = true;	
	char someName[200];
	char profName[200];

	if (numTS==1)
		cut_toUse = PEDCUT_1TS;
	else
		cut_toUse = PEDCUT_2TS;
	cout << "CUT=" << cut_toUse << endl;;

	//
	//	Read in GCF and Tubes Map
	//
	TubesMap map;
	readTubesMap(hfSide, mapFileName, map, gcfFile, 0);
	int checkSkips = 0;

	vector<string> vStr;
	//
	//	Main Loop for analysis
	//
	int counter_Gain = 0;
	char histName[200];
	bool checkSwap = false;

	//	Sorry for no indentation here 
	for (int itube=0; itube<NUMTUBETYPES; itube++)
	for (int iphi=0; iphi<NUMPHIS; iphi++)
		for (int ieta=0; ieta<NUMETAS; ieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				checkSwap = false;
				cout << "### Processing: " << 2*iphi+1 << "  " << ieta+29
					<< idepth+1 << endl;

				if (swaps.numSrcErrors>0)
				{
					int iswap = isSrcError(swaps, 2*iphi+1, ieta+29, itube, 
							numTS);
					if (iswap==1)
						continue;
				}

				//
				//	Get a Sig histo
				//
				sprintf(histName, "1D_F3B3_wSrc_%d_PHI%d_ETA%d_D%d",
						itube, 2*iphi+1, ieta+29, idepth+1);
				in->cd("SRC/1D_F3B3");
				TH1D *hRAW_S = (TH1D*)gDirectory->Get(histName);
				if (hRAW_S->GetEntries()<100)
					continue;

				//
				//	Get a Sig_F histo
				//
				sprintf(histName, "1D_F_wSrc_%d_PHI%d_ETA%d_D%d",
						itube, 2*iphi+1, ieta+29, idepth+1);
				in->cd("SRC/1D_F");
				TH1D *hRAW_S_F = (TH1D*)gDirectory->Get(histName);
				if (hRAW_S_F->GetEntries()<100)
					continue;

				//
				//	Get a Sig_B histo
				//
				sprintf(histName, "1D_B_wSrc_%d_PHI%d_ETA%d_D%d",
						itube, 2*iphi+1, ieta+29, idepth+1);
				in->cd("SRC/1D_B");
				TH1D *hRAW_S_B = (TH1D*)gDirectory->Get(histName);
				if (hRAW_S_B->GetEntries()<100)
					continue;

				//
				//	Get a BG histo
				//	
				in->cd("NOSRC/1D");
				sprintf(histName, "1D_NoSrc_PHI%d_ETA%d_D%d",
						2*iphi+1, ieta+29, idepth+1);
				TH1D *hRAW_B = (TH1D*)gDirectory->Get(histName);
				if (hRAW_B->GetEntries()<100)
					continue;

				//
				//	Get a Sig_B23 histo
				//
				sprintf(histName, "1D_B23_wSrc_%d_PHI%d_ETA%d_D%d",
						itube, 2*iphi+1, ieta+29, idepth+1);
				in->cd("SRC/1D_B23");
				TH1D *hRAW_S_B23 = (TH1D*)gDirectory->Get(histName);
				if (hRAW_S_B23->GetEntries()<100)
					continue;

				//
				//	Get the GCF
				//
				double gcf = 1;
				if (idepth==0 &&map.tubes[iphi][ieta][itube].groove=='H'
						|| idepth==1 && map.tubes[iphi][ieta][itube].groove
						=='E')
					gcf = map.tubes[iphi][ieta][itube].em;
				else if (idepth==0 && map.tubes[iphi][ieta][itube].groove
						=='E' || idepth==1 
						&& map.tubes[iphi][ieta][itube].groove=='H')
					gcf = map.tubes[iphi][ieta][itube].had;
				cout << gcf << endl;

				//
				//	Fit PED Region of S/B Histos,
				//	Extract the Means of Pedestals
				//	Compute the (Qs - Ps) - (Qb - Pb)
				//
				fit(hRAW_S);
				vOut.rootFile->cd("Histos/1D/SRC");
				hRAW_S->Write();
				fit(hRAW_B);
				vOut.rootFile->cd("Histos/1D/BKG");
				hRAW_B->Write();
				fit(hRAW_S_F);
				vOut.rootFile->cd("Histos/1D_F");
				hRAW_S_F->Write();
				fit(hRAW_S_B);
				vOut.rootFile->cd("Histos/1D_B");
				hRAW_S_B->Write();
				fit(hRAW_S_B23);
				vOut.rootFile->cd("Histos/1D_B23");
				hRAW_S_B23->Write();

				TF1 *fit_S = hRAW_S->GetFunction("myGaus");
				TF1 *fit_S_F = hRAW_S_F->GetFunction("myGaus");
				TF1 *fit_S_B = hRAW_S_B->GetFunction("myGaus");
				TF1 *fit_S_B23 = hRAW_S_B23->GetFunction("myGaus");
				TF1 *fit_B = hRAW_B->GetFunction("myGaus");
				Double_t qie_mean_S = fit_S->GetParameter(1);
				Double_t qie_mean_B = fit_B->GetParameter(1);
				Double_t qie_sigma_S = fit_S->GetParameter(2);
				Double_t qie_sigma_B = fit_B->GetParameter(2);

				//
				//	Fill the 2D maps of PED Mean/Sigmas
				//
				vOut.hPedMean_S[itube][idepth]->Fill(2*iphi+1, 
						ieta+29, qie_mean_S);
				vOut.hPedMean_B[itube][idepth]->Fill(2*iphi+1, 
						ieta+29, qie_mean_B);
				vOut.hPedSigma_S[itube][idepth]->Fill(2*iphi+1, 
						ieta+29, qie_sigma_S);
				vOut.hPedSigma_B[itube][idepth]->Fill(2*iphi+1, 
						ieta+29, qie_sigma_B);

				double ratio = getSwapRatio(hRAW_S, numTS);
				vOut.hSwapRatio[itube][idepth]->Fill(2*iphi+1, ieta+29, ratio);

				Double_t qie_mean_S_F = fit_S_F->GetParameter(1);
				Double_t qie_mean_S_B = fit_S_B->GetParameter(1);
				Double_t qie_mean_S_B23 = fit_S_B23->GetParameter(1);

				Double_t totQSum_S_wOF = 0;
				Double_t totQSum_B_wOF = 0;
				Double_t totNSum_S_wOF = 0;
				Double_t totNSum_B_wOF = 0;
				Double_t totQSum_S_F_wOF = 0;
				Double_t totNSum_S_F_wOF = 0;
				Double_t totQSum_S_B_wOF = 0;
				Double_t totNSum_S_B_wOF = 0;
				Double_t totQSum_S_B23_wOF = 0;
				Double_t totNSum_S_B23_wOF = 0;

				Double_t totQSum_S_NoOF = 0;
				Double_t totNSum_S_NoOF = 0;
				Double_t totQSum_B_NoOF = 0;
				Double_t totNSum_B_NoOF = 0;
				Double_t totQSum_S_F_NoOF = 0;
				Double_t totNSum_S_F_NoOF = 0;
				Double_t totQSum_S_B_NoOF = 0;
				Double_t totNSum_S_B_NoOF = 0;
				Double_t totQSum_S_B23_NoOF = 0;
				Double_t totNSum_S_B23_NoOF = 0;


				for (int iBin=0; iBin<NUMBINS; iBin++)
				{
					totNSum_S_wOF += hRAW_S->GetBinContent(iBin+1)*
						hRAW_S->GetBinWidth(iBin+1);
					totNSum_B_wOF += hRAW_B->GetBinContent(iBin+1)*
						hRAW_B->GetBinWidth(iBin+1);
					totNSum_S_F_wOF += hRAW_S_F->GetBinContent(iBin+1)*
						hRAW_S_F->GetBinWidth(iBin+1);
					totNSum_S_B_wOF += hRAW_S_B->GetBinContent(iBin+1)*
						hRAW_S_B->GetBinWidth(iBin+1);
					totNSum_S_B23_wOF += hRAW_S_B23->GetBinContent(iBin+1)*
						hRAW_S_B23->GetBinWidth(iBin+1);

					totQSum_S_wOF += hRAW_S->GetBinCenter(iBin+1)*
						hRAW_S->GetBinContent(iBin+1)*
						hRAW_S->GetBinWidth(iBin+1);
					totQSum_B_wOF += hRAW_B->GetBinCenter(iBin+1)*
						hRAW_B->GetBinContent(iBin+1)*
						hRAW_B->GetBinWidth(iBin+1);
					totQSum_S_F_wOF +=hRAW_S_F->GetBinCenter(iBin+1)*
						hRAW_S_F->GetBinContent(iBin+1)*
						hRAW_S_F->GetBinWidth(iBin+1);
					totQSum_S_B_wOF +=hRAW_S_B->GetBinCenter(iBin+1)*
						hRAW_S_B->GetBinContent(iBin+1)*
						hRAW_S_B->GetBinWidth(iBin+1);
					totQSum_S_B23_wOF += hRAW_S_B23->GetBinCenter(iBin+1)*
						hRAW_S_B23->GetBinContent(iBin+1)*
						hRAW_S_B23->GetBinWidth(iBin+1);

					if (iBin<31)
					{
						totNSum_S_NoOF += hRAW_S->GetBinContent(iBin+1)*
							hRAW_S->GetBinWidth(iBin+1);
						totNSum_B_NoOF += hRAW_B->GetBinContent(iBin+1)*
							hRAW_B->GetBinWidth(iBin+1);
						totNSum_S_F_NoOF += hRAW_S_F->GetBinContent(iBin+1)*
							hRAW_S_F->GetBinWidth(iBin+1);
						totNSum_S_B_NoOF += hRAW_S_B->GetBinContent(iBin+1)*
							hRAW_S_B->GetBinWidth(iBin+1);
						totNSum_S_B23_NoOF += hRAW_S_B23->GetBinContent(iBin+1)*
							hRAW_S_B23->GetBinWidth(iBin+1);

						totQSum_S_NoOF += hRAW_S->GetBinCenter(iBin+1)*
							hRAW_S->GetBinContent(iBin+1)*
							hRAW_S->GetBinWidth(iBin+1);
						totQSum_B_NoOF += hRAW_B->GetBinCenter(iBin+1)*
							hRAW_B->GetBinContent(iBin+1)*
							hRAW_B->GetBinWidth(iBin+1);
						totQSum_S_F_NoOF +=hRAW_S_F->GetBinCenter(iBin+1)*
							hRAW_S_F->GetBinContent(iBin+1)*
							hRAW_S_F->GetBinWidth(iBin+1);
						totQSum_S_B_NoOF +=hRAW_S_B->GetBinCenter(iBin+1)*
							hRAW_S_B->GetBinContent(iBin+1)*
							hRAW_S_B->GetBinWidth(iBin+1);
						totQSum_S_B23_NoOF += hRAW_S_B23->GetBinCenter(iBin+1)*
							hRAW_S_B23->GetBinContent(iBin+1)*
							hRAW_S_B23->GetBinWidth(iBin+1);
					}
				}

				Double_t signal_S_wOF = totQSum_S_wOF/totNSum_S_wOF - 
					qie_mean_S;
				Double_t signal_B_wOF = totQSum_B_wOF/totNSum_B_wOF - 
					qie_mean_B;
				Double_t signal_S_F_wOF = totQSum_S_F_wOF/totNSum_S_F_wOF - 
					qie_mean_S_F;
				Double_t signal_S_B_wOF = totQSum_S_B_wOF/totNSum_S_B_wOF - 
					qie_mean_S_B;
				Double_t signal_S_B23_wOF = totQSum_S_B23_wOF/totNSum_S_B23_wOF- 
					qie_mean_S_B23;

				Double_t signal_S_NoOF = totQSum_S_NoOF/totNSum_S_NoOF - 
					qie_mean_S;
				Double_t signal_B_NoOF = totQSum_B_NoOF/totNSum_B_NoOF - 
					qie_mean_B;
				Double_t signal_S_F_NoOF = totQSum_S_F_NoOF/totNSum_S_F_NoOF - 
					qie_mean_S_F;
				Double_t signal_S_B_NoOF = totQSum_S_B_NoOF/totNSum_S_B_NoOF - 
					qie_mean_S_B;
				Double_t signal_S_B23_NoOF =
					totQSum_S_B23_NoOF/totNSum_S_B23_NoOF - qie_mean_S_B23;

				Double_t signal_wOF_UNCORR = signal_S_wOF - signal_B_wOF;
				Double_t signal_F_wOF_UNCORR = signal_S_F_wOF - signal_B_wOF;
				Double_t signal_B_wOF_UNCORR = signal_S_B_wOF - signal_B_wOF;
				Double_t signal_B23_wOF_UNCORR = signal_S_B23_wOF - signal_B_wOF;

				Double_t signal_NoOF_UNCORR = signal_S_NoOF - signal_B_NoOF;
				Double_t signal_F_NoOF_UNCORR = signal_S_F_NoOF - signal_B_NoOF;
				Double_t signal_B_NoOF_UNCORR = signal_S_B_NoOF - signal_B_NoOF;
				Double_t signal_B23_NoOF_UNCORR = signal_S_B23_NoOF - 
					signal_B_NoOF;
				
				Double_t signal_wOF_CORR = signal_wOF_UNCORR/gcf;
				Double_t signal_F_wOF_CORR = signal_F_wOF_UNCORR/gcf;
				Double_t signal_B_wOF_CORR = signal_B_wOF_UNCORR/gcf;
				Double_t signal_B23_wOF_CORR = signal_B23_wOF_UNCORR/gcf;
				
				Double_t signal_NoOF_CORR = signal_NoOF_UNCORR/gcf;
				Double_t signal_F_NoOF_CORR = signal_F_NoOF_UNCORR/gcf;
				Double_t signal_B_NoOF_CORR = signal_B_NoOF_UNCORR/gcf;
				Double_t signal_B23_NoOF_CORR = signal_B23_NoOF_UNCORR/gcf;

				Double_t ratioSFSB_wOF = signal_F_wOF_CORR/signal_B_wOF_CORR;
				Double_t ratioB23toS_wOF = signal_B23_wOF_CORR/signal_wOF_CORR;
				
				Double_t ratioSFSB_NoOF = signal_F_NoOF_CORR/signal_B_NoOF_CORR;
				Double_t ratioB23toS_NoOF = signal_B23_NoOF_CORR/
					signal_NoOF_CORR;

				Double_t ratio_NoOF2wOF = signal_NoOF_CORR/signal_wOF_CORR;
				vOut.hRatio_NoOF2wOF->Fill(ratio_NoOF2wOF);

				//
				//	Fill the Histos, Profiles, Graphs...
				//
				vOut.pSignalvsETA_UNCORR[itube][idepth][0]->Fill(ieta+29, 
						signal_wOF_UNCORR);
				vOut.pSignalvsETA_CORR[itube][idepth][0]->Fill(ieta+29, 
						signal_wOF_CORR);
				vOut.pSFSBvsETA[itube][idepth][0]->Fill(ieta+29, 
						ratioSFSB_wOF);
				vOut.pSB23toSvsETA[itube][idepth][0]->Fill(ieta+29,
						ratioB23toS_wOF);

				vOut.pSignalvsETA_UNCORR[itube][idepth][1]->Fill(ieta+29, 
						signal_NoOF_UNCORR);
				vOut.pSignalvsETA_CORR[itube][idepth][1]->Fill(ieta+29, 
						signal_NoOF_CORR);
				vOut.pSFSBvsETA[itube][idepth][1]->Fill(ieta+29, 
						ratioSFSB_NoOF);
				vOut.pSB23toSvsETA[itube][idepth][1]->Fill(ieta+29,
						ratioB23toS_NoOF);

				vOut.hSFSB[itube][idepth][0]->Fill(ratioSFSB_wOF);
				vOut.hSB23toS[itube][idepth][0]->Fill(ratioB23toS_wOF);
				
				vOut.hSFSB[itube][idepth][1]->Fill(ratioSFSB_NoOF);
				vOut.hSB23toS[itube][idepth][1]->Fill(ratioB23toS_NoOF);

				double gain_OV1 = gainMap[iphi][ieta][idepth].gain_OV1;
				double gain_OV1P100 = gainMap[iphi][ieta][idepth].gain_OV1P100;
				double gain_OV2 = gainMap[iphi][ieta][idepth].gain_OV2;
				double gainToUse, timeFactor;

				//
				//	Compute and Fill the Gain Ratio.
				//	OV2/OV1
				//
				double gainRatio(0);
				if (itube==0)
				{
					gainRatio = gain_OV2/gain_OV1;
					vOut.hGainsRatio->Fill(gainRatio);
				}

				//	set TS-dependent parameters
				if (numTS==1)
				{
					timeFactor = 1;
					gainToUse = gain_OV1;
				}
				else if (numTS==2)
				{
					timeFactor = 0.5;
					gainToUse = gain_OV1P100;
				}

				//
				//	Compute:
				//
				//
				Double_t signal_OV2_wOF = timeFactor*signal_wOF_CORR*
					gain_OV2/gainToUse;
				Double_t adc2GeV_OV2_wOF = GEVPER25NS[idepth]/signal_OV2_wOF*
					SOURCEACTIVITYCORRECTION[hfSide];

				Double_t signal_OV2_NoOF = timeFactor*signal_NoOF_CORR*
					gain_OV2/gainToUse;
				Double_t adc2GeV_OV2_NoOF = GEVPER25NS[idepth]/signal_OV2_NoOF*
					SOURCEACTIVITYCORRECTION[hfSide];

				vOut.hSignal_UNCORR[itube][idepth][0]->Fill(signal_wOF_UNCORR);
				vOut.hSignal_UNCORR[itube][idepth][1]->Fill(signal_NoOF_UNCORR);
				vOut.hSignal_CORR[itube][idepth][0]->Fill(signal_wOF_CORR);
				vOut.hSignal_CORR[itube][idepth][1]->Fill(signal_NoOF_CORR);
				vOut.hSignal_OV2[itube][idepth][0]->Fill(signal_OV2_wOF);
				vOut.hSignal_OV2[itube][idepth][1]->Fill(signal_OV2_NoOF);
				vOut.hADC2GeV_OV2_CORR[itube][idepth][0]->Fill(adc2GeV_OV2_wOF);
				vOut.hADC2GeV_OV2_CORR[itube][idepth][1]->Fill(adc2GeV_OV2_NoOF);

				//
				//	Print the Results to txt File
				//
				vOut.txtFile
					<< itube << " " << 2*iphi+1 << " " << ieta+29 << " " 
					<< idepth << " " << gcf << " " 
					<< signal_S_NoOF << " " << signal_B_NoOF << " "
					<< qie_mean_S << " " << qie_mean_B << " "
					<< signal_wOF_CORR << " "
					<< signal_NoOF_CORR << " " << signal_OV2_wOF << " " 
					<< signal_OV2_NoOF << " " << adc2GeV_OV2_wOF << " "
					<< adc2GeV_OV2_NoOF << " " << gain_OV2 << " "
					<< gain_OV1 << " " << gain_OV1P100
					<< endl;
			}

	//
	//	Save and Close
	//
	for (itube=0; itube<NUMTUBETYPES; itube++)
	for (int i=1; i<=2; i++)
	{
		vOut.pSignalvsETA_UNCORR[itube][i-1][0]->GetXaxis()->SetTitle("ieta");
		vOut.pSignalvsETA_UNCORR[itube][i-1][1]->GetXaxis()->SetTitle("ieta");
		vOut.pSignalvsETA_UNCORR[itube][i-1][0]->GetYaxis()->SetTitle("Signal");
		vOut.pSignalvsETA_UNCORR[itube][i-1][1]->GetYaxis()->SetTitle("Signal");
		vOut.pSignalvsETA_UNCORR[itube][i-1][0]->GetYaxis()->SetTitleOffset(1.2);
		vOut.pSignalvsETA_UNCORR[itube][i-1][1]->GetYaxis()->SetTitleOffset(1.2);
		vOut.pSignalvsETA_CORR[itube][i-1][0]->GetXaxis()->SetTitle("ieta");
		vOut.pSignalvsETA_CORR[itube][i-1][1]->GetXaxis()->SetTitle("ieta");
		vOut.pSignalvsETA_CORR[itube][i-1][1]->GetYaxis()
			->SetTitle("Signal(Corr. by GCF)");
		vOut.pSignalvsETA_CORR[itube][i-1][1]->GetYaxis()
			->SetTitle("Signal(Corr. by GCF)");
		vOut.pSignalvsETA_CORR[itube][i-1][0]->GetYaxis()->SetTitleOffset(1.2);
		vOut.pSignalvsETA_CORR[itube][i-1][1]->GetYaxis()->SetTitleOffset(1.2);
		if (numTS==1)
		{
			vOut.pSignalvsETA_UNCORR[itube][i-1][0]->GetYaxis()
				->SetRangeUser(0.001, 0.01);
			vOut.pSignalvsETA_UNCORR[itube][i-1][1]->GetYaxis()
				->SetRangeUser(0.001, 0.01);
			vOut.pSignalvsETA_CORR[itube][i-1][0]->GetYaxis()
				->SetRangeUser(0.001, 0.01);
			vOut.pSignalvsETA_CORR[itube][i-1][1]->GetYaxis()
				->SetRangeUser(0.001, 0.01);
		}
		else if (numTS==2)
		{
			vOut.pSignalvsETA_UNCORR[itube][i-1][0]->GetYaxis()
				->SetRangeUser(0.01, 0.1);
			vOut.pSignalvsETA_UNCORR[itube][i-1][1]->GetYaxis()
				->SetRangeUser(0.01, 0.1);
			vOut.pSignalvsETA_CORR[itube][i-1][0]->GetYaxis()
				->SetRangeUser(0.01, 0.1);
			vOut.pSignalvsETA_CORR[itube][i-1][1]->GetYaxis()
				->SetRangeUser(0.01, 0.1);
		}

		vOut.pSFSBvsETA[itube][i-1][0]->GetXaxis()->SetTitle("ieta");
		vOut.pSFSBvsETA[itube][i-1][1]->GetXaxis()->SetTitle("ieta");
		vOut.pSFSBvsETA[itube][i-1][0]->GetYaxis()->SetTitle("SF/SB");
		vOut.pSFSBvsETA[itube][i-1][1]->GetYaxis()->SetTitle("SF/SB");
		vOut.pSFSBvsETA[itube][i-1][0]->GetYaxis()->SetTitleOffset(1.2);
		vOut.pSFSBvsETA[itube][i-1][1]->GetYaxis()->SetTitleOffset(1.2);
		vOut.pSFSBvsETA[itube][i-1][0]->GetYaxis()->SetRangeUser(0.8, 1.2);
		vOut.pSFSBvsETA[itube][i-1][1]->GetYaxis()->SetRangeUser(0.8, 1.2);

		vOut.pSB23toSvsETA[itube][i-1][0]->GetXaxis()->SetTitle("ieta");
		vOut.pSB23toSvsETA[itube][i-1][1]->GetXaxis()->SetTitle("ieta");
		vOut.pSB23toSvsETA[itube][i-1][0]->GetYaxis()->SetTitle("SB2_3/SAll");
		vOut.pSB23toSvsETA[itube][i-1][1]->GetYaxis()->SetTitle("SB2_3/SAll");
		vOut.pSB23toSvsETA[itube][i-1][0]->GetYaxis()->SetTitleOffset(1.2);
		vOut.pSB23toSvsETA[itube][i-1][1]->GetYaxis()->SetTitleOffset(1.2);
		vOut.pSB23toSvsETA[itube][i-1][0]->GetYaxis()->SetRangeUser(0.8, 1.2);
		vOut.pSB23toSvsETA[itube][i-1][1]->GetYaxis()->SetRangeUser(0.8, 1.2);

		vOut.hSFSB[itube][i-1][0]->GetXaxis()->SetTitle("SF/SB");
		vOut.hSFSB[itube][i-1][1]->GetXaxis()->SetTitle("SF/SB");
		vOut.hSB23toS[itube][i-1][0]->GetXaxis()->SetTitle("SB2_3/SAll");
		vOut.hSB23toS[itube][i-1][1]->GetXaxis()->SetTitle("SB2_3/SAll");

		vOut.hSignal_UNCORR[itube][i-1][0]->GetXaxis()
			->SetTitle("Signal_UNCORR");
		vOut.hSignal_UNCORR[itube][i-1][1]->GetXaxis()
			->SetTitle("Signal_UNCORR");
		vOut.hSignal_CORR[itube][i-1][0]->GetXaxis()
			->SetTitle("Signal_CORR");
		vOut.hSignal_CORR[itube][i-1][1]->GetXaxis()
			->SetTitle("Signal_CORR");
		vOut.hSignal_OV2[itube][i-1][0]->GetXaxis()
			->SetTitle("Signal_OV2");
		vOut.hSignal_OV2[itube][i-1][1]->GetXaxis()
			->SetTitle("Signal_OV2");
		vOut.hADC2GeV_OV2_CORR[itube][i-1][0]->GetXaxis()
			->SetTitle("ADC2GeV @OV2");
		vOut.hADC2GeV_OV2_CORR[itube][i-1][1]->GetXaxis()
			->SetTitle("ADC2GeV @OV2");
	}

	vOut.hGainsRatio->GetXaxis()->SetTitle("Gain Ratios OV2/OV1");
	vOut.rootFile->Write();

	return;
}

//
//	Init all the swapped Channels
//
int initSwaps()
{
	Coordinates coord;

	coord.iphi = 57;
	coord.ieta = 29;
	coord.idepth=0;
	swappedChs[0] = coord;
	coord.iphi = 57;
	coord.ieta = 30;
	coord.idepth=0;
	swappedChs[1] = coord;
	coord.iphi = 57;
	coord.ieta = 31;
	coord.idepth=0;
	swappedChs[2] = coord;
	coord.iphi = 57;
	coord.ieta = 32;
	coord.idepth=0;
	swappedChs[3] = coord;
	coord.iphi = 57;
	coord.ieta = 33;
	coord.idepth=0;
	swappedChs[4] = coord;
	coord.iphi = 57;
	coord.ieta = 34;
	coord.idepth=0;
	swappedChs[5] = coord;


	coord.iphi = 61;
	coord.ieta = 35;
	coord.idepth=0;
	swappedChs[6] = coord;
	coord.iphi = 61;
	coord.ieta = 36;
	coord.idepth=0;
	swappedChs[7] = coord;
	coord.iphi = 61;
	coord.ieta = 37;
	coord.idepth=0;
	swappedChs[8] = coord;
	coord.iphi = 61;
	coord.ieta = 38;
	coord.idepth=0;
	swappedChs[9] = coord;
	coord.iphi = 61;
	coord.ieta = 39;
	coord.idepth=0;
	swappedChs[10] = coord;
	coord.iphi = 59;
	coord.ieta = 40;
	coord.idepth=0;
	swappedChs[11] = coord;

	coord.iphi = 53;
	coord.ieta = 29;
	coord.idepth=1;
	swappedChs[12] = coord;
	coord.iphi = 53;
	coord.ieta = 30;
	coord.idepth=1;
	swappedChs[13] = coord;
	coord.iphi = 53;
	coord.ieta = 31;
	coord.idepth=1;
	swappedChs[14] = coord;
	coord.iphi = 53;
	coord.ieta = 32;
	coord.idepth=1;
	swappedChs[15] = coord;
	coord.iphi = 53;
	coord.ieta = 33;
	coord.idepth=1;
	swappedChs[16] = coord;
	coord.iphi = 53;
	coord.ieta = 34;
	coord.idepth=1;
	swappedChs[17] = coord;

	coord.iphi = 57;
	coord.ieta = 35;
	coord.idepth=1;
	swappedChs[18] = coord;
	coord.iphi = 57;
	coord.ieta = 36;
	coord.idepth=1;
	swappedChs[19] = coord;
	coord.iphi = 57;
	coord.ieta = 37;
	coord.idepth=1;
	swappedChs[20] = coord;
	coord.iphi = 57;
	coord.ieta = 38;
	coord.idepth=1;
	swappedChs[21] = coord;
	coord.iphi = 57;
	coord.ieta = 39;
	coord.idepth=1;
	swappedChs[22] = coord;
	coord.iphi = 55;
	coord.ieta = 40;
	coord.idepth=1;
	swappedChs[23] = coord;

	coord.iphi = 25;
	coord.ieta = 35;
	coord.idepth=1;
	swappedChs[24] = coord;
	coord.iphi = 25;
	coord.ieta = 36;
	coord.idepth=1;
	swappedChs[25] = coord;
	coord.iphi = 25;
	coord.ieta = 37;
	coord.idepth=1;
	swappedChs[26] = coord;
	coord.iphi = 25;
	coord.ieta = 38;
	coord.idepth=1;
	swappedChs[27] = coord;
	coord.iphi = 25;
	coord.ieta = 39;
	coord.idepth=1;
	swappedChs[28] = coord;
	coord.iphi = 23;
	coord.ieta = 40;
	coord.idepth=1;
	swappedChs[29] = coord;

	coord.iphi = 23;
	coord.ieta = 34;
	coord.idepth=0;
	swappedChs[30] = coord;
	coord.iphi = 23;
	coord.ieta = 34;
	coord.idepth=1;
	swappedChs[31] = coord;
	coord.iphi = 23;
	coord.ieta = 35;
	coord.idepth=0;
	swappedChs[32] = coord;
	coord.iphi = 23;
	coord.ieta = 35;
	coord.idepth=1;
	swappedChs[33] = coord;

	/*
	coord.iphi = 27;
	coord.ieta = 12;
	coord.idepth = 0;
	swappedChs[0] = coord;

	coord.iphi = 28;
	coord.ieta = 10;
	coord.idepth = 0;
	swappedChs[1] = coord;

	coord.iphi = 28;
	coord.ieta = 6;
	coord.idepth = 0;
	swappedChs[2] = coord;

	coord.iphi = 28;
	coord.ieta = 7;
	coord.idepth = 0;
	swappedChs[3] = coord;

	coord.iphi = 28;
	coord.ieta = 8;
	coord.idepth = 0;
	swappedChs[4] = coord;

	coord.iphi = 28;
	coord.ieta = 9;
	coord.idepth = 0;
	swappedChs[5] = coord;

	coord.iphi = 29;
	coord.ieta = 6;
	coord.idepth = 0;
	swappedChs[6] = coord;

	coord.iphi = 29;
	coord.ieta = 7;
	coord.idepth = 0;
	swappedChs[7] = coord;

	coord.iphi = 29;
	coord.ieta = 8;
	coord.idepth = 0;
	swappedChs[8] = coord;

	coord.iphi = 29;
	coord.ieta = 9;
	coord.idepth = 0;
	swappedChs[9] = coord;

	coord.iphi = 29;
	coord.ieta = 10;
	coord.idepth = 0;
	swappedChs[10] = coord;

	coord.iphi = 29;
	coord.ieta = 11;
	coord.idepth = 0;
	swappedChs[11] = coord;

	coord.iphi = 5;
	coord.ieta = 6;
	coord.idepth = 1;
	swappedChs[12] = coord;

	coord.iphi = 5;
	coord.ieta = 7;
	coord.idepth = 1;
	swappedChs[13] = coord;
	
	coord.iphi = 5;
	coord.ieta = 8;
	coord.idepth = 1;
	swappedChs[14] = coord;

	coord.iphi = 5;
	coord.ieta = 9;
	coord.idepth = 1;
	swappedChs[15] = coord;

	coord.iphi = 5;
	coord.ieta = 10;
	coord.idepth = 1;
	swappedChs[16] = coord;

	coord.iphi = 5;
	coord.ieta = 11;
	coord.idepth = 1;
	swappedChs[17] = coord;

	coord.iphi = 33;
	coord.ieta = 6;
	coord.idepth = 1;
	swappedChs[18] = coord;

	coord.iphi = 33;
	coord.ieta = 7;
	coord.idepth = 1;
	swappedChs[19] = coord;

	coord.iphi = 33;
	coord.ieta = 8;
	coord.idepth = 1;
	swappedChs[20] = coord;

	coord.iphi = 33;
	coord.ieta = 9;
	coord.idepth = 1;
	swappedChs[21] = coord;

	coord.iphi = 33;
	coord.ieta = 10;
	coord.idepth = 1;
	swappedChs[22] = coord;

	coord.iphi = 33;
	coord.ieta = 11;
	coord.idepth = 1;
	swappedChs[23] = coord;

	coord.iphi = 0;
	coord.ieta = 7;
	coord.idepth = 0;
	swappedChs[24] = coord;

	coord.iphi = 0;
	coord.ieta = 8;
	coord.idepth = 0;
	swappedChs[25] = coord;

	coord.iphi = 0;
	coord.ieta = 7;
	coord.idepth = 1;
	swappedChs[26] = coord;

	coord.iphi = 0;
	coord.ieta = 8;
	coord.idepth = 1;
	swappedChs[27] = coord;

	coord.iphi = 26;
	coord.ieta = 0;
	coord.idepth = 0;
	swappedChs[28] = coord;

	coord.iphi = 26;
	coord.ieta = 0;
	coord.idepth = 1;
	swappedChs[29] = coord;
	coord.iphi = 26;
	coord.ieta = 1;
	coord.idepth = 0;
	swappedChs[30] = coord;

	coord.iphi = 26;
	coord.ieta = 1;
	coord.idepth = 1;
	swappedChs[31] = coord;
*/
	
}

int initOut(TOut &out, int iTS, string rootFileName, string txtFileName)
{
	out.rootFile = new TFile(rootFileName.c_str(), "recreate");
	out.txtFile.open(txtFileName.c_str());

	out.rootFile->mkdir("Profiles");
	out.rootFile->mkdir("Histos");
	out.rootFile->mkdir("Graphs");

	out.rootFile->cd("Histos");
	gDirectory->mkdir("2D");
	gDirectory->mkdir("1D");
	gDirectory->mkdir("1D_F");
	gDirectory->mkdir("1D_B");
	gDirectory->mkdir("1D_B23");
	gDirectory->cd("1D");
	gDirectory->mkdir("SRC");
	gDirectory->mkdir("BKG");

	char histName[200];
	out.hGainsRatio = new TH1D("GainsRatio", "GainsRatio", 100, 0, 1);
	out.hRatio_NoOF2wOF = new TH1D("Ratio_NoOF2wOF", 
			"Ratio NoOF2wOF", 100, 0, 2);
	for (int itube=0; itube<NUMTUBETYPES; itube++)
	{
		for (int id=0; id<NUMDEPTHS; id++)
		{
			out.rootFile->cd("Profiles");

			for (int iOF=0; iOF<2; iOF++)
			{
				sprintf(histName, "SignalvsETA_UNCORR_TT%d_D%d_%s", itube, id+1,
						(iOF==0) ? "withOF" : "withoutOF");
				out.pSignalvsETA_UNCORR[itube][id][iOF] = new TProfile(histName,
						histName, 13, 29, 42);
				sprintf(histName, "SignalvsETA_CORR_TT%d_D%d_%s", itube, id+1,
						(iOF==0) ? "withOF" : "withoutOF");
				out.pSignalvsETA_CORR[itube][id][iOF] = new TProfile(histName, 
						histName, 13, 29, 42);
				sprintf(histName, "SFSBvsETA_TT%d_D%d_%s", itube, id+1,
						iOF==0 ? "withOF" : "withoutOF");
				out.pSFSBvsETA[itube][id][iOF] = new TProfile(histName, histName,
						13, 29, 42);
				sprintf(histName, "SB23toSvsETA_TT%d_D%d_%s", itube, id+1,
						iOF==0 ? "withOF" : "withoutOF");
				out.pSB23toSvsETA[itube][id][iOF] = new TProfile(histName, 
						histName, 13, 29, 42);
			}

			out.rootFile->cd("Histos/2D");

			sprintf(histName, "PedMean_S_TT%d_D%d", itube, id+1);
			out.hPedMean_S[itube][id] = new TH2D(histName, histName,
					72, 0, 72, 13, 29, 42);
			sprintf(histName, "PedMean_B_TT%d_D%d", itube, id+1);
			out.hPedMean_B[itube][id] = new TH2D(histName, histName,
					72, 0, 72, 13, 29, 42);
			sprintf(histName, "PedSigma_S_TT%d_D%d", itube, id+1);
			out.hPedSigma_S[itube][id] = new TH2D(histName, histName,
					72, 0, 72, 13, 29, 42);
			sprintf(histName, "PedSigma_B_TT%d_D%d", itube, id+1);
			out.hPedSigma_B[itube][id] = new TH2D(histName, histName,
					72, 0, 72, 13, 29, 42);
			sprintf(histName, "SwapRatio_TT%d_D%d", itube, id+1);
			out.hSwapRatio[itube][id] = new TH2D(histName, histName,
					72, 0, 72, 13, 29, 42);
//			sprintf(histName, "ADC2GeV_OV2_TT%d_D%d", itube. id+1);
//			out.hADC2GeV_OV2[itube][id] = new TH2D(histName ,histName, 
//					72, 0, 72, 13, 29, 42);

			out.rootFile->cd("Histos/1D");

			for (int iOF=0; iOF<2; iOF++)
			{
				if (iTS==1)
				{
					sprintf(histName, "Signal_UNCORR_TT%d_D%d_%s", itube, id+1,
							iOF==0 ? "withOF" : "withoutOF");
					out.hSignal_UNCORR[itube][id][iOF] = new TH1D(histName, 
							histName, 100, 0, 0.02);
					sprintf(histName, "Signal_CORR_TT%d_D%d_%s", itube, id+1,
							iOF==0 ? "withOF" : "withoutOF");
					out.hSignal_CORR[itube][id][iOF] = new TH1D(histName, 
							histName, 100, 0, 0.02);
				}
				else 
				{
					sprintf(histName, "Signal_UNCORR_TT%d_D%d_%s", itube, id+1,
							iOF==0 ? "withOF" : "withoutOF");
					out.hSignal_UNCORR[itube][id][iOF] = new TH1D(histName, 
							histName, 1000, 0, 0.2);
					sprintf(histName, "Signal_CORR_TT%d_D%d_%s", itube, id+1,
							iOF==0 ? "withOF" : "withoutOF");
					out.hSignal_CORR[itube][id][iOF] = new TH1D(histName, 
							histName, 1000, 0, 0.2);
				}

				sprintf(histName, "Signal_OV2_TT%d_D%d_%s", itube, id+1,
						iOF==0 ? "withOF" : "withoutOF");
				out.hSignal_OV2[itube][id][iOF] = new TH1D(histName, histName,
						1000, 0, 0.02);
				sprintf(histName, "ADC2GeV_OV2_UNCORR_TT%d_D%d_%s", itube, id+1,
						iOF==0 ? "withOF" : "withoutOF");
				out.hADC2GeV_OV2_UNCORR[itube][id][iOF] = new TH1D(histName,
						histName, 1000, 0, 2);
				sprintf(histName, "ADC2GeV_OV2_CORR_TT%d_D%d_%s", itube, id+1,
						iOF==0 ? "withOF" : "withoutOF");
				out.hADC2GeV_OV2_CORR[itube][id][iOF] = new TH1D(histName, 
						histName, 1000, 0, 2);
				sprintf(histName, "SFSB_TT%d_D%d_%s", itube, id+1, 
						iOF==0 ? "withOF" : "withoutOF");
				out.hSFSB[itube][id][iOF] = new TH1D(histName, histName,
						100, 0, 2);
				sprintf(histName, "SB23toS_TT%d_D%d_%s", itube, id+1,
						iOF==0 ? "withOF" : "withoutOF");
				out.hSB23toS[itube][id][iOF] = new TH1D(histName, histName,
						100, 0, 2);
			}
		}
	}

	return 0;
}

//
//	Fit the Histo
//
void fit(TH1D *hist)
{
	Double_t qie_mean_raw = hist->GetMean();
	int qie_mean_bin = floor(qie_mean_raw)+1;
	Double_t qie_mean = hist->GetBinCenter(qie_mean_bin);
	Double_t qie_rms_raw = hist->GetRMS();
	Double_t min = ((qie_mean-4*qie_rms_raw)>0) ? (qie_mean-4*qie_rms_raw) :
		0;
	Double_t max = qie_mean+4*qie_rms_raw;
	TF1 *myGaus = new TF1("myGaus", "gaus", min, max);
	hist->SetLineColor(kBlack);
	myGaus->SetLineColor(kRed);
	hist->Fit("myGaus", "QR");
}

//
//	Get the Swap Ratio
//
double getSwapRatio(TH1D *h, int iTS)
{
	double ratio = 0;
	double qie_mean = h->GetFunction("myGaus")->GetParameter(1);
	int qie_mean_Bin = floor(qie_mean)+1;
	int qie_src_Bin = 23;

	double eventsAtMean = h->GetBinContent(qie_mean_Bin)*
		h->GetBinWidth(qie_mean_Bin);
	double eventsAtSrc = h->GetBinContent(qie_src_Bin)*
		h->GetBinWidth(qie_src_Bin);

	ratio = eventsAtSrc/eventsAtMean;

	return ratio;
}

//
//	Reads GCF and Tubes Map
//
int readTubesMap(int hfSide, string fileNameMap, TubesMap &map, string gcfFileName, 
		int verb)
{
	cout << "### Reading in Tubes " << endl;

	ifstream mapFile(fileNameMap.c_str());
	if (!mapFile)
	{
		cout << "### ERROR no Tubes Map File" << endl;
		return -1;
	}

	ifstream gcfFile(gcfFileName.c_str());
	double em[31];
	double had[31];
	double uncert[31];
	if (!gcfFile)
	{
		cout << "### ERROR No GCF File" << endl;
		return -1;
	}

	int counter = 0;
	while (gcfFile >> em[counter])
	{
		gcfFile >> had[counter] >> uncert[counter];
		counter++;
	}

	string tubeName, wedgeName, swaste;
	int eta, phi, tubeStart, tubeEnd, iwaste;
	char groove;

	counter = 0;
	while (mapFile >> tubeName)
	{
		if (hfSide==0)
			mapFile >> wedgeName >> iwaste >> swaste
				>> eta >> eta >> phi >> iwaste >> iwaste 
				>> tubeStart >> tubeEnd >> groove;
		else if (hfSide==1)
			mapFile >> wedgeName
				>> eta >> eta >> phi
				>> tubeStart >> tubeEnd >> groove;

		int iphi = (phi-1)/2;
		int ieta = abs(eta) - 29;
		int checkRBX, checkPhi, checkEta, checkTubeNum;
		char tubeType;

		if (hfSide==0)
			sscanf(tubeName.c_str(), "HFP%d_ETA%d_PHI%d_T%d%c%*[^]",
					&checkRBX, &checkEta, &checkPhi, &checkTubeNum,
					&tubeType);
		else if (hfSide==1)
			sscanf(tubeName.c_str(), "HFM%d_ETA%d_PHI%d_T%d%c%*[^]",
					&checkRBX, &checkEta, &checkPhi, &checkTubeNum,
					&tubeType);

		int iTubeType = 0;
		if (tubeType=='A' || tubeType=='_')
			iTubeType = 0;
		else if (tubeType=='B')
			iTubeType = 1;
		else 
			cout << "### ERROR: Wrong Tube Type" << endl;
		map.tubes[iphi][ieta][iTubeType].tubeStart = tubeStart;
		map.tubes[iphi][ieta][iTubeType].tubeEnd = tubeEnd;
		map.tubes[iphi][ieta][iTubeType].tubeName = tubeName;
		map.tubes[iphi][ieta][iTubeType].em = em[counter%31];
		map.tubes[iphi][ieta][iTubeType].had = had[counter%31];
		map.tubes[iphi][ieta][iTubeType].uncert = uncert[counter%31];
		map.tubes[iphi][ieta][iTubeType].groove = groove;

		counter++;
	}

	return 0;
}

int main(int arc, char **argv)
{
	int hfSide = atoi(argv[1]);
	string inFileName = argv[2];
	string outFileName = argv[3];
	string outTxtFileName = argv[4];
	string mapFileName = argv[5];
	string gcfFile = argv[6];
	string eMapFileName = argv[7];
	string gainMapFileName = argv[8];
	string swapsFileName = argv[9];
	int numTS = atoi(argv[10]);

	analyze(hfSide, inFileName, outFileName, outTxtFileName,
			mapFileName, gcfFile, eMapFileName, gainMapFileName,
			swapsFileName, numTS);

	return 0;
}

