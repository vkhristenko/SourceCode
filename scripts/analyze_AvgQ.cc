


#define NUMPHIS 36
#define NUMETAS 13
#define NUMTUBETYPES 2
#define NUMDEPTHS 2
#define NUMBINS 32
#define PEDCUT_1TS 10
#define PEDCUT_2TS 18

//
//	Adapt this for either HFP or HFM
//	#swaps in HFP is 34
//	#swaps in HFM is 32
//
//#define NUMSWAPS 32
#define NUMSWAPS 34

#include "map.cc"

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

void analyze_AvgQ(string inFileName, string outFileName, string outTxtFileName, 
		string mapFileName, string gcfFile, string eMapFileName,
		string gainMapFileName, int numTS)
{
	analyze(inFileName, outFileName, outTxtFileName,
			mapFileName, gcfFile, eMapFileName, gainMapFileName, numTS);
}

double getSwapRatio(TH1D *h, int iTS);

void analyze(string inFileName, string outFileName, string outTxtFileName, 
		string mapFileName, string gcfFile, 
		string eMapFileName, string gainMapFileName, int numTS)
{
	
	//	Open up the In File
	//
	TFile *in = new TFile(inFileName.c_str());

	//
	//	Get Electronics and Gain Maps
	//
	TEMap eMap;
	getMap(eMapFileName, eMap);
	TGainMap gainMap;
	getGainMap(gainMapFileName, gainMap, eMap);	
/*	for (int iphi=0; iphi<NUMPHIS; iphi++)
		for (int ieta=0; ieta<NUMETAS; ieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				cout << iphi << "  " << ieta << "  " << idepth 
					<< "  " << gainMap[iphi][ieta][idepth].gain_OV1
					<< "  " << gainMap[iphi][ieta][idepth].gain_OV2 << endl;
			}
*/			

	//
	//	Create the Output File, initialzie whatever
	//
	ofstream outTxt(outTxtFileName.c_str());
	TFile *out = new TFile(outFileName.c_str(), "recreate");
	out->mkdir("Profiles");
	out->mkdir("Histos");
	out->mkdir("Graphs");

	out->cd("Profiles");
	TProfile *pSPEvsETA_UNCORR[NUMDEPTHS];
	TProfile *pSPEvsETA_CORR[NUMDEPTHS];
	TProfile *pSPEvsPHI_Avg[NUMDEPTHS];
	TProfile *pSPEvsPHI[NUMETAS][NUMDEPTHS];
	TProfile *pSFSBvsPHI[NUMETAS][NUMDEPTHS];
	TProfile *pSFoverSB[NUMDEPTHS];
	TProfile *pSFoverSB_D12;
	TProfile *pSFSBvsPHI_Avg[NUMDEPTHS];
	pSPEvsETA_UNCORR[0] = new TProfile("SignalvsETA_D1_UNCORR", 
			"SignalvsETA_D1_UNCORR",
			13, 29, 42);
	pSPEvsETA_UNCORR[1] = new TProfile("SignalvsETA_D2_UNCORR", 
			"SignalvsETA_D2_UNCORR",
			13, 29, 42);
	pSPEvsETA_CORR[0] = new TProfile("SignalvsETA_D1_CORR", 
			"SignalvsETA_D1_CORR", 
			13, 29, 42);
	pSPEvsETA_CORR[1] = new TProfile("SignalvsETA_D2_CORR", 
			"SignalvsETA_D2_CORR", 
			13, 29, 42);
	pSPEvsPHI_Avg[0] = new TProfile("SignalvsPHI_D1_Avg", "SignalvsPHI_D1_Avg",
			72, 0, 72);
	pSPEvsPHI_Avg[1] = new TProfile("SignalvsPHI_D2_Avg", "SignalvsPHI_D2_Avg",
			72, 0, 72);
	pSFoverSB[0] = new TProfile("SFoverSB_D1", "SFoverSB_D1",
			13, 29, 42);
	pSFoverSB[1] = new TProfile("SFoverSB_D2", "SFoverSB_D2",
			13, 29, 42);
	pSFoverSB_D12 = new TProfile("SFoverSB_D12", "SFoverSB_D12", 
			13, 29, 42);
	pSFSBvsPHI_Avg[0] = new TProfile("SFoverSBvsPHI_D1_Avg", "SFoverSBvsPHI_D1_Avg",
			72, 0, 72);
	pSFSBvsPHI_Avg[1] = new TProfile("SFoverSBvsPHI_D2_Avg", "SFoverSBvsPHI_D2_Avg",
			72, 0, 72);

	gDirectory->mkdir("SPEvsPHI");
	gDirectory->mkdir("SFSBvsPHI");
	char profName[200];
	for (int ieta=0; ieta<NUMETAS; ieta++)
	{
		for (int idepth=0; idepth<NUMDEPTHS; idepth++)
		{
			out->cd("Profiles/SPEvsPHI");
			sprintf(profName, "SignalvsPHI_ETA%d_D%d", ieta+29, idepth+1);
			pSPEvsPHI[ieta][idepth] = new TProfile(profName, profName,
					72, 0, 72);

			out->cd("Profiles/SFSBvsPHI");
			sprintf(profName, "SFoverSBvsPHI_ETA%d_D%d", ieta+29, idepth+1);
			pSFSBvsPHI[ieta][idepth] = new TProfile(profName, profName,
					72, 0, 72);
		}
	}

	out->cd("Histos");
	gDirectory->mkdir("2D");
	gDirectory->mkdir("1D");
	gDirectory->mkdir("1D_F");
	gDirectory->mkdir("1D_B");
	gDirectory->cd("1D");
	gDirectory->mkdir("SIG");
	gDirectory->mkdir("BACK");
	out->cd("Histos/2D");
	TH2D *hPedMean_S[NUMDEPTHS];
	TH2D *hPedMean_B[NUMDEPTHS];
	TH2D *hPedSigma_S[NUMDEPTHS];
	TH2D *hPedSigma_B[NUMDEPTHS];
	TH2D *hSwapRatio[NUMDEPTHS];
	char someName[200];
	for (int i=0; i<2; i++)
	{
		sprintf(someName, "PedMean_S_D%d", i+1);
		hPedMean_S[i] = new TH2D(someName, someName, 72, 0, 72, 13, 29, 42);
		
		sprintf(someName, "PedMean_B_D%d", i+1);
		hPedMean_B[i] = new TH2D(someName, someName, 72, 0, 72, 13, 29, 42);

		sprintf(someName, "PedSigma_S_D%d", i+1);
		hPedSigma_S[i] = new TH2D(someName, someName, 72, 0, 72, 13, 29, 42);

		sprintf(someName, "PedSigma_B_D%d", i+1);
		hPedSigma_B[i] = new TH2D(someName, someName, 72, 0, 72, 13, 29, 42);

		sprintf(someName, "SwapRatio_D%d", i+1);
		hSwapRatio[i] = new TH2D(someName, someName, 72, 0, 72, 13, 29, 42);
	}
	out->cd("Histos");
	TH1D *hSPE_UNCORR;
	TH1D *hSPE_CORR;
	if(numTS==1) 
	{
		hSPE_UNCORR = new TH1D("Signals_UNCORR", "Signals_UNCORR",100, 0, 0.02);
		hSPE_CORR = new TH1D("Signals_CORR", "Signals_CORR", 100, 0, 0.02);
	}
	else
	{
		hSPE_UNCORR = new TH1D("Signals_UNCORR", "Signals_UNCORR", 1000, 0, 0.2);
		hSPE_CORR = new TH1D("Signals_CORR", "Signals_CORR", 1000, 0, 0.2);
	}
	TH1D *hSFoverSB[NUMDEPTHS];
	hSFoverSB[0] = new TH1D("SFoverSB_D1", "SFoverSB_D1", 100, 0, 2);
	hSFoverSB[1] = new TH1D("SFoverSB_D2", "SFoverSB_D2", 100, 0, 2);
	TH1D *hSFoverSB_D12 = new TH1D("SFoverSB_D12", "SFoverSB_D12",
			100, 0, 2);
	TGraph *gSPEvsGain = new TGraph();
	gSPEvsGain->SetTitle("SignalvsGain");
	gSPEvsGain->SetName("SignalvsGain");

	if (numTS==1)
		cut_toUse = PEDCUT_1TS;
	else
		cut_toUse = PEDCUT_2TS;
	cout << "CUT=" << cut_toUse << endl;;

	//
	//	Init Swaps
	//
	initSwaps();

	//
	//	Read in GCF and Tubes Map
	//
	TubesMap map;
	readTubesMap(mapFileName, map, gcfFile, 0);
	int iTubeType = 0;
	int checkSkips = 0;

	vector<string> vStr;
	//
	//	Main Loop for analysis
	//
	int counter_Gain = 0;
	char histName[200];
	bool checkSwap = false;
	for (int iphi=0; iphi<NUMPHIS; iphi++)
		for (int ieta=0; ieta<NUMETAS; ieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				checkSwap = false;
				cout << "### Processing: " << 2*iphi+1 << "  " << ieta+29
					<< idepth+1 << endl;

				//
				//	Check if this channel is among swaps...
				//
				for (int iswap=0; iswap<NUMSWAPS; iswap++)
					if (swappedChs[iswap].iphi==(2*iphi+1) && 
							swappedChs[iswap].ieta==(ieta+29) &&
							swappedChs[iswap].idepth==idepth)
						checkSwap = true;
				if (checkSwap)
				{
					checkSkips++;
					continue;
				}

//				if (numTS==2)
//					if ((2*iphi+1)==51 && (ieta+29)==31)
//						continue;

				//
				//	Get a Sig histo
				//
				sprintf(histName, "1D_wSrc_0_PHI%d_ETA%d_D%d",
						2*iphi+1, ieta+29, idepth+1);
				in->cd("SRC/1D");
				TH1D *hRAW_S = (TH1D*)gDirectory->Get(histName);
				if (hRAW_S->GetEntries()<100)
					continue;

				//
				//	Get a Sig_F histo
				//
				sprintf(histName, "1D_F_wSrc_0_PHI%d_ETA%d_D%d",
						2*iphi+1, ieta+29, idepth+1);
				in->cd("SRC/1D_F");
				TH1D *hRAW_S_F = (TH1D*)gDirectory->Get(histName);
				if (hRAW_S_F->GetEntries()<100)
					continue;

				//
				//	Get a Sig_B histo
				//
				sprintf(histName, "1D_B_wSrc_0_PHI%d_ETA%d_D%d",
						2*iphi+1, ieta+29, idepth+1);
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
				//	Get the GCF
				//
				double gcf = 1;
				if (idepth==0 &&map.tubes[iphi][ieta][iTubeType].groove=='H'
						|| idepth==1 && map.tubes[iphi][ieta][iTubeType].groove
						=='E')
					gcf = map.tubes[iphi][ieta][iTubeType].em;
				else if (idepth==0 && map.tubes[iphi][ieta][iTubeType].groove
						=='E' || idepth==1 
						&& map.tubes[iphi][ieta][iTubeType].groove=='H')
					gcf = map.tubes[iphi][ieta][iTubeType].had;
				cout << gcf << endl;

				//
				//	Fit PED Region of S/B Histos,
				//	Extract the Means of Pedestals
				//	Compute the (Qs - Ps) - (Qb - Pb)
				//
				fit(hRAW_S);
				out->cd("Histos/1D/SIG");
				hRAW_S->Write();
				fit(hRAW_B);
				out->cd("Histos/1D/BACK");
				hRAW_B->Write();
				fit(hRAW_S_F);
				out->cd("Histos/1D_F");
				hRAW_S_F->Write();
				fit(hRAW_S_B);
				out->cd("Histos/1D_B");
				hRAW_S_B->Write();

				TF1 *fit_S = hRAW_S->GetFunction("myGaus");
				TF1 *fit_S_F = hRAW_S_F->GetFunction("myGaus");
				TF1 *fit_S_B = hRAW_S_B->GetFunction("myGaus");
				TF1 *fit_B = hRAW_B->GetFunction("myGaus");
				Double_t qie_mean_S = fit_S->GetParameter(1);
				Double_t qie_mean_B = fit_B->GetParameter(1);
				Double_t qie_sigma_S = fit_S->GetParameter(2);
				Double_t qie_sigma_B = fit_B->GetParameter(2);
				
				//
				//	Fill the 2D maps of PED Mean/Sigmas
				//
				hPedMean_S[idepth]->Fill(2*iphi+1, ieta+29, qie_mean_S);
				hPedMean_B[idepth]->Fill(2*iphi+1, ieta+29, qie_mean_B);
				hPedSigma_S[idepth]->Fill(2*iphi+1, ieta+29, qie_sigma_S);
				hPedSigma_B[idepth]->Fill(2*iphi+1, ieta+29, qie_sigma_B);

				double ratio = getSwapRatio(hRAW_S, numTS);
				hSwapRatio[idepth]->Fill(2*iphi+1, ieta+29, ratio);

				Double_t qie_mean_S_F = fit_S_F->GetParameter(1);
				Double_t qie_mean_S_B = fit_S_B->GetParameter(1);
				Double_t totQSum_S = 0;
				Double_t totQSum_B = 0;
				Double_t totNSum_S = 0;
				Double_t totNSum_B = 0;
				Double_t totQSum_S_F = 0;
				Double_t totNSum_S_F = 0;
				Double_t totQSum_S_B = 0;
				Double_t totNSum_S_B = 0;
				for (int iBin=0; iBin<NUMBINS; iBin++)
				{
					totNSum_S += hRAW_S->GetBinContent(iBin+1)*
						hRAW_S->GetBinWidth(iBin+1);
					totNSum_B += hRAW_B->GetBinContent(iBin+1)*
						hRAW_B->GetBinWidth(iBin+1);
					totNSum_S_F += hRAW_S_F->GetBinContent(iBin+1)*
						hRAW_S_F->GetBinWidth(iBin+1);
					totNSum_S_B += hRAW_S_B->GetBinContent(iBin+1)*
						hRAW_S_B->GetBinWidth(iBin+1);
					totQSum_S+= hRAW_S->GetBinCenter(iBin+1)*
						hRAW_S->GetBinContent(iBin+1)*
						hRAW_S->GetBinWidth(iBin+1);
					totQSum_B += hRAW_B->GetBinCenter(iBin+1)*
						hRAW_B->GetBinContent(iBin+1)*
						hRAW_B->GetBinWidth(iBin+1);
					totQSum_S_F+=hRAW_S_F->GetBinCenter(iBin+1)*
						hRAW_S_F->GetBinContent(iBin+1)*
						hRAW_S_F->GetBinWidth(iBin+1);
					totQSum_S_B+=hRAW_S_B->GetBinCenter(iBin+1)*
						hRAW_S_B->GetBinContent(iBin+1)*
						hRAW_S_B->GetBinWidth(iBin+1);
				}

				Double_t speQ_S = totQSum_S/totNSum_S - qie_mean_S;
				Double_t speQ_B = totQSum_B/totNSum_B - qie_mean_B;
				Double_t speQ_UNCORR = speQ_S - speQ_B;
				Double_t speQ_CORR = speQ_UNCORR/gcf;
				Double_t speQ_S_F = totQSum_S_F/totNSum_S_F - qie_mean_S_F;
				Double_t speQ_S_B = totQSum_S_B/totNSum_S_B - qie_mean_S_B;
				Double_t speQ_S_F_mBG = speQ_S_F - speQ_B;
				Double_t speQ_S_B_mBG = speQ_S_B - speQ_B;
				Double_t ratioSFSB = speQ_S_F_mBG/speQ_S_B_mBG;

				//
				//	Fill the Histos, Profiles, Graphs...
				//
				cout << speQ_S << "  " << speQ_B << endl;
				cout << speQ_UNCORR << "  " << speQ_CORR << endl;
				cout << speQ_S_F << "  " << speQ_S_B << endl;
				cout << speQ_S_F_mBG << "  " << speQ_S_B_mBG << endl;
				cout << totNSum_S_F << "  " << totNSum_S_B << endl;
				cout << totQSum_S_F << "  " << totQSum_S_B << endl;
				pSPEvsETA_UNCORR[idepth]->Fill(ieta+29, speQ_UNCORR);
				pSPEvsETA_CORR[idepth]->Fill(ieta+29, speQ_CORR);
				pSFoverSB[idepth]->Fill(ieta+29, ratioSFSB);
				pSFoverSB_D12->Fill(ieta+29, ratioSFSB);
				pSPEvsPHI_Avg[idepth]->Fill(2*iphi+1, speQ_CORR);
				pSFSBvsPHI_Avg[idepth]->Fill(2*iphi+1, ratioSFSB);
				pSPEvsPHI[ieta][idepth]->Fill(2*iphi+1, speQ_CORR);
				pSPEvsPHI[ieta][idepth]->Fill(2*iphi+1, speQ_CORR+0.0000001);
				pSFSBvsPHI[ieta][idepth]->Fill(2*iphi+1, ratioSFSB);
				pSFSBvsPHI[ieta][idepth]->Fill(2*iphi+1, ratioSFSB+0.000001);


				hSPE_UNCORR->Fill(speQ_UNCORR);
				hSPE_CORR->Fill(speQ_CORR);
				hSFoverSB[idepth]->Fill(ratioSFSB);
				hSFoverSB_D12->Fill(ratioSFSB);

				if (numTS==1)
					outTxt << 2*iphi+1 << "  " << ieta+29 << "  " << idepth+1
						<< "  " << speQ_CORR << "  " 
						<< gainMap[iphi][ieta][idepth].gain_OV1 << "  "
						<< qie_mean_S << "  " 
						<< gainMap[iphi][ieta][idepth].gain_OV2
						<< endl;
				else if (numTS==2)
					outTxt << 2*iphi+1 << "  " << ieta+29 << "  " << idepth+1
						<< "  " << speQ_CORR << "  " 
						<< gainMap[iphi][ieta][idepth].gain_OV1P100 << "  "
						<< qie_mean_S << "  "
						<< gainMap[iphi][ieta][idepth].gain_OV2
						<< endl;

				//
				//	Plot SPE vs Gain
				//
				if (numTS==1)
					gSPEvsGain->SetPoint(counter_Gain, 
							gainMap[iphi][ieta][idepth].gain_OV1, speQ_CORR);
				else if (numTS==2)
					gSPEvsGain->SetPoint(counter_Gain, 
							gainMap[iphi][ieta][idepth].gain_OV1P100, 
							speQ_CORR);


				counter_Gain++;
			}

	for (int ieta=0; ieta<NUMETAS; ieta++)
	{
		for (int idepth=0; idepth<NUMDEPTHS; idepth++)
		{
			pSPEvsPHI[ieta][idepth]->GetXaxis()->SetTitle("iphi");
			pSPEvsPHI[ieta][idepth]->GetYaxis()->SetTitle("Signal(Corr. by GCF)");
			pSPEvsPHI[ieta][idepth]->GetYaxis()->SetTitleOffset(1.2);
			pSFSBvsPHI[ieta][idepth]->GetXaxis()->SetTitle("ieta");
			pSFSBvsPHI[ieta][idepth]->GetYaxis()->SetTitle("SF/SB");
			pSFSBvsPHI[ieta][idepth]->GetYaxis()->SetTitleOffset(1.2);
			pSFSBvsPHI[ieta][idepth]->GetYaxis()->SetRangeUser(0.8, 1.2);
			if (numTS==1)
				pSPEvsPHI[ieta][idepth]->GetYaxis()->SetRangeUser(0.001, 0.01);
			else if (numTS==2)
				pSPEvsPHI[ieta][idepth]->GetYaxis()->SetRangeUser(0.01, 0.1);
		}
	}

	//
	//	Save and Close
	//
	for (int i=1; i<=2; i++)
	{
		pSPEvsETA_UNCORR[i-1]->GetXaxis()->SetTitle("ieta");
		pSPEvsETA_UNCORR[i-1]->GetYaxis()->SetTitle("Signal");
		pSPEvsETA_UNCORR[i-1]->GetYaxis()->SetTitleOffset(1.2);
		pSPEvsETA_CORR[i-1]->GetXaxis()->SetTitle("ieta");
		pSPEvsETA_CORR[i-1]->GetYaxis()->SetTitle("Signal(Corr. by GCF)");
		pSPEvsETA_CORR[i-1]->GetYaxis()->SetTitleOffset(1.2);
		pSPEvsPHI_Avg[i-1]->GetXaxis()->SetTitle("iphi");
		pSPEvsPHI_Avg[i-1]->GetYaxis()->SetTitle("Signal(Corr. by GCF)");
		pSPEvsPHI_Avg[i-1]->GetYaxis()->SetTitleOffset(1.2);
		if (numTS==1)
		{
			pSPEvsETA_UNCORR[i-1]->GetYaxis()->SetRangeUser(0.001, 0.01);
			pSPEvsETA_CORR[i-1]->GetYaxis()->SetRangeUser(0.001, 0.01);
			pSPEvsPHI_Avg[i-1]->GetYaxis()->SetRangeUser(0.001, 0.01);
		}
		else if (numTS==2)
		{
			pSPEvsETA_UNCORR[i-1]->GetYaxis()->SetRangeUser(0.01, 0.1);
			pSPEvsETA_CORR[i-1]->GetYaxis()->SetRangeUser(0.01, 0.1);
			pSPEvsPHI_Avg[i-1]->GetYaxis()->SetRangeUser(0.01, 0.1);
		}

		pSFoverSB[i-1]->GetXaxis()->SetTitle("ieta");
		pSFoverSB[i-1]->GetYaxis()->SetTitle("SF/SB");
		pSFoverSB[i-1]->GetYaxis()->SetTitleOffset(1.2);
		pSFoverSB[i-1]->GetYaxis()->SetRangeUser(0.8, 1.2);

		pSFSBvsPHI_Avg[i-1]->GetXaxis()->SetTitle("iphi");
		pSFSBvsPHI_Avg[i-1]->GetYaxis()->SetTitle("SF/SB");
		pSFSBvsPHI_Avg[i-1]->GetYaxis()->SetTitleOffset(1.2);
		pSFSBvsPHI_Avg[i-1]->GetYaxis()->SetRangeUser(0.8, 1.2);

		hSFoverSB[i-1]->GetXaxis()->SetTitle("SF/SB");
	}
	pSFoverSB_D12->GetXaxis()->SetTitle("ieta");
	pSFoverSB_D12->GetYaxis()->SetTitle("SF/SB");
	pSFoverSB_D12->GetYaxis()->SetTitleOffset(1.1);
	pSFoverSB_D12->GetYaxis()->SetRangeUser(0.8, 1.2);
	hSFoverSB_D12->GetXaxis()->SetTitle("SF/SB");
	hSPE_UNCORR->GetXaxis()->SetTitle("Signal");
	hSPE_CORR->GetXaxis()->SetTitle("Signal(Corr. by GCF)");

	out->cd("Graphs");
	gSPEvsGain->Write();

	out->Write();

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
	hist->Fit("myGaus", "R");
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
int readTubesMap(string fileNameMap, TubesMap &map, string gcfFileName, 
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
		mapFile >> wedgeName >> iwaste >> swaste
			>> eta >> eta >> phi >> iwaste >> iwaste 
			>> tubeStart >> tubeEnd >> groove;

		int iphi = (phi-1)/2;
		int ieta = abs(eta) - 29;
		int checkRBX, checkPhi, checkEta, checkTubeNum;
		char tubeType;
		sscanf(tubeName.c_str(), "HFP%d_ETA%d_PHI%d_T%d%c%*[^]",
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
