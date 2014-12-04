
















struct TResult 
{
	double signal_S_NoOF;
	double signal_B_NoOF;
	double qie_mean_S;
	double qie_mean_B;
	double signal_wOF;
	double signal_NoOF;
	double signal_OV2_wOF;
	double signal_OV2_NoOF;
	double ratioSB23toS;
	double adc2GeV_OV2_wOF;
	double adc2GeV_OV2_NoOF;
	double gain_OV1;
	double gain_OV2;
	double gain_OV1P100;
	double gcf;
};

#define ADC2fC 0.3593

#include "analyze_AvgQ.cc"

typedef TResult TResults[NUMTUBETYPES][NUMPHIS][NUMETAS][NUMDEPTHS];

void readTxt(TResults r, string txt);
void check(TResults r);

void analyze_Final(int hfSide, string rootFileName_1TS, string txtFileName_1TS,
		string rootFileName_2TS, string txtFileName_2TS,
		string swapsFileName, string outFileName, string outTxtFileName)
{
	TFile *rootIn_1TS = new TFile(rootFileName_1TS.c_str());
	TFile *rootIn_2TS = new TFile(rootFileName_2TS.c_str());
	TSwaps swaps;
	readSwaps(swapsFileName, swaps, 0);
	TResults results_1TS;
	init(results_1TS);
	TResults results_2TS;
	init(results_2TS);
	readTxt(results_1TS, txtFileName_1TS);
	readTxt(results_2TS, txtFileName_2TS);
	ofstream outTxt(outTxtFileName.c_str());
//	check(results_1TS);
//	check(results_2TS);

	TFile *out = new TFile(outFileName.c_str(), "recreate");
	TH1D *hADC2GeV_OV2_1TS[2][2];
	TH1D *hADC2GeV_OV2_2TS[2][2];
	TH1D *hADC2GeV_OV2_Final[2][2];
	TH1D *hADC2GeV_OV2_1TS2TSAvg[2][2];
	TH1D *hfC2GeV_OV2_Final[2][2];
	TH1D *hSignal_OV2_1TS[2][2]; 
	TH1D *hSignal_OV2_2TS[2][2];
	TH1D *hSignal_OV2_Final[2][2];
	TH1D *hSignal_OV2_1TS2TSAvg[2][2];
	TH1D *h1TSto2TS[2][2];
	TH1D *hNoOF2wOF[2];
	TH1D *hA2B[2];

	char histName[200];
	char titleName[200];
	for (int iOF=0; iOF<2; iOF++)
	{
		sprintf(histName, "A2B_%s", iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "A2B %s", iOF==0 ? "withOF" : "withoutOF");
		hA2B[iOF] = new TH1D(histName, titleName, 1000, 0, 2);
	for (int itube=0; itube<NUMTUBETYPES; itube++)
	{
		sprintf(histName, "fC2GeV_OV2_Final_TT%d_%s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "fC2GeV @OV2 Final TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hfC2GeV_OV2_Final[itube][iOF] = new TH1D(histName, titleName, 1000, 0, 2);

		sprintf(histName, "ADC2GeV_OV2_1TS_TT%d_%s", itube, 
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "ADC2GeV @OV2 for 1TS TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hADC2GeV_OV2_1TS[itube][iOF] = new TH1D(histName, titleName, 1000, 0, 2);

		sprintf(histName, "ADC2GeV_OV2_2TS_TT%d_%s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "ADC2GeV @OV2 for 2TS TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hADC2GeV_OV2_2TS[itube][iOF] = new TH1D(histName, titleName, 1000, 0, 2);

		sprintf(histName, "ADC2GeV_OV2_Final_TT%d_%s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "ADC2GeV @OV2 TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hADC2GeV_OV2_Final[itube][iOF] = new TH1D(histName, titleName, 1000, 0, 2);

		sprintf(histName, "Signal_OV2_1TS_TT%d_%s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "Signal @OV2 for 1TS TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hSignal_OV2_1TS[itube][iOF] = new TH1D(histName, titleName, 1000, 0, 
				0.02);

		sprintf(histName, "Signal_OV2_2TS_TT%d_%s", itube, 
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "Signal @OV2 for 2TS TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hSignal_OV2_2TS[itube][iOF] = new TH1D(histName, titleName, 1000, 0, 
				0.02);

		sprintf(histName, "Signal_OV2_Final_TT%d_%s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "Signal @OV2 TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hSignal_OV2_Final[itube][iOF] = new TH1D(histName, titleName, 1000, 0, 
				0.02);

		sprintf(histName, "1TSto2TS_TT%d_%s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "Ratio of 1TS to 2TS Results TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		h1TSto2TS[itube][iOF] = new TH1D(histName, titleName, 1000, 0, 2);
/*
		sprintf(histName, "Signal_OV2_1TS2TSAvg_TT%d_%s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "Signal @OV2 1TS2TSAvg TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hSignal_OV2_1TS2TS[itube][iOF] = new TH1D(histName, titleName,
				1000, 0, 0.02);

		sprintf(histName, "ADC2GeV_OV2_1TS2TSAvg_TT%d_%s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		sprintf(titleName, "ADC2GeV @OV2 1TS2TSAvg TT%d %s", itube,
				iOF==0 ? "withOF" : "withoutOF");
		hADC2GeV_OV2_1TS2TS[itube][iOF] = new TH1D(histName, titleName,
				1000, 0, 2);
*/
		if (iOF==0)
		{
			sprintf(histName, "NoOF2wOF_TT%d", itube);
			sprintf(titleName, "NoOF 2 wOF TT%d", itube);
			hNoOF2wOF[itube] = new TH1D(histName, titleName, 1000, 0, 2);
		}
	}
	}

	for (int itube=0; itube<NUMTUBETYPES; itube++)
	for (int iiphi=0; iiphi<NUMPHIS; iiphi++)
		for (int iieta=0; iieta<NUMETAS; iieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				//	Whenever we have results for 1TS 
				//	0 for wOF
				//	1 for NoOF
				if (results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF!=0)
				{	
					double signal_OV2_wOF = 
						results_1TS[itube][iiphi][iieta][idepth].signal_OV2_wOF;
					double signal_OV2_NoOF = 
						results_1TS[itube][iiphi][iieta][idepth].signal_OV2_NoOF;
					double adc2GeV_OV2_wOF = 
						results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF;
					double adc2GeV_OV2_NoOF = 
						results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_NoOF;
					double rNoOF2wOF = adc2GeV_OV2_NoOF/adc2GeV_OV2_wOF;
					double fC2GeV_OV2_wOF = adc2GeV_OV2_wOF*ADC2fC;
					double fC2GeV_OV2_NoOF = adc2GeV_OV2_NoOF*ADC2fC;

					hSignal_OV2_1TS[itube][0]->Fill(signal_OV2_wOF);
					hSignal_OV2_1TS[itube][1]->Fill(signal_OV2_NoOF);
					hSignal_OV2_Final[itube][0]->Fill(signal_OV2_wOF);
					hSignal_OV2_Final[itube][1]->Fill(signal_OV2_NoOF);
					hADC2GeV_OV2_1TS[itube][0]->Fill(adc2GeV_OV2_wOF);
					hADC2GeV_OV2_1TS[itube][1]->Fill(adc2GeV_OV2_NoOF);
					hADC2GeV_OV2_Final[itube][0]->Fill(adc2GeV_OV2_wOF);
					hADC2GeV_OV2_Final[itube][1]->Fill(adc2GeV_OV2_NoOF);
					hfC2GeV_OV2_Final[itube][0]->Fill(fC2GeV_OV2_wOF);
					hfC2GeV_OV2_Final[itube][1]->Fill(fC2GeV_OV2_NoOF);
					hNoOF2wOF[itube]->Fill(rNoOF2wOF);

					if (itube==1 &&
							results_1TS[0][iiphi][iieta][idepth].adc2GeV_OV2_wOF
							!=0)
					{
						double adc2GeV_OV2_wOF_A = 
							results_1TS[0][iiphi][iieta][idepth].adc2GeV_OV2_wOF;
						double adc2GeV_OV2_NoOF_A = 
							results_1TS[0][iiphi][iieta][idepth].adc2GeV_OV2_NoOF;

						double rA2B_wOF = adc2GeV_OV2_wOF_A/adc2GeV_OV2_wOF;
						double rA2B_NoOF = adc2GeV_OV2_NoOF_A/adc2GeV_OV2_NoOF;
						hA2B[0]->Fill(rA2B_wOF);
						hA2B[1]->Fill(rA2B_NoOF);
					}
				}

				if (results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF!=0 
						&& 
						results_2TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF
						!=0)
		i		{
					double adc2GeV_OV2_1TS_wOF = 
						results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF;
					double adc2GeV_OV2_2TS_wOF = 
						results_2TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF;
					double adc2GeV_OV2_1TS_NoOF = 
						results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_NoOF;
					double adc2GeV_OV2_2TS_NoOF =
						results_2TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_NoOF;

					double r1TSto2TS_wOF = adc2GeV_OV2_1TS_wOF/
						adc2GeV_OV2_2TS_wOF;
					double r1TSto2TS_NoOF = adc2GeV_OV2_1TS_NoOF/
						adc2GeV_OV2_2TS_NoOF;
					h1TSto2TS[itube][0]->Fill(r1TSto2TS_wOF);
					h1TSto2TS[itube][1]->Fill(r1TSto2TS_NoOF);
				}

				if (results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF==0 
						&& isSrcError(swaps, 2*iiphi+1, iieta+29, itube, 1)==1)
				{
					double signal_OV2_wOF = 
						results_2TS[itube][iiphi][iieta][idepth].signal_OV2_wOF;
					double signal_OV2_NoOF = 
						results_2TS[itube][iiphi][iieta][idepth].signal_OV2_NoOF;
					double adc2GeV_OV2_wOF = 
						results_2TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF;
					double adc2GeV_OV2_NoOF =
						results_2TS[itube][iiphi][iieta][idepth].adc2GeV_OV2_NoOF;
					double rNoOF2wOF = signal_OV2_NoOF/signal_OV2_wOF;
					double fC2GeV_OV2_wOF = adc2GeV_OV2_wOF*ADC2fC;
					double fC2GeV_OV2_NoOF = adc2GeV_OV2_NoOF*ADC2fC;

					hSignal_OV2_Final[itube][0]->Fill(signal_OV2_wOF);
					hSignal_OV2_Final[itube][1]->Fill(signal_OV2_NoOF);
					hADC2GeV_OV2_Final[itube][0]->Fill(adc2GeV_OV2_wOF);
					hADC2GeV_OV2_Final[itube][1]->Fill(adc2GeV_OV2_NoOF);
					hfC2GeV_OV2_Final[itube][0]->Fill(fC2GeV_OV2_wOF);
					hfC2GeV_OV2_Final[itube][1]->Fill(fC2GeV_OV2_NoOF);
					hNoOF2wOF[itube]->Fill(rNoOF2wOF);
				}
			}

	out->Write();
	return;
}

void readTxt(TResults r, string txtName)
{
	ifstream txt(txtName.c_str());
	int iphi, ieta, depth;
	int iiphi, iieta, idepth;
	int itube;
	double gcf;
	double signal_wOF, signal_NoOF, signal_OV2_wOF, signal_OV2_NoOF, 
		   adc2GeV_OV2_wOF, adc2GeV_OV2_NoOF, 
		   gain_OV1, gain_OV2, mean, gain_OV1P100,
		   signal_S_NoOF, signal_B_NoOF, qie_mean_S, qie_mean_B;
	double ratioSB23toS;
	while(txt >> itube)
	{
		txt >> iphi >> ieta >> depth >> gcf >> signal_S_NoOF >> signal_B_NoOF
			>> qie_mean_S >> qie_mean_B
			>> signal_wOF >> signal_NoOF >> signal_OV2_wOF >> signal_OV2_NoOF 
			>> adc2GeV_OV2_wOF >> adc2GeV_OV2_NoOF >> gain_OV1
			>> gain_OV2 >> gain_OV1P100;

		iiphi = (iphi-1)/2;
		iieta = ieta-29;
		idepth = depth;
		r[itube][iiphi][iieta][idepth].signal_wOF = signal_wOF;
		r[itube][iiphi][iieta][idepth].signal_NoOF = signal_NoOF;
		r[itube][iiphi][iieta][idepth].signal_OV2_wOF = signal_OV2_wOF;
		r[itube][iiphi][iieta][idepth].signal_OV2_NoOF = signal_OV2_NoOF;
//		r[itube][iiphi][iieta][idepth].ratioSB23toS = ratioSB23toS;
		r[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF = adc2GeV_OV2_wOF;
		r[itube][iiphi][iieta][idepth].adc2GeV_OV2_NoOF = adc2GeV_OV2_NoOF;
		r[itube][iiphi][iieta][idepth].gain_OV1 = gain_OV1;
		r[itube][iiphi][iieta][idepth].gain_OV2 = gain_OV2;
		r[itube][iiphi][iieta][idepth].gain_OV1P100 = gain_OV1P100;
		r[itube][iiphi][iieta][idepth].gcf = gcf;
	}

	txt.close();
	return;
}

init(TResults r)
{
	for (int itube=0; itube<NUMTUBETYPES; itube++)
	for (int iiphi=0; iiphi<NUMPHIS; iiphi++)
		for (int iieta=0; iieta<NUMETAS; iieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				r[itube][iiphi][iieta][idepth].signal_wOF = 0;
				r[itube][iiphi][iieta][idepth].signal_NoOF = 0;
				r[itube][iiphi][iieta][idepth].signal_OV2_wOF = 0;
				r[itube][iiphi][iieta][idepth].signal_OV2_NoOF = 0;
				r[itube][iiphi][iieta][idepth].adc2GeV_OV2_wOF = 0;
				r[itube][iiphi][iieta][idepth].adc2GeV_OV2_NoOF = 0;
				r[itube][iiphi][iieta][idepth].gain_OV1 = 0;
				r[itube][iiphi][iieta][idepth].gain_OV2 = 0;
				r[itube][iiphi][iieta][idepth].gain_OV1P100 = 0;
				r[itube][iiphi][iieta][idepth].gcf = 0;
			}
}

void check(TResults r)
{
	for (int itube=0; itube<NUMTUBETYPES; itube++)
	for (int iiphi=0; iiphi<NUMPHIS; iiphi++)
		for (int iieta=0; iieta<NUMETAS; iieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
				cout << r[itube][iiphi][iieta][idepth].signal_NoOF << "  " 
					<< r[itube][iiphi][iieta][idepth].adc2GeV_OV2_NoOF << endl;
}

