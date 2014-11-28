
















struct TResult 
{
	double signal;
	double signal_OV2;
	double ratioSB23toS;
	double adc2GeV_OV2;
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
	TH1D *hADC2GeV_OV2_1TS = new TH1D("ADC2GeV_OV2_1TS", "ADC2GeV @OV2 1TS", 
			10000, 0, 2);
	TH1D *hADC2GeV_OV2_2TS = new TH1D("ADC2GeV_OV2_2TS", "ADC2GeV_OV2_2TS",
			10000, 0, 2);
	TH1D *hADC2GeV_OV2_Final = new TH1D("ADC2GeV_OV2_Final", "ADC2GeV_OV2_Final",
			10000, 0, 2);
	TH1D *hfC2GeV_OV2_Final = new TH1D("fC2GeV_OV2_Final", "fC2GeV_OV2_Final",
			10000, 0, 2);
	TH1D *hADC2GeV_OV2_Recovered[2];
	hADC2GeV_OV2_Recovered[0] = new TH1D("ADC2GeV_OV2_Recovered_D1",
			"ADC2GeV_OV2_Recovered_D1", 10000, 0, 2);
	hADC2GeV_OV2_Recovered[1] = new TH1D("ADC2GeV_OV2_Recovered_D2",
			"ADC2GeV_OV2_Recovered_D2", 10000, 0, 2);
	TH1D *hGainsRatio = new TH1D("GainsRatio", "Gains Ratio", 10000, 0, 1);
	TH1D *hSignal_1TS = new TH1D("Signal_1TS", "Signal 1TS", 10000, 0, 0.02);
	TH1D *hSignal_1TS_Recovered[2];
	hSignal_1TS_Recovered[0] = new TH1D("Signal_1TS_Recovered_D1", 
			"Signal_1TS_Recovered_D1",10000,  0, 0.02);
	hSignal_1TS_Recovered[1] = new TH1D("Signal_1TS_Recovered_D2", 
			"Signal_1TS_Recovered_D2", 10000, 0, 0.02);
	TH1D *hSignal_1TS_Final = new TH1D("Signal_1TS_Final", "Signal_1TS_Final",
			10000, 0, 0.02);
	TH1D *hSignal_2TS = new TH1D("Signal_2TS", "Signal 2TS", 10000, 0, 0.2);
	TH1D *hSignal_OV2_1TS = new TH1D("Signal_OV2_1TS", "Signal @OV2 from 1TS",
			10000, 0, 0.02);
	TH1D *hSignal_OV2_Final = new TH1D("Signal_OV2_Final", "Signal @OV2 Final",
			10000, 0, 0.02);
	TH1D *hSignal_OV2_2TS = new TH1D("Signal_OV2_2TS", "Signal @OV2 from 2TS",
			10000, 0, 0.02);
	TH1D *h1TSto2TS = new TH1D("1TS/2TS", "1TS/2TS", 1000, 0.6, 1.4);
	TProfile* pSigVSieta[2];
	pSigVSieta[0] = new TProfile("Signal_OV1_VSieta_D1", "Signal @OV1 VS ieta D1",
			13, 29, 42);
	pSigVSieta[1] = new TProfile("Signal_OV1_VSieta_D2", "Signal @OV1 VS ieta D2",
			13, 29, 42);
	TH1D *hSB23toS[2];
	hSB23toS[0] = new TH1D("SB23toS_D1", "SB23toS_D1", 100, 0, 2);
	hSB23toS[1] = new TH1D("SB23toS_D2", "SB23toS_D2", 100, 0, 2);
	TProfile *pSB23toS[2];
	pSB23toS[0] = new TProfile("SB23toSvsieta_D1", "SB23toS vs ieta D1",
			13, 29, 42);
	pSB23toS[1] = new TProfile("SB23toSvsieta_D2", "SB23toS vs ieta D2",
			13, 29, 42);

	for (int itube=0; itube<NUMTUBETYPES; itube++)
	for (int iiphi=0; iiphi<NUMPHIS; iiphi++)
		for (int iieta=0; iieta<NUMETAS; iieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				if (results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2!=0)
//						results_2TS[iiphi][iieta][idepth].adc2GeV_OV2!=0)
				{
					double adc2GeV_1TS=
						results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2;
					double fC2GeV_OV2 = adc2GeV_1TS*ADC2fC;
					hADC2GeV_OV2_1TS->Fill(adc2GeV_1TS);
					hADC2GeV_OV2_Final->Fill(adc2GeV_1TS);
					hfC2GeV_OV2_Final->Fill(fC2GeV_OV2);
					double gain_OV1 = 
						results_1TS[itube][iiphi][iieta][idepth].gain_OV1;
					double gain_OV2 = 
						results_1TS[itube][iiphi][iieta][idepth].gain_OV2;
					double gainRatio = gain_OV2/gain_OV1;
					hGainsRatio->Fill(gainRatio);
					double signal_1TS = 
						results_1TS[itube][iiphi][iieta][idepth].signal;
					hSignal_1TS->Fill(signal_1TS);
					hSignal_1TS_Final->Fill(signal_1TS);
					double signal_OV2_1TS = 
						results_1TS[itube][iiphi][iieta][idepth].signal_OV2;
					hSignal_OV2_1TS->Fill(signal_OV2_1TS);
					hSignal_OV2_Final->Fill(signal_OV2_1TS);
					pSigVSieta[idepth]->Fill(iieta+29, signal_1TS);

					hSB23toS[idepth]->Fill(
						results_1TS[itube][iiphi][iieta][idepth].ratioSB23toS);
					pSB23toS[idepth]->Fill(iieta+29,
						results_1TS[itube][iiphi][iieta][idepth].ratioSB23toS);

					if (isTubeSwap(swaps, 2*iiphi+1, iieta+29)>=0 || 
							isChSwap(swaps, 2*iiphi+1, iieta+29, idepth+1)>=0 ||
							isSrcError(swaps, 2*iiphi+1, iieta+29, 2)==1)
					{
						hSignal_1TS_Recovered[idepth]->Fill(signal_1TS);
						hADC2GeV_OV2_Recovered[idepth]->Fill(adc2GeV_1TS);
					}

					double gcf = results_1TS[itube][iiphi][iieta][idepth].gcf;

					outTxt << itube << "  " << 2*iiphi+1 << "  " 
						<< iieta+29 << "  "	<< idepth+1 << "  " 
						<< gcf << "  " << signal_1TS <<  "  "
						<< adc2GeV_1TS << "  " << gain_OV1 << "  "
						<< gain_OV2
						<< endl;
				}

				if (results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2!=0 && 
						results_2TS[itube][iiphi][iieta][idepth].adc2GeV_OV2!=0)
				{
					double r1TSto2TS = 
						results_1TS[itube][iiphi][iieta][idepth].signal_OV2/
						results_2TS[itube][iiphi][iieta][idepth].signal_OV2;
					h1TSto2TS->Fill(r1TSto2TS);
				}

				if (results_1TS[itube][iiphi][iieta][idepth].adc2GeV_OV2==0 && 
						isSrcError(swaps, 2*iiphi+1, iieta+29, 1)==1)
				{
					double adc2GeV_OV2 = 
						results_2TS[itube][iiphi][iieta][idepth].adc2GeV_OV2;
					double fC2GeV_OV2 = adc2GeV_OV2*ADC2fC;
					hADC2GeV_OV2_Final->Fill(
							results_2TS[itube][iiphi][iieta][idepth].adc2GeV_OV2);
					hfC2GeV_OV2_Final->Fill(fC2GeV_OV2);
					double signal_OV2_2TS = 
						results_2TS[itube][iiphi][iieta][idepth].signal_OV2;
					double signal_OV1_2TS = 
						results_2TS[itube][iiphi][iieta][idepth].signal_OV2/
						results_2TS[itube][iiphi][iieta][idepth].gain_OV2*
						results_2TS[itube][iiphi][iieta][idepth].gain_OV1;
					hSignal_OV2_Final->Fill(signal_OV2_2TS);
					hSignal_1TS_Final->Fill(signal_OV1_2TS);
					pSigVSieta[idepth]->Fill(iieta+29, signal_OV1_2TS);

					double gain_OV1 = 
						results_2TS[itube][iiphi][iieta][idepth].gain_OV1;
					double gain_OV2 = 
						results_2TS[itube][iiphi][iieta][idepth].gain_OV2;
					hGainsRatio->Fill(gain_OV2/gain_OV1);

					hSB23toS[idepth]->Fill(
						results_2TS[itube][iiphi][iieta][idepth].ratioSB23toS);
					pSB23toS[idepth]->Fill(iieta+29,
						results_2TS[itube][iiphi][iieta][idepth].ratioSB23toS);
				
					hSignal_1TS_Recovered[idepth]->Fill(signal_OV1_2TS);
					hADC2GeV_OV2_Recovered[idepth]->Fill(adc2GeV_OV2);

					double gcf = results_1TS[itube][iiphi][iieta][idepth].gcf;

					outTxt << itube << "  " << 2*iiphi+1 << "  " 
						<< iieta+29 << "  "	<< idepth+1 << "  " 
						<< gcf << "  " << signal_OV1_2TS <<  "  "
						<< adc2GeV_OV2 << "  " << gain_OV1 << "  "
						<< gain_OV2
						<< endl;
				}
			}

	pSigVSieta[0]->GetYaxis()->SetRangeUser(0.001, 0.01);
	pSigVSieta[1]->GetYaxis()->SetRangeUser(0.001, 0.01);



//	pSB23toS[0]->GetYaxis()->SetRangeUser(0.8, 1.2);
//	pSB23toS[1]->GetYaxis()->SetRangeUser(0.8, 1.2);
	out->Write();
//	out->Close();
	return;
}

void readTxt(TResults r, string txtName)
{
	ifstream txt(txtName.c_str());
	int iphi, ieta, depth;
	int iiphi, iieta, idepth;
	int iTubeType;
	double gcf;
	double signal, signal_OV2, adc2GeV, gain_OV1, gain_OV2, mean, gain_OV1P100;
	double ratioSB23toS;
	while(txt >> iTubeType)
	{
		txt >> gcf >> iphi >> ieta >> depth 
			>> signal >> signal_OV2 >> ratioSB23toS
			>> adc2GeV >> gain_OV1
			>> mean >> gain_OV2 >> gain_OV1P100;

		iiphi = (iphi-1)/2;
		iieta = ieta-29;
		idepth = depth-1;
		r[iTubeType][iiphi][iieta][idepth].signal = signal;
		r[iTubeType][iiphi][iieta][idepth].signal_OV2 = signal_OV2;
		r[iTubeType][iiphi][iieta][idepth].ratioSB23toS = ratioSB23toS;
		r[iTubeType][iiphi][iieta][idepth].adc2GeV_OV2 = adc2GeV;
		r[iTubeType][iiphi][iieta][idepth].gain_OV1 = gain_OV1;
		r[iTubeType][iiphi][iieta][idepth].gain_OV2 = gain_OV2;
		r[iTubeType][iiphi][iieta][idepth].gain_OV1P100 = gain_OV1P100;
		r[iTubeType][iiphi][iieta][idepth].gcf = gcf;
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
				r[itube][iiphi][iieta][idepth].signal = 0;
				r[itube][iiphi][iieta][idepth].signal_OV2 = 0;
				r[itube][iiphi][iieta][idepth].adc2GeV_OV2 = 0;
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
				cout << r[itube][iiphi][iieta][idepth].signal << "  " 
					<< r[itube][iiphi][iieta][idepth].adc2GeV_OV2 << endl;
}

