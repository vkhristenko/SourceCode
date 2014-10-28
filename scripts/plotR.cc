














#define GEVPER25NS 744.6E-6

//
//	Define Source Activity Correction since November 2013 till 
//	the month when the Sourcing has been done
//
//#define SOURCEACTIVITYCORRECTION 0.88117
#define SOURCEACTIVITYCORRECTION 0.92398

void plotR(string results_AvgQ_1TS, string results_AvgQ_2TS,
		string results_Fixed_1TS, string results_Fixed_2TS,
		string results_Extr_1TS, string results_Extr_2TS)
{
	ifstream in_AvgQ_1TS(results_AvgQ_1TS.c_str());
	ifstream in_AvgQ_2TS(results_AvgQ_2TS.c_str());
	ifstream in_Fixed_1TS(results_Fixed_1TS.c_str());
	ifstream in_Fixed_2TS(results_Fixed_2TS.c_str());
	ifstream in_Extr_1TS(results_Extr_1TS.c_str());
	ifstream in_Extr_2TS(results_Extr_2TS.c_str());

	TFile *out_root = new TFile("out_plotR.root", "recreate");

	TH1D *hR_AvgQ = new TH1D("R_AvgQ", "R_AvgQ", 1000, 0, 10);
	TH1D *hR_Fixed = new TH1D("R_Fixed", "R_Fixed", 1000, 0, 10);
	TH1D *hR_Extr = new TH1D("R_Extr", "R_Extr", 1000, 0, 10);

	TH1D *hComp_AvgQ = new TH1D("1TSover2TS_AvgQ", "1TS Over 2TS", 1000, 0, 2);
	TH1D *hSignal_1TS_AvgQ = new TH1D("SignalOverVoltage_1TS_Avg",
			"Signal Over OV1 1TS Avg", 1000, 0, 0.0000002);
	TH1D *hSignal_2TS_AvgQ = new TH1D("SignalOverVoltage_2TS_Avg",
			"Signal Over OV1P100 2TS Avg", 1000, 0, 0.0000002);
	TH1D *hComp_Fixed = new TH1D("1TSover2TS_Fixed", "1TS Over 2TS", 1000, 0, 2);
	TH1D *hComp_Extr = new TH1D("1TSover2TS_Extr", "1TS Over 2TS", 1000, 0, 2);
	TH1D *hSignal_OV2_AvgQ = new TH1D("Signal_OV2_AvgQ", "Signal for OV2",
			10000, 0, 0.02);
	TH1D *hADC2GeV_OV2_AvgQ = new TH1D("ADC2GeV_AvgQ_OV2", "ADC2GeV OV2 AvgQ",
			10000, 0, 2);
	TH1D *hADC2GeV_OV2_AvgQ_WSAD = new TH1D("ADC2GeV_OV2_AvgQ_WSAD",
			"ADC2GeV OV2 AvgQ Correcting for Source Activity Decrease",
			10000, 0, 2);
	TH2D *hADC2GeV_OV2_Map_D1 = new TH2D("ADC2GeV_OV2_Map_D1",
			"ADC2GeV OV2 Map D1", 72, 0, 72, 13, 29, 42);
	TH2D *hADC2GeV_OV2_Map_D2 = new TH2D("ADC2GeV_OV2_Map_D2",
			"ADC2GeV OV2 Map D2", 72, 0, 72, 13, 29, 42);
	TH2D *hADC2GeV_OV2_Map_WSAD_D1 = new TH2D("ADC2GeV_OV2_Map_WSAD_D1",
			"ADC2GeV OV2 Map D1 Correcting for Source Activity Decrease", 
			72, 0, 72, 13, 29, 42);
	TH2D *hADC2GeV_OV2_Map_WSAD_D2 = new TH2D("ADC2GeV_OV2_Map_WSAD_D2",
			"ADC2GeV OV2 Map D2 Correcting for Source Activity Decrease", 
			72, 0, 72, 13, 29, 42);

	TH1D *hGainRatio = new TH1D("GainRatios", "GainRatios", 10000, 0, 1);
	ofstream out("ADC2GeV.txt");

	int iphi, ieta, depth;
	Double_t signal_AvgQ_1TS, signal_AvgQ_2TS,
			 signal_Fixed_1TS, signal_Fixed_2TS,
			 signal_Extr_1TS, signal_Extr_2TS;
	Double_t gain_OV1, gain_OV1P100, gain_OV2;
	Double_t mean_qie_1TS, mean_qie_2TS;
	while(in_AvgQ_1TS >> iphi)
	{
		in_AvgQ_1TS >> ieta >> depth >> signal_AvgQ_1TS
			>> gain_OV1 >> mean_qie_1TS >> gain_OV2;
		in_AvgQ_2TS >> iphi >> ieta >> depth >> signal_AvgQ_2TS
			>> gain_OV1P100 >> mean_qie_2TS >> gain_OV2;
		Double_t r = (10 - mean_qie_1TS)/(18 - mean_qie_2TS)*
			(gain_OV1P100/gain_OV1);
		hR_AvgQ->Fill(r);
		
		in_Fixed_1TS >> iphi >> ieta >> depth >> signal_Fixed_1TS
			>> gain_OV1 >> mean_qie_1TS >> gain_OV2;
		in_Fixed_2TS >> iphi >> ieta >> depth >> signal_Fixed_2TS
			>> gain_OV1P100 >> mean_qie_2TS >> gain_OV2;
		Double_t r = (10 - mean_qie_1TS)/(18 - mean_qie_2TS)*
			(gain_OV1P100/gain_OV1);
		hR_Fixed->Fill(r);


		in_Extr_1TS >> iphi >> ieta >> depth >> signal_Extr_1TS
			>> gain_OV1 >> mean_qie_1TS >> gain_OV2;
		in_Extr_2TS >> iphi >> ieta >> depth >> signal_Extr_2TS
			>> gain_OV1P100 >> mean_qie_2TS >> gain_OV2;
		Double_t r = (10 - mean_qie_1TS)/(18 - mean_qie_2TS)*
			(gain_OV1P100/gain_OV1);
		hR_Extr->Fill(r);

		Double_t S1G1 = signal_AvgQ_1TS/gain_OV1;
		Double_t S2G2 = signal_AvgQ_2TS/2./gain_OV1P100;
		Double_t comp = S1G1/S2G2;
		Double_t signal_OV2 = signal_AvgQ_1TS/gain_OV1*gain_OV2;
		Double_t adc2GeV = GEVPER25NS/signal_OV2; 
		Double_t adc2GeV_Corr = adc2GeV*SOURCEACTIVITYCORRECTION;
		hSignal_OV2_AvgQ->Fill(signal_OV2);
		hComp_AvgQ->Fill(comp);
		hSignal_1TS_AvgQ->Fill(S1G1);
		hSignal_2TS_AvgQ->Fill(S2G2);
		hADC2GeV_OV2_AvgQ->Fill(adc2GeV);
		hADC2GeV_OV2_AvgQ_WSAD->Fill(adc2GeV_Corr);

		if (depth==1)
		{
			hADC2GeV_OV2_Map_D1->Fill(iphi, ieta, adc2GeV);
			hADC2GeV_OV2_Map_WSAD_D1->Fill(iphi, ieta, adc2GeV_Corr);
		}
		else
		{
			hADC2GeV_OV2_Map_D2->Fill(iphi, ieta, adc2GeV);
			hADC2GeV_OV2_Map_WSAD_D2->Fill(iphi, ieta, adc2GeV_Corr);
		}

		out << iphi << "  " << ieta << "  " << depth << "  "
			<< signal_AvgQ_1TS << "  " << gain_OV1 << "  "
			<< signal_OV2 << "  " << gain_OV2 << "  "
			<< adc2GeV << endl;

		S1G1 = signal_Fixed_1TS/gain_OV1;
		S2G2 = signal_Fixed_2TS/2./gain_OV1P100;
		comp = S1G1/S2G2;
		hComp_Fixed->Fill(comp);

		S1G1 = signal_Extr_1TS/gain_OV1;
		S2G2 = signal_Extr_2TS/2./gain_OV1P100;
		comp = S1G1/S2G2;
		hComp_Extr->Fill(comp);

		Double_t gainRatio = gain_OV1/gain_OV1P100;
		hGainRatio->Fill(gainRatio);
	}

	TCanvas *cMaps = new TCanvas("cMaps", "cMaps", 200, 10, 800, 600);
	cMaps->Divide(2,2);
	cMaps->cd(1);
	hADC2GeV_OV2_Map_D1->GetZaxis()->SetRangeUser(0, 1.0);
	hADC2GeV_OV2_Map_D1->Draw("colz");
	cMaps->cd(2);
	hADC2GeV_OV2_Map_D2->GetZaxis()->SetRangeUser(0, 1.0);
	hADC2GeV_OV2_Map_D2->Draw("colz");
	cMaps->cd(3);
	hADC2GeV_OV2_Map_WSAD_D1->GetZaxis()->SetRangeUser(0, 1.0);
	hADC2GeV_OV2_Map_WSAD_D1->Draw("colz");
	cMaps->cd(4);
	hADC2GeV_OV2_Map_WSAD_D2->GetZaxis()->SetRangeUser(0, 1.0);
	hADC2GeV_OV2_Map_WSAD_D2->Draw("colz");

	TCanvas *cADC2GeV = new TCanvas("cADC2GeV_OV2", "ADC2GeV_OV2",
			200, 10, 800, 600);
	cADC2GeV->Divide(1,2);
	cADC2GeV->cd(1);
	hADC2GeV_OV2_AvgQ->Draw();
	cADC2GeV->cd(2);
	hADC2GeV_OV2_AvgQ_WSAD->Draw();

	TCanvas *cSignal_OV2 = new TCanvas("Signal_OV2", "Signal OV2",
			200, 10, 800, 600);
	cSignal_OV2->cd();
	hSignal_OV2_AvgQ->Draw();

	TCanvas *cSignals = new TCanvas("Signals", "Signals", 200, 10, 800, 600);
	hSignal_1TS_AvgQ->SetLineColor(kRed);
	hSignal_2TS_AvgQ->SetLineColor(kBlack);
	hSignal_1TS_AvgQ->Draw();
	hSignal_2TS_AvgQ->Draw("sames");

	char cName[200];
	string name("SGRatio");
	sprintf(cName, "c%s", name.c_str());
	TCanvas *cSGRatio = new TCanvas(cName, cName, 200, 10, 800, 600);
	cSGRatio->cd();

	hComp_AvgQ->SetLineColor(kRed);
	hComp_Fixed->SetLineColor(kBlue);
	hComp_Extr->SetLineColor(kBlack);

	hComp_AvgQ->SetLineWidth(2);
	hComp_Fixed->SetLineWidth(2);
	hComp_Extr->SetLineWidth(2);

	TLegend *legend = new TLegend(0.6, 0.8, 0.8, 0.9);
	legend->AddEntry(hComp_AvgQ, "1TS Over 2TS AvgQ");
	legend->AddEntry(hComp_Fixed, "1TS Over 2TS Fixed Interval");
	legend->AddEntry(hComp_Extr, "1TS Over 2TS Extrapolation");


	hComp_AvgQ->Draw();
//	hComp_Fixed->Draw("sames");
//	hComp_Extr->Draw("sames");
//	legend->Draw();
	

//	hR_AvgQ->Draw();
//	hR_Fixed->Draw("sames");
//	hR_Extr->Draw("sames");
//	legend->Draw();
	gPad->Update();
	gPad->Modified();

	TCanvas *cGains = new TCanvas("GainRatios", "GainRatios", 200, 10, 800, 600);
	cGains->cd();
	hGainRatio->Rebin(10);
	hGainRatio->Draw();

//	out_root->Write();
//	out_root->Close();

//	TPaveStats *st_AvgQ
}










