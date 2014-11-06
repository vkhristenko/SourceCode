










#include "analyze_Fixed.cc"

void genPlots(string results_AvgQ, string results_Fixed,
		string results_Extr)
{
/*
	plotDistrs(results_AvgQ, results_Fixed, results_Extr, string("UNCORR"),
			string("Signals"));
	plotDistrs(results_AvgQ, results_Fixed, results_Extr, string("CORR"),
			string("Signals"));
	plotDistrs(results_AvgQ, results_Fixed, results_Extr, string("D1"),
			string("SFoverSB"));
	plotDistrs(results_AvgQ, results_Fixed, results_Extr, string("D2"),
			string("SFoverSB"));


	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D1_UNCORR"),
			string("SignalvsETA"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D2_UNCORR"),
			string("SignalvsETA"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D1_CORR"),
			string("SignalvsETA"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D2_CORR"),
			string("SignalvsETA"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D1_Avg"),
			string("SignalvsPHI"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D2_Avg"),
			string("SignalvsPHI"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D1"),
			string("SFoverSB"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D2"),
			string("SFoverSB"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D1_Avg"),
			string("SFoverSBvsPHI"));
	plotProfiles(results_AvgQ, results_Fixed, results_Extr, string("D2_Avg"),
			string("SFoverSBvsPHI"));
*/	
	plotGraphs(results_AvgQ, results_Fixed, results_Extr, string(""), 
			string("SignalvsGain"));
//	plotSignalRatiosForDiffMethods(results_AvgQ, results_Fixed,
//			results_Extr, string(""), string("SignalRatios"));
//	plot2D(results_AvgQ);
}

void plot2D(string data_AvgQ)
{
	TFile *in = new TFile(data_AvgQ.c_str());
//	TFile *in_results = new TFile(results_AvgQ.c_str());

	TH2D *hPED_RMSOfMeans[2];
	TH2D *hPED_RMSOfRMS[2];
	hPED_RMSOfMeans[0] = new TH2D("PED_RMSOfMeans_D1", "RMS of PED Means D1",
			72, 0, 72, 13, 29, 42);
	hPED_RMSOfMeans[1] = new TH2D("PED_RMSOfMeans_D2", "RMS of PED Means D2",
			72, 0, 72, 13, 29, 42);
	hPED_RMSOfRMS[0] = new TH2D("PED_RMSOfRMS_D1", "RMS of RMS of PED D1",
			72, 0, 72, 13, 29, 42);
	hPED_RMSOfRMS[1] = new TH2D("PED_RMSOfRMS_D2", "RMS of RMS of PED D2",
			72, 0, 72, 13, 29, 42);

	char name[200];
	for (int iiphi=0; iiphi<NUMPHIS; iiphi++)
	{
		for (int iieta=0; iieta<NUMETAS; iieta++)
		{
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				int iphi = 2*iiphi+1;
				int ieta = iieta+29;
				int depth = idepth+1;

				in->cd("NOSRC/PEDMEANS");
				sprintf(name, "PEDMeans_PHI%d_ETA%d_D%d",
						iphi, ieta, depth);
				TH1D *h = (TH1D*)gDirectory->Get(name);
				double rmsOfMean = h->GetRMS();
				hPED_RMSOfMeans[idepth]->Fill(iphi, ieta, rmsOfMean);

				in->cd("NOSRC/PEDRMSS");
				sprintf(name, "PEDRMSs_PHI%d_ETA%d_D%d",
						iphi, ieta, depth);
				h = (TH1D*)gDirectory->Get(name);
				double rmsOfrms = h->GetRMS();
				hPED_RMSOfRMS[idepth]->Fill(iphi, ieta, rmsOfrms);

			}
		}
	}

	TCanvas *cPED = new TCanvas("cBGStability", "cBGStability",
			200, 10, 800, 600);
	cPED->Divide(2, 2);

	cPED->cd(1);
	hPED_RMSOfMeans[0]->Draw("colz");
	cPED->cd(2);
	hPED_RMSOfMeans[1]->Draw("colz");
	cPED->cd(3);
	hPED_RMSOfRMS[0]->Draw("colz");
	cPED->cd(4);
	hPED_RMSOfRMS[1]->Draw("colz");

	return;
}

void plotSignalRatiosForDiffMethods(string results_AvgQ, 
		string results_Fixed,
		string results_Extr, 
		string type, string histBaseName)
{
	ifstream in_AvgQ(results_AvgQ.c_str());
	ifstream in_Fixed(results_Fixed.c_str());
	ifstream in_Extr(results_Extr.c_str());

	TH1D *hRatioFixed2Average = new TH1D("Fixed2Average", 
			"Fixed2Average", 500, 0.8, 1.2);
	TH1D *hRatioExtr2Average = new TH1D("Extr2Average",
			"Extr2Average", 500, 0.8, 1.2);
	int iphi, ieta, depth;
	Double_t signal_AvgQ, signal_Fixed, signal_Extr;
	Double_t gain, qie_mean_ped; 
	while (in_AvgQ >> iphi)
	{
		in_AvgQ >> ieta >> depth >> signal_AvgQ >> gain >> qie_mean_ped;
		in_Fixed >> iphi >> ieta >> depth >> signal_Fixed >> gain
			>> qie_mean_ped;
		in_Extr >> iphi >> ieta >> depth >> signal_Extr >> gain 
			>> qie_mean_ped;

		Double_t fixed2Average = signal_Fixed/signal_AvgQ;
		Double_t extr2Average = signal_Extr/signal_AvgQ;
		hRatioFixed2Average->Fill(fixed2Average);
		hRatioExtr2Average->Fill(extr2Average);
	}

	char cName[200];
	sprintf(cName, "c%s", histBaseName.c_str());
	TCanvas *cSignals = new TCanvas(cName, cName, 200, 10, 800, 600);
	cSignals->cd();

	hRatioFixed2Average->SetLineColor(kRed);
	hRatioExtr2Average->SetLineColor(kBlue);
	hRatioFixed2Average->SetLineWidth(2);
	hRatioExtr2Average->SetLineWidth(2);

	TLegend * legend = new TLegend(0.6, 0.8, 0.8, 0.9);
	legend->AddEntry(hRatioFixed2Average);
	legend->AddEntry(hRatioExtr2Average);

	hRatioFixed2Average->Draw();
	hRatioExtr2Average->Draw("sames");
	legend->Draw();
	gPad->Update();
	
	TPaveStats *st_Fixed = 
		(TPaveStats*)hRatioFixed2Average->FindObject("stats");
	st_Fixed->SetY1NDC(0.8);
	st_Fixed->SetY2NDC(1);
	st_Fixed->SetX1NDC(0.9);
	st_Fixed->SetX2NDC(1);
	
	TPaveStats *st_Extr = 
		(TPaveStats*)hRatioExtr2Average->FindObject("stats");
	st_Extr->SetY1NDC(0.6);
	st_Extr->SetY2NDC(0.8);
	st_Extr->SetX1NDC(0.9);
	st_Extr->SetX2NDC(1.);
	
	gPad->Modified();
	
}

void plotGraphs(string results_AvgQ, string results_Fixed,
		string results_Extr, string type, string histBaseName)
{
	TFile *in_AvgQ = new TFile(results_AvgQ.c_str());
	TFile *in_Fixed = new TFile(results_Fixed.c_str());
	TFile *in_Extr = new TFile(results_Extr.c_str());

	char grName[200];
	char grTitle[200];
	sprintf(grName, "%s", histBaseName.c_str());
	in_AvgQ->cd("Graphs");
	TGraph *gSPEvsGain_AvgQ = (TGraph*)gDirectory->Get(grName);
	sprintf(grTitle, "%s_AvgQ", histBaseName.c_str());
	gSPEvsGain_AvgQ->SetNameTitle(grTitle, grTitle);

	in_Fixed->cd("Graphs");
	TGraph *gSPEvsGain_Fixed = (TGraph*)gDirectory->Get(grName);
	sprintf(grTitle, "%s_Fixed", histBaseName.c_str());
	gSPEvsGain_Fixed->SetNameTitle(grTitle, grTitle);

	in_Extr->cd("Graphs");
	TGraph *gSPEvsGain_Extr = (TGraph*)gDirectory->Get(grName);
	sprintf(grTitle, "%s_Extr", histBaseName.c_str());
	gSPEvsGain_Extr->SetNameTitle(grTitle, grTitle);

	gSPEvsGain_AvgQ->SetMarkerStyle(20);
	gSPEvsGain_Fixed->SetMarkerStyle(21);
	gSPEvsGain_Extr->SetMarkerStyle(22);

	gSPEvsGain_AvgQ->SetMarkerColor(kRed);
	gSPEvsGain_Fixed->SetMarkerColor(kBlue);
	gSPEvsGain_Extr->SetMarkerColor(kGreen);

	TLegend *legend = new TLegend(0.6, 0.8, 0.8, 0.9);
	legend->AddEntry(gSPEvsGain_AvgQ);
	legend->AddEntry(gSPEvsGain_Fixed);
	legend->AddEntry(gSPEvsGain_Extr);

	char cName[200];
	sprintf(cName, "c%s", histBaseName.c_str());
	TCanvas *cSPEvsGain = new TCanvas(cName, cName, 200, 10, 800, 600);
//	gSPEvsGain_AvgQ->Draw("ap");
//	gSPEvsGain_Fixed->Draw("ap");
//	gSPEvsGain_Extr->Draw("sameap");

	TMultiGraph *mg = new TMultiGraph();
	mg->Add(gSPEvsGain_AvgQ);
//	mg->Add(gSPEvsGain_Fixed);
//	mg->Add(gSPEvsGain_Extr);
	mg->Draw("ap");
	legend->Draw();
	cSPEvsGain->Update();
	cSPEvsGain->Modified();
}

void plotProfiles(string results_AvgQ, string results_Fixed,
		string results_Extr, string type, string histBaseName)
{		
	TFile *in_AvgQ = new TFile(results_AvgQ.c_str());
	TFile *in_Fixed = new TFile(results_Fixed.c_str());
	TFile *in_Extr = new TFile(results_Extr.c_str());

	char histName[200];
	char histTitle[200];
	//
	//	Plot Uncorrected SPE Distrs overlaid
	//
	sprintf(histName, "%s_%s", histBaseName.c_str(), type.c_str());
	in_AvgQ->cd("Profiles");
	TProfile *hSPEs_UNCORR_AvgQ = (TProfile*)gDirectory->Get(histName);
	sprintf(histTitle, "%s_%s_AvgQ", histBaseName.c_str(), type.c_str());
	hSPEs_UNCORR_AvgQ->SetNameTitle(histTitle, histTitle);
	in_Fixed->cd("Profiles");
	TProfile *hSPEs_UNCORR_Fixed = (TProfile*)gDirectory->Get(histName);
	sprintf(histTitle, "%s_%s_Fixed", histBaseName.c_str(), type.c_str());
	hSPEs_UNCORR_Fixed->SetNameTitle(histTitle, histTitle);
	in_Extr->cd("Profiles");
	TProfile *hSPEs_UNCORR_Extr = (TProfile*)gDirectory->Get(histName);
	sprintf(histTitle, "%s_%s_Extr", histBaseName.c_str(), type.c_str());
	hSPEs_UNCORR_Extr->SetNameTitle(histTitle, histTitle);
	
	TLegend *legend = new TLegend(0.6, 0.8, 0.8, 0.9);

	char cName[200];
	sprintf(cName, "c%s_%s", histBaseName.c_str(), type.c_str());
	TCanvas *cSPEs_UNCORR = new TCanvas(cName, cName, 200, 10, 800, 600);
	hSPEs_UNCORR_AvgQ->SetLineColor(kRed);
	hSPEs_UNCORR_Fixed->SetLineColor(kBlue);
	hSPEs_UNCORR_Extr->SetLineColor(kBlack);

	hSPEs_UNCORR_AvgQ->SetLineWidth(2);
	hSPEs_UNCORR_Fixed->SetLineWidth(2);
	hSPEs_UNCORR_Extr->SetLineWidth(2);

	legend->AddEntry(hSPEs_UNCORR_AvgQ);
	legend->AddEntry(hSPEs_UNCORR_Fixed);
	legend->AddEntry(hSPEs_UNCORR_Extr);

	hSPEs_UNCORR_AvgQ->Draw();
	hSPEs_UNCORR_Fixed->Draw("sames");
	hSPEs_UNCORR_Extr->Draw("sames");
	legend->Draw();

	gPad->Update();
	TPaveStats *st_AvgQ = (TPaveStats*)hSPEs_UNCORR_AvgQ->FindObject("stats");
	st_AvgQ->SetY1NDC(0.8);
	st_AvgQ->SetY2NDC(1);
	st_AvgQ->SetX1NDC(0.9);
	st_AvgQ->SetX2NDC(1);
	TPaveStats *st_Fixed = (TPaveStats*)hSPEs_UNCORR_Fixed->FindObject("stats");
	TPaveStats *st_Fixed = (TPaveStats*)hSPEs_UNCORR_Fixed->FindObject("stats");
	st_Fixed->SetY1NDC(0.6);
	st_Fixed->SetY2NDC(0.8);
	st_Fixed->SetX1NDC(0.9);
	st_Fixed->SetX2NDC(1);
	TPaveStats *st_Fixed = (TPaveStats*)hSPEs_UNCORR_Fixed->FindObject("stats");
	TPaveStats *st_Extr = (TPaveStats*)hSPEs_UNCORR_Extr->FindObject("stats");
	st_Extr->SetY1NDC(0.6);
	st_Extr->SetY2NDC(0.4);
	st_Extr->SetX1NDC(0.9);
	st_Extr->SetX2NDC(1);
	TPaveStats *st_Fixed = (TPaveStats*)hSPEs_UNCORR_Fixed->FindObject("stats");
	gPad->Modified();

	char fileName[200];
	sprintf(fileName, "pics_2TS_02102014/%s.png", cSPEs_UNCORR->GetName());
	cSPEs_UNCORR->SaveAs(fileName);
}

void plotDistrs(string results_AvgQ, string results_Fixed,
		string results_Extr, string type, string histBaseName)
{
	TFile *in_AvgQ = new TFile(results_AvgQ.c_str());
	TFile *in_Fixed = new TFile(results_Fixed.c_str());
	TFile *in_Extr = new TFile(results_Extr.c_str());

	char histName[200];
	char histTitle[200];
	//
	//	Plot Uncorrected SPE Distrs overlaid
	//
	sprintf(histName, "%s_%s", histBaseName.c_str(), type.c_str());
	in_AvgQ->cd("Histos");
	TH1D *hSPEs_UNCORR_AvgQ = (TH1D*)gDirectory->Get(histName);
	sprintf(histTitle, "%s_%s_AvgQ", histBaseName.c_str(), type.c_str());
	hSPEs_UNCORR_AvgQ->SetNameTitle(histTitle, histTitle);
	hSPEs_UNCORR_AvgQ->Rebin(2);
	in_Fixed->cd("Histos");
	TH1D *hSPEs_UNCORR_Fixed = (TH1D*)gDirectory->Get(histName);
	sprintf(histTitle, "%s_%s_Fixed", histBaseName.c_str(), type.c_str());
	hSPEs_UNCORR_Fixed->SetNameTitle(histTitle, histTitle);
	hSPEs_UNCORR_Fixed->Rebin(2);
	in_Extr->cd("Histos");
	TH1D *hSPEs_UNCORR_Extr = (TH1D*)gDirectory->Get(histName);
	sprintf(histTitle, "%s_%s_Extr", histBaseName.c_str(), type.c_str());
	hSPEs_UNCORR_Extr->SetNameTitle(histTitle, histTitle);
	hSPEs_UNCORR_Extr->Rebin(2);
	
	TLegend *legend = new TLegend(0.8, 0.2, 1, 0.4);

	char cName[200];
	sprintf(cName, "c%s_%s", histBaseName.c_str(), type.c_str());
	TCanvas *cSPEs_UNCORR = new TCanvas(cName, cName, 200, 10, 800, 600);
	hSPEs_UNCORR_AvgQ->SetLineColor(kRed);
	hSPEs_UNCORR_Fixed->SetLineColor(kBlue);
	hSPEs_UNCORR_Extr->SetLineColor(kBlack);

	hSPEs_UNCORR_AvgQ->SetLineWidth(2);
	hSPEs_UNCORR_Fixed->SetLineWidth(2);
	hSPEs_UNCORR_Extr->SetLineWidth(2);

	legend->AddEntry(hSPEs_UNCORR_AvgQ);
	legend->AddEntry(hSPEs_UNCORR_Fixed);
	legend->AddEntry(hSPEs_UNCORR_Extr);

	hSPEs_UNCORR_AvgQ->Draw();
	hSPEs_UNCORR_Fixed->Draw("sames");
	hSPEs_UNCORR_Extr->Draw("sames");
	legend->Draw();

	gPad->Update();
	TPaveStats *st_AvgQ = (TPaveStats*)hSPEs_UNCORR_AvgQ->FindObject("stats");
	st_AvgQ->SetY1NDC(0.8);
	st_AvgQ->SetY2NDC(1);
	st_AvgQ->SetX1NDC(0.8);
	st_AvgQ->SetX2NDC(1);
	TPaveStats *st_Fixed = (TPaveStats*)hSPEs_UNCORR_Fixed->FindObject("stats");
	TPaveStats *st_Fixed = (TPaveStats*)hSPEs_UNCORR_Fixed->FindObject("stats");
	st_Fixed->SetY1NDC(0.6);
	st_Fixed->SetY2NDC(0.8);
	st_Fixed->SetX1NDC(0.8);
	st_Fixed->SetX2NDC(1);
	TPaveStats *st_Fixed = (TPaveStats*)hSPEs_UNCORR_Fixed->FindObject("stats");
	TPaveStats *st_Extr = (TPaveStats*)hSPEs_UNCORR_Extr->FindObject("stats");
	st_Extr->SetY1NDC(0.6);
	st_Extr->SetY2NDC(0.4);
	st_Extr->SetX1NDC(0.8);
	st_Extr->SetX2NDC(1);
	TPaveStats *st_Fixed = (TPaveStats*)hSPEs_UNCORR_Fixed->FindObject("stats");
	gPad->Modified();

	char fileName[200];
	sprintf(fileName, "pics_2TS_02102014/%s.png", cSPEs_UNCORR->GetName());
	cSPEs_UNCORR->SaveAs(fileName);
}







