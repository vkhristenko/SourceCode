





#define NUMPHIS 36
#define NUMETAS 13
#define NUMDEPTHS 2



struct TCoords
{
	double signal;
	double gain_OV1;
	double gain_OV1P100;
	double gain_OV2;
};


#include "map.cc"

void plot2TSvs1TS(string fileName_1TS, string fileName_2TS)
{
	ifstream in_1TS(fileName_1TS.c_str());
	TCoords spe_1TS[NUMPHIS][NUMETAS][NUMDEPTHS];
	int phi, eta, depth;
	int iphi, ieta, idepth;
	double dwaste;
	double gain_OV1, gain_OV1P100, gain_OV2;
	Double_t spe;
	while (in_1TS >> phi)
	{
		in_1TS >> eta >> depth >> spe >> gain_OV1 >> dwaste >> gain_OV2;
		iphi = (phi-1)/2;
		ieta = eta-29;
		idepth = depth-1;
		spe_1TS[iphi][ieta][idepth].signal = spe;
		spe_1TS[iphi][ieta][idepth].gain_OV1 = gain_OV1;
		spe_1TS[iphi][ieta][idepth].gain_OV2 = gain_OV2;
	}

	ifstream in_2TS(fileName_2TS.c_str());
	TCoords spe_2TS[NUMPHIS][NUMETAS][NUMDEPTHS];
	while (in_2TS >> phi)
	{
		in_2TS >> eta >> depth >> spe >> gain_OV1P100 >> dwaste >> gain_OV2;
		iphi = (phi-1)/2;
		ieta = eta-29;
		idepth = depth-1;
		spe_2TS[iphi][ieta][idepth].signal = spe;
		spe_2TS[iphi][ieta][idepth].gain_OV1P100 = gain_OV1P100;
		spe_2TS[iphi][ieta][idepth].gain_OV2 = gain_OV2;
	}

	TFile *out = new TFile("out_scatter.root", "recreate");
	TGraph *gScatter = new TGraph();
	TH1D *h1TSover2TS = new TH1D("1TSover2TS", "1TS over 2TS", 1000, 0, 2);
	gScatter->SetName("1TSvs2TS");
	gScatter->SetTitle("1TSvs2TS");

	int counter=0;
	for (iphi=0; iphi<NUMPHIS; iphi++)
		for (ieta=0; ieta<NUMETAS; ieta++)
			for (idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				cout << "Running " << 2*iphi+1 << "  " << ieta+29 << "  "
					<< idepth+1 << endl;
				if (spe_1TS[iphi][ieta][idepth].signal!=0 && 
						spe_2TS[iphi][ieta][idepth].signal!=0)
				{
					cout << spe_1TS[iphi][ieta][idepth].signal << "  "
						<< spe_2TS[iphi][ieta][idepth].signal << endl;
					gScatter->SetPoint(counter, spe_1TS[iphi][ieta][idepth].signal,
							spe_2TS[iphi][ieta][idepth].signal);

					double S1G1 = spe_1TS[iphi][ieta][idepth].signal/
						spe_1TS[iphi][ieta][idepth].gain_OV1;
					double S2G2 = spe_2TS[iphi][ieta][idepth].signal/
						spe_2TS[iphi][ieta][idepth].gain_OV1P100/2.;
					double r1TSover2TS = S1G1/S2G2;
					h1TSover2TS->Fill(r1TSover2TS);
					counter++;
				}
			}

	gScatter->SetMarkerStyle(20);
	gScatter->SetMarkerColor(kBlack);
	gScatter->GetXaxis()->SetTitle("1TS");
	gScatter->GetXaxis()->SetTitleSize(0.06);
	gScatter->GetXaxis()->SetTitleOffset(0.8);
	gScatter->GetYaxis()->SetTitle("2TS");
	gScatter->GetYaxis()->SetTitleSize(0.06);
	gScatter->GetYaxis()->SetTitleOffset(0.8);
	gScatter->Draw("AP");

	out->Write();
	out->Close();

	return;
}
