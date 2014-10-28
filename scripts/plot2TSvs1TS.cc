





#define NUMPHIS 36
#define NUMETAS 13
#define NUMDEPTHS 2






#include "map.cc"

void plot2TSvs1TS(string fileName_1TS, string fileName_2TS)
{
	ifstream in_1TS(fileName_1TS.c_str());
	Double_t spe_1TS[NUMPHIS][NUMETAS][NUMDEPTHS];
	int phi, eta, depth;
	int iphi, ieta, idepth;
	Double_t spe;
	while (in_1TS >> phi)
	{
		in_1TS >> eta >> depth >> spe;
		iphi = (phi-1)/2;
		ieta = eta-29;
		idepth = depth-1;
		spe_1TS[iphi][ieta][idepth] = spe;
	}

	ifstream in_2TS(fileName_2TS.c_str());
	Double_t spe_2TS[NUMPHIS][NUMETAS][NUMDEPTHS];
	while (in_2TS >> phi)
	{
		in_2TS >> eta >> depth >> spe;
		iphi = (phi-1)/2;
		ieta = eta-29;
		idepth = depth-1;
		spe_2TS[iphi][ieta][idepth] = spe;
	}

	TFile *out = new TFile("out_scatter.root", "recreate");
	TGraph *gScatter = new TGraph();
	gScatter->SetName("1TSvs2TS");
	gScatter->SetTitle("1TSvs2TS");

	int counter=0;
	for (iphi=0; iphi<NUMPHIS; iphi++)
		for (ieta=0; ieta<NUMETAS; ieta++)
			for (idepth=0; idepth<NUMDEPTHS; idepth++)
			{
				cout << "Running " << 2*iphi+1 << "  " << ieta+29 << "  "
					<< idepth+1 << endl;
				if (spe_1TS[iphi][ieta][idepth]!=0 && 
						spe_2TS[iphi][ieta][idepth]!=0)
				{
					cout << spe_1TS[iphi][ieta][idepth] << "  "
						<< spe_2TS[iphi][ieta][idepth] << endl;
					gScatter->SetPoint(counter, spe_1TS[iphi][ieta][idepth],
							spe_2TS[iphi][ieta][idepth]);
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

	return;
}
