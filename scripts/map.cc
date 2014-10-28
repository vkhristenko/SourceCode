




#define NUMPHIS 36
#define NUMETAS 13
#define NUMDEPTHS 2
#define NUMRBXS 36
#define NUMBBS 3
#define NUMPMTS 8




struct TEtaPhiDepth
{
	int ieta;
	int iphi;
	int depth;
};

struct TRbxBBPmt
{
	int rbx;
	int bb;
	int pmt;
};

struct TEMap
{
	TRbxBBPmt eta2rbxSpace[NUMPHIS][NUMETAS][NUMDEPTHS];
	TEtaPhiDepth rbx2etaSpace[NUMRBXS][NUMBBS][NUMPMTS];
};

typedef struct TEtaPhiDepthGain
{
	Double_t gain_OV1;
	Double_t gain_OV2;
	Double_t gain_OV1P100;
};

typedef TEtaPhiDepthGain TGainMap[NUMPHIS][NUMETAS][NUMDEPTHS];

//
//	Map rbx Space to Eta/Phi Space
//
TEtaPhiDepth rbx2eta(TRbxBBPmt rbxBBPmt, TEMap& map)
{
	TEtaPhiDepth etaPhiDepth;
	int iirbx = rbxBBPmt.rbx-1;
	int iibb = rbxBBPmt.bb-1;
	int iipmt = rbxBBPmt.pmt-1;
	etaPhiDepth.ieta = map.rbx2etaSpace[iirbx][iibb][iipmt].ieta;
	etaPhiDepth.iphi = map.rbx2etaSpace[iirbx][iibb][iipmt].iphi;
	etaPhiDepth.depth = map.rbx2etaSpace[iirbx][iibb][iipmt].depth;

	return etaPhiDepth;
}

//
//	Map eta/phi space to rbx space
//
TRbxBBPmt eta2rbx(TEtaPhiDepth etaPhiDepth, TEMap &map)
{
	TRbxBBPmt rbxBBPmt;
	int iiphi = (etaPhiDepth.iphi-1)/2;
	int iieta = etaPhiDepth.ieta-29;
	int iidepth = etaPhiDepth.depth-1;
	rbxBBPmt.rbx = map.eta2rbxSpace[iiphi][iieta][iidepth].rbx;
	rbxBBPmt.bb = map.eta2rbxSpace[iiphi][iieta][iidepth].bb;
	rbxBBPmt.pmt = map.eta2rbxSpace[iiphi][iieta][iidepth].pmt;

	return rbxBBPmt;
}

//
//	Get Gain Map
//
void getGainMap(string gainMapFileName, TGainMap gainMap, TEMap eMap)
{
	cout << "### Reading Gain Map..." << endl;

	ifstream in(gainMapFileName.c_str());
	string shf, sbb, spmt;
	Double_t gain_OV1, gain_OV2, gain_OV1P100;
	int rbx, bb, pmt;
	char hfSide;
	while (in >> shf)
	{
		in >> sbb >> spmt >> gain_OV1 >> gain_OV2 >> gain_OV1P100;
		sscanf(shf.c_str(), "HF%c%d", &hfSide, &rbx);

//		cout << shf << " " << sbb << " " << spmt << " " << endl;

		//
		//	MOdify for HFP
		//
		if (hfSide=='M')
			continue;
		sscanf(sbb.c_str(), "BB%d", &bb);
		sscanf(spmt.c_str(), "PMT%d", &pmt);

//		cout << rbx << "  " << bb << "  " << pmt << endl;
//			<< "  " << gain_OV2 << endl;;

		TRbxBBPmt rbxBBPmt;
		rbxBBPmt.rbx = rbx;
		rbxBBPmt.bb = bb;
		rbxBBPmt.pmt = pmt;

		TEtaPhiDepth etaPhiDepth = rbx2eta(rbxBBPmt, eMap);
		int iiphi = (etaPhiDepth.iphi-1)/2;
		int iieta = (etaPhiDepth.ieta-29);
		int iidepth = etaPhiDepth.depth-1;

//		cout << etaPhiDepth.iphi << "  " << etaPhiDepth.ieta << "  "
//			<< etaPhiDepth.depth << endl;

//	cout << iiphi << "  " << iieta << "  " << iidepth << endl;

		gainMap[iiphi][iieta][iidepth].gain_OV1 = gain_OV1;
		gainMap[iiphi][iieta][iidepth].gain_OV2 = gain_OV2;
		gainMap[iiphi][iieta][iidepth].gain_OV1P100 = gain_OV1P100;
	}
	cout << "### Done..." << endl;

	/*
	for (int iphi=0; iphi<NUMPHIS; iphi++)
		for (int ieta=0; ieta<NUMETAS; ieta++)
			for (int idepth=0; idepth<NUMDEPTHS; idepth++)
				cout << iphi << "  " << ieta << "  " << idepth << "  "
					<< gainMap[iphi][ieta][idepth].gain_OV1 << endl;
					*/

	return;
}

//
//	Get the map
//
void getMap(string mapFileName, TEMap &map)
{
	cout << "### Reading Electronics Map..." << endl;

	ifstream in(mapFileName.c_str());
	int iphi, ieta, depth, rbx, bb, pmt;
	while(in >> iphi)
	{
		in >> ieta >> depth >> rbx >> bb >> pmt;
//		rbx -= 36;
		int iiphi = (iphi-1)/2;
		int iieta = abs(ieta)-29;
		int iidepth = depth-1;
		int iirbx = rbx-1;
		int iibb = bb-1;
		int iipmt = pmt-1;
		map.eta2rbxSpace[iiphi][iieta][iidepth].rbx = rbx;
		map.eta2rbxSpace[iiphi][iieta][iidepth].bb = bb;
		map.eta2rbxSpace[iiphi][iieta][iidepth].pmt = pmt;
		map.rbx2etaSpace[iirbx][iibb][iipmt].iphi = iphi;
		map.rbx2etaSpace[iirbx][iibb][iipmt].ieta = abs(ieta);
		map.rbx2etaSpace[iirbx][iibb][iipmt].depth = depth;
	}
	cout << "### Done..." << endl;
	
	return;
}










