







#define MAXSWAPS 50

#include <iostream>
#include <string>

using namespace std;


//
//	Define structs to simplify the workflow with Swaps/Src Errors 
//
struct TChSwap
{
	int o_iphi, o_ieta, o_depth;
	int f_iphi, f_ieta, f_depth;
};

struct TSrcError
{
	int iphi, ieta;
	int q1TS;
	int q2TS;
};

typedef TChSwap TChSwaps[MAXSWAPS];
typedef TSrcError TSrcErrors[MAXSWAPS];

struct TSwaps
{
	int numChSwaps;
	int numTubeSwaps;
	int numSrcErrors;
	TChSwaps chSwaps;
	TChSwaps tubeSwaps;
	TSrcErrors srcErrors;
};

int isTubeSwap(TSwaps &swaps, int iphi, int ieta)
{
	int iswap = -1;
	for (int i=0; i<swaps.numTubeSwaps; i++)
		if (iphi==swaps.tubeSwaps[i].o_iphi && ieta==swaps.tubeSwaps[i].o_ieta)
		{
			iswap = i;	
			break;
		}

	return iswap;
}

int isChSwap(TSwaps &swaps, int iphi, int ieta, int depth)
{
	int iswap = -1;
	for (int i=0; i<swaps.numChSwaps; i++)
		if (iphi==swaps.chSwaps[i].o_iphi && ieta==swaps.chSwaps[i].o_ieta
				&& depth==swaps.chSwaps[i].o_depth)
		{
			iswap = i;
			break;
		}

	return iswap;
}

int isSrcError(TSwaps &swaps, int iphi, int ieta, int iTS)
{
	int iswap = 0;
	for (int i=0; i<swaps.numSrcErrors; i++)
		if (iphi==swaps.srcErrors[i].iphi && ieta==swaps.srcErrors[i].ieta)
		{
			if (iTS==1)
				iswap = swaps.srcErrors[i].q1TS;
			else if (iTS==2)
				iswap = swaps.srcErrors[i].q2TS;
		}

	return iswap;
}

void readSwaps(string inSwapsFileName, TSwaps &swaps, int verbosity=0)
{
	ifstream in(inSwapsFileName.c_str());
	int numChSwaps, numTubeSwaps, numSrcErrors;
	int o_iphi, o_ieta, o_depth, f_iphi, f_ieta, f_depth;
	int q1TS, q2TS;

	//
	//	Parse the Swaps file
	//
	in >> numChSwaps;
	swaps.numChSwaps = numChSwaps;
	for (int iChSwap=0; iChSwap<numChSwaps; iChSwap++)
	{
		in >> o_iphi >> o_ieta >> o_depth >> f_iphi >> f_ieta >> f_depth;
		swaps.chSwaps[iChSwap].o_iphi = o_iphi;
		swaps.chSwaps[iChSwap].o_ieta = o_ieta;
		swaps.chSwaps[iChSwap].o_depth = o_depth;
		swaps.chSwaps[iChSwap].f_iphi = f_iphi;
		swaps.chSwaps[iChSwap].f_ieta = f_ieta;
		swaps.chSwaps[iChSwap].f_depth = f_depth;
	}
	in >> numTubeSwaps;
	swaps.numTubeSwaps = numTubeSwaps;
	for (int iTubeSwap=0; iTubeSwap<numTubeSwaps; iTubeSwap++)
	{
		in >> o_iphi >> o_ieta >> f_iphi >> f_ieta;
		swaps.tubeSwaps[iTubeSwap].o_iphi = o_iphi;
		swaps.tubeSwaps[iTubeSwap].o_ieta = o_ieta;
		swaps.tubeSwaps[iTubeSwap].f_iphi = f_iphi;
		swaps.tubeSwaps[iTubeSwap].f_ieta = f_ieta;
	}
	in >> numSrcErrors;
	swaps.numSrcErrors = numSrcErrors;
	for (int iErr=0; iErr<numSrcErrors; iErr++)
	{
		in >>  o_iphi >> o_ieta >> q1TS >> q2TS;
		swaps.srcErrors[iErr].iphi = o_iphi;
		swaps.srcErrors[iErr].ieta = o_ieta;
		swaps.srcErrors[iErr].q1TS = q1TS;
		swaps.srcErrors[iErr].q2TS = q2TS;
	}

	if (verbosity>1)
	{
		cout << "### Channel Swaps: " << numChSwaps<< endl;
		for (int i=0; i<numChSwaps; i++)
			cout << swaps.chSwaps[i].o_iphi << "  " << swaps.chSwaps[i].o_ieta
				<< "  " << swaps.chSwaps[i].o_depth << "  ->  " 
				<< swaps.chSwaps[i].f_iphi
				<< "  " << swaps.chSwaps[i].f_ieta << "  " 
				<< swaps.chSwaps[i].f_depth 
				<< endl;
		cout << "### Tube Swaps: " << numTubeSwaps << endl;
		for (int i=0; i<numTubeSwaps; i++)
			cout << swaps.tubeSwaps[i].o_iphi << "  " << swaps.tubeSwaps[i].o_ieta
				<< "  ->  " << swaps.tubeSwaps[i].f_iphi << "  " 
				<< swaps.tubeSwaps[i].f_ieta
				<< endl;
		cout << "### Src Errors: " << numSrcErrors << endl;
		for (int i=0; i<numSrcErrors; i++)
			cout << swaps.srcErrors[i].iphi << "  " << swaps.srcErrors[i].ieta 
				<< "  "
				<< swaps.srcErrors[i].q1TS << "  " << swaps.srcErrors[i].q2TS
				<< endl;
	}

	return;
}




































