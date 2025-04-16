#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

#include <bitset>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;

#define MAXMID 5
#define MAXEVT 100000
#define bDEBUG false 

#define bNKFADC2CH false 

map<int, int> daq_gatewidth = {{31, 2}, {41, 2}, {42, 2}};
map<int, int> daq_nch = {{1, 4}, {2, 4}, {31, 32}, {41, 32}, {42, 32}}; 
map<int, int> daq_datalength = {{1, 8192/4}, {2, 512/4}, {31, 256}, {41, 512}, {42, 512}};

//convert DAQ channel index to Geo channel index
map<int, int> GetHodoChMap(void);

//convert DAQ channel idnex to module index
map<std::pair<int,int>, vector<int>> GetCaloChMap(void);

std::map<std::pair<int,int>, vector<int>> chMapCalo;
std::map<int, int> chMapHodo;

void decode_for_ALL(int RunNo = 2080, int nEvtToRead = 10000, const char* inPath = "./25KEKDATA"){

	gInterpreter->GenerateDictionary("vector<vector<short>>", "vector");

	//check mid
	//1: NKFADC500, Gas chamber
	//2: NKFADC500, Trigger
	//3: JBNU DAQ,  
	vector<int> mid = {31, 42, 41, 1, 2};
	const int nmid = mid.size();

	chMapCalo = GetCaloChMap();
	chMapHodo = GetHodoChMap();

	cout << "*****calo map check*****" << endl;
	for (int imid=0; imid<2; imid++){
		for (int ich=0; ich<32; ich++){
			vector a = chMapCalo[std::make_pair(imid+41, ich+1)];
			cout << Form("MID %d CH %02d L/R %d MODID %02d COL %d ROW %d",imid+41, ich+1, a[0], a[1], a[2], a[2]) << endl;
		}
	}
	cout << "*****calo map check*****" << endl;

	if ( bNKFADC2CH ){
		cout << "*****decode NKFADC ONLY 2 channels, CH2 and CH3 *****" << endl;
	}

	//Output file

	int t_run_number;
	int t_tcb_trigger_number;
	ULong64_t t_tcb_trigger_time;

	int t_mid1_nch;
	int t_mid1_wlength;
	vector<vector <short>> t_mid1_adc;
	int t_mid2_nch;
	int t_mid2_wlength;
	vector<vector <short>> t_mid2_adc;

	int t_mid31_nch;
	int t_mid31_wlength;
	vector<vector <short>> t_mid31_adc;
	int t_mid41_nch;
	int t_mid41_wlength;
	vector<vector <short>> t_mid41_adc;
	int t_mid42_nch;
	int t_mid42_wlength;
	vector<vector <short>> t_mid42_adc;

	TFile* F = new TFile(Form("Run_%d_Waveform.root",RunNo), "recreate");
	TTree* T = new TTree("T", "T");
	T->Branch("run_number", &t_run_number);
	T->Branch("tcb_trigger_number", &t_tcb_trigger_number);
	T->Branch("tcb_trigger_time", &t_tcb_trigger_time);
	if ( find(mid.begin(), mid.end(), 1)!=mid.end() ){
		T->Branch("mid1_nch", &t_mid1_nch);
		T->Branch("mid1_wlength", &t_mid1_wlength);
		T->Branch("mid1_adc", &t_mid1_adc);
	}
	if ( find(mid.begin(), mid.end(), 2)!=mid.end() ){
		T->Branch("mid2_nch", &t_mid2_nch);
		T->Branch("mid2_wlength", &t_mid2_wlength);
		T->Branch("mid2_adc", &t_mid2_adc);
	}
	if ( find(mid.begin(), mid.end(), 31)!=mid.end() ){
		T->Branch("mid31_nch", &t_mid31_nch);
		T->Branch("mid31_wlength", &t_mid31_wlength);
		T->Branch("mid31_adc", &t_mid31_adc);
	}
	if ( find(mid.begin(), mid.end(), 41)!=mid.end() ){
		T->Branch("mid41_nch", &t_mid41_nch);
		T->Branch("mid41_wlength", &t_mid41_wlength);
		T->Branch("mid41_adc", &t_mid41_adc);
	}
	if ( find(mid.begin(), mid.end(), 42)!=mid.end() ){
		T->Branch("mid42_nch", &t_mid42_nch);
		T->Branch("mid42_wlength", &t_mid42_wlength);
		T->Branch("mid42_adc", &t_mid42_adc);
	}

	T->Print();

	//Input file
	char *inFile[MAXMID];
	FILE *fp[MAXMID];
	unsigned int file_size[MAXMID];
	TH1C *hEvtCount[MAXMID];
	TH1C *hEvtCount_ERR[MAXMID];
	TH1C *hEvtCount_ALL = new TH1C("hEvtCount_ALL", "hEvtCount_ALL", MAXEVT, 0, MAXEVT);

	char header[32];
	char data[10000];

	for (int imid=0; imid<nmid; imid++){

		hEvtCount[imid] = new TH1C(Form("hEvtCount_MID%d",mid[imid]), Form("hEvtCount_MID%d",mid[imid]), MAXEVT, 0, MAXEVT);
		hEvtCount_ERR[imid] = new TH1C(Form("hEvtCount_ERR_MID%d",mid[imid]), Form("hEvtCount_ERR_MID%d",mid[imid]), MAXEVT, 0, MAXEVT);

		if ( mid[imid]==1 || mid[imid]==2 ){
			inFile[imid] = Form("%s/Run_%d/Run_%d_MID_%d/FADCData_%d_%i.dat", inPath, RunNo, RunNo, mid[imid], mid[imid], RunNo);
		}else if ( mid[imid]==31 ){
			inFile[imid] = Form("%s/Run_%d/Run_%d_MID_%d/jbnu_daq_%d_%i.dat", inPath, RunNo, RunNo, mid[imid], mid[imid], RunNo);
		}else{
			inFile[imid] = Form("%s/Run_%d/Run_%d_MID_%d/bic_daq_%d_%i.dat", inPath, RunNo, RunNo, mid[imid], mid[imid], RunNo);
		}

		cout <<Form("File check for decoding %s...\n", inFile[imid]);

		//Get data file size (the size should be < 2 GB)
		fp[imid] = fopen(inFile[imid], "rb");
		fseek(fp[imid], 0L, SEEK_END);
		file_size[imid] = ftell(fp[imid]);
		//fclose(fp[imid]);
		fseek(fp[imid], 0L, SEEK_SET);

		fread(header, 1, 32, fp[imid]);
		int data_length = 0;
		if ( mid[imid]==1 || mid[imid]==2 ){
			for (int a=0; a<4; a++) data_length += ((int)(header[4*a] & 0xFF) << 8*a);
		}else{
			for (int a=0; a<4; a++) data_length += ((int)(header[a] & 0xFF) << 8*a);
		}
		cout << "File size: " << file_size[imid] << " Bytes, Data length: " << data_length << endl;

		fseek(fp[imid], 0L, SEEK_SET);

	}//imid

	map<pair<int, int>, int> map_ltn_wlength;
	map<pair<int, int>, unsigned long> map_ltn_tcbtrigtime;
	map<pair<int, int>, vector<vector<short>>> map_ltn_waveform;

	//Variables
	bool bSTOP_ALL = 0;
	bool bSTOP[MAXMID] = {0};
	unsigned int data_read[MAXMID] = {0};
	int nPacketProcessed[MAXMID] = {0};
	int trigN[MAXMID] = {0};
	unsigned long long trigT[MAXMID] = {0};

	vector<int> filled_tcb_trignum; 

	while ( !bSTOP_ALL ){

		//mid: 31, 41, 42
		for (int imid=0; imid<nmid; imid++){

			if ( bSTOP[imid] ) continue;

			const int thismid = mid[imid];
			const int thisDL = daq_datalength[thismid];
			const int thisnch = daq_nch[thismid];

			//Read 100*daq_nch packets 
			for (int ii=0; ii<(thisnch+1)*10; ii++){

				if ( thismid==1 || thismid==2 ){
/**********************************************************************************************************************/

					const int nADC = (thisDL - 32)/2; //# of ADC samples per event
					const int nTDC = nADC/4;

					fread(data, 1, thisDL*thisnch, fp[imid]); data_read[imid] += thisDL*thisnch;
					nPacketProcessed[imid] += thisnch;

					char dataChop[thisnch][thisDL];

					for (int jj=0; jj<thisnch; jj++){

						//Skip channel 1 and 4 for the 2CH mode for NKFADC500
						if ( bNKFADC2CH && (jj==0 || jj==3) ) continue;

						for (int kk=0; kk<thisDL; kk++) dataChop[jj][kk] = data[jj + 4*kk];

						int data_length = 0;
						for (int a=0; a<4; a++) data_length += ((int)(dataChop[jj][a] & 0xFF) << 8*a);

						int trigger_type = ((int)dataChop[jj][6] & 0x0F);

						int tcb_trigger_number = 0;
						for (int a=0; a<4; a++) tcb_trigger_number += ((int)(dataChop[jj][7+a] & 0xFF) << 8*a);

						int tcb_trigger_fine_time = (int)(dataChop[jj][11] & 0xFF);
						int tcb_trigger_coarse_time = 0;
						for (int a=0; a<3; a++) tcb_trigger_coarse_time += ((int)(dataChop[jj][a+12] & 0xFF) << 8*a);
						unsigned long long tcb_trigger_time = (tcb_trigger_fine_time * 8) + (tcb_trigger_coarse_time * 1000);

						int channel = ((int)dataChop[jj][16] & 0xFF);

						if ( bDEBUG ){
							cout << "NKFADC INFO: mid " << thismid 
								<< ", trigger number " << tcb_trigger_number 
								<< ", trigger time " << tcb_trigger_time 
								<< ", channel "  << channel 
								<< ", data length " << data_length 
								<< endl;
						}

						//Skip bad channel number
						if ( channel==0 || channel>thisnch ){
							cout << "NKFADC WARNNING! suspicious channel number! MID: " << thismid << ", CH:" << channel << endl; 
							if ((nPacketProcessed[imid]/(thisnch)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
								bSTOP[imid] = true;
								break;
							}
							continue;
						}

						//Skip bad trigger type 
						if ( trigger_type!=3 ){
							cout << "NKFADC WARNNING! suspicious trigger type! MID: " << thismid << ", Trigger Type:" << trigger_type << endl; 
							if ((nPacketProcessed[imid]/(thisnch)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
								bSTOP[imid] = true;
								break;
							}
							continue;
						}

						if ( fabs(tcb_trigger_number-trigN[imid])>10000 ){
							cout << "NFADC WARNNING! suspicious tcb trigger number! MID: " << thismid <<", TCB TrigN: " << tcb_trigger_number << " " << trigN[imid] 
								<< ", CH: " << channel << endl;
							if ((nPacketProcessed[imid]/(thisnch)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
								bSTOP[imid] = true;
								break;
							}
							continue;
						}

						if ( trigN[imid]!=tcb_trigger_number ){
							if ( trigT[imid]!=tcb_trigger_time ){
								trigN[imid] = tcb_trigger_number;
								trigT[imid] = tcb_trigger_time;
							}else{
								//cout << "WARNNING! different trigger number but same trigger time!" << endl;
							}
						}

						if ( hEvtCount_ALL->GetBinContent(tcb_trigger_number+1)==mid.size() ){
							cout << "NKFADC WARNNING! already filled trigger number! MID: " << thismid 
								<< ", Trigger Number: " << tcb_trigger_number 
								<< ", CH: " << channel 
								<< endl;
							if ((nPacketProcessed[imid]/(thisnch)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
								bSTOP[imid] = true;
								break;
							}
							continue;
						}

						bool bGOOD = true;
						const int wlength = (data_length - 32) / 2;

						if ( hEvtCount[imid]->GetBinContent(tcb_trigger_number+1)==0 ){

							//Initialization
							map_ltn_wlength.insert(make_pair(make_pair(thismid, tcb_trigger_number), wlength));
							map_ltn_tcbtrigtime.insert(make_pair(make_pair(thismid, tcb_trigger_number), tcb_trigger_time));
							map_ltn_waveform[make_pair(thismid, tcb_trigger_number)] = vector<vector<short>>(thisnch, vector <short>(wlength,0));

							if ( bDEBUG ){
								cout <<Form("NKFADC INFO: Initial trigger entry, MID %d, Trigger Number %d, Time %llu, ADC matrix %lu x %lu\n", 
										thismid,
										tcb_trigger_number,
										tcb_trigger_time,
										map_ltn_waveform[make_pair(thismid, tcb_trigger_number)].size(),
										map_ltn_waveform[make_pair(thismid, tcb_trigger_number)][0].size()
										);
							}
						}else{

							//Trigger time check!
							if ( tcb_trigger_time!=map_ltn_tcbtrigtime[make_pair(thismid, tcb_trigger_number)] ){
								cout << "NKFADC WARNNING! different trigger number but same trigger time!" << endl;
								cout << "MID: " << thismid 
									<< ", Trigger Number: " << tcb_trigger_number 
									<< ", Time1: " << tcb_trigger_time 
									<< ", Time2: " << map_ltn_tcbtrigtime[make_pair(thismid, tcb_trigger_number)] 
									<< endl;
								bGOOD = false;
							}
						}

						if ( bGOOD ){
							hEvtCount[imid]->Fill(tcb_trigger_number);
							//For the 2CH mode, fill the trigger number again to satisfy the check of filled channel number  
							if ( bNKFADC2CH ) hEvtCount[imid]->Fill(tcb_trigger_number);

							for (int a=0; a<wlength; a++)
							{
								const int iSmp = 32 + 2*a;
								int t_adc = (dataChop[jj][iSmp] & 0xFF) + ((dataChop[jj][iSmp + 1] & 0xF) << 8);
								map_ltn_waveform[make_pair(thismid, tcb_trigger_number)][channel-1][a] = (short)(t_adc);

								if ( bDEBUG ){
									cout << map_ltn_waveform[make_pair(thismid, tcb_trigger_number)][channel-1][a] << " ";
								}
							}

							if ( bDEBUG )	cout << endl;
						}//bGOOD

						//Check whether all channels are collected for a given trigger number
						if ( bGOOD && hEvtCount[imid]->GetBinContent(tcb_trigger_number+1)==thisnch )
						{
							hEvtCount_ALL->Fill(tcb_trigger_number);

							if ( bDEBUG ){
								cout <<Form("NKFADC INFO: Collect all channels, MID %d, Trigger Number %d\n", 
										thismid,
										tcb_trigger_number 
										);
							}

							if ( hEvtCount_ALL->GetBinContent(tcb_trigger_number+1)==mid.size() ){
								filled_tcb_trignum.push_back(tcb_trigger_number);
							}
						}//

					}//jj - ch


					if (nPacketProcessed[imid]%10000 == 0) cout << "MID" << thismid << " Processed eventNum = " << nPacketProcessed[imid]/(thisnch+1) << endl;
					if ((nPacketProcessed[imid]/(thisnch)) >= nEvtToRead
							|| data_read[imid] >= file_size[imid]
						 ){
						bSTOP[imid] = true;
						break;
					}

/**********************************************************************************************************************/
				}else{
/**********************************************************************************************************************/

					//Read header
					//++++++++++++++++++++++++++++++++++++++++++++
					fread(header, 1, 32, fp[imid]);

					int data_length = 0;
					for (int a=0; a<4; a++) data_length += ((int)(header[a] & 0xFF) << 8*a);

					/*
						 int run_number = 0;
						 for (int a=0; a<2; a++) run_number += ((int)(header[a+4] & 0xFF) << 8*a);
						 int trigger_type = ((int)header[6] & 0xFF);
					 */

					int tcb_trigger_number = 0;
					for (int a=0; a<4; a++) tcb_trigger_number += ((int)(header[a+7] & 0xFF) << 8*a);

					int tcb_trigger_fine_time = ((int)header[11] & 0xFF);
					int tcb_trigger_coarse_time = 0;
					for (int a=0; a<3; a++) tcb_trigger_coarse_time += ((int)(header[a+12] & 0xFF) << 8*a);
					unsigned long long tcb_trigger_time = (tcb_trigger_fine_time * 8) + (tcb_trigger_coarse_time * 1000);

					//int mid = ((int)header[15] & 0xFF);
					int channel = ((int)header[16] & 0xFF);

					if ( bDEBUG ){
						cout << "BICDAQ INFO: mid " << thismid 
							<< ", trigger number " << tcb_trigger_number 
							<< ", trigger time " << tcb_trigger_time 
							<< ", channel "  << channel 
							<< ", data length " << data_length 
							<< endl;
					}

					//Skip bad channel number
					if ( channel>thisnch ){
						cout << "BICDAQ WARNNING! suspicious channel number! MID: " << thismid << ", CH:" << channel << endl; 
						fread(data, 1, thisDL - 32, fp[imid]); data_read[imid] += thisDL;
						nPacketProcessed[imid]++;
						if ((nPacketProcessed[imid]/(thisnch+1)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
							bSTOP[imid] = true;
							break;
						}
						continue;
					}

					//Skip channel 0 (fast analysis data) for JBNU and BIC DAQ
					if ( channel==0 ){
						fread(data, 1, 256 - 32, fp[imid]); data_read[imid] += 256;
						nPacketProcessed[imid]++;
						if ((nPacketProcessed[imid]/(thisnch+1)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
							bSTOP[imid] = true;
							break;
						}
						continue;
					}

					if ( fabs(tcb_trigger_number-trigN[imid])>10000 ){
						cout << "BICDAQ WARNNING! suspicious tcb trigger number! MID: " << thismid <<", TCB TrigN: " << tcb_trigger_number << " " << trigN[imid] 
							<< ", CH: " << channel << endl;
						fread(data, 1, thisDL - 32, fp[imid]); data_read[imid] += thisDL;
						nPacketProcessed[imid]++;
						if ((nPacketProcessed[imid]/(thisnch+1)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
							bSTOP[imid] = true;
							break;
						}
						continue;
					}

					if ( trigN[imid]!=tcb_trigger_number ){
						if ( trigT[imid]!=tcb_trigger_time ){
							trigN[imid] = tcb_trigger_number;
							trigT[imid] = tcb_trigger_time;
						}else{
							cout << "WARNNING! different trigger number but same trigger time!" << endl;
						}
					}

					if ( data_length!=thisDL ){
						fread(data, 1, thisDL - 32, fp[imid]); data_read[imid] += thisDL;
						cout << "BICDAQ WARNNING! suspicious data length! MID: " << thismid << ", DLength: " << data_length << ", CH: " << channel << endl;
						nPacketProcessed[imid]++;
						if ((nPacketProcessed[imid]/(thisnch+1)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
							bSTOP[imid] = true;
							break;
						}
						continue;
					}

					if ( hEvtCount_ALL->GetBinContent(tcb_trigger_number+1)==mid.size() ){
						cout << "BICDAQ WARNNING! already filled trigger number! MID: " 
							<< thismid << ", Trigger Number: " 
							<< tcb_trigger_number << ", CH: " 
							<< channel 
							<< endl;
						fread(data, 1, data_length - 32, fp[imid]); data_read[imid] += data_length;
						nPacketProcessed[imid]++;
						hEvtCount_ERR[imid]->Fill(tcb_trigger_number);
						if ((nPacketProcessed[imid]/(thisnch+1)) >= nEvtToRead || data_read[imid] >= file_size[imid]){
							bSTOP[imid] = true;
							break;
						}
						continue;
					}

					//If you want to debug from a certain trigger number
#if 0
					if ( tcb_trigger_number<9980 ){
						fread(data, 1, data_length - 32, fp[imid]); data_read[imid] += data_length;
						nPacketProcessed[imid]++;
						if ((nPacketProcessed[imid]/(thisnch+1)) == nEvtToRead || data_read[imid] >= file_size[imid]){
							bSTOP[imid] = true;
							break;
						}
						continue;
					}
#endif

					//Read body, data_length - 32 bytes (header)
					//++++++++++++++++++++++++++++++++++++++++++++
					fread(data, 1, data_length - 32, fp[imid]);	data_read[imid] += data_length;
					nPacketProcessed[imid]++;

					bool bGOOD = true;
					const int wlength = (thismid==31) ? (data_length - 32) / 2 : (data_length - 32) / 4;

					if ( hEvtCount[imid]->GetBinContent(tcb_trigger_number+1)==0 ){

						//Initialization
						map_ltn_wlength.insert(make_pair(make_pair(thismid, tcb_trigger_number), wlength));
						map_ltn_tcbtrigtime.insert(make_pair(make_pair(thismid, tcb_trigger_number), tcb_trigger_time));
						map_ltn_waveform[make_pair(thismid, tcb_trigger_number)] = vector<vector<short>>(thisnch, vector <short>(wlength,0));

						if ( bDEBUG ){
							cout <<Form("BICDAQ INFO: Initial trigger entry, MID %d, Trigger Number %d, Time %llu, ADC matrix %lu x %lu\n", 
									thismid,
									tcb_trigger_number,
									tcb_trigger_time,
									map_ltn_waveform[make_pair(thismid, tcb_trigger_number)].size(),
									map_ltn_waveform[make_pair(thismid, tcb_trigger_number)][0].size()
									);
						}
					}else{

						//Trigger time check!
						if ( tcb_trigger_time!=map_ltn_tcbtrigtime[make_pair(thismid, tcb_trigger_number)] ){
							cout << "BICDAQ WARNNING! different trigger number but same trigger time!";
							cout << " MID: " << thismid 
								<< ", CH: " << channel
								<< ", Trigger Number: " << tcb_trigger_number 
								<< ", Time1: " << tcb_trigger_time 
								<< ", Time2: " << map_ltn_tcbtrigtime[make_pair(thismid, tcb_trigger_number)] 
								<< endl;
							bGOOD = false;
							hEvtCount_ERR[imid]->Fill(tcb_trigger_number);
						}
					}

					if ( bGOOD ){
						hEvtCount[imid]->Fill(tcb_trigger_number);

						if ( thismid==31 ){
							for (int a=0; a<wlength; a++)
							{
								int t_adc1 = (data[2 * a + 0] & 0xFF);
								int t_adc2 = (data[2 * a + 1] & 0xFF) << 8;
								map_ltn_waveform[make_pair(thismid, tcb_trigger_number)][channel-1][a] = (short)(t_adc1 + t_adc2);
								if ( bDEBUG ){
									cout << map_ltn_waveform[make_pair(thismid, tcb_trigger_number)][channel-1][a] << " ";
								}
							}//a
						}else if ( thismid==41 || thismid==42 ){
							for (int a=0; a<wlength; a++)
							{
								int t_adc1 = (data[4 * a + 0] & 0xFF);
								int t_adc2 = (data[4 * a + 1] & 0xFF) << 8;
								map_ltn_waveform[make_pair(thismid, tcb_trigger_number)][channel-1][a] = (short)(t_adc1 + t_adc2);
								if ( bDEBUG ){
									cout << map_ltn_waveform[make_pair(thismid, tcb_trigger_number)][channel-1][a] << " ";
								}
							}
						}

						if ( bDEBUG )	cout << endl;
					}//bGOOD


					//Check whether all channels are collected for a given trigger number
					if ( bGOOD && hEvtCount[imid]->GetBinContent(tcb_trigger_number+1)==thisnch )
					{
						hEvtCount_ALL->Fill(tcb_trigger_number);

						if ( bDEBUG ){
							cout <<Form("BICDAQ INFO: Collect all channels, MID %d, Trigger Number %d\n", 
									thismid,
									tcb_trigger_number 
									);
						}

						if ( hEvtCount_ALL->GetBinContent(tcb_trigger_number+1)==mid.size() ){
							filled_tcb_trignum.push_back(tcb_trigger_number);
						}
					}//

					//++++++++++++++++++++++++++++++++++++++++++++

					if (nPacketProcessed[imid]%10000 == 0) cout << "MID" << thismid << " Processed eventNum = " << nPacketProcessed[imid]/(thisnch+1) << endl;
					if ((nPacketProcessed[imid]/(thisnch+1)) == nEvtToRead
							|| data_read[imid] >= file_size[imid]
						 ){
						bSTOP[imid] = true;
						break;
					}
/**********************************************************************************************************************/
				}//mid check


			}//ii

		}//imid

		//Fill TTree
		for (int ii=0; ii<filled_tcb_trignum.size(); ii++){
			t_run_number = RunNo; 
			t_tcb_trigger_number = filled_tcb_trignum[ii];
			t_tcb_trigger_time = map_ltn_tcbtrigtime[make_pair(mid[0], t_tcb_trigger_number)]; 

			bool bGOODTT = true;

			for (int imid=0; imid<nmid; imid++){

				if ( t_tcb_trigger_time!=map_ltn_tcbtrigtime[make_pair(mid[imid], t_tcb_trigger_number)] ){
					cout << "WARNNING AT FILL! same trigger number but differnt trigger time!" 
						<< " Trigger Number: " << t_tcb_trigger_number 
						<< ", MID: " << mid[0] 
						<< ", Time: " << t_tcb_trigger_time 
						<< ", MID: " << mid[imid] 
						<< ", Time: " << map_ltn_tcbtrigtime[make_pair(mid[imid], t_tcb_trigger_number)] 
						<< endl;
					bGOODTT = false;
				}

				if ( mid[imid]==1 ){
					t_mid1_nch = daq_nch[mid[imid]];
					t_mid1_wlength = map_ltn_wlength[make_pair(mid[imid], t_tcb_trigger_number)];
					t_mid1_adc = map_ltn_waveform[make_pair(mid[imid], t_tcb_trigger_number)];
				}else if ( mid[imid]==2 ){
					t_mid2_nch = daq_nch[mid[imid]];
					t_mid2_wlength = map_ltn_wlength[make_pair(mid[imid], t_tcb_trigger_number)];
					t_mid2_adc = map_ltn_waveform[make_pair(mid[imid], t_tcb_trigger_number)];
				}else if ( mid[imid]==31 ){
					t_mid31_nch = daq_nch[mid[imid]];
					t_mid31_wlength = map_ltn_wlength[make_pair(mid[imid], t_tcb_trigger_number)];
					t_mid31_adc = map_ltn_waveform[make_pair(mid[imid], t_tcb_trigger_number)];
				}else if ( mid[imid]==41 ){
					t_mid41_nch = daq_nch[mid[imid]];
					t_mid41_wlength = map_ltn_wlength[make_pair(mid[imid], t_tcb_trigger_number)];
					t_mid41_adc = map_ltn_waveform[make_pair(mid[imid], t_tcb_trigger_number)];
				}else if ( mid[imid]==42 ){
					t_mid42_nch = daq_nch[mid[imid]];
					t_mid42_wlength = map_ltn_wlength[make_pair(mid[imid], t_tcb_trigger_number)];
					t_mid42_adc = map_ltn_waveform[make_pair(mid[imid], t_tcb_trigger_number)];
				}

				//map_ltn_wlength.erase(make_pair(mid[imid], t_tcb_trigger_number));
				//map_ltn_tcbtrigtime.erase(make_pair(mid[imid], t_tcb_trigger_number));
				map_ltn_waveform.erase(make_pair(mid[imid], t_tcb_trigger_number));

			}//imid

			if ( bGOODTT ){
				T->Fill();
			}

			if ( bDEBUG ){
				cout << "Fill TTree, Trigger Number " << t_tcb_trigger_number << ", Time " << t_tcb_trigger_time << endl;
			}

		}//ii

		filled_tcb_trignum.clear();

		int bSTOP_SUM = 0;
		for (int ii=0; ii<nmid; ii++) bSTOP_SUM += bSTOP[ii];
		
		if ( bSTOP_SUM==nmid ) bSTOP_ALL = true;

	}//while

	F->cd();
	hEvtCount_ALL->Write();
	for (int imid=0; imid<nmid; imid++){
		hEvtCount[imid]->Write();
		hEvtCount_ERR[imid]->Write();
	}
	T->Write();

	F->Close();
	
	return;
}

//___________________________________________________________________________________________//
//convert DAQ channel index to geometry index
map<int, int> GetHodoChMap(void)
{
    std::map<int, int> chMap;

    chMap.insert( std::pair<int, int> ( 1,  1) );
    chMap.insert( std::pair<int, int> ( 2,  2) );
    chMap.insert( std::pair<int, int> ( 3,  3) );
    chMap.insert( std::pair<int, int> ( 4,  4) );
    chMap.insert( std::pair<int, int> ( 5,  8) );
    chMap.insert( std::pair<int, int> ( 6,  7) );
    chMap.insert( std::pair<int, int> ( 7,  6) );
    chMap.insert( std::pair<int, int> ( 8,  5) );
    chMap.insert( std::pair<int, int> ( 9, 12) );
    chMap.insert( std::pair<int, int> (10, 11) );
    chMap.insert( std::pair<int, int> (11, 10) );
    chMap.insert( std::pair<int, int> (12,  9) );
    chMap.insert( std::pair<int, int> (13, 13) );
    chMap.insert( std::pair<int, int> (14, 14) );
    chMap.insert( std::pair<int, int> (15, 15) );
    chMap.insert( std::pair<int, int> (16, 16) );

    chMap.insert( std::pair<int, int> ( 1+16,  1+16) );
    chMap.insert( std::pair<int, int> ( 2+16,  2+16) );
    chMap.insert( std::pair<int, int> ( 3+16,  3+16) );
    chMap.insert( std::pair<int, int> ( 4+16,  4+16) );
    chMap.insert( std::pair<int, int> ( 5+16,  8+16) );
    chMap.insert( std::pair<int, int> ( 6+16,  7+16) );
    chMap.insert( std::pair<int, int> ( 7+16,  6+16) );
    chMap.insert( std::pair<int, int> ( 8+16,  5+16) );
    chMap.insert( std::pair<int, int> ( 9+16, 12+16) );
    chMap.insert( std::pair<int, int> (10+16, 11+16) );
    chMap.insert( std::pair<int, int> (11+16, 10+16) );
    chMap.insert( std::pair<int, int> (12+16,  9+16) );
    chMap.insert( std::pair<int, int> (13+16, 13+16) );
    chMap.insert( std::pair<int, int> (14+16, 14+16) );
    chMap.insert( std::pair<int, int> (15+16, 15+16) );
    chMap.insert( std::pair<int, int> (16+16, 16+16) );

    return chMap;
}//map

//___________________________________________________________________________________________//
//convert DAQ channel index to module/geometry index 
map<pair<int,int>, vector<int>> GetCaloChMap(void)
{

	std::map<std::pair<int, int>, vector<int>> chMap;

	//1st value: left 0, right 1 
	//2nd value: module ID 
	//3rd value: column number 0-7 //0:upstream, 7:downstream
	//4th value: row number 0-3 //0:1st floor, 3:4th floor
	vector<int> mid41ch01 = {0,  1, 7, 0}; //M1L
	vector<int> mid41ch02 = {1,  1, 7, 0}; //M1R
	vector<int> mid41ch03 = {0,  6, 7, 1}; //M6L
	vector<int> mid41ch04 = {1,  6, 7, 1}; //M6R
	vector<int> mid41ch05 = {1,  5, 4, 0}; //M5R
	vector<int> mid41ch06 = {1, 26, 4, 1}; //M26R
	vector<int> mid41ch07 = {1, 30, 4, 2}; //M30R
	vector<int> mid41ch08 = {1, 29, 4, 3}; //M29R
	vector<int> mid41ch09 = {1,  3, 5, 0}; //M3R
	vector<int> mid41ch10 = {1, 18, 5, 1}; //M18R
	vector<int> mid41ch11 = {1, 17, 5, 2}; //M17R
	vector<int> mid41ch12 = {1, 31, 5, 3}; //M31R
	vector<int> mid41ch13 = {1,  2, 6, 0}; //M2R
	vector<int> mid41ch14 = {1, 27, 6, 1}; //M27R
	vector<int> mid41ch15 = {1, 22, 6, 2}; //M22R
	vector<int> mid41ch16 = {1, 28, 6, 3}; //M28R
	vector<int> mid41ch17 = {0,  7, 7, 2}; //M7L
	vector<int> mid41ch18 = {1,  7, 7, 2}; //M7R
	vector<int> mid41ch19 = {0, 15, 7, 3}; //M15L
	vector<int> mid41ch20 = {1, 15, 7, 3}; //M15R
	vector<int> mid41ch21 = {0,  5, 4, 0}; //M5L
	vector<int> mid41ch22 = {0, 26, 4, 1}; //M26L
	vector<int> mid41ch23 = {0, 30, 4, 2}; //M30L
	vector<int> mid41ch24 = {0, 29, 4, 3}; //M29L
	vector<int> mid41ch25 = {0,  3, 5, 0}; //M3L
	vector<int> mid41ch26 = {0, 18, 5, 1}; //M18L
	vector<int> mid41ch27 = {0, 17, 5, 2}; //M17L
	vector<int> mid41ch28 = {0, 31, 5, 3}; //M31L
	vector<int> mid41ch29 = {0,  2, 6, 0}; //M2L
	vector<int> mid41ch30 = {0, 27, 6, 1}; //M27L
	vector<int> mid41ch31 = {0, 22, 6, 2}; //M22L
	vector<int> mid41ch32 = {0, 28, 6, 3}; //M28L

	vector<int> mid42ch01 = {1,  8, 0, 0}; //M8R
	vector<int> mid42ch02 = {1,  9, 0, 1}; //M9R
	vector<int> mid42ch03 = {1, 24, 0, 2}; //M24R
	vector<int> mid42ch04 = {1, 12, 0, 3}; //M12R
	vector<int> mid42ch05 = {1, 14, 1, 0}; //M14R
	vector<int> mid42ch06 = {1, 13, 1, 1}; //M13R
	vector<int> mid42ch07 = {1, 21, 1, 2}; //M21R
	vector<int> mid42ch08 = {1, 10, 1, 3}; //M10R
	vector<int> mid42ch09 = {1,  4, 2, 0}; //M4R
	vector<int> mid42ch10 = {1, 25, 2, 1}; //M25R
	vector<int> mid42ch11 = {1, 23, 2, 2}; //M23R
	vector<int> mid42ch12 = {1, 32, 2, 3}; //M32R
	vector<int> mid42ch13 = {1, 11, 3, 0}; //M11R
	vector<int> mid42ch14 = {1, 33, 3, 1}; //M33R
	vector<int> mid42ch15 = {1, 20, 3, 2}; //M20R
	vector<int> mid42ch16 = {1, 19, 3, 3}; //M19R
	vector<int> mid42ch17 = {0,  8, 0, 0}; //M8L
	vector<int> mid42ch18 = {0,  9, 0, 1}; //M9L
	vector<int> mid42ch19 = {0, 24, 0, 2}; //M24L
	vector<int> mid42ch20 = {0, 12, 0, 3}; //M12L
	vector<int> mid42ch21 = {0, 14, 1, 0}; //M14L
	vector<int> mid42ch22 = {0, 13, 1, 1}; //M13L
	vector<int> mid42ch23 = {0, 21, 1, 2}; //M21L
	vector<int> mid42ch24 = {0, 10, 1, 3}; //M10L
	vector<int> mid42ch25 = {0,  4, 2, 0}; //M4L
	vector<int> mid42ch26 = {0, 25, 2, 1}; //M25L
	vector<int> mid42ch27 = {0, 23, 2, 2}; //M23L
	vector<int> mid42ch28 = {0, 32, 2, 3}; //M32L
	vector<int> mid42ch29 = {0, 11, 3, 0}; //M11L
	vector<int> mid42ch30 = {0, 33, 3, 1}; //M33L
	vector<int> mid42ch31 = {0, 20, 3, 2}; //M20L
	vector<int> mid42ch32 = {0, 19, 3, 3}; //M19L

	chMap.insert( std::make_pair(std::make_pair(41, 1), mid41ch01) );
	chMap.insert( std::make_pair(std::make_pair(41, 2), mid41ch02) );
	chMap.insert( std::make_pair(std::make_pair(41, 3), mid41ch03) );
	chMap.insert( std::make_pair(std::make_pair(41, 4), mid41ch04) );
	chMap.insert( std::make_pair(std::make_pair(41, 5), mid41ch05) );
	chMap.insert( std::make_pair(std::make_pair(41, 6), mid41ch06) );
	chMap.insert( std::make_pair(std::make_pair(41, 7), mid41ch07) );
	chMap.insert( std::make_pair(std::make_pair(41, 8), mid41ch08) );
	chMap.insert( std::make_pair(std::make_pair(41, 9), mid41ch09) );
	chMap.insert( std::make_pair(std::make_pair(41,10), mid41ch10) );
	chMap.insert( std::make_pair(std::make_pair(41,11), mid41ch11) );
	chMap.insert( std::make_pair(std::make_pair(41,12), mid41ch12) );
	chMap.insert( std::make_pair(std::make_pair(41,13), mid41ch13) );
	chMap.insert( std::make_pair(std::make_pair(41,14), mid41ch14) );
	chMap.insert( std::make_pair(std::make_pair(41,15), mid41ch15) );
	chMap.insert( std::make_pair(std::make_pair(41,16), mid41ch16) );
	chMap.insert( std::make_pair(std::make_pair(41,17), mid41ch17) );
	chMap.insert( std::make_pair(std::make_pair(41,18), mid41ch18) );
	chMap.insert( std::make_pair(std::make_pair(41,19), mid41ch19) );
	chMap.insert( std::make_pair(std::make_pair(41,20), mid41ch20) );
	chMap.insert( std::make_pair(std::make_pair(41,21), mid41ch21) );
	chMap.insert( std::make_pair(std::make_pair(41,22), mid41ch22) );
	chMap.insert( std::make_pair(std::make_pair(41,23), mid41ch23) );
	chMap.insert( std::make_pair(std::make_pair(41,24), mid41ch24) );
	chMap.insert( std::make_pair(std::make_pair(41,25), mid41ch25) );
	chMap.insert( std::make_pair(std::make_pair(41,26), mid41ch26) );
	chMap.insert( std::make_pair(std::make_pair(41,27), mid41ch27) );
	chMap.insert( std::make_pair(std::make_pair(41,28), mid41ch28) );
	chMap.insert( std::make_pair(std::make_pair(41,29), mid41ch29) );
	chMap.insert( std::make_pair(std::make_pair(41,30), mid41ch30) );
	chMap.insert( std::make_pair(std::make_pair(41,31), mid41ch31) );
	chMap.insert( std::make_pair(std::make_pair(41,32), mid41ch32) );

	chMap.insert( std::make_pair(std::make_pair(42, 1), mid42ch01) );
	chMap.insert( std::make_pair(std::make_pair(42, 2), mid42ch02) );
	chMap.insert( std::make_pair(std::make_pair(42, 3), mid42ch03) );
	chMap.insert( std::make_pair(std::make_pair(42, 4), mid42ch04) );
	chMap.insert( std::make_pair(std::make_pair(42, 5), mid42ch05) );
	chMap.insert( std::make_pair(std::make_pair(42, 6), mid42ch06) );
	chMap.insert( std::make_pair(std::make_pair(42, 7), mid42ch07) );
	chMap.insert( std::make_pair(std::make_pair(42, 8), mid42ch08) );
	chMap.insert( std::make_pair(std::make_pair(42, 9), mid42ch09) );
	chMap.insert( std::make_pair(std::make_pair(42,10), mid42ch10) );
	chMap.insert( std::make_pair(std::make_pair(42,11), mid42ch11) );
	chMap.insert( std::make_pair(std::make_pair(42,12), mid42ch12) );
	chMap.insert( std::make_pair(std::make_pair(42,13), mid42ch13) );
	chMap.insert( std::make_pair(std::make_pair(42,14), mid42ch14) );
	chMap.insert( std::make_pair(std::make_pair(42,15), mid42ch15) );
	chMap.insert( std::make_pair(std::make_pair(42,16), mid42ch16) );
	chMap.insert( std::make_pair(std::make_pair(42,17), mid42ch17) );
	chMap.insert( std::make_pair(std::make_pair(42,18), mid42ch18) );
	chMap.insert( std::make_pair(std::make_pair(42,19), mid42ch19) );
	chMap.insert( std::make_pair(std::make_pair(42,20), mid42ch20) );
	chMap.insert( std::make_pair(std::make_pair(42,21), mid42ch21) );
	chMap.insert( std::make_pair(std::make_pair(42,22), mid42ch22) );
	chMap.insert( std::make_pair(std::make_pair(42,23), mid42ch23) );
	chMap.insert( std::make_pair(std::make_pair(42,24), mid42ch24) );
	chMap.insert( std::make_pair(std::make_pair(42,25), mid42ch25) );
	chMap.insert( std::make_pair(std::make_pair(42,26), mid42ch26) );
	chMap.insert( std::make_pair(std::make_pair(42,27), mid42ch27) );
	chMap.insert( std::make_pair(std::make_pair(42,28), mid42ch28) );
	chMap.insert( std::make_pair(std::make_pair(42,29), mid42ch29) );
	chMap.insert( std::make_pair(std::make_pair(42,30), mid42ch30) );
	chMap.insert( std::make_pair(std::make_pair(42,31), mid42ch31) );
	chMap.insert( std::make_pair(std::make_pair(42,32), mid42ch32) );

	return chMap;
}

