#include <algorithm>

using namespace std;

map<int, int> wform_factor = {{31, 150}, {41, 1800}, {42, 1800}};
map<int, int> cutPed = {{1, 20}, {2, 20}, {31, 20}, {41, 20}, {42, 20}};
const float xCVS = 1.0; //Adjust canvas size

//convert DAQ channel index to Geo channel index
map<int, int> GetHodoChMap(void);

//convert DAQ channel idnex to module index
map<std::pair<int,int>, vector<int>> GetCaloChMap(void);

std::map<std::pair<int,int>, vector<int>> chMapCalo;
std::map<int, int> chMapHodo;

void readback_for_ALL(int RunNo = 61130, int nEvtToRead=1000){

	gInterpreter->GenerateDictionary("vector<vector<short>>", "vector");

	chMapCalo = GetCaloChMap();
	chMapHodo = GetHodoChMap();

	//check mid
	//1: NKFADC500, Gas chamber
	//2: NKFADC500, Trigger
	//3: JBNU DAQ,  
	vector<int> mid = {1, 2, 31, 41, 42};

	int t_run_number;
	int t_tcb_trigger_number;
	ULong64_t t_tcb_trigger_time;
	int t_mid1_nch;
	int t_mid1_wlength = 1;
	vector<vector <short>> *t_mid1_adc = 0;
	int t_mid2_nch;
	int t_mid2_wlength = 1;
	vector<vector <short>> *t_mid2_adc = 0;
	int t_mid31_nch;
	int t_mid31_wlength = 1;
	vector<vector <short>> *t_mid31_adc = 0;
	int t_mid41_nch;
	int t_mid41_wlength = 1;
	vector<vector <short>> *t_mid41_adc = 0;
	int t_mid42_nch;
	int t_mid42_wlength = 1;
	vector<vector <short>> *t_mid42_adc = 0;

	TFile* F = new TFile(Form("Run_%d_Waveform.root",RunNo), "read");
	TTree* T = (TTree*)F->Get("T"); 

	T->SetBranchAddress("run_number", &t_run_number);
	T->SetBranchAddress("tcb_trigger_number", &t_tcb_trigger_number);
	T->SetBranchAddress("tcb_trigger_time", &t_tcb_trigger_time);
	if ( find(mid.begin(), mid.end(), 1)!=mid.end() ){
		T->SetBranchAddress("mid1_nch", &t_mid1_nch);
		T->SetBranchAddress("mid1_wlength", &t_mid1_wlength);
		T->SetBranchAddress("mid1_adc", &t_mid1_adc);
	}
	if ( find(mid.begin(), mid.end(), 2)!=mid.end() ){
		T->SetBranchAddress("mid2_nch", &t_mid2_nch);
		T->SetBranchAddress("mid2_wlength", &t_mid2_wlength);
		T->SetBranchAddress("mid2_adc", &t_mid2_adc);
	}
	if ( find(mid.begin(), mid.end(), 31)!=mid.end() ){
		T->SetBranchAddress("mid31_nch", &t_mid31_nch);
		T->SetBranchAddress("mid31_wlength", &t_mid31_wlength);
		T->SetBranchAddress("mid31_adc", &t_mid31_adc);
	}
	if ( find(mid.begin(), mid.end(), 41)!=mid.end() ){
		T->SetBranchAddress("mid41_nch", &t_mid41_nch);
		T->SetBranchAddress("mid41_wlength", &t_mid41_wlength);
		T->SetBranchAddress("mid41_adc", &t_mid41_adc);
	}
	if ( find(mid.begin(), mid.end(), 42)!=mid.end() ){
		T->SetBranchAddress("mid42_nch", &t_mid42_nch);
		T->SetBranchAddress("mid42_wlength", &t_mid42_wlength);
		T->SetBranchAddress("mid42_adc", &t_mid42_adc);
	}

	int nentries = T->GetEntries();

	T->GetEntry(nentries-1);
	int maxtrignum = int(1.05*t_tcb_trigger_number);
	cout << "Maximum Trigger Number: " << t_tcb_trigger_number << endl;

	TH2D *hTrigNumTime = new TH2D("hTrigNumTime", Form("; Trigger Number; Trigger Time"), maxtrignum,0,maxtrignum, 200,0,0.5*1.E10);
	TH2D *hPulse_mid1[4];
	TH2D *hPulse_mid2[4];
	TH2D *hPulse_mid31[32];
	TH2D *hPulse_mid41[32];
	TH2D *hPulse_mid42[32];

	for (int ich=0; ich<4; ich++){
		hPulse_mid1[ich] = new TH2D(Form("hPulse_mid1_ch%d",1+ich),"",t_mid1_wlength,0,t_mid1_wlength,3000,0,6000);
		hPulse_mid2[ich] = new TH2D(Form("hPulse_mid2_ch%d",1+ich),"",t_mid2_wlength,0,t_mid2_wlength,3000,0,6000);
	}

	for (int ich=0; ich<32; ich++){
		hPulse_mid31[ich] = new TH2D(Form("hPulse_mid31_ch%d",1+ich),"",t_mid31_wlength,0,t_mid31_wlength,3000,-100,5900);
		hPulse_mid41[ich] = new TH2D(Form("hPulse_mid41_ch%d",1+ich),"",t_mid41_wlength,0,t_mid41_wlength,3000,-100,14900);
		hPulse_mid42[ich] = new TH2D(Form("hPulse_mid42_ch%d",1+ich),"",t_mid42_wlength,0,t_mid42_wlength,3000,-100,14900);
	}

	//Canvas for NKFADC500
	TCanvas *c1a = new TCanvas(Form("c1a"), Form("RUN %d, NKFADC500 QA", RunNo), -1, 0, 200*1.2*4*xCVS, 200*2*xCVS); 
	gPad->SetMargin(0,0,0,0);
	c1a->Divide(4,2,0,0);

	TH1F* H1Temp = new TH1F("", "", 1024, 0, 1024);
	H1Temp->GetXaxis()->SetRangeUser(0, 124);
	H1Temp->GetXaxis()->SetTitle("Time bin");
	H1Temp->GetXaxis()->SetTitleSize(0.055);
	H1Temp->GetXaxis()->SetLabelSize(0.050);
	H1Temp->GetYaxis()->SetTitle("ADC");
	H1Temp->GetYaxis()->SetTitleSize(0.055);
	H1Temp->GetYaxis()->SetLabelSize(0.050);
	H1Temp->GetYaxis()->SetTitleOffset(1.4);
	H1Temp->GetZaxis()->SetTitleSize(0.055);
	H1Temp->GetZaxis()->SetLabelSize(0.050);

	for (int ii=0; ii<2; ii++){
		for (int jj=0; jj<4; jj++){
			c1a->cd(4*ii+jj+1);
			gPad->SetTicks();
			gPad->SetLogz();
			gPad->SetMargin(0.15,0.12,0.14,0.01);
			if ( ii==0 ){
				H1Temp->GetXaxis()->SetRangeUser(0, t_mid1_wlength);
			}else if ( ii==1 ){
				H1Temp->GetXaxis()->SetRangeUser(0, t_mid2_wlength);
			}
			H1Temp->GetYaxis()->SetRangeUser(0, 6000);
			TH1F *htmp = (TH1F*)H1Temp->DrawCopy();
			htmp->SetStats(0);
		}
	}

	//Canvas for JBNUDAQ
	TCanvas *c1b = new TCanvas(Form("c1b"), Form("RUN %d, JBNUDAQ QA", RunNo), -1, 0, 200*1.2*8*xCVS, 200*4*xCVS); 
	gPad->SetMargin(0,0,0,0);
	c1b->Divide(8,4,0,0);

	for (int ii=0; ii<4; ii++){
		for (int jj=0; jj<8; jj++){
			c1b->cd(8*ii+jj+1);
			gPad->SetTicks();
			gPad->SetLogz();
			gPad->SetMargin(0.15,0.12,0.14,0.01);
			H1Temp->GetYaxis()->SetRangeUser(-1*wform_factor[31], 10*wform_factor[31]);
			H1Temp->GetXaxis()->SetRangeUser(0, t_mid31_wlength);
			TH1F *htmp = (TH1F*)H1Temp->DrawCopy();
			htmp->SetStats(0);
		}
	}

	//Canvas for BICDAQ 
	TCanvas *c1c[2];
	for (int kk=0; kk<2; kk++){
		c1c[kk] = new TCanvas(Form("c1c_%d",kk), Form("RUN %d, BICDAQ QA, L/R %d", RunNo, kk), -1, 0, 200*1.2*8*xCVS, 200*4*xCVS); 
		gPad->SetMargin(0,0,0,0);
		c1c[kk]->Divide(8,4,0,0);

		for (int ii=0; ii<4; ii++){
			for (int jj=0; jj<8; jj++){
				c1c[kk]->cd(8*ii+jj+1);
				gPad->SetTicks();
				gPad->SetLogz();
				gPad->SetMargin(0.15,0.12,0.14,0.01);
				H1Temp->GetYaxis()->SetRangeUser(-1*wform_factor[41], 10*wform_factor[41]);
				H1Temp->GetXaxis()->SetRangeUser(0, t_mid41_wlength);
				TH1F *htmp = (TH1F*)H1Temp->DrawCopy();
				htmp->SetStats(0);
			}
		}
	}


	/*
	for (int ii=0; ii<7; ii++){
		for (int jj=0; jj<16; jj++){
			c1->cd(16*ii + jj + 1);
			if ( ii==0 ){
			}else if ( ii==1 || ii==2 ){
				H1Temp->GetYaxis()->SetRangeUser(-1*wform_factor[41], 10*wform_factor[41]);
				H1Temp->GetXaxis()->SetRangeUser(0, t_mid41_wlength);
			}else if ( ii==3 || ii==4 ){
				H1Temp->GetYaxis()->SetRangeUser(-1*wform_factor[42], 10*wform_factor[42]);
				H1Temp->GetXaxis()->SetRangeUser(0, t_mid42_wlength);
			}else if ( ii==5 || ii==6 ){
			}else{
				H1Temp->GetYaxis()->SetRangeUser(0, 6000);
			}

			gPad->SetTicks();
			gPad->SetMargin(0.15,0.05,0.1,0.01);
			TH1F *htmp = (TH1F*)H1Temp->DrawCopy();
			htmp->SetStats(0);
		}
	}
	c1->cd();
	c1->Update();
	*/

	//return;

	//Analysis
	for (unsigned int ien=0; ien<std::min(nentries, nEvtToRead); ien++){

		T->GetEntry(ien);

		hTrigNumTime->Fill(t_tcb_trigger_number, t_tcb_trigger_time);

		for (int imid=0; imid<mid.size(); imid++){

			if ( mid[imid]==1 ){
				for (int ich=0; ich<t_mid1_nch; ich++){
					for (int iw=0; iw<t_mid1_wlength; iw++){
						int adc = (*t_mid1_adc)[ich][iw];
						hPulse_mid1[ich]->Fill(iw+0.5, adc);
					}//iw
				}//ich
			}else if ( mid[imid]==2 ){
				for (int ich=0; ich<t_mid2_nch; ich++){
					for (int iw=0; iw<t_mid2_wlength; iw++){
						int adc = (*t_mid2_adc)[ich][iw];
						hPulse_mid2[ich]->Fill(iw+0.5, adc);
					}//iw
				}//ich
			}else if ( mid[imid]==31 ){
				for (int ich=0; ich<t_mid31_nch; ich++){
					for (int iw=0; iw<t_mid31_wlength; iw++){
						int adc = (*t_mid31_adc)[ich][iw];
						hPulse_mid31[ich]->Fill(iw+0.5, adc);
					}//iw

				}//ich
			}else if ( mid[imid]==41 ){
				for (int ich=0; ich<t_mid41_nch; ich++){
					for (int iw=0; iw<t_mid41_wlength; iw++){
						int adc = (*t_mid41_adc)[ich][iw];
						hPulse_mid41[ich]->Fill(iw+0.5, adc);
					}//iw

				}//ich
			}else if ( mid[imid]==42 ){
				for (int ich=0; ich<t_mid42_nch; ich++){
					for (int iw=0; iw<t_mid42_wlength; iw++){
						int adc = (*t_mid42_adc)[ich][iw];
						hPulse_mid42[ich]->Fill(iw+0.5, adc);
					}//iw

				}//ich
			}//

		}//imid

	}//ien


	TCanvas *c2 = new TCanvas("c2", "c2", -1, 0, 1.2*400, 400);
	gPad->SetMargin(0.15,0.12,0.12,0.05);
	gPad->SetTicks();
	hTrigNumTime->GetXaxis()->SetTitleSize(0.05);
	hTrigNumTime->GetXaxis()->SetLabelSize(0.045);
	hTrigNumTime->GetYaxis()->SetTitleSize(0.05);
	hTrigNumTime->GetYaxis()->SetLabelSize(0.045);
	hTrigNumTime->GetZaxis()->SetTitleSize(0.05);
	hTrigNumTime->GetZaxis()->SetLabelSize(0.045);
	hTrigNumTime->GetZaxis()->SetRangeUser(0, 1.5);
	TH2D *hTrigNumTime_cp = (TH2D*)hTrigNumTime->DrawCopy("colz");
	hTrigNumTime_cp->SetStats(0);

	if ( find(mid.begin(), mid.end(), 1)!=mid.end() ){
		for (int ich=0; ich<4; ich++){
			c1a->cd(ich+1);
			hPulse_mid1[ich]->Draw("colz same");

			TLegend *leg = new TLegend(0.2, 0.95-0.07*2, 0.5, 0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->AddEntry("","MID 1","h"); 
			leg->AddEntry("",Form("CH %d", ich+1),"h"); 
			leg->Draw();
		}
	}

	if ( find(mid.begin(), mid.end(), 2)!=mid.end() ){
		for (int ich=0; ich<4; ich++){
			c1a->cd(ich+1+4);
			hPulse_mid2[ich]->Draw("colz same");

			TLegend *leg = new TLegend(0.2, 0.95-0.07*2, 0.5, 0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->AddEntry("","MID 2","h"); 
			leg->AddEntry("",Form("CH %d", ich+1),"h"); 
			leg->Draw();
		}
	}

	if ( find(mid.begin(), mid.end(), 31)!=mid.end() ){
		for (int ii=0; ii<2; ii++){
			for (int ich=0; ich<16; ich++){

				int geoch = chMapHodo[16*ii + ich + 1];
				c1b->cd(geoch);
				hPulse_mid31[16*ii + ich]->Draw("colz same");

				TLegend *leg = new TLegend(0.2, 0.95-0.07*3, 0.5, 0.95);
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->SetTextSize(0.06);
				leg->AddEntry("","MID 31","h"); 
				leg->AddEntry("",Form("DAQ CH %d", ich+1),"h"); 
				leg->AddEntry("",Form("GEO CH %d", (geoch<=16) ? geoch : geoch-16),"h"); 
				leg->Draw();
			}
		}
	}

	if ( find(mid.begin(), mid.end(), 41)!=mid.end() ){
		for (int ich=0; ich<32; ich++){

			if ( ich<=3 || (ich>=17 && ich<=20) ) continue;

			vector a = chMapCalo[std::make_pair(41, ich+1)];
			int lrid = a[0];
			int modid = a[1];
			int colid = a[2]; 
			int rowid = a[3];

			c1c[lrid]->cd( colid + 1 + 8*(3 - rowid));
			hPulse_mid41[ich]->Draw("colz same");

			TLegend *leg = new TLegend(0.2, 0.95-0.07*3, 0.5, 0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->AddEntry("","MID 41","h"); 
			leg->AddEntry("",Form("DAQ CH %d", ich+1),"h"); 
			leg->AddEntry("",Form("MOD ID %d L/R %d", modid, lrid),"h"); 
			leg->Draw();

		}
	}

	if ( find(mid.begin(), mid.end(), 42)!=mid.end() ){
		for (int ich=0; ich<32; ich++){

			vector a = chMapCalo[std::make_pair(42, ich+1)];
			int lrid = a[0];
			int modid = a[1];
			int colid = a[2]; 
			int rowid = a[3];

			c1c[lrid]->cd( colid + 1 + 8*(3 - rowid));
			hPulse_mid42[ich]->Draw("colz same");

			TLegend *leg = new TLegend(0.2, 0.95-0.07*3, 0.5, 0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextSize(0.06);
			leg->AddEntry("","MID 42","h"); 
			leg->AddEntry("",Form("DAQ CH %d", ich+1),"h"); 
			leg->AddEntry("",Form("MOD ID %d L/R %d", modid, lrid),"h"); 
			leg->Draw();

		}
	}


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
