#define ReadHallCData_cxx
#include "ReadHallCData.h"

// To write to VCS database
#include "Database.h"
#include "nlohmann/json.hpp" 
using nlohmann::json;

using namespace std;



void ReadHallCData::Loop (vector <string> vector_name, int runID, string  process, string what, string target){

	iChain=new TChain("T");
	for (unsigned int i=0;i<vector_name.size();i++){
		cout<<"including file: "<<vector_name[i].c_str()<<endl;
		iChain->Add(Form("%s",vector_name[i].c_str())); 
	}
	test=Init(process);
	if (test==0 || iChain == 0){
		cout<<"problem init files, tree not found"<<endl;
		return;
	}
	cout<<"number of entries in Chain: "<<iChain->GetEntries()<<endl;	
  
	cout<<"Read cuts"<<endl;
	test=ReadHeader(inc);
	if (!test){ cout<<"don't find option / cut file FATAL ERROR"<<endl; return; }
	LL = inc[0]; trig_time_cut_low=inc[12]; trig_time_cut_up=inc[13]; roc_time_cut_low=inc[9];
	roc_time_cut_up=inc[10]; time_shift=inc[11];
	hdp_low=inc[1]; hdp_up=inc[2]; sdp_low=inc[3]; sdp_up=inc[4];
	th_hms_up=inc[5]; phi_hms_up=inc[6]; th_shms_up=inc[7]; phi_shms_up=inc[8];
	hms_cal_en_cut=inc[14]; shms_cal_en_cut = inc[15];
	hms_cer_npe_cut=inc[16]; shms_aero_npe_cut=inc[17]; shms_hgcer_npe_cut=inc[18];
	hms_beta_min=inc[19]; hms_beta_max=inc[20];
	shms_beta_min=inc[21]; shms_beta_max=inc[22];

	cout<<"user options: "<<endl;
	cout<<"target lenght = "<<LL<<", trigger time cut low/up = "<<trig_time_cut_low<<" "<<trig_time_cut_up<<endl;
	cout<<"roc time cut low/up = "<<roc_time_cut_low<<" "<<roc_time_cut_up<<" time shift = "<<time_shift<<endl;
	cout<<"dp cut low/up HMS = "<<hdp_low<<" "<<hdp_up<<" SHMS = "<<sdp_low<<" "<<sdp_up<<endl;
	cout<<"theta & phi cuts up HMS = "<<th_hms_up<<" "<<phi_hms_up<<" SHMS = "<<th_shms_up<<" "<<phi_shms_up<<endl;
	cout<<"HMS cal energy electron = "<<hms_cal_en_cut<<" SHMS cal energy electron = "<<shms_cal_en_cut<<endl;
	cout<<" HMS Ch npe limit = "<<hms_cer_npe_cut<<endl;
	cout<<"SHMS aerogel npe limit = "<<shms_aero_npe_cut<<" SHMS HGC npe limit = "<<shms_hgcer_npe_cut<<endl;
	cout<<"HMS beta min, max = "<<hms_beta_min<<" "<<hms_beta_max<<endl;
	cout<<"SHMS beta min, max = "<<shms_beta_min<<" "<<shms_beta_max<<endl;
	cout<<"Start producing the new Tree"<<endl;
	
  //InitBins(process);
  //printf("min=%f %f %f %f %f %f %f %f %f %f %f ",M2miss_min, Mmiss_min, Emiss_min, miss_mom_min, PTmiss_min, PT2miss_min, M2miss_max, Mmiss_max, miss_mom_max, Emiss_max, PTmiss_max, PT2miss_max);
  //printf("%d %d %d %d %d %d %d %d %d",M2miss_bins, Mmiss_bins, miss_mom_bins, Emiss_bins, PT2miss_bins, PTmiss_bins, abs_time_bins, rel_time_bins, trig_time_bins);
cout<< hms_beta_min<<" "<< hms_beta_max<<" "<<hms_beta_bins<<" "<<shms_beta_min<<" "<< shms_beta_max<<" "<<shms_beta_bins<<" "<<endl;

	runindex=runID;
	test = FillSingleRunRunInfos(process, runID, Eb, targetmass, HMS_p_central, SHMS_p_central, HMS_th_central, SHMS_th_central, HMS_run_l, SHMS_run_l, 
				     HMS_B2_cur, SHMS_B2_cur, HMS_B4_cut, SHMS_B4_cut, HMS_act_time, SHMS_act_time, HMS_B2_cur_cut, SHMS_B2_cur_cut, HMS_live, SHMS_live);
	if (!test){ cout<<"ERROR: cannot read run information from replay"<<endl; return ; }
	
	cout<<"\n*** database run information:"<<endl;
	cout<<"run ID: "<<runID<<endl;
	cout<<"beam E= "<<Eb<<" target mass = "<<targetmass<<endl;
	cout<<"HMS p central "<<HMS_p_central<<" SHMS p central "<<SHMS_p_central<<" HMS th central: "<<HMS_th_central<<" SHMS th central "<<SHMS_th_central<<endl;
	cout<<"HMS / SHMS run lenght = "<<HMS_run_l<<" "<<SHMS_run_l<<endl;
	cout<<"HMS / SHMS B2 current = "<<HMS_B2_cur<<" "<<SHMS_B2_cur<<" HMS / SHMS B2 current after B4 cut :"<<HMS_B2_cur_cut<<" "<<SHMS_B2_cur_cut<<endl;
	cout<<"cut on B4 HMS/SHMS = "<<HMS_B4_cut<<" "<<SHMS_B4_cut<<endl;
	cout<<"Actual time HMS / SHMS above cut= "<<HMS_act_time<<" "<< SHMS_act_time<<endl;
	cout<<"HMS live time = "<<HMS_live <<" SHMS live time = "<<SHMS_live<<endl;
	
	test = FillSingleRunTriggerInfos(process, runID, HMS34rates, SHMS34rates, HMS_ST, SHMS_ST, C_T, All_T);
	if (!test){ cout<<"ERROR: cannot read trigger information from replay"<<endl; return ; }
	
	cout<<"\n*** database trigger info:"<<endl;
	cout<<"HMS / SHMS 3/4 trigger rates = "<<HMS34rates<<" "<< SHMS34rates<<endl;
        cout<<"HMS / SHMS single trigger rates = "<<HMS_ST<<" "<<SHMS_ST<<endl;
        cout<<"coincidence trigger rates = "<<C_T<<" All triggers "<<All_T<<endl;
	
	test = FillSingleRunEffInfos(process, runID, HMS_E_eff, HMS_H_eff, SHMS_E_eff, SHMS_H_eff);
	if (!test){ cout<<"ERROR: cannot read efficiency information from replay"<<endl; return ; }
	
	cout<<"\n*** database efficiency information from replay"<<endl;
	cout<<"HMS / SHMS electron efficiency= "<<HMS_E_eff<<" "<<SHMS_E_eff<<endl;
	cout<<"HMS / SHMS hadron efficiency= "<<HMS_H_eff<<" "<<SHMS_H_eff<<"\n"<<endl;
	cout<<"Start loop over entries"<<endl;
	if (target.compare("LH2")==0){
		lumiexp_HMS = HMS_B2_cur_cut * 1e-6 / Ce *LL * rhoH * NA/MH * 1e-36 * HMS_act_time * 1000;  
		lumiexp_SHMS = SHMS_B2_cur_cut * 1e-6 / Ce *LL * rhoH * NA/MH * 1e-36 * SHMS_act_time * 1000;
	} else if (target.compare("dummy")==0){
		lumiexp_HMS = HMS_B2_cur_cut * 1e-6 / Ce * 0.1816 * NA/26.9815 * 1e-36 * HMS_act_time * 1000; // check thickness
		lumiexp_SHMS = SHMS_B2_cur_cut * 1e-6 / Ce * 0.1816 * NA/26.9815 * 1e-36 * SHMS_act_time * 1000;
	}
	cout<<"\n*** Measured luminosity HMS / SHMS = "<<lumiexp_HMS<< " "<<lumiexp_SHMS<< " nb"<<endl;  

  test = FillSingleRunMissRefInfos(process, misstimeref, runID);
  if (!test){ cout<<"WARNING: cannot read missing time information"<<endl;  }

	TFile *file = new TFile(Form(VCS_REPLAY_PATH "/Ana/files/HallCData_%d.root",runID),"RECREATE"); 
        if ( file->IsOpen() ) printf("File opened successfully\n");
	TTree *HallCTree = new TTree("HallCTree","HallC tree + internal format");
	gDirectory->ls("-m");
	
	HallCTree->Branch("runindex",&runindex, "runindex/I");
        HallCTree->Branch("ALV_el_out_data",&ALV_el_out_data,"ALV_el_out_data[4]/F");
        HallCTree->Branch("ALV_el_in_data",&ALV_el_in_data,"ALV_el_in_data[4]/F");
        if (process.compare("vcs")==0 || process.compare("vcsLT")==0 
	    || process.compare("elastic")==0 || process.compare("elasticLT")==0){
		HallCTree->Branch("ALV_gamma_out_data",&ALV_gamma_out_data,"ALV_gamma_out_data[4]/F");
		HallCTree->Branch("ALV_proton_out_data",&ALV_proton_out_data,"ALV_proton_out_data[4]/F");
	} else if (process.compare("pi0")==0 || process.compare("pi0LT")==0
		   ){
		HallCTree->Branch("ALV_pi_out_data",&ALV_pi_out_data,"ALV_pi_out_data[4]/F");
		HallCTree->Branch("ALV_proton_out_data",&ALV_proton_out_data,"ALV_proton_out_data[4]/F");
	} else if (process.compare("pi+")==0 || process.compare("pi+LT")==0
		   || process.compare("K+")==0 || process.compare("K+LT")==0
		   ){
		HallCTree->Branch("ALV_pi_out_data",&ALV_pi_out_data,"ALV_pi_out_data[4]/F");
		HallCTree->Branch("ALV_neutron_out_data",&ALV_neutron_out_data,"ALV_neutron_out_data[4]/F");
	}
	HallCTree->Branch("M2miss_data",&M2miss_data,"M2miss_data/F");
	HallCTree->Branch("Mmiss_data",&Mmiss_data,"Mmiss_data/F");
	HallCTree->Branch("Emiss_data",&Emiss_data,"Emiss_data/F");
	HallCTree->Branch("PTmiss_data",&PTmiss_data,"PTmiss_data/F");
	HallCTree->Branch("PT2miss_data",&PT2miss_data,"PT2miss_data/F");
	HallCTree->Branch("Q2_kinmod", &Q2_kinmod,"Q2_kinmod/F");
	HallCTree->Branch("W_kinmod", &W_kinmod,"W_kinmod/F");
	HallCTree->Branch("nu_kinmod", &nu_kinmod,"nu_kinmod/F");
	HallCTree->Branch("mt_kinmod", &mt_kinmod,"mt_kinmod/F");
	HallCTree->Branch("epsilon_kinmod", &epsilon_kinmod,"epsilon_kinmod/F");
	HallCTree->Branch("Xbj_kinmod", &Xbj_kinmod,"Xbj_kinmod/F");
	HallCTree->Branch("Q2_data", &Q2_data,"Q2_data/F");
	HallCTree->Branch("W_data", &W_data,"W_data/F");
	HallCTree->Branch("nu_data", &nu_data,"nu_data/F");
	HallCTree->Branch("mt_data", &mt_data,"mt_data/F");
	HallCTree->Branch("epsilon_data", &epsilon_data,"epsilon_data/F");
	HallCTree->Branch("Xbj_data", &Xbj_data,"Xbj_data/F");
	HallCTree->Branch("CosThCM_data", &CosThCM_data,"CosThCM_data/F");
	HallCTree->Branch("ThCM_data", &ThCM_data,"ThCM_data/F");
	HallCTree->Branch("Phi_data", &Phi_data,"Phi_data/F");
	HallCTree->Branch("beta_SHMS", &beta_SHMS, "beta_SHMS/F");
	HallCTree->Branch("beta_HMS", &beta_HMS, "beta_HMS/F");
	HallCTree->Branch("Eb", &Eb, "Eb/F");
	HallCTree->Branch("targetmass", &targetmass, "targetmass/F");
	HallCTree->Branch("HMS_p_central", &HMS_p_central, "HMS_p_central/F");
	HallCTree->Branch("SHMS_p_central", &SHMS_p_central, "SHMS_p_central/F");
	HallCTree->Branch("HMS_th_central", &HMS_th_central, "HMS_th_central/F");
	HallCTree->Branch("SHMS_th_central", &SHMS_th_central, "SHMS_th_central/F");
	HallCTree->Branch("HMS_run_l", &HMS_run_l, "HMS_run_l/F");
	HallCTree->Branch("SHMS_run_l", &SHMS_run_l, "SHMS_run_l/F");
	HallCTree->Branch("HMS_B2_cur", &HMS_B2_cur, "HMS_B2_cur/F");
	HallCTree->Branch("SHMS_B2_cur", &SHMS_B2_cur, "SHMS_B2_cur/F");
	HallCTree->Branch("HMS_B4_cut", &HMS_B4_cut, "HMS_B4_cut/F");
	HallCTree->Branch("SHMS_B4_cut", &SHMS_B4_cut, "SHMS_B4_cut/F");
	HallCTree->Branch("HMS_B2_cur_cut", &HMS_B2_cur_cut, "HMS_B2_cur_cut/F");
	HallCTree->Branch("SHMS_B2_cur_cut", &SHMS_B2_cur_cut, "SHMS_B2_cur_cut/F");
	HallCTree->Branch("HMS_live", &HMS_live, "HMS_live/F");
	HallCTree->Branch("SHMS_live", &SHMS_live, "SHMS_live/F");
	HallCTree->Branch("HMS34rates", &HMS34rates, "HMS34rates/F");
	HallCTree->Branch("SHMS34rates", &SHMS34rates, "SHMS34rates/F");
	HallCTree->Branch("HMS_ST", &HMS_ST, "HMS_ST/F");
	HallCTree->Branch("SHMS_ST", &SHMS_ST, "SHMS_ST/F");
	HallCTree->Branch("C_T", &C_T, "C_T/F");
	HallCTree->Branch("All_T", &All_T, "All_T/F");
	HallCTree->Branch("HMS_E_eff", &HMS_E_eff, "HMS_E_eff/F");
	HallCTree->Branch("HMS_H_eff", &HMS_H_eff, "HMS_H_eff/F");
	HallCTree->Branch("SHMS_E_eff", &SHMS_E_eff, "SHMS_E_eff/F");
	HallCTree->Branch("SHMS_H_eff", &SHMS_H_eff, "SHMS_H_eff/F");
	HallCTree->Branch("HMS_act_time", &HMS_act_time, "HMS_act_time/F");
	HallCTree->Branch("SHMS_act_time", &SHMS_act_time, "SHMS_act_time/F");
	HallCTree->Branch("lumiexp_HMS", &lumiexp_HMS, "lumiexp_HMS/F");
	HallCTree->Branch("lumiexp_SHMS", &lumiexp_SHMS, "lumiexp_SHMS/F");
	HallCTree->Branch("time_shift", &time_shift, "time_shift/F");
	HallCTree->Branch("time_roc1", &time_roc1, "time_roc1/F");
	HallCTree->Branch("time_roc2", &time_roc2, "time_roc2/F");
	HallCTree->Branch("time_trig1", &time_trig1, "time_trig1/F");
	HallCTree->Branch("time_trig2", &time_trig2, "time_trig2/F");
  HallCTree->Branch("misstimeref", &misstimeref, "misstimeref[22]/I");

	cout<<"out tree init "<<HallCTree->GetEntries()<<endl;
        for (Int_t i=0; i<iChain -> GetEntries(); i++) {

                if (i%10000==0) cout<<"# "<<i<<endl;
                entry = iChain->GetEntry(i);
                if (i==0) cerr<<"enter loop"<<endl;
		
		h_CTime_epCoinTime_ROC1[0]->Fill(CTime_epCoinTime_ROC1-time_shift);
		h_CTime_epCoinTime_ROC2[0]->Fill(CTime_epCoinTime_ROC2-time_shift);
		h_CTime_epCoinTime_ROC2_large[0]->Fill(CTime_epCoinTime_ROC2);
		h_CTime_epCoinTime_TRIG1[0]->Fill(CTime_epCoinTime_TRIG1-time_shift);
		h_CTime_epCoinTime_TRIG2[0]->Fill(CTime_epCoinTime_TRIG4-time_shift);
	
		h_encal_hms->Fill(H_cal_etotnorm);
		h_encal_shms->Fill(P_cal_etotnorm);
		h_cer_npe_hms->Fill(H_cer_npeSum);
		h_aero_npe_shms->Fill(P_aero_npeSum);
		h_hgcer_npe_shms->Fill(P_hgcer_npeSum);
		h_delta_hms->Fill(H_gtr_dp);
		h_delta_shms->Fill(P_gtr_dp);
		h_phi_hms->Fill(H_gtr_ph);
		h_theta_hms->Fill(P_gtr_th);
		h_phi_shms->Fill(P_gtr_ph);
		h_theta_shms->Fill(P_gtr_th);
		h2_encal_cer->Fill(H_cal_etotnorm,H_cer_npeSum);
		h2_aero_hgcer_npe->Fill(P_aero_npeSum,P_hgcer_npeSum);

		// acceptance fiducial cuts
		if (H_gtr_dp<hdp_low || H_gtr_dp>hdp_up) continue; // hms fiducial cut
		if (P_gtr_dp<sdp_low || P_gtr_dp>sdp_up) continue; // shms fiducial cut
		if (fabs(H_gtr_th)>th_hms_up || fabs(H_gtr_ph)>phi_hms_up || fabs(P_gtr_th)>th_shms_up || fabs(P_gtr_ph)>phi_shms_up) continue; // angular fiducial cuts; process dependent adjust
		
		h2_Pebeta[1]->Fill(H_gtr_p, H_gtr_beta);
		h2_elebeta[1]->Fill(P_gtr_p, P_gtr_beta);
		h2_COIN_H_beta[1]->Fill(CTime_epCoinTime_ROC1-time_shift,H_gtr_beta);		
		h2_COIN_P_beta[1]->Fill(CTime_epCoinTime_ROC1-time_shift,P_gtr_beta);		
		
		h_CTime_epCoinTime_ROC1[1]->Fill(CTime_epCoinTime_ROC1-time_shift);
		h_CTime_epCoinTime_ROC2[1]->Fill(CTime_epCoinTime_ROC2-time_shift);
		h_CTime_epCoinTime_ROC2_large[1]->Fill(CTime_epCoinTime_ROC2);
		h_CTime_epCoinTime_TRIG1[1]->Fill(CTime_epCoinTime_TRIG1-time_shift);
		h_CTime_epCoinTime_TRIG2[1]->Fill(CTime_epCoinTime_TRIG4-time_shift);

		// PID cuts
		if (process.compare("elastic")==0 || process.compare("vcs")==0 || process.compare("pi0")==0){
			if ((hms_cal_en_cut>0 && H_cal_etotnorm>hms_cal_en_cut) || (hms_cer_npe_cut>0 && H_cer_npeSum>hms_cer_npe_cut)
			   || (hms_beta_min>0 && H_gtr_beta<hms_beta_min) || (hms_beta_max>0 && H_gtr_beta>hms_beta_max)
			    ) continue; // proton in HMS
			if ((shms_aero_npe_cut>0 && P_aero_npeSum<shms_aero_npe_cut) || (shms_hgcer_npe_cut>0 && P_hgcer_npeSum<shms_hgcer_npe_cut)
			    || (shms_cal_en_cut>0 && P_cal_etotnorm>shms_cal_en_cut) 
			    || (shms_beta_min>0 && P_gtr_beta<shms_beta_min) || (shms_beta_max>0 && P_gtr_beta>shms_beta_max)
			    ) continue; // electron in SHMS 
		} else if (process.compare("elasticLT")==0 || process.compare("vcsLT")==0 || process.compare("pi0LT")==0){
			if ((hms_cal_en_cut>0 && H_cal_etotnorm<hms_cal_en_cut) || (hms_cer_npe_cut>0 && H_cer_npeSum<hms_cer_npe_cut)
			   || (hms_beta_min>0 && H_gtr_beta<hms_beta_min) || (hms_beta_max>0 && H_gtr_beta>hms_beta_max)
			    ) continue; // electron in HMS
			if ((shms_aero_npe_cut>0 && P_aero_npeSum>shms_aero_npe_cut) || (shms_hgcer_npe_cut>0 && P_hgcer_npeSum>shms_hgcer_npe_cut)
			    || (shms_cal_en_cut>0 && P_cal_etotnorm>shms_cal_en_cut) 
			    || (shms_beta_min>0 && P_gtr_beta<shms_beta_min) || (shms_beta_max>0 && P_gtr_beta>shms_beta_max)
			    ) continue;  // proton in SHMS
		} else if (process.compare("pi+LT")==0){
			if ((hms_cal_en_cut>0 && H_cal_etotnorm<hms_cal_en_cut) || (hms_cer_npe_cut>0 && H_cer_npeSum<hms_cer_npe_cut)
			   || (hms_beta_min>0 && H_gtr_beta<hms_beta_min) || (hms_beta_max>0 && H_gtr_beta>hms_beta_max)
			    ) continue; // electron in HMS
			if ((shms_aero_npe_cut>0 && P_aero_npeSum<shms_aero_npe_cut) 
			    || (shms_cal_en_cut>0 && P_cal_etotnorm>shms_cal_en_cut) 
			    || (shms_beta_min>0 && P_gtr_beta<shms_beta_min) || (shms_beta_max>0 && P_gtr_beta>shms_beta_max)
			    ) continue;  // pion in SHMS
		} else if (process.compare("K+LT")==0){
			if ((shms_aero_npe_cut>0 && P_aero_npeSum<shms_aero_npe_cut) || (shms_hgcer_npe_cut>0 && P_hgcer_npeSum>shms_hgcer_npe_cut)
			    || (shms_cal_en_cut>0 && P_cal_etotnorm>shms_cal_en_cut) 
			    || (shms_beta_min>0 && P_gtr_beta<shms_beta_min) || (shms_beta_max>0 && P_gtr_beta>shms_beta_max)
			    ) continue;  // kaon in SHMS
			if ((hms_cal_en_cut>0 && H_cal_etotnorm<hms_cal_en_cut) || (hms_cer_npe_cut>0 && H_cer_npeSum<hms_cer_npe_cut)
			   || (hms_beta_min>0 && H_gtr_beta<hms_beta_min) || (hms_beta_max>0 && H_gtr_beta>hms_beta_max)
			    ) continue; // electron in HMS
		} else if (process.compare("pi+")==0){
			if ((shms_aero_npe_cut>0 && P_aero_npeSum<shms_aero_npe_cut) || (shms_hgcer_npe_cut>0 && P_hgcer_npeSum<shms_hgcer_npe_cut)
			    || (shms_cal_en_cut>0 && P_cal_etotnorm>shms_cal_en_cut) 
			    || (shms_beta_min>0 && P_gtr_beta<shms_beta_min) || (shms_beta_max>0 && P_gtr_beta>shms_beta_max)
			    ) continue; // electron in SHMS 
			if ((hms_cal_en_cut>0 && H_cal_etotnorm<hms_cal_en_cut) || (hms_cer_npe_cut>0 && H_cer_npeSum<hms_cer_npe_cut)
			   || (hms_beta_min>0 && H_gtr_beta<hms_beta_min) || (hms_beta_max>0 && H_gtr_beta>hms_beta_max)
			    ) continue; // pi+ in HMS
		} else if (process.compare("K+")==0){
			if ((shms_aero_npe_cut>0 && P_aero_npeSum<shms_aero_npe_cut) || (shms_hgcer_npe_cut>0 && P_hgcer_npeSum<shms_hgcer_npe_cut)
			    || (shms_cal_en_cut>0 && P_cal_etotnorm>shms_cal_en_cut) 
			    || (shms_beta_min>0 && P_gtr_beta<shms_beta_min) || (shms_beta_max>0 && P_gtr_beta>shms_beta_max)
			    ) continue; // electron in SHMS 
			if ((hms_cal_en_cut>0 && H_cal_etotnorm<hms_cal_en_cut) || (hms_cer_npe_cut>0 && H_cer_npeSum>hms_cer_npe_cut)
			   || (hms_beta_min>0 && H_gtr_beta<hms_beta_min) || (hms_beta_max>0 && H_gtr_beta>hms_beta_max)
			    ) continue; // kaon+ in HMS
		}

		h_CTime_epCoinTime_ROC1[2]->Fill(CTime_epCoinTime_ROC1-time_shift);
		h_CTime_epCoinTime_ROC2[2]->Fill(CTime_epCoinTime_ROC2-time_shift);
		h_CTime_epCoinTime_ROC2_large[2]->Fill(CTime_epCoinTime_ROC2);
		h_CTime_epCoinTime_TRIG1[2]->Fill(CTime_epCoinTime_TRIG1-time_shift);
		h_CTime_epCoinTime_TRIG2[2]->Fill(CTime_epCoinTime_TRIG4-time_shift);

		h2_Pebeta[2]->Fill(H_gtr_p, H_gtr_beta);
		h2_elebeta[2]->Fill(P_gtr_p, P_gtr_beta);
		h2_COIN_H_beta[2]->Fill(CTime_epCoinTime_ROC1-time_shift,H_gtr_beta);		
		h2_COIN_P_beta[2]->Fill(CTime_epCoinTime_ROC1-time_shift,P_gtr_beta);		
		
		// allocations to new tree
		beta_HMS = H_gtr_beta; 
		beta_SHMS = P_gtr_beta;
   	Q2_kinmod = H_kin_primary_Q2;
    W_kinmod = H_kin_primary_W;
	  epsilon_kinmod = H_kin_primary_epsilon;
	  nu_kinmod = H_kin_primary_nu;
		Xbj_kinmod = H_kin_primary_x_bj;
		mt_kinmod = -P_kin_secondary_MandelT; 
		time_trig1 = CTime_epCoinTime_TRIG1; 
		time_trig2 = CTime_epCoinTime_TRIG4; 
		time_roc1 = CTime_epCoinTime_ROC1; 
		time_roc2 = CTime_epCoinTime_ROC2;

		// calculation of kinematics and "missing"
		if (process.compare("elastic")==0 || process.compare("vcs")==0 || process.compare("pi0")==0){
			LV_el_out_data.SetPxPyPzE(P_gtr_px, P_gtr_py, P_gtr_pz, sqrt(pow(P_gtr_px,2)+pow(P_gtr_py,2)+pow(P_gtr_pz,2)+pow(m_el,2)));  // if e in SHMS
			LV_proton_out_data.SetPxPyPzE(H_gtr_px, H_gtr_py, H_gtr_pz, sqrt(pow(H_gtr_px,2)+pow(H_gtr_py,2)+pow(H_gtr_pz,2)+pow(M_Proton,2)));  // if p in HMS
		} else if (process.compare("K+")==0){
			LV_el_out_data.SetPxPyPzE(P_gtr_px, P_gtr_py, P_gtr_pz, sqrt(pow(P_gtr_px,2)+pow(P_gtr_py,2)+pow(P_gtr_pz,2)+pow(m_el,2))); 
			LV_proton_out_data.SetPxPyPzE(H_gtr_px, H_gtr_py, H_gtr_pz, sqrt(pow(H_gtr_px,2)+pow(H_gtr_py,2)+pow(H_gtr_pz,2)+pow(M_kp,2)));  // if k+ in HMS
		} else if (process.compare("pi+")==0){
			LV_el_out_data.SetPxPyPzE(P_gtr_px, P_gtr_py, P_gtr_pz, sqrt(pow(P_gtr_px,2)+pow(P_gtr_py,2)+pow(P_gtr_pz,2)+pow(m_el,2))); 
			LV_proton_out_data.SetPxPyPzE(H_gtr_px, H_gtr_py, H_gtr_pz, sqrt(pow(H_gtr_px,2)+pow(H_gtr_py,2)+pow(H_gtr_pz,2)+pow(M_pip,2)));  // if pi+ in HMS
		} else if (process.compare("pi+LT")==0){
			LV_el_out_data.SetPxPyPzE(H_gtr_px, H_gtr_py, H_gtr_pz, sqrt(pow(H_gtr_px,2)+pow(H_gtr_py,2)+pow(H_gtr_pz,2)+pow(m_el,2))); // if e in HMS 
			LV_proton_out_data.SetPxPyPzE(P_gtr_px, P_gtr_py, P_gtr_pz, sqrt(pow(P_gtr_px,2)+pow(P_gtr_py,2)+pow(P_gtr_pz,2)+pow(M_pip,2))); // if pi+ in SHMS
		} else if (process.compare("K+LT")==0){
			LV_el_out_data.SetPxPyPzE(H_gtr_px, H_gtr_py, H_gtr_pz, sqrt(pow(H_gtr_px,2)+pow(H_gtr_py,2)+pow(H_gtr_pz,2)+pow(m_el,2))); // if e in HMS 
			LV_proton_out_data.SetPxPyPzE(P_gtr_px, P_gtr_py, P_gtr_pz, sqrt(pow(P_gtr_px,2)+pow(P_gtr_py,2)+pow(P_gtr_pz,2)+pow(M_kp,2))); // if k+ in SHMS
		} else if (process.compare("elasticLT")==0 || process.compare("pi0LT")==0 || process.compare("vcsLT")) {
			LV_el_out_data.SetPxPyPzE(H_gtr_px, H_gtr_py, H_gtr_pz, sqrt(pow(H_gtr_px,2)+pow(H_gtr_py,2)+pow(H_gtr_pz,2)+pow(m_el,2))); // if e in HMS 
			LV_proton_out_data.SetPxPyPzE(P_gtr_px, P_gtr_py, P_gtr_pz, sqrt(pow(P_gtr_px,2)+pow(P_gtr_py,2)+pow(P_gtr_pz,2)+pow(M_Proton,2))); // if P in SHMS
		}
		LV_proton_in_data.SetPxPyPzE(0,0,0,M_Proton);
		LV_el_in_data.SetPxPyPzE(0,0,sqrt(pow(Eb,2)-pow(m_el,2)),Eb);
		LV_virtual.SetPxPyPzE(-LV_el_out_data.Px()+LV_el_in_data.Px(), -LV_el_out_data.Py()+LV_el_in_data.Py(), 
				     -LV_el_out_data.Pz()+LV_el_in_data.Pz(), -LV_el_out_data.E()+LV_el_in_data.E());	
		LV_gamma_out_data.SetPxPyPzE(LV_virtual.Px()-LV_proton_out_data.Px(),LV_virtual.Py()-LV_proton_out_data.Py(),
					 LV_virtual.Pz()-LV_proton_out_data.Pz(),LV_proton_in_data.E()-LV_proton_out_data.E()+LV_virtual.E() );

		Q2_data = -LV_virtual.M2();
		W_data = (LV_proton_in_data+LV_virtual).M();
		Xbj_data = Q2_data/(2.*M_Proton*LV_virtual.E());
		nu_data = LV_virtual.E();
		if ( process.compare("pi+")==0 || process.compare("pi+LT")==0
		    || process.compare("K+")==0 || process.compare("K+LT")==0
		    ) mt_data = -(LV_proton_in_data-LV_gamma_out_data).M2(); // for pi+ data
		else mt_data = -(LV_proton_in_data-LV_proton_out_data).M2(); // for VCS and pi0 data
		epsilon_data = pow(1.+2*(pow(LV_virtual.E() ,2)+Q2_data)/Q2_data*pow(tan(LV_el_out_data.Theta()/2.),2),-1);

		M2miss_data= LV_gamma_out_data.M2();
		Mmiss_data= LV_gamma_out_data.M();
		PTmiss_data = LV_gamma_out_data.Perp();
		PT2miss_data = LV_gamma_out_data.Perp2();
		Emiss_data = LV_gamma_out_data.E(); // put vs ref
		
		h2_M2miss_CT1[0]->Fill(CTime_epCoinTime_ROC1,M2miss_data);
		h2_M2miss_CT2[0]->Fill(CTime_epCoinTime_ROC2,M2miss_data);


		// timing cuts
		if ((CTime_epCoinTime_ROC2-time_shift)<roc_time_cut_low) continue;
		if ((CTime_epCoinTime_ROC2-time_shift)>roc_time_cut_up) continue;
		if ((CTime_epCoinTime_TRIG1-time_shift)<trig_time_cut_low) continue;
		if ((CTime_epCoinTime_TRIG1-time_shift)>trig_time_cut_up) continue;

		// fill array
		FillArray_LV(LV_el_out_data, ALV_el_out_data);
		FillArray_LV(LV_el_in_data, ALV_el_in_data);
		if (process.compare("elasticLT")==0 || process.compare("elastic")==0
		    || process.compare("vcs")==0 || process.compare("vcsLT")==0) {
			FillArray_LV(LV_gamma_out_data, ALV_gamma_out_data);
			FillArray_LV(LV_proton_out_data, ALV_proton_out_data);
		} else if (process.compare("pi0LT")==0 || process.compare("pi0")==0){
			FillArray_LV(LV_proton_out_data, ALV_proton_out_data);
			FillArray_LV(LV_gamma_out_data, ALV_pi_out_data);
		} else if (process.compare("pi+LT")==0 || process.compare("pi+")==0
			   || process.compare("K+")==0 || process.compare("K+LT")==0
			   ){
			FillArray_LV(LV_proton_out_data, ALV_pi_out_data);
			FillArray_LV(LV_gamma_out_data, ALV_neutron_out_data);
		}

		HallCTree->Fill();
		if (what.compare("reduce")==0) continue;

		if (process.compare("elastic")==0 || process.compare("elasticLT")==0){
			Phi_data = 0; CosThCM_data=0; ThCM_data=0;
		} else if (process.compare("vcs")==0 || process.compare("pi0")==0  || process.compare("pi+")==0
			   || process.compare("vcsLT")==0 || process.compare("pi0LT")==0  || process.compare("pi+LT")==0
			   || process.compare("K+")==0 || process.compare("K+LT")==0
			   ){ 
			V_Normal_Mu = (LV_el_in_data.Vect().Cross(LV_el_out_data.Vect())).Unit();
			V_Normal_Final = (LV_virtual.Vect().Cross(LV_proton_out_data.Vect())).Unit();
			Phi_data = TMath::ACos(V_Normal_Mu * V_Normal_Final);
			Phi_data = 360.-Phi_data*180./PI;
			if ( (V_Normal_Mu *LV_gamma_out_data.Vect()) < 0.)  Phi_data = 360.-Phi_data;
			PinCM = M_Proton*LV_virtual.P()/W_data;
			if (process.compare("vcs")==0 || process.compare("vcsLT")==0 || process.compare("pi0")==0 || process.compare("pi0LT")==0){
				if (process.compare("pi0")==0 || process.compare("pi0LT")==0) {
					Egout_CMeP = (pow(W_data, 2)-pow(M_Proton,2))/(2.*W_data);
					PoutCM = sqrt(pow(Egout_CMeP,2)-pow(M_piz,2));
					CosThCM_data = (-mt_data +Q2_data -pow(M_piz,2) + 2*Egout_CMeP*sqrt(pow(PinCM,2)-Q2_data))/ (2.*PinCM*PoutCM);
				} else { 
					Egout_CMeP = (pow(W_data, 2)-pow(M_Proton,2) + pow(M_piz,2))/(2.*W_data);
					PoutCM =Egout_CMeP;
					CosThCM_data = (-mt_data +Q2_data + 2*Egout_CMeP*sqrt(pow(PinCM,2)-Q2_data)) / (2.*PinCM*Egout_CMeP);
				}
			} else if (process.compare("pi+")==0 || process.compare("pi+LT")==0 || process.compare("K+")==0 || process.compare("K+LT")==0) {
				if (process.compare("pi+")==0 || process.compare("pi+LT")==0){
					PoutCM=sqrt( (pow((pow(W_data,2)- pow(M_pip,2) -M_Neutron*M_Neutron),2 )-4.* pow(M_pip,2) *pow(M_Neutron,2))/(4.*pow(W_data,2)));
					Egout_CMeP=(pow(W_data,2)-pow(M_Neutron,2)+ pow(M_pip,2))/(2.*W_data) ;
					CosThCM_data = (-mt_data +Q2_data -pow(M_pip,2) + 2*Egout_CMeP*sqrt(pow(PinCM,2)-Q2_data))/ (2.*PinCM*PoutCM);
				} else {
					PoutCM=sqrt( (pow((pow(W_data,2)- pow(M_kp,2) -M_Neutron*M_Neutron),2 )-4.* pow(M_kp,2) *pow(M_Neutron,2))/(4.*pow(W_data,2)));
					Egout_CMeP=(pow(W_data,2)-pow(M_Neutron,2)+ pow(M_kp,2))/(2.*W_data) ;
					CosThCM_data = (-mt_data +Q2_data -pow(M_kp,2) + 2*Egout_CMeP*sqrt(pow(PinCM,2)-Q2_data))/ (2.*PinCM*PoutCM);
				}
			}
			ThCM_data = acos(CosThCM_data)*180./PI;

		}

		// no cut, from kinematic module
		h_Q2[0]->Fill(Q2_kinmod);
		h_epsilon[0]->Fill(epsilon_kinmod);
		h_Xbj[0]->Fill(Xbj_kinmod);
		h_W[0]->Fill(W_kinmod);
		h_nu[0]->Fill(nu_kinmod);
		h_mt[0]->Fill(mt_kinmod);

		h2_Pemom[0]->Fill(LV_proton_out_data.E(), LV_proton_out_data.P());
		h2_elemom[0]->Fill(LV_el_out_data.E(), LV_el_out_data.P());
		h2_gemom[0]->Fill(LV_gamma_out_data.E(), LV_gamma_out_data.P());
		h_Pmom[0]->Fill(LV_proton_out_data.P());
		h_elmom[0]->Fill(LV_el_out_data.P());
		h_gmom[0]->Fill(LV_gamma_out_data.P());

		h2_Q2W[0]->Fill(Q2_kinmod, W_kinmod);
		h2_XQ2[0]->Fill(Xbj_kinmod, Q2_kinmod);
		h2_Q2mt[0]->Fill(Q2_kinmod, mt_kinmod);
		h2_nuep[0]->Fill(nu_kinmod,epsilon_kinmod);

		// no cut, from calculations
		h_Q2[1]->Fill(Q2_data);
		h_epsilon[1]->Fill(epsilon_data);
		h_Xbj[1]->Fill(Xbj_data);
		h_W[1]->Fill(W_data);
		h_nu[1]->Fill(nu_data);
		h_mt[1]->Fill(mt_data);

		h2_Q2W[1]->Fill(Q2_data, W_data);
		h2_XQ2[1]->Fill(Xbj_data, Q2_data);
		h2_Q2mt[1]->Fill(Q2_data, mt_data);
		h2_nuep[1]->Fill(nu_data,epsilon_data);

		h_Mmiss[1]->Fill(Mmiss_data);
		h_M2miss[1]->Fill(M2miss_data);
		h_Emiss[1]->Fill(Emiss_data);
		h_PTmiss[1]->Fill(PTmiss_data);
		h_PT2miss[1]->Fill(PT2miss_data);
		
		if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		    ){
			h2_Q2Th[0]->Fill(Q2_kinmod, ThCM_data);
			h2_WTh[0]->Fill(W_kinmod, ThCM_data);
			h2_ThCMPhi[0]->Fill(Phi_data, ThCM_data);
	        	h2_mtTh[0]->Fill(mt_kinmod, ThCM_data);
		
			h2_Q2Th[1]->Fill(Q2_data, ThCM_data);
			h2_WTh[1]->Fill(W_data, ThCM_data);
			h2_ThCMPhi[1]->Fill(Phi_data, ThCM_data);
	        	h2_mtTh[1]->Fill(mt_data, ThCM_data);

			h_Phi[1]->Fill(Phi_data);
			h_ThCM[1]->Fill(ThCM_data);
			h_CosThCM[1]->Fill(CosThCM_data);
		}

		if (fabs(CTime_epCoinTime_ROC2-time_shift)<1.0){

			h_Pmom[2]->Fill(LV_proton_out_data.P());
			h_elmom[2]->Fill(LV_el_out_data.P());
			h_Q2[2]->Fill(Q2_data);
			h_epsilon[2]->Fill(epsilon_data);
			h_Xbj[2]->Fill(Xbj_data);
			h_W[2]->Fill(W_data);
			h_nu[2]->Fill(nu_data);
			h_mt[2]->Fill(mt_data);
			
			h2_Q2W[2]->Fill(Q2_data, W_data);
			h2_XQ2[2]->Fill(Xbj_data, Q2_data);
			h2_Q2mt[2]->Fill(Q2_data, mt_data);
			h2_nuep[2]->Fill(nu_data,epsilon_data);
			h2_Pemom[2]->Fill(LV_proton_out_data.E(), LV_proton_out_data.P());
			h2_elemom[2]->Fill(LV_el_out_data.E(), LV_el_out_data.P());
			h2_gemom[2]->Fill(LV_gamma_out_data.E(), LV_gamma_out_data.P());

			h_Mmiss[2]->Fill(Mmiss_data);
			h_M2miss[2]->Fill(M2miss_data);
			h_Emiss[2]->Fill(Emiss_data);
			h_PTmiss[2]->Fill(PTmiss_data);
			h_PT2miss[2]->Fill(PT2miss_data);
			h_gmom[2]->Fill(LV_gamma_out_data.P());
		
			h_Mmiss[4]->Fill(Mmiss_data);
			h_M2miss[4]->Fill(M2miss_data);
			h_Emiss[4]->Fill(Emiss_data);
			h_PTmiss[4]->Fill(PTmiss_data);
			h_PT2miss[4]->Fill(PT2miss_data);
			h_gmom[4]->Fill(LV_gamma_out_data.P());
		
			h_Q2[4]->Fill(Q2_data);
			h_epsilon[4]->Fill(epsilon_data);
			h_Xbj[4]->Fill(Xbj_data);
			h_W[4]->Fill(W_data);
			h_nu[4]->Fill(nu_data);
			h_mt[4]->Fill(mt_data);
			h_Pmom[4]->Fill(LV_proton_out_data.P());
			h_elmom[4]->Fill(LV_el_out_data.P());

			h2_Q2W[4]->Fill(Q2_data, W_data);
			h2_XQ2[4]->Fill(Xbj_data, Q2_data);
			h2_Q2mt[4]->Fill(Q2_data, mt_data);
			h2_nuep[4]->Fill(nu_data,epsilon_data);
			h2_Pemom[4]->Fill(LV_proton_out_data.E(), LV_proton_out_data.P());
			h2_elemom[4]->Fill(LV_el_out_data.E(), LV_el_out_data.P());
			h2_gemom[4]->Fill(LV_gamma_out_data.E(), LV_gamma_out_data.P());
		
				
			if (process.compare("vcs")==0 || process.compare("pi0")==0
			    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
			    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
			    || process.compare("K+")==0 || process.compare("pi+")==0
				){
				h_Phi[2]->Fill(Phi_data);
				h_ThCM[2]->Fill(ThCM_data);
				h_CosThCM[2]->Fill(CosThCM_data);
				h2_Q2Th[2]->Fill(Q2_data, ThCM_data);
				h2_WTh[2]->Fill(W_data, ThCM_data);
				h2_ThCMPhi[2]->Fill(Phi_data, ThCM_data);
				h2_mtTh[2]->Fill(mt_data, ThCM_data);
			
				h_Phi[4]->Fill(Phi_data);
				h_ThCM[4]->Fill(ThCM_data);
				h_CosThCM[4]->Fill(CosThCM_data);
				h2_Q2Th[4]->Fill(Q2_data, ThCM_data);
				h2_WTh[4]->Fill(W_data, ThCM_data);
				h2_ThCMPhi[4]->Fill(Phi_data, ThCM_data);
				h2_mtTh[4]->Fill(mt_data, ThCM_data);
			}

		} else if (fabs(CTime_epCoinTime_ROC1-time_shift)>7 && fabs(CTime_epCoinTime_ROC1-time_shift)<13){
		
			h_Mmiss[3]->Fill(Mmiss_data, 1./6.);
			h_M2miss[3]->Fill(M2miss_data, 1./6.);
			h_Emiss[3]->Fill(Emiss_data, 1./6.);
			h_PTmiss[3]->Fill(PTmiss_data, 1./6.);
			h_PT2miss[3]->Fill(PT2miss_data, 1./6.);
			h_gmom[3]->Fill(LV_gamma_out_data.P(), 1./6.);
			
			h_Q2[3]->Fill(Q2_data, 1./6.);
			h_epsilon[3]->Fill(epsilon_data, 1./6.);
			h_Xbj[3]->Fill(Xbj_data, 1./6.);
			h_W[3]->Fill(W_data, 1./6.);
			h_nu[3]->Fill(nu_data, 1./6.);
			h_mt[3]->Fill(mt_data, 1./6.);
			h_Pmom[3]->Fill(LV_proton_out_data.P(), 1./6.);
			h_elmom[3]->Fill(LV_el_out_data.P(), 1./6.);

			h2_Q2W[3]->Fill(Q2_data, W_data, 1./6.);
			h2_XQ2[3]->Fill(Xbj_data, Q2_data, 1./6.);
			h2_Q2mt[3]->Fill(Q2_data, mt_data, 1./6.);
			h2_nuep[3]->Fill(nu_data,epsilon_data, 1./6.);
			h2_Pemom[3]->Fill(LV_proton_out_data.E(), LV_proton_out_data.P(), 1./6.);
			h2_elemom[3]->Fill(LV_el_out_data.E(), LV_el_out_data.P(), 1./6.);
			h2_gemom[3]->Fill(LV_gamma_out_data.E(), LV_gamma_out_data.P(), 1./6.);
			
			if (process.compare("vcs")==0 || process.compare("pi0")==0){
				h_Phi[3]->Fill(Phi_data, 1./6.);
				h_ThCM[3]->Fill(ThCM_data, 1./6.);
				h_CosThCM[3]->Fill(CosThCM_data, 1./6.);
				h2_Q2Th[3]->Fill(Q2_data, ThCM_data, 1./6.);
				h2_WTh[3]->Fill(W_data, ThCM_data, 1./6.);
				h2_ThCMPhi[3]->Fill(Phi_data, ThCM_data, 1./6.);
				h2_mtTh[3]->Fill(mt_data, ThCM_data, 1./6.);
			}
		}

	}
	//gDirectory->ls("-m");
	//HallCTree->Print();
	file->Write();
	file->Close();
	
	out2.open(Form(VCS_REPLAY_PATH "/Ana/files/lumi_%d.dat", runID));
	out2 << runID << " "<< lumiexp_HMS << " "<<lumiexp_SHMS <<" "<< HMS_act_time << " "<< SHMS_act_time << " "<<HMS_B2_cur_cut<<" "<< SHMS_B2_cur_cut<<endl;
	out2.close();
  
  // Also write out luminosity information to the VCS database
  {
    const std::string db_fname = VCS_REPLAY_PATH "/Online/database/VCSdb.json";
    std::cout << "Writing luminosity to VCS database: " << db_fname << "\n";
    // step 1. --> load database to memory
    auto db = ReadJSON(db_fname);
    // step 2. --> Add info to database
    const std::string run_str = std::to_string(runID);
    db[run_str]["luminosity"]["HMS"]["luminosity"] = lumiexp_HMS;
    db[run_str]["luminosity"]["HMS"]["active_time"] = HMS_act_time;
    db[run_str]["luminosity"]["HMS"]["B2_current_cut"] = HMS_B2_cur_cut;
    db[run_str]["luminosity"]["SHMS"]["luminosity"] = lumiexp_SHMS;
    db[run_str]["luminosity"]["SHMS"]["active_time"] = SHMS_act_time;
    db[run_str]["luminosity"]["SHMS"]["B2_current_cut"] = SHMS_B2_cur_cut;
    // step 3. --> write database
    WriteJSON(db_fname, db);
  }
			
	return ;

}

////////////////////////////////////////////////////////////////////////

int ReadHallCData::InitHist(){

	h_encal_hms = new TH1F("h_encal_hms","energy cal HMS;E (GeV);events", hms_en_bins, hms_en_min_hist,hms_en_max_hist);
	h_encal_shms = new TH1F("h_encal_shms","energy cal SHMS;E (GeV);events", shms_en_bins, shms_en_min_hist, shms_en_max_hist);
	h_cer_npe_hms= new TH1F("h_cer_npe_hms","npe cer HMS;npe;events", (int)(hms_cer_npe_max_hist-hms_cer_npe_min_hist),hms_cer_npe_min_hist,hms_cer_npe_max_hist);
	h_aero_npe_shms= new TH1F("h_aero_npe_shms","npe aerogel HMS;npe;events", (int)(shms_aero_npe_max_hist-shms_aero_npe_min_hist), shms_aero_npe_min_hist, shms_aero_npe_max_hist);
	h_hgcer_npe_shms= new TH1F("h_hgcer_npe_shms",";npe hgcer;events", (int)(shms_hg_npe_max_hist-shms_hg_npe_min_hist), shms_hg_npe_min_hist, shms_hg_npe_max_hist);
	h_delta_hms= new TH1F("h_delta_hms","#Delta p HMS;#Delta p HMS;events",hms_delta_bins, hms_delta_min_hist, hms_delta_max_hist);
	h_delta_shms= new TH1F("h_delta_shms","#Delta p SHMS;# p SHMS;events", shms_delta_bins, shms_delta_min_hist, shms_delta_max_hist);
	h_phi_hms= new TH1F("h_phi_hms","#phi HMS;#phi;events", hms_phi_bins, hms_phi_min_hist, hms_phi_max_hist);
	h_theta_hms= new TH1F("h_theta_hms","#theta HMS;#theta;events", hms_th_bins, hms_th_min_hist, hms_th_max_hist);
	h_phi_shms= new TH1F("h_phi_shms","#phi SHMS;E (GeV);events", shms_phi_bins, shms_phi_min_hist, shms_phi_max_hist);
	h_theta_shms= new TH1F("h_theta_shms","#theta SHMS;E (GeV);events", shms_th_bins, shms_th_min_hist, shms_th_max_hist);
	h2_encal_cer= new TH2F("h_encal_cer","HMS npe cer vs energy cal;E (GeV);npe", hms_en_bins/2, hms_en_min_hist,hms_en_max_hist, (int)(hms_cer_npe_max_hist-hms_cer_npe_min_hist),hms_cer_npe_min_hist,hms_cer_npe_max_hist  );
	h2_aero_hgcer_npe= new TH2F("h_aero_hgcer_npe","SHMS npe hgcer vs aero;npe;npe", (int)(shms_hg_npe_max_hist-shms_hg_npe_min_hist), shms_hg_npe_min_hist, shms_hg_npe_max_hist,(int)(shms_aero_npe_max_hist-shms_aero_npe_min_hist), shms_aero_npe_min_hist, shms_aero_npe_max_hist);

	for (int i=0; i<5; i++){
	
		h_CTime_epCoinTime_ROC1[i] = new TH1F(Form("h_CTime_epCoinTime_ROC1[%d]",i),"coin time roc 1 - delay;time (ns);events", rel_time_bins, rel_time_min_hist,rel_time_max_hist);
		h_CTime_epCoinTime_ROC2[i] = new TH1F(Form("h_CTime_epCoinTime_ROC2[%d]",i),"coin time roc 2 - delay;time (ns);events", rel_time_bins, rel_time_min_hist,rel_time_max_hist);
		h_CTime_epCoinTime_ROC2_large[i] = new TH1F(Form("h_CTime_epCoinTime_ROC2_large[%d]",i),"coin time roc 2 absolute;time (ns);events", abs_time_bins, abs_time_min_hist,abs_time_max_hist);
		h_CTime_epCoinTime_TRIG1[i] = new TH1F(Form("h_CTime_epCoinTime_TRIG1[%d]",i),"trigger 1 time - delay;time (ns);events",trig_time_bins, trig_time_min_hist,trig_time_max_hist);
		h_CTime_epCoinTime_TRIG2[i] = new TH1F(Form("h_CTime_epCoinTime_TRIG2[%d]",i),"trigger 4 time - delay;time (ns);events", trig_time_bins, trig_time_min_hist,trig_time_max_hist);
		h_Pmom[i] = new TH1F(Form("h_Pmom[%d]",i),"HMS mom.;p_{P} (GeV);events",hms_mom_bins,hms_mom_min_hist, hms_mom_max_hist);	
		h_elmom[i] = new TH1F(Form("h_elmom[%d]",i),"SHMS mom.;p_{e} (GeV);events", shms_mom_bins, shms_mom_min_hist, shms_mom_max_hist);	
		h_gmom[i] = new TH1F(Form("h_gmom[%d]",i),"missing mom. (g);p_{miss} (GeV);events",miss_mom_bins, miss_mom_min_hist, miss_mom_max_hist);	
		h_Emiss[i] = new TH1F(Form("h_Emiss[%d]",i),"missing energy;E_{miss} (GeV);events",Emiss_bins, Emiss_min_hist, Emiss_max_hist);	
		h_Mmiss[i] = new TH1F(Form("h_Mmiss[%d]",i),"missing mass;M_{miss} (GeV);events",Mmiss_bins,Mmiss_min_hist,Mmiss_max_hist);	
		h_M2miss[i] = new TH1F(Form("h_M2miss[%d]",i),"missing mass sq;M^{2}_{miss} (GeV^{2});events",M2miss_bins,M2miss_min_hist,M2miss_max_hist);	
		h_PTmiss[i] = new TH1F(Form("h_PTmiss[%d]",i),"missing PT;P_{T,miss} (GeV);events",PTmiss_bins, PTmiss_min_hist, PTmiss_max_hist);	
		h_PT2miss[i] = new TH1F(Form("h_PT2miss[%d]",i),"missing PT2;P_{T,miss}^{2} (GeV^{2});events",PT2miss_bins, PT2miss_min_hist, PT2miss_max_hist);	
		h_Q2[i] = new TH1F(Form("h_Q2[%d]",i),"Q^{2};Q^{2} (GeV^{2});events", Q2_bins, Q2_min_hist, Q2_max_hist);
		h_epsilon[i] = new TH1F(Form("h_epsilon[%d]",i),"#epsilon;#epsilon;events", eps_bins, eps_min_hist, eps_max_hist);
		h_Xbj[i] = new TH1F(Form("h_Xbj[%d]",i),"x_{bj};x_{bj};events", xbj_bins, xbj_min_hist, xbj_max_hist);
		h_CosThCM[i] = new TH1F(Form("h_CosThCM[%d]",i),"cos(#Theta_{CM});cos(#Theta_{CM});events", costh_bins, costh_min_hist, costh_max_hist);
		h_ThCM[i] = new TH1F(Form("h_ThCM[%d]",i),"#Theta_{CM};#Theta_{CM} (deg.);events", th_bins, th_min_hist, th_max_hist);
		h_mt[i] = new TH1F(Form("h_mt[%d]",i),"-t;-t (GeV^{2});events", mt_bins, mt_min_hist, mt_max_hist);
		h_W[i] = new TH1F(Form("h_W[%d]",i),"W;W (GeV);events", W_bins, W_min_hist, W_max_hist);
		h_nu[i] = new TH1F(Form("h_nu[%d]",i),"#nu;#nu (GeV);events", nu_bins, nu_min_hist, nu_max_hist);
		h_Phi[i] = new TH1F(Form("h_Phi[%d]",i),"#Phi;#Phi (deg.);events", phi_bins, phi_min_hist, phi_max_hist);

		h2_Q2W[i] = new TH2F(Form("h2_Q2W[%d]",i),"W vs Q^{2};Q^{2} (GeV^{2}); W (GeV)",Q2_bins/2, Q2_min_hist, Q2_max_hist, W_bins/2, W_min_hist, W_max_hist);
		h2_XQ2[i] = new TH2F(Form("h2_XQ2[%d]",i),"Q^{2} vs x_{bj};x_{bj}; Q^{2} (GeV^{2})", xbj_bins/2, xbj_min_hist, xbj_max_hist, Q2_bins/2, Q2_min_hist, Q2_max_hist);
		h2_Q2mt[i] = new TH2F(Form("h2_Q2mt[%d]",i),"-t vs Q^{2};Q^{2} (GeV^{2}); -t (GeV^{2})",  Q2_bins/2, Q2_min_hist, Q2_max_hist, mt_bins/2, mt_min_hist, mt_max_hist);
		h2_ThCMPhi[i] = new TH2F(Form("h2_ThCMPhi[%d]",i),"#Theta_{CM} vs #Phi;#Phi (deg.);#Theta_{CM} (deg.)", phi_bins/2, phi_min_hist, phi_max_hist, th_bins/2, th_min_hist, th_max_hist);
		h2_nuep[i] = new TH2F(Form("h2_nuep[%d]",i),"#epsilon vs #nu;#nu (GeV);#epsilon",nu_bins/2, nu_min_hist, nu_max_hist, eps_bins/2, eps_min_hist, eps_max_hist);
		h2_Q2Th[i] = new TH2F(Form("h2_Q2Th[%d]",i),"#Theta_{CM} vs Q^{2};Q^{2} (GeV^{2});#Theta_{CM} (deg.)",Q2_bins/2, Q2_min_hist, Q2_max_hist,th_bins/2, th_min_hist, th_max_hist);
		h2_WTh[i] = new TH2F(Form("h2_WTh[%d]",i),"#Theta_{CM} vs W;W (GeV);#Theta_{CM} (deg.)", W_bins/2, W_min_hist, W_max_hist, th_bins/2, th_min_hist, th_max_hist);
		h2_mtTh[i] = new TH2F(Form("h2_mtTh[%d]",i),"#Theta_{CM} vs -t;-t (GeV2);#Theta_{CM} (deg.)",mt_bins/2, mt_min_hist, mt_max_hist, th_bins/2, th_min_hist, th_max_hist);

		h2_Pebeta[i] = new TH2F(Form("h2_Pebeta[%d]",i),"HMS #beta vs mom.;E (GeV); #beta", hms_mom_bins/2,hms_mom_min_hist, hms_mom_max_hist , hms_beta_bins/2, hms_beta_min_hist, hms_beta_max_hist);
		h2_elebeta[i] = new TH2F(Form("h2_elebeta[%d]",i),"SHMS #beta vs mom.;E (GeV); #beta", shms_mom_bins/2, shms_mom_min_hist, shms_mom_max_hist , shms_beta_bins/2, shms_beta_min_hist, shms_beta_max_hist);
		h2_COIN_P_beta[i] = new TH2F(Form("h2_COIN_P_beta[%d]",i),"SHMS: #beta vs ROC2 coin time - delay;time (ns); #beta", rel_time_bins, rel_time_min_hist,rel_time_max_hist, hms_beta_bins/2, hms_beta_min_hist, hms_beta_max_hist);
		h2_COIN_H_beta[i] = new TH2F(Form("h2_COIN_H_beta[%d]",i),"HMS: #beta vs ROC2 coin time - delay;time (ns); #beta", rel_time_bins, rel_time_min_hist,rel_time_max_hist, shms_beta_bins/2, shms_beta_min_hist, shms_beta_max_hist);
		h2_Pemom[i] = new TH2F(Form("h2_Pemom[%d]",i),"HMS p vs E;E (GeV); p (GeV)", hms_mom_bins/2,hms_mom_min_hist, hms_mom_max_hist ,hms_mom_bins/2,hms_mom_min_hist, hms_mom_max_hist );
		h2_elemom[i] = new TH2F(Form("h2_elemom[%d]",i),"SHMS p vs E;E (GeV); p (GeV)", shms_mom_bins/2, shms_mom_min_hist, shms_mom_max_hist ,shms_mom_bins/2, shms_mom_min_hist, shms_mom_max_hist);
		h2_gemom[i] = new TH2F(Form("h2_gemom[%d]",i),"missing p vs E;E_{miss} (GeV); p_{miss} (GeV)",miss_mom_bins/2, miss_mom_min_hist, miss_mom_max_hist,miss_mom_bins/2, miss_mom_min_hist, miss_mom_max_hist);
		h2_M2miss_CT1[i] = new TH2F(Form("h2_M2miss_CT1[%d]",i),"missing mass sq vs ROC1 time - delay;time (ns);M^{2}_{miss} (GeV^{2})", rel_time_bins, rel_time_min_hist,rel_time_max_hist, M2miss_bins,M2miss_min_hist,M2miss_max_hist);
		h2_M2miss_CT2[i] = new TH2F(Form("h2_M2miss_CT2[%d]",i),"missing mass sq vs ROC2 time - delay;time (ns);M^{2}_{miss} (GeV^{2})",rel_time_bins, rel_time_min_hist,rel_time_max_hist, M2miss_bins,M2miss_min_hist,M2miss_max_hist);

	}

	return 1;

}

/////////////////////////////////////////////////////////////////////////////

int ReadHallCData::DeleteHist(){

	h_encal_shms->Delete(); h_encal_hms->Delete(); 
	h_cer_npe_hms->Delete();
	h_aero_npe_shms->Delete();
	h_hgcer_npe_shms->Delete();
	h_delta_hms->Delete();
	h_delta_shms->Delete();
	h_phi_hms->Delete();
	h_theta_hms->Delete();
	h_phi_shms->Delete();
	h_theta_shms->Delete();
	
	for (int i=0; i<5; i++){
		h2_M2miss_CT1[i]->Delete();h2_M2miss_CT2[i]->Delete();
		h2_Q2W[i]->Delete(); h2_XQ2[i]->Delete(); h2_Q2mt[i]->Delete();  h2_ThCMPhi[i]->Delete();  
		h2_nuep[i]->Delete(); h2_Q2Th[i]->Delete(); h2_WTh[i]->Delete();  h2_mtTh[i]->Delete();
		h2_elebeta[i]->Delete(); h2_Pebeta[i]->Delete(); h2_COIN_P_beta[i]->Delete(); h2_COIN_H_beta[i]->Delete();
		h_CTime_epCoinTime_ROC1[i]->Delete();h_CTime_epCoinTime_ROC2[i]->Delete(); h_CTime_epCoinTime_TRIG1[i]->Delete(); h_CTime_epCoinTime_TRIG2[i]->Delete();
		h_Pmom[i]->Delete();  h_elmom[i]->Delete();  h_gmom[i]->Delete();  h_nu[i]->Delete(); h_Emiss[i]->Delete();  h_Mmiss[i]->Delete();  
		h_M2miss[i]->Delete();  h_PTmiss[i]->Delete();  h_PT2miss[i]->Delete();   
		h_Q2[i]->Delete();   h_epsilon[i]->Delete();  h_Xbj[i]->Delete();   h_CTime_epCoinTime_ROC2_large[i]->Delete(); 
		h_CosThCM[i]->Delete();  h_ThCM[i]->Delete();    h_mt[i]->Delete();    h_W[i]->Delete();   h_Phi[i]->Delete(); 
		h2_Pemom[i]->Delete(); h2_elemom[i]->Delete(); h2_gemom[i]->Delete();
	}
	return 1;
}


////////////////////////////////////////////

int ReadHallCData::DrawHist(string process, int run){

	TCanvas *c1 = new TCanvas();
	TSpectrum *spec = new TSpectrum(100);
	TSpectrum *smass = new TSpectrum(10);
	int nfounds =0, bin=0, npeaks=0, bmin, bmax;
	float xp, yp, integral;
	double gf_par[9], mass_integral[3], mass_sum1s[3];
	float pos_mem[3]={-1}, int_mem[3]={0}, max_mem[3]={0};
	int itag[3]={-1};
	float sum_bkg_up=0, sum_bkg_low=0, sum_bkg_av=0;
	vector <float> peak_pos, peak_int, peak_max;
	vector <float> peak_pos_mass, peak_int_mass, peak_max_mass;
	outfile.open(Form(VCS_REPLAY_PATH "/Ana/files/timing_%d.dat", run));

	// timing
	c1->cd(1);
	h_CTime_epCoinTime_ROC2[0]->Draw(); gPad->SetLogy(1); 
	h_CTime_epCoinTime_ROC2[0]->SetMinimum(1);
	h_CTime_epCoinTime_ROC2[1]->SetLineColor(8); h_CTime_epCoinTime_ROC2[1]->Draw("same");
	h_CTime_epCoinTime_ROC2[2]->SetLineColor(2); h_CTime_epCoinTime_ROC2[2]->Draw("same");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf(",run)); 
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf(",run)); 
	cout<<"start fitting timing peaks, can take time"<<endl;
	nfounds = spec->Search(h_CTime_epCoinTime_ROC2[2],1,"new");
	cout<<"found "<<nfounds<<" peaks "<<endl; c1->Update();
	
	double *xpeaks=spec->GetPositionX();
	TH1F *h2=(TH1F*) h_CTime_epCoinTime_ROC2[2]->Clone("h2");
	c1->Clear(); gPad->SetLogy(0);
	c1->cd(1); h2->Draw(""); gPad->SetLogy(0);
	cout<<"peaks positions (ns) and counts: "<<endl; 
	for (int i=0; i<nfounds; i++){
		xp=xpeaks[i];
		bin=h2->GetXaxis()->FindBin(xp);
		yp=h2->GetBinContent(bin);
		peak_max.push_back(yp);
		peak_pos.push_back(xpeaks[i]);
  	bmin = h2->GetXaxis()->FindBin(xpeaks[i]-1); 
    bmax = h2->GetXaxis()->FindBin(xpeaks[i]+1); 
    if (bmin>0 && bmax<h2->GetSize()-2) integral = h2->Integral(bmin,bmax);
		else integral = 0;
		peak_int.push_back(integral);
		cout<<" "<<xpeaks[i]<<" "<<integral<<" "<<yp<<endl;
		if (peak_int[i] > int_mem[0]) {
			max_mem[0] = peak_max[i];
			int_mem[0] = peak_int[i];
			pos_mem[0] = xpeaks[i];
			itag[0]=i;
		}
	} 
	for (int i=0; i<nfounds; i++){
		if (i==itag[0]) continue;
		if (peak_int[i] > int_mem[1]) {
			max_mem[1] = peak_max[i];
			int_mem[1] = peak_int[i];
			pos_mem[1] = xpeaks[i];
			itag[1]=i;
		}
	}
	for (int i=0; i<nfounds; i++){
		if (i==itag[0] || i==itag[1]) continue;
		if (peak_int[i] > int_mem[2]) {
			max_mem[2] = peak_max[i];
			int_mem[2] = peak_int[i];
			pos_mem[2] = xpeaks[i];
			itag[2]=i;
		}
	}
  bmin = h2->GetXaxis()->FindBin(xpeaks[itag[0]]+11); 
  bmax = h2->GetXaxis()->FindBin(xpeaks[itag[0]]+19); 
	if (xpeaks[itag[0]]+19<28) sum_bkg_up = (float) h2->Integral(bmin, bmax) / 4.;
  bmin = h2->GetXaxis()->FindBin(xpeaks[itag[0]]-19); 
  bmax = h2->GetXaxis()->FindBin(xpeaks[itag[0]]-11); 
	if (xpeaks[itag[0]]-19 > -28) sum_bkg_low = (float) h2->Integral(bmin, bmax) / 4.;
	if (sum_bkg_up>0 && sum_bkg_low>0) sum_bkg_av=(sum_bkg_up+sum_bkg_low)*0.5;
	else if (sum_bkg_up>0) sum_bkg_av=sum_bkg_up;
	else sum_bkg_av=sum_bkg_low;

	cout<<"main timing peaks:"<<endl;
	cout<<"1) delay: "<<pos_mem[0]<<" N events: "<<int_mem[0]<<" peak index "<< itag[0] << " peak max "<<max_mem[0]<<endl;
	cout<<"2) delay: "<<pos_mem[1]<<" N events: "<<int_mem[1]<<" peak index "<< itag[1] << " peak max "<<max_mem[1] <<endl;
	cout<<"3) delay: "<<pos_mem[2]<<" N events: "<<int_mem[2]<<" peak index "<< itag[2] << " peak max "<<max_mem[2] <<endl;
	cout<<"background up: "<<sum_bkg_up<<" background low: "<<sum_bkg_low<<endl;
	cout<<"N events in each peak: "<<int_mem[0]-sum_bkg_av<<" "<<int_mem[1]-sum_bkg_av<<" "<<int_mem[2]-sum_bkg_av<<endl;
  cout<<"absolute time positions: "<<pos_mem[0]+time_shift<<" "<<pos_mem[1]+time_shift<<" "<<pos_mem[2]+time_shift<<endl;
	outfile<< run<<" "<<pos_mem[0]<<" "<<int_mem[0]-sum_bkg_av<<" "<<int_mem[0]<<" "<<sum_bkg_up<<" "<<sum_bkg_low<<" "<<max_mem[0]<<" "<<(int_mem[0]-sum_bkg_av)/lumiexp_HMS<<endl;
	outfile<< run<<" "<<pos_mem[1]<<" "<<int_mem[1]-sum_bkg_av<<" "<<int_mem[1]<<" "<<sum_bkg_up<<" "<<sum_bkg_low<<" "<<max_mem[1]<<" "<<(int_mem[1]-sum_bkg_av)/lumiexp_HMS<<endl;
	outfile<< run<<" "<<pos_mem[2]<<" "<<int_mem[2]-sum_bkg_av<<" "<<int_mem[2]<<" "<<sum_bkg_up<<" "<<sum_bkg_low<<" "<<max_mem[2]<<" "<<(int_mem[2]-sum_bkg_av)/lumiexp_HMS<<endl;
	outfile.close();

  // Also write out timing information to the VCS database
  {
    const std::string db_fname = VCS_REPLAY_PATH  "/Online/database/VCSdb.json";
    std::cout << "Writing timing info to VCS database: " << db_fname << "\n";
    // step 1. --> load database to memory
    auto db = ReadJSON(db_fname);
    // step 2. --> Add info to database
    const std::string run_str = std::to_string(run);
    std::vector<std::string> peak_names = {"peak1", "peak2", "peak3"};
    for (int i = 0; i < 3; ++i) {
      db[run_str]["timing"][peak_names[i]]["delay"] = pos_mem[i];
      db[run_str]["timing"][peak_names[i]]["events"] = int_mem[i] - sum_bkg_av;
      db[run_str]["timing"][peak_names[i]]["raw_events"] = int_mem[i];
      db[run_str]["timing"][peak_names[i]]["background_avg"] = sum_bkg_av;
      db[run_str]["timing"][peak_names[i]]["background_low"] = sum_bkg_low;
      db[run_str]["timing"][peak_names[i]]["background_high"] = sum_bkg_up;
      db[run_str]["timing"][peak_names[i]]["max"] = max_mem[i];
      db[run_str]["timing"][peak_names[i]]["N_over_L"] = (int_mem[i]-sum_bkg_av)/lumiexp_HMS;
    }
    // step 3. --> write database
    WriteJSON(db_fname, db);
  }

	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf",run)); c1->Clear(); c1->cd(1);
	h_CTime_epCoinTime_ROC1[0]->Draw(); gPad->SetLogy(1);
	h_CTime_epCoinTime_ROC1[0]->SetMinimum(1);
	h_CTime_epCoinTime_ROC1[1]->SetLineColor(8); h_CTime_epCoinTime_ROC1[1]->Draw("same");
	h_CTime_epCoinTime_ROC1[2]->SetLineColor(2); h_CTime_epCoinTime_ROC1[2]->Draw("same");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf",run)); c1->Clear(); c1->cd(1);
	h_CTime_epCoinTime_ROC2_large[0]->Draw(); gPad->SetLogy(0); 
	h_CTime_epCoinTime_ROC2_large[0]->SetMinimum(0);
  c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf",run)); 
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->Clear(); c1->cd(1);
	h_CTime_epCoinTime_TRIG1[0]->Draw(); gPad->SetLogy(1); h_CTime_epCoinTime_TRIG1[0]->SetMinimum(1);
	h_CTime_epCoinTime_TRIG1[0]->SetMinimum(1);
	h_CTime_epCoinTime_TRIG1[1]->SetLineColor(8); h_CTime_epCoinTime_TRIG1[1]->Draw("same");
	h_CTime_epCoinTime_TRIG1[2]->SetLineColor(2); h_CTime_epCoinTime_TRIG1[2]->Draw("same");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf",run)); c1->Clear(); c1->cd(1);
	h_CTime_epCoinTime_TRIG2[0]->Draw(); gPad->SetLogy(1); h_CTime_epCoinTime_TRIG2[0]->SetMinimum(1);
	h_CTime_epCoinTime_TRIG2[0]->SetMinimum(1);
	h_CTime_epCoinTime_TRIG2[1]->SetLineColor(8); h_CTime_epCoinTime_TRIG2[1]->Draw("same");
	h_CTime_epCoinTime_TRIG2[2]->SetLineColor(2); h_CTime_epCoinTime_TRIG2[2]->Draw("same");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf",run)); c1->Clear(); 
	c1->Divide(1,2); 
	c1->cd(1); h2_COIN_H_beta[1]->Draw("colz"); gPad->SetLogy(0);
	c1->cd(2); h2_COIN_H_beta[2]->Draw("colz"); gPad->SetLogy(0);
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf",run)); c1->Clear();
	c1->Divide(1,2); 
	c1->cd(1); h2_COIN_P_beta[1]->Draw("colz"); gPad->SetLogy(0);
	c1->cd(2); h2_COIN_P_beta[2]->Draw("col"); gPad->SetLogy(0);
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf",run)); 
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->Clear(); c1->Divide(1,2); 
	c1->cd(1); h2_M2miss_CT1[0]->Draw("colz"); gPad->SetLogy(0);
	c1->cd(2); h2_M2miss_CT2[0]->Draw("colz"); gPad->SetLogy(0);
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/cointime_%d.pdf)",run)); 
	c1->Clear(); 

	// positions
	c1->Divide(2,2);
	c1->cd(1); h_encal_hms->Draw();
	c1->cd(2); h_cer_npe_hms->Draw();
	c1->cd(3); h2_encal_cer->Draw("colz");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/positions_%d.pdf(",run));
	c1->Clear(); c1->Divide(2,2);
	c1->cd(1); h_encal_shms->Draw();
	c1->cd(2); h_aero_npe_shms->Draw();
	c1->cd(3); h_hgcer_npe_shms->Draw();
	c1->cd(4); h2_aero_hgcer_npe->Draw("colz");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/positions_%d.pdf",run));
	c1->Clear(); c1->Divide(2,2);
	c1->cd(1); h_delta_hms->Draw();
	c1->cd(2); h_phi_hms->Draw();
	c1->cd(3); h_theta_hms->Draw();
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/positions_%d.pdf",run));
	c1->Clear(); c1->Divide(2,2);
	c1->cd(1); h_delta_shms->Draw();
	c1->cd(2); h_phi_shms->Draw();
	c1->cd(3); h_theta_shms->Draw();
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/positions_%d.pdf)",run));

	// kinematics
	for (int i=2;i<=4;i++){
		if (!h_Mmiss[i]->GetSumw2N()) h_Mmiss[i]->Sumw2();
		if (!h_M2miss[i]->GetSumw2N()) h_M2miss[i]->Sumw2();
		if (!h_Emiss[i]->GetSumw2N())	h_Emiss[i]->Sumw2();
		if (!h_PTmiss[i]->GetSumw2N())	h_PTmiss[i]->Sumw2();
		if (!h_PT2miss[i]->GetSumw2N())	h_PT2miss[i]->Sumw2();
		if (!h_gmom[i]->GetSumw2N())	h_gmom[i]->Sumw2();
		if (!h_Q2[i]->GetSumw2N())	h_Q2[i]->Sumw2();
		if (!h_epsilon[i]->GetSumw2N())	h_epsilon[i]->Sumw2();
		if (!h_Xbj[i]->GetSumw2N())	h_Xbj[i]->Sumw2();
		if (!h_mt[i]->GetSumw2N())	h_mt[i]->Sumw2();
		if (!h_nu[i]->GetSumw2N())	h_nu[i]->Sumw2();
		if (!h_W[i]->GetSumw2N())	h_W[i]->Sumw2();
		if (!h_Pmom[i]->GetSumw2N())	h_Pmom[i]->Sumw2();
		if (!h_elmom[i]->GetSumw2N())	h_elmom[i]->Sumw2();
		//if (!h2_Q2W[i]->GetSumw2N())	h2_Q2W[i]->Sumw2();
		//if (!h2_Q2mt[i]->GetSumw2N())	h2_Q2mt[i]->Sumw2();
		//if (!h2_XQ2[i]->GetSumw2N())	h2_XQ2[i]->Sumw2();
		//if (!h2_nuep[i]->GetSumw2N())	h2_nuep[i]->Sumw2();
		//if (!h2_Pemom[i]->GetSumw2N())	h2_Pemom[i]->Sumw2();
		//if (!h2_gemom[i]->GetSumw2N())	h2_gemom[i]->Sumw2();
		//if (!h2_elemom[i]->GetSumw2N())	h2_elemom[i]->Sumw2();
		if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		){	
	    	    //if (!h2_ThCMPhi[i]->GetSumw2N()) h2_ThCMPhi[i]->Sumw2();
			//if (!h2_mtTh[i]->GetSumw2N()) h2_mtTh[i]->Sumw2();
			if (!h_Phi[i]->GetSumw2N()) h_Phi[i]->Sumw2();
			if (!h_ThCM[i]->GetSumw2N()) h_ThCM[i]->Sumw2();
			if (!h_CosThCM[i]->GetSumw2N()) h_CosThCM[i]->Sumw2();
			//if (!h2_Q2Th[i]->GetSumw2N()) h2_Q2Th[i]->Sumw2();
			//if (!h2_WTh[i]->GetSumw2N()) h2_WTh[i]->Sumw2();
		}
	}

	h_Mmiss[4]->Add(h_Mmiss[3],-1);
	h_M2miss[4]->Add(h_M2miss[3],-1);
	h_Emiss[4]->Add(h_Emiss[3],-1);
	h_PTmiss[4]->Add(h_PTmiss[3],-1);
	h_PT2miss[4]->Add(h_PT2miss[3],-1);
	h_gmom[4]->Add(h_gmom[3],-1);
	h_Q2[4]->Add(h_Q2[3],-1);	
	h_epsilon[4]->Add(h_epsilon[3],-1);	
	h_Xbj[4]->Add(h_Xbj[3],-1);	
	h_mt[4]->Add(h_mt[3],-1);	
	h_nu[4]->Add(h_nu[3],-1);	
	h_W[4]->Add(h_W[3],-1);	
	h_Pmom[4]->Add(h_Pmom[3],-1);	
	h_elmom[4]->Add(h_elmom[3],-1);	
	h2_Q2W[4]->Add(h2_Q2W[3],-1);	
	h2_Q2mt[4]->Add(h2_Q2mt[3],-1);	
	h2_XQ2[4]->Add(h2_XQ2[3],-1);	
	h2_nuep[4]->Add(h2_nuep[3],-1);	
	h2_Pemom[4]->Add(h2_Pemom[3],-1);	
	h2_gemom[4]->Add(h2_gemom[3],-1);	
	h2_elemom[4]->Add(h2_elemom[3],-1);	
	if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		){	
		h2_Q2Th[4]->Add(h2_Q2Th[3],-1);	
		h2_WTh[4]->Add(h2_WTh[3],-1);	
		h2_ThCMPhi[4]->Add(h2_ThCMPhi[3],-1);	
		h2_mtTh[4]->Add(h2_mtTh[3],-1);	
		h_Phi[4]->Add(h_Phi[3],-1);	
		h_ThCM[4]->Add(h_ThCM[3],-1);	
		h_CosThCM[4]->Add(h_CosThCM[3],-1);	
	}
	for (int i=2;i<=4;i++){
		h2_Q2W[i]->SetMinimum(1);
		h2_Q2mt[i]->SetMinimum(1);
		h2_XQ2[i]->SetMinimum(1);
		h2_nuep[i]->SetMinimum(1);
		h2_Pemom[i]->SetMinimum(1);
		h2_gemom[i]->SetMinimum(1);
		h2_elemom[i]->SetMinimum(1);
		if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		){	
			h2_ThCMPhi[i]->SetMinimum(1);
			h2_mtTh[i]->SetMinimum(1);
			h2_Q2Th[i]->SetMinimum(1);
			h2_WTh[i]->SetMinimum(1);
		}
	}

	c1->Clear(); 
	c1->Divide(2,2); 
	c1->cd(1); 
	h_Q2[2]->SetLineColor(9); h_Q2[2]->Draw("HIST"); gPad->SetLogy(0); h_Q2[2]->SetMinimum(0); 
	h_Q2[4]->SetLineColor(2); h_Q2[4]->Draw("HISTsame");
	h_Q2[3]->SetLineColor(8); h_Q2[3]->Draw("HISTsame");
	c1->cd(2); 
	h_Xbj[2]->SetLineColor(9); h_Xbj[2]->Draw("HIST"); gPad->SetLogy(0); h_Xbj[2]->SetMinimum(0);
	h_Xbj[4]->SetLineColor(2); h_Xbj[4]->Draw("HISTsame");
	h_Xbj[3]->SetLineColor(8); h_Xbj[3]->Draw("HISTsame");
	c1->cd(3); 
	h_mt[2]->SetLineColor(9); h_mt[2]->Draw("HIST"); gPad->SetLogy(0); h_mt[2]->SetMinimum(0);
	h_mt[4]->SetLineColor(2); h_mt[4]->Draw("HISTsame");
	h_mt[3]->SetLineColor(8); h_mt[3]->Draw("HISTsame");
	c1->cd(4); 
	h_W[2]->SetLineColor(9); h_W[2]->Draw("HIST"); gPad->SetLogy(0); h_W[2]->SetMinimum(0);
	h_W[4]->SetLineColor(2); h_W[4]->Draw("HISTsame");
	h_W[3]->SetLineColor(8); h_W[3]->Draw("HISTsame");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics_%d.pdf(",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->Clear(); 
	c1->Divide(2,2); 
	c1->cd(1); 
	h_nu[2]->SetLineColor(9); h_nu[2]->Draw("HIST"); gPad->SetLogy(0); h_nu[2]->SetMinimum(0);
	h_nu[4]->SetLineColor(2); h_nu[4]->Draw("HISTsame");
	h_nu[3]->SetLineColor(8); h_nu[3]->Draw("HISTsame");
	c1->cd(2); 
	h_epsilon[2]->SetLineColor(9); h_epsilon[2]->Draw("HIST"); gPad->SetLogy(0); h_epsilon[2]->SetMinimum(0);
	h_epsilon[4]->SetLineColor(2); h_epsilon[4]->Draw("HISTsame");
	h_epsilon[3]->SetLineColor(8); h_epsilon[3]->Draw("HISTsame");
	c1->cd(3); 
	h_Pmom[2]->SetLineColor(9); h_Pmom[2]->Draw("HIST"); gPad->SetLogy(0); h_Pmom[2]->SetMinimum(0);
	h_Pmom[4]->SetLineColor(2); h_Pmom[4]->Draw("HISTsame");
	h_Pmom[3]->SetLineColor(8); h_Pmom[3]->Draw("HISTsame");
	c1->cd(4); 
	h_elmom[2]->SetLineColor(9); h_elmom[2]->Draw("HIST"); gPad->SetLogy(0); h_elmom[2]->SetMinimum(0);
	h_elmom[4]->SetLineColor(2); h_elmom[4]->Draw("HISTsame");
	h_elmom[3]->SetLineColor(8); h_elmom[3]->Draw("HISTsame");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics_%d.pdf",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->Clear(); 
	if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		){	
		c1->Divide(2,2); 
		c1->cd(1); 
		h_Phi[2]->SetLineColor(9); h_Phi[2]->Draw("HIST"); gPad->SetLogy(0); h_Phi[2]->SetMinimum(0);
		h_Phi[4]->SetLineColor(2); h_Phi[4]->Draw("HISTsame");
		h_Phi[3]->SetLineColor(8); h_Phi[3]->Draw("HISTsame");
		c1->cd(2); 
		h_ThCM[2]->SetLineColor(9); h_ThCM[2]->Draw("HIST"); gPad->SetLogy(0); h_ThCM[2]->SetMinimum(0);
		h_ThCM[4]->SetLineColor(2); h_ThCM[4]->Draw("HISTsame");
		h_ThCM[3]->SetLineColor(8); h_ThCM[3]->Draw("HISTsame");
		c1->cd(3); 
		h_CosThCM[2]->SetLineColor(9); h_CosThCM[2]->Draw("HIST"); gPad->SetLogy(0); h_CosThCM[2]->SetMinimum(0);
		h_CosThCM[4]->SetLineColor(2); h_CosThCM[4]->Draw("HISTsame");
		h_CosThCM[3]->SetLineColor(8); h_CosThCM[3]->Draw("HISTsame");
		c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics_%d.pdf",run));
		c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
		c1->Clear(); 
	}
	c1->Divide(2,2); 
	c1->cd(1); h_Q2[0]->Draw("HIST");  h_Q2[1]->Draw("HISTsame"); 
	c1->cd(2); h_Xbj[0]->Draw("HIST"); h_Xbj[1]->Draw("HISTsame");
	c1->cd(3); h_mt[0]->Draw("HIST"); h_mt[1]->Draw("HISTsame");
	c1->cd(4); h_W[0]->Draw("HIST"); h_W[1]->Draw("HISTsame");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics_%d.pdf",run));
	c1->Clear(); 
	c1->Divide(2,2); 
	c1->cd(1); h_nu[0]->Draw("HIST"); h_nu[1]->Draw("HISTsame");
	c1->cd(2); h_epsilon[0]->Draw("HIST"); h_epsilon[1]->Draw("HISTsame");
	c1->cd(3); h_Pmom[0]->Draw("HIST");
	c1->cd(4); h_elmom[0]->Draw("HIST");
	if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		){	
		c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics_%d.pdf",run));
		c1->Clear(); 
		c1->Divide(2,2); 
		c1->cd(1); h_Phi[1]->Draw("HIST");
		c1->cd(2); h_ThCM[1]->Draw("HIST"); 
		c1->cd(3); h_CosThCM[1]->Draw("HIST");
	}
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics_%d.pdf)",run));
	c1->Clear();

	// missing mass...
	c1->Clear();c1->Divide(2,2);
	c1->cd(1); 
	h_M2miss[4]->SetLineColor(2); h_M2miss[4]->Draw("E0"); gPad->SetLogy(0); h_M2miss[4]->SetMinimum(0); 
	c1->cd(2); 
	h_Mmiss[4]->SetLineColor(2); h_Mmiss[4]->Draw("E0");gPad->SetLogy(0); h_Mmiss[4]->SetMinimum(0);
	c1->cd(3); 
	h_Emiss[4]->SetLineColor(2); h_Emiss[4]->Draw("E0"); gPad->SetLogy(0); h_Emiss[4]->SetMinimum(0); 
	c1->cd(4); 
	h_PTmiss[4]->SetLineColor(2); h_PTmiss[4]->Draw("E0"); gPad->SetLogy(0); h_PTmiss[4]->SetMinimum(0);
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/exclusivity_%d.pdf(",run));
	c1->Clear();c1->Divide(2,2);
	c1->cd(1); 
	h_M2miss[2]->SetLineColor(6); h_M2miss[2]->Draw("HIST"); gPad->SetLogy(0); h_M2miss[2]->SetMinimum(0);
	h_M2miss[4]->Draw("HISTsame"); h_M2miss[3]->SetLineColor(8); h_M2miss[3]->Draw("HISTsame");
	c1->cd(2); 
	h_Mmiss[2]->SetLineColor(6); h_Mmiss[2]->Draw("HIST"); gPad->SetLogy(0); h_Mmiss[2]->SetMinimum(0);
	h_Mmiss[4]->Draw("HISTsame"); h_Mmiss[3]->SetLineColor(8); h_Mmiss[3]->Draw("HISTsame"); 
	c1->cd(3); 
	h_Emiss[2]->SetLineColor(6); h_Emiss[2]->Draw("HIST"); gPad->SetLogy(0); h_Emiss[2]->SetMinimum(0);
	h_Emiss[4]->Draw("HISTsame"); h_Emiss[3]->SetLineColor(8); h_Emiss[3]->Draw("HISTsame"); 
	c1->cd(4); 
	h_PTmiss[2]->SetLineColor(6); h_PTmiss[2]->Draw("HIST"); gPad->SetLogy(0); h_PTmiss[2]->SetMinimum(0);
	h_PTmiss[4]->Draw("HISTsame"); h_PTmiss[3]->SetLineColor(8); h_PTmiss[3]->Draw("HISTsame"); 
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/exclusivity_%d.pdf",run));	
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/exclusivity_%d.root",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->Clear(); c1->Divide(2,2);
	c1->cd(1); 
	h_PT2miss[4]->SetLineColor(2); h_PT2miss[4]->Draw("E0"); gPad->SetLogy(0); h_PT2miss[4]->SetMinimum(0);
	c1->cd(2); 
	h_gmom[4]->SetLineColor(2); h_gmom[4]->Draw("E0"); gPad->SetLogy(0); h_gmom[4]->SetMinimum(0);
	c1->cd(3); 
	h_PT2miss[2]->SetLineColor(6); h_PT2miss[2]->Draw("HIST");  gPad->SetLogy(0);
	h_PT2miss[4]->Draw("HISTsame"); h_PT2miss[3]->SetLineColor(8); h_PT2miss[3]->Draw("HISTsame"); 
	c1->cd(4); 
	h_gmom[2]->SetLineColor(6); h_gmom[2]->Draw("HIST"); gPad->SetLogy(0); h_gmom[2]->SetMinimum(0);
	h_gmom[4]->Draw("HISTsame"); h_gmom[3]->SetLineColor(8); h_gmom[3]->Draw("HISTsame"); 
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/exclusivity_%d.pdf",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	// fit of exclusivity
	cout<<"fit of missing mass"<<endl;
	nfounds = 0; npeaks = 0; 
	for (int i=0;i<3;i++){
		pos_mem[i]=-1; int_mem[i]=0; max_mem[i]=0;	
	}
	nfounds = smass->Search(h_M2miss[4],1,"new");
	cout<<"found "<<nfounds<<" peaks in missing mass distribution "<<endl;
	double *xpM=smass->GetPositionX();
	TH1F *hM= (TH1F*) h_M2miss[4]->Clone("hM");
	c1->Clear(); c1->Divide(1,1); c1->cd(1); hM->SetLineColor(1); hM->Draw("HIST");gPad->SetLogy(0);
	for (int i=0; i<nfounds; i++){
		xp=xpM[i];
	        bin=hM->GetXaxis()->FindBin(xp);
		yp=hM->GetBinContent(bin);
		peak_max_mass.push_back(yp);
		peak_pos_mass.push_back(xpM[i]);
		bmin = hM->GetXaxis()->FindBin(xpM[i]-0.01);
		bmax = hM->GetXaxis()->FindBin(xpM[i]+0.01);
		if (bmin>0 && bmax<hM->GetSize()-2) integral = hM->Integral(bmin,bmax);
		peak_int_mass.push_back(integral);
		cout<<" "<<xpM[i]<<" "<<integral<<" "<<yp<<endl;
		if (peak_int_mass[i] > int_mem[0]) {
			max_mem[0] = peak_max_mass[i];
			int_mem[0] = peak_int_mass[i];
			pos_mem[0] = xpM[i];
			itag[0]=i;
		}
	}
	if (nfounds>1){
	    for (int i=0; i<nfounds; i++){
		if (i==itag[0]) continue;
		if (peak_int_mass[i] > int_mem[1]) {
			max_mem[1] = peak_max_mass[i];
			int_mem[1] = peak_int_mass[i];
			pos_mem[1] = xpM[i];
			itag[1]=i;
		}
	    }
	}
	if (nfounds>2){
	   for (int i=0; i<nfounds; i++){
		if (i==itag[0] || i==itag[1]) continue;
		if (peak_int_mass[i] > int_mem[2]) {
			max_mem[2] = peak_max_mass[i];
			int_mem[2] = peak_int_mass[i];
			pos_mem[2] = xpM[i];
			itag[2]=i;
		}
	   }
	}
	
	TF1 *fgaus =new TF1("fgaus","[0] * exp(-(x-[1])*(x-[1])/2./[2]/[2])", pos_mem[0]-0.1,pos_mem[0]+0.1);
	TF1 *fgaus2 =new TF1("fgaus2","[0] * exp(-(x-[1])*(x-[1])/2./[2]/[2])", pos_mem[1]-0.1,pos_mem[1]+0.1);
	TF1 *fgaus3 =new TF1("fgaus3","[0] * exp(-(x-[1])*(x-[1])/2./[2]/[2])", pos_mem[2]-0.1,pos_mem[2]+0.1);
	fgaus->SetParameter(0,max_mem[0]); 
	fgaus->SetParameter(1,pos_mem[0]); 
	fgaus->SetParameter(2,0.005); 
	fgaus->SetParLimits(1,pos_mem[0]-0.015,pos_mem[0]+0.015);
	fgaus->SetParLimits(0,max_mem[0]/3.,max_mem[0]*2);
	fgaus->SetParLimits(2,0.0001, 0.1);
	if (nfounds>1){
		fgaus2->SetParameter(0,max_mem[1]); 
		fgaus2->SetParameter(1,pos_mem[1]); 
		fgaus2->SetParameter(2,0.005); 
		fgaus2->SetParLimits(1,pos_mem[1]-0.015,pos_mem[1]+0.015);
		fgaus2->SetParLimits(0,max_mem[1]/3.,max_mem[1]*2);
		fgaus2->SetParLimits(2,0.0001, 0.1);
	}	
	if (nfounds>2){
		fgaus3->SetParameter(0,max_mem[2]); 
		fgaus3->SetParameter(1,pos_mem[2]); 
		fgaus3->SetParameter(2,0.005); 
		fgaus3->SetParLimits(1,pos_mem[2]-0.015,pos_mem[2]+0.015);
		fgaus3->SetParLimits(0,max_mem[2]/3.,max_mem[2]*2);
		fgaus3->SetParLimits(2,0.0001, 0.1);
	}
		
	hM->Fit("fgaus", "ER");
	fgaus->Draw("same");
	gf_par[0] = (float) fgaus->GetParameter(0);
	gf_par[1] = (float) fgaus->GetParameter(1);
	gf_par[2] = (float) fgaus->GetParameter(2);
	//integral =(float) fgaus->Integral((double) (gf_par[1]-3*gf_par[2]), (double) (gf_par[1]+3*gf_par[2]), 1e-6);
	bmin = hM->GetXaxis()->FindBin((gf_par[1]-gf_par[2]));
	bmax = hM->GetXaxis()->FindBin((gf_par[1]+gf_par[2]));
	mass_integral[0] = gf_par[0]*sqrt(PI*(2.*pow(gf_par[2],2))) *(M2miss_max_hist-M2miss_min_hist)/((float) M2miss_bins); // adjust bin width of hM
	if (bmin>0 && bmax<hM->GetSize()-2) mass_sum1s[0] = hM->Integral(bmin,bmax)/0.68;
	cout<<" gaus mmass fit par "<<gf_par[0]<<" "<<gf_par[1]<<" "<<gf_par[2]<< " 3sig integral "<<mass_integral[0]<<" sum 1 sigm "<<mass_sum1s[0]<<endl;
	
	if (nfounds>1){
		hM->Fit("fgaus2", "ER");
		fgaus3->Draw("same");
		gf_par[3] = (float) fgaus2->GetParameter(0);
		gf_par[4] = (float) fgaus2->GetParameter(1);
		gf_par[5] = (float) fgaus2->GetParameter(2);
		bmin = hM->GetXaxis()->FindBin((gf_par[4]-gf_par[5]));
		bmax = hM->GetXaxis()->FindBin((gf_par[4]+gf_par[5]));
		mass_integral[1] = gf_par[3]*sqrt(PI*(2.*pow(gf_par[5],2)))  *(M2miss_max_hist-M2miss_min_hist)/((float) M2miss_bins);
		if (bmin>0 && bmax<hM->GetSize()-2) mass_sum1s[1] = hM->Integral(bmin,bmax)/0.68;
		cout<<" gaus mmass fit par "<<gf_par[3]<<" "<<gf_par[4]<<" "<<gf_par[5]<< " 3sig integral "<<mass_integral[1]<<" sum 1 sigm "<<mass_sum1s[1]<<endl;
	}
	if (nfounds>2){
		hM->Fit("fgaus3", "ER");
		fgaus3->Draw("same");
		gf_par[6] = (float) fgaus2->GetParameter(0);
		gf_par[7] = (float) fgaus2->GetParameter(1);
		gf_par[8] = (float) fgaus2->GetParameter(2);
		bmin = hM->GetXaxis()->FindBin((gf_par[7]-gf_par[8]));
		bmax = hM->GetXaxis()->FindBin((gf_par[7]+gf_par[8]));
		mass_integral[2] = gf_par[6]*sqrt(PI*(2.*pow(gf_par[8],2))) *(M2miss_max_hist-M2miss_min_hist)/((float) M2miss_bins);  
		if (bmin>0 && bmax<hM->GetSize()-2) mass_sum1s[2] = hM->Integral(bmin,bmax)/0.68;
		cout<<" gaus mmass fit par "<<gf_par[6]<<" "<<gf_par[7]<<" "<<gf_par[8]<< " integral "<<mass_integral[2]<<" sum with bkg "<<mass_sum1s[2]<<endl;
	}
	
	out3.open(Form(VCS_REPLAY_PATH "/Ana/files/missmass2_%d.dat", run));
	cout<<"main mass peaks"<<endl;
	cout<<"1) mass : "<<pos_mem[0]<<" "<<" peak sum (tot)= "<<int_mem[0]<<" "<<mass_integral[0]<<" "<<mass_sum1s[0]<<" peak index= "<<itag[0]<<" peak max= "<<max_mem[0]<<endl;
	cout<<"2) mass : "<<pos_mem[1]<<" "<<" peak sum (tot)= "<<int_mem[1]<<" "<<mass_integral[1]<<" "<<mass_sum1s[1]<<" peak index= "<<itag[1]<<" peak max= "<<max_mem[1]<<endl;
	cout<<"3) mass : "<<pos_mem[2]<<" "<<" peak sum (tot)= "<<int_mem[2]<<" "<<mass_integral[2]<<" "<<mass_sum1s[2]<<" peak index= "<<itag[2]<<" peak max= "<<max_mem[2]<<endl;
	out3<<run<<" "<<pos_mem[0]<<" "<<mass_integral[0]<<" "<<mass_sum1s[0]<<" "<< mass_integral[0]/lumiexp_HMS <<endl;
	out3<<run<<" "<<pos_mem[1]<<" "<<mass_integral[1]<<" "<<mass_sum1s[1]<<" "<< mass_integral[1]/lumiexp_HMS <<endl;
	out3<<run<<" "<<pos_mem[2]<<" "<<mass_integral[2]<<" "<<mass_sum1s[2]<<" "<< mass_integral[2]/lumiexp_HMS <<endl;
	out3.close();

  // Also write out timing missing mass info to the VCS database
  {
    const std::string db_fname = VCS_REPLAY_PATH "/Online/database/VCSdb.json";
    std::cout << "Writing MM info to VCS database: " << db_fname << "\n";
    // step 1. --> load database to memory
    auto db = ReadJSON(db_fname);
    // step 2. --> Add info to database
    const std::string run_str = std::to_string(run);
    std::vector<std::string> peak_names = {"peak1", "peak2", "peak3"};
    for (int i = 0; i < 3; ++i) {
      db[run_str]["missing_mass"][peak_names[i]]["mass"] = pos_mem[i];
      db[run_str]["missing_mass"][peak_names[i]]["integral"] = mass_integral[i];
      db[run_str]["missing_mass"][peak_names[i]]["sum_withbkg"] = mass_sum1s[i];
      db[run_str]["missing_mass"][peak_names[i]]["sum_arbcut"] = int_mem[i];
      db[run_str]["missing_mass"][peak_names[i]]["max"] = max_mem[i];
    }
    // step 3. --> write database
    WriteJSON(db_fname, db);
  }

	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/exclusivity_%d.pdf)",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 


	// kinematics 2D plots
	c1->Clear();
	c1->Divide(2,2); 
	c1->cd(1); h2_Q2W[4]->Draw("HISTcolz"); 
	c1->cd(2); h2_XQ2[4]->Draw("HISTcolz");
	c1->cd(3); h2_Q2mt[4]->Draw("HISTcolz");
	c1->cd(4); h2_nuep[4]->Draw("HISTcolz");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf(",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->Clear();
	if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		){	
		c1->Divide(2,2); 
		c1->cd(1); h2_Q2Th[4]->Draw("HISTcolz"); 
		c1->cd(2); h2_WTh[4]->Draw("HISTcolz");
		c1->cd(3); h2_ThCMPhi[4]->Draw("HISTcolz");
		c1->cd(4); h2_mtTh[4]->Draw("HISTcolz");
		c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf",run));
		c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
		c1->Clear();
	}
	c1->Divide(2,1); 
	c1->cd(1); h2_elemom[4]->Draw("HISTcolz");
	c1->cd(2); h2_Pemom[4]->Draw("HISTcolz");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->Clear();
	c1->Divide(2,2); 
	c1->cd(1); h2_Q2W[2]->Draw("HISTcolz"); 
	c1->cd(2); h2_XQ2[2]->Draw("HISTcolz");
	c1->cd(3); h2_Q2mt[2]->Draw("HISTcolz");
	c1->cd(4); h2_nuep[2]->Draw("HISTcolz");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
	c1->Clear();
	if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		){	
		c1->Divide(2,2); 
		c1->cd(1); h2_Q2Th[2]->Draw("HISTcolz"); 
		c1->cd(2); h2_WTh[2]->Draw("HISTcolz");
		c1->cd(3); h2_ThCMPhi[2]->Draw("HISTcolz");
		c1->cd(4); h2_mtTh[2]->Draw("HISTcolz");
		c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf",run));
		c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf",run)); 
		c1->Clear();
	}
	c1->Divide(2,2); 
	c1->cd(1); h2_elebeta[2]->Draw("HISTcolz"); 
	c1->cd(2); h2_Pebeta[2]->Draw("HISTcolz"); 
	c1->cd(3); h2_elemom[2]->Draw("HISTcolz");
	c1->cd(4); h2_Pemom[2]->Draw("HISTcolz");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf",run));
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/ana_monitor_%d.pdf)",run)); 
	c1->Clear();
	c1->Divide(2,2); 
	c1->cd(1); h2_Q2W[0]->Draw("HISTcolz"); 
	c1->cd(2); h2_XQ2[0]->Draw("HISTcolz");
	c1->cd(3); h2_Q2mt[0]->Draw("HISTcolz");
	c1->cd(4); h2_nuep[0]->Draw("HISTcolz");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf",run));
	c1->Clear();
	if (process.compare("vcs")==0 || process.compare("pi0")==0
		    || process.compare("vcsLT")==0 || process.compare("pi0LT")==0
		    || process.compare("pi+LT")==0 || process.compare("K+LT")==0
		    || process.compare("K+")==0 || process.compare("pi+")==0
		){	
		c1->Divide(2,2); 
		c1->cd(1); h2_Q2Th[0]->Draw("HISTcolz"); 
		c1->cd(2); h2_WTh[0]->Draw("HISTcolz");
		c1->cd(3); h2_ThCMPhi[0]->Draw("HISTcolz");
		c1->cd(4); h2_mtTh[0]->Draw("HISTcolz");
		c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf",run));
		c1->Clear();
	}
	c1->Divide(2,2); 
	c1->cd(1); h2_elebeta[1]->Draw("HISTcolz"); 
	c1->cd(2); h2_Pebeta[1]->Draw("HISTcolz"); 
	c1->cd(3); h2_elemom[0]->Draw("HISTcolz");
	c1->cd(4); h2_Pemom[0]->Draw("HISTcolz");
	c1->SaveAs(Form(VCS_REPLAY_PATH "/Ana/Results/kinematics2D_%d.pdf)",run));

	return 1;

}

//////////////////////////////////////////////////////////////////////////
bool ReadHallCData::Notify()
{
        return kTRUE;

}


//////////////////////////////////////////////////////////////////////////

ReadHallCData::~ReadHallCData()
{
   	cerr<<"call destructor"<<endl;
	if (!iChain) return;
	delete iChain->GetCurrentFile();
   // if (iChain) delete iChain->  GetCurrentFile()->GetName();
}

//////////////////////////////////////////////////////////////////////////

ReadHallCData::ReadHallCData () : iChain(0) 
{
	cout<<"init class"<< endl;
}
//////////////////////////////////////////////////////////////////////////

ReadHallCData::ReadHallCData ( vector <string> files ) : iChain(0) 
{

        iChain=new TChain("T");
        for (unsigned int i=0;i<files.size();i++){
                iChain->Add(Form("%s",files[i].c_str()));

        }
        //int test=Init();

        if (test==0) {
                cout<<"Initialization issues, check it!"<<endl;
        }


}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


Int_t ReadHallCData::GetEntry(long long entry)
{
   if (!iChain) {
        cout<<"cannot open the tree"<<endl;
        return 0;
   }
   return iChain->GetEntry(entry);
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

long long ReadHallCData::LoadTree(long long entry)
{
// Set the environment to read one entry
   if (!iChain) return -5;
   long long centry = iChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!iChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)iChain;
   if (chain->GetTreeNumber() != fCurrent) {
        fCurrent = iChain->GetTreeNumber();
        Notify();
   }
   return centry;
}

////////////////////////////////////////

int ReadHallCData::Init(string process){

	InitValues();
	InitBranch(process);
	return 1;
}


//////////////////////////////////////////
//

Int_t ReadHallCData::InitValues(){

	//Event_Branch=0;
   g_evtyp =0;
   IPM3H07A_XPOS =0;
   IPM3H07A_YPOS =0;
   IPM3H07B_XPOS =0;
   IPM3H07B_YPOS =0;
   IPM3H07C_XPOS =0;
   IPM3H07C_YPOS =0;
   IPM3H07A_XRAW =0;
   IPM3H07A_YRAW =0;
   IPM3H07B_XRAW =0;
   IPM3H07B_YRAW =0;
   IPM3H07C_XRAW =0;
   IPM3H07C_YRAW =0;
   IPM3H07A_XSOF =0;
   IPM3H07A_YSOF =0;
   IPM3H07B_XSOF =0;
   IPM3H07B_YSOF =0;
   IPM3H07C_XSOF =0;
   IPM3H07C_YSOF=0;
   fEvtHdr_fEvtTime=0;
   fEvtHdr_fEvtNum=0;
   fEvtHdr_fEvtType=0;
   fEvtHdr_fEvtLen=0;
   fEvtHdr_fHelicity=0;
   fEvtHdr_fTargetPol=0;
   fEvtHdr_fRun=0;
   IBC3H00CRCUR4=0;
   IGL1I00OD16_16=0;
   HELFREQ=0;
   HELDELAYs=0;
   hac_bcm_average=0;
   ibcm1=0;
   ibcm2=0;
   iunser=0;
   itov3out=0;
   itov4out=0;
   H_gtr_xptar=0;
    H_gtr_yptar=0;
    P_gtr_xptar=0;
    P_gtr_yptar=0;
    mmiss=0;
    CTime_epCoinTime_ROC1=0;
    CTime_epCoinTime_ROC2=0;
    CTime_epCoinTime_TRIG1=0;
    CTime_epCoinTime_TRIG4=0;
    H_cal_eprtrack=0;
    H_cal_eprtracknorm=0;
    H_cal_etot=0;
    H_cal_etotnorm=0;
    H_cal_etottracknorm=0;
    H_cal_etrack=0;
    H_cal_etracknorm=0;
    H_cer_npeSum=0;
    H_dc_ntrack=0;
    H_dc_x_fp=0;
    H_dc_xp_fp=0;
    H_dc_y_fp=0;
    H_dc_yp_fp=0;
    H_extcor_xsieve=0;
    H_extcor_ysieve=0;
    H_gtr_beta=0;
    H_gtr_dp=0;
    H_gtr_index=0;
    H_gtr_ok=0;
    H_gtr_p=0;
    H_gtr_ph=0;
    H_gtr_px=0;
    H_gtr_py=0;
    H_gtr_pz=0;
    H_gtr_th=0;
    H_gtr_x=0;
    H_gtr_y=0;
    H_hod_beta=0;
    H_hod_betachisqnotrack=0;
    H_hod_betanotrack=0;
    H_hod_fpHitsTime=0;
    H_hod_starttime=0;
    H_kin_primary_Q2=0;
    H_kin_primary_W=0;
    H_kin_primary_W2=0;
    H_kin_primary_epsilon=0;
    H_kin_primary_nu=0;
    H_kin_primary_omega=0;
    H_kin_primary_ph_q=0;
    H_kin_primary_q3m=0;
    H_kin_primary_q_x=0;
    H_kin_primary_q_y=0;
    H_kin_primary_q_z=0;
    H_kin_primary_scat_ang_deg=0;
    H_kin_primary_scat_ang_rad=0;
    H_kin_primary_th_q=0;
    H_kin_primary_x_bj=0;
    H_rb_raster_fr_xbpm_tar=0;
    H_rb_raster_fr_ybpm_tar=0;
    H_react_x=0;
    H_react_y=0;
    H_react_z=0;
    P_aero_npeSum=0;
    P_cal_eprtrack=0;
    P_cal_eprtracknorm=0;
    P_cal_etot=0;
    P_cal_etotnorm=0;
    P_cal_etottracknorm=0;
    P_cal_etracknorm=0;
    P_cal_pr_eplane=0;
    P_dc_ntrack=0;
    P_dc_x_fp=0;
    P_dc_xp_fp=0;
    P_dc_y_fp=0;
    P_dc_yp_fp=0;
    P_extcor_xsieve=0;
    P_extcor_ysieve=0;
    P_gtr_beta=0;
    P_gtr_dp=0;
    P_gtr_index=0;
    P_gtr_ok=0;
    P_gtr_p=0;
    P_gtr_ph=0;
    P_gtr_px=0;
    P_gtr_py=0;
    P_gtr_pz=0;
    P_gtr_th=0;
    P_gtr_x=0;
    P_gtr_y=0;
    P_hgcer_npeSum=0;
    P_hod_beta=0;
    P_hod_betachisqnotrack=0;
    P_hod_betanotrack=0;
    P_hod_fpHitsTime=0;
    P_hod_goodscinhit=0;
    P_hod_starttime=0;
    P_kin_secondary_Erecoil=0;
    P_kin_secondary_MandelS=0;
    P_kin_secondary_MandelT=0;
    P_kin_secondary_MandelU=0;
    P_kin_secondary_Mrecoil=0;
    P_kin_secondary_Prec_x=0;
    P_kin_secondary_Prec_y=0;
    P_kin_secondary_Prec_z=0;
    P_kin_secondary_emiss=0;
    P_kin_secondary_emiss_nuc=0;
    P_kin_secondary_ph_bq=0;
    P_kin_secondary_ph_xq=0;
    P_kin_secondary_phb_cm=0;
    P_kin_secondary_phx_cm=0;
    P_kin_secondary_pmiss=0;
    P_kin_secondary_pmiss_x=0;
    P_kin_secondary_pmiss_y=0;
    P_kin_secondary_pmiss_z=0;
    P_kin_secondary_px_cm=0;
    P_kin_secondary_t_tot_cm=0;
    P_kin_secondary_tb=0;
    P_kin_secondary_tb_cm=0;
    P_kin_secondary_th_bq=0;
    P_kin_secondary_th_xq=0;
    P_kin_secondary_thb_cm=0;
    P_kin_secondary_thx_cm=0;
    P_kin_secondary_tx=0;
    P_kin_secondary_tx_cm=0;
    P_kin_secondary_xangle=0;
    P_ngcer_npeSum=0;
    P_react_z=0;
    T_coin_h1T_tdcMultiplicity=0;
    T_coin_h1T_tdcTime=0;
    T_coin_h1T_tdcTimeRaw=0;
    T_coin_h1X_tdcMultiplicity=0;
    T_coin_h1X_tdcTime=0;
    T_coin_h1X_tdcTimeRaw=0;
    T_coin_h1Y_tdcMultiplicity=0;
    T_coin_h1Y_tdcTime=0;
    T_coin_h1Y_tdcTimeRaw=0;
    T_coin_h2T_tdcMultiplicity=0;
    T_coin_h2T_tdcTime=0;
    T_coin_h2T_tdcTimeRaw=0;
    T_coin_h2X_tdcMultiplicity=0;
    T_coin_h2X_tdcTime=0;
    T_coin_h2X_tdcTimeRaw=0;
    T_coin_h2Y_tdcMultiplicity=0;
    T_coin_h2Y_tdcTime=0;
    T_coin_h2Y_tdcTimeRaw=0;
    T_coin_hAER_adcMultiplicity=0;
    T_coin_hAER_adcPed=0;
    T_coin_hAER_adcPedRaw=0;
    T_coin_hAER_adcPulseAmp=0;
    T_coin_hAER_adcPulseAmpRaw=0;
    T_coin_hAER_adcPulseInt=0;
    T_coin_hAER_adcPulseIntRaw=0;
    T_coin_hAER_adcPulseTime=0;
    T_coin_hAER_adcPulseTimeRaw=0;
    T_coin_hASUM_adcMultiplicity=0;
    T_coin_hASUM_adcPed=0;
    T_coin_hASUM_adcPedRaw=0;
    T_coin_hASUM_adcPulseAmp=0;
    T_coin_hASUM_adcPulseAmpRaw=0;
    T_coin_hASUM_adcPulseInt=0;
    T_coin_hASUM_adcPulseIntRaw=0;
    T_coin_hASUM_adcPulseTime=0;
    T_coin_hASUM_adcPulseTimeRaw=0;
    T_coin_hASUM_tdcMultiplicity=0;
    T_coin_hASUM_tdcTime=0;
    T_coin_hASUM_tdcTimeRaw=0;
    T_coin_hBSUM_adcMultiplicity=0;
    T_coin_hBSUM_adcPed=0;
    T_coin_hBSUM_adcPedRaw=0;
    T_coin_hBSUM_adcPulseAmp=0;
    T_coin_hBSUM_adcPulseAmpRaw=0;
    T_coin_hBSUM_adcPulseInt=0;
    T_coin_hBSUM_adcPulseIntRaw=0;
    T_coin_hBSUM_adcPulseTime=0;
    T_coin_hBSUM_adcPulseTimeRaw=0;
    T_coin_hBSUM_tdcMultiplicity=0;
    T_coin_hBSUM_tdcTime=0;
    T_coin_hBSUM_tdcTimeRaw=0;
    T_coin_hCER_adcMultiplicity=0;
    T_coin_hCER_adcPed=0;
    T_coin_hCER_adcPedRaw=0;
    T_coin_hCER_adcPulseAmp=0;
    T_coin_hCER_adcPulseAmpRaw=0;
    T_coin_hCER_adcPulseInt=0;
    T_coin_hCER_adcPulseIntRaw=0;
    T_coin_hCER_adcPulseTime=0;
    T_coin_hCER_adcPulseTimeRaw=0;
    T_coin_hCER_tdcMultiplicity=0;
    T_coin_hCER_tdcTime=0;
    T_coin_hCER_tdcTimeRaw=0;
    T_coin_hCSUM_adcMultiplicity=0;
    T_coin_hCSUM_adcPed=0;
    T_coin_hCSUM_adcPedRaw=0;
    T_coin_hCSUM_adcPulseAmp=0;
    T_coin_hCSUM_adcPulseAmpRaw=0;
    T_coin_hCSUM_adcPulseInt=0;
    T_coin_hCSUM_adcPulseIntRaw=0;
    T_coin_hCSUM_adcPulseTime=0;
    T_coin_hCSUM_adcPulseTimeRaw=0;
    T_coin_hCSUM_tdcMultiplicity=0;
    T_coin_hCSUM_tdcTime=0;
    T_coin_hCSUM_tdcTimeRaw=0;
    T_coin_hDCREF1_tdcMultiplicity=0;
    T_coin_hDCREF1_tdcTime=0;
    T_coin_hDCREF1_tdcTimeRaw=0;
    T_coin_hDCREF2_tdcMultiplicity=0;
    T_coin_hDCREF2_tdcTime=0;
    T_coin_hDCREF2_tdcTimeRaw=0;
    T_coin_hDCREF3_tdcMultiplicity=0;
    T_coin_hDCREF3_tdcTime=0;
    T_coin_hDCREF3_tdcTimeRaw=0;
    T_coin_hDCREF4_tdcMultiplicity=0;
    T_coin_hDCREF4_tdcTime=0;
    T_coin_hDCREF4_tdcTimeRaw=0;
    T_coin_hDCREF5_tdcMultiplicity=0;
    T_coin_hDCREF5_tdcTime=0;
    T_coin_hDCREF5_tdcTimeRaw=0;
    T_coin_hDSUM_adcMultiplicity=0;
    T_coin_hDSUM_adcPed=0;
    T_coin_hDSUM_adcPedRaw=0;
    T_coin_hDSUM_adcPulseAmp=0;
    T_coin_hDSUM_adcPulseAmpRaw=0;
    T_coin_hDSUM_adcPulseInt=0;
    T_coin_hDSUM_adcPulseIntRaw=0;
    T_coin_hDSUM_adcPulseTime=0;
    T_coin_hDSUM_adcPulseTimeRaw=0;
    T_coin_hDSUM_tdcMultiplicity=0;
    T_coin_hDSUM_tdcTime=0;
    T_coin_hDSUM_tdcTimeRaw=0;
    T_coin_hEDTM_tdcMultiplicity=0;
    T_coin_hEDTM_tdcTime=0;
    T_coin_hEDTM_tdcTimeRaw=0;
    T_coin_hEL_CLEAN_ROC1_tdcMultiplicity=0;
    T_coin_hEL_CLEAN_ROC1_tdcTime=0;
    T_coin_hEL_CLEAN_ROC1_tdcTimeRaw=0;
    T_coin_hEL_CLEAN_ROC2_tdcMultiplicity=0;
    T_coin_hEL_CLEAN_ROC2_tdcTime=0;
    T_coin_hEL_CLEAN_ROC2_tdcTimeRaw=0;
    T_coin_hEL_HI_ROC1_tdcMultiplicity=0;
    T_coin_hEL_CLEAN_ROC1_tdcTime=0;
    T_coin_hEL_CLEAN_ROC1_tdcTimeRaw=0;
    T_coin_hEL_CLEAN_ROC2_tdcMultiplicity=0;
    T_coin_hEL_CLEAN_ROC2_tdcTime=0;
    T_coin_hEL_CLEAN_ROC2_tdcTimeRaw=0;
    T_coin_hEL_HI_ROC1_tdcMultiplicity=0;
    T_coin_hEL_HI_ROC1_tdcTime=0;
    T_coin_hEL_HI_ROC1_tdcTimeRaw=0;
    T_coin_hEL_HI_ROC2_tdcMultiplicity=0;
    T_coin_hEL_HI_ROC2_tdcTime=0;
    T_coin_hEL_HI_ROC2_tdcTimeRaw=0;
    T_coin_hEL_LO_LO_ROC1_tdcMultiplicity=0;
    T_coin_hEL_LO_LO_ROC1_tdcTime=0;
    T_coin_hEL_LO_LO_ROC1_tdcTimeRaw=0;
    T_coin_hEL_LO_LO_ROC2_tdcMultiplicity=0;
    T_coin_hEL_LO_LO_ROC2_tdcTime=0;
    T_coin_hEL_LO_LO_ROC2_tdcTimeRaw=0;
    T_coin_hEL_LO_ROC1_tdcMultiplicity=0;
    T_coin_hEL_LO_ROC1_tdcTime=0;
    T_coin_hEL_LO_ROC1_tdcTimeRaw=0;
    T_coin_hEL_LO_ROC2_tdcMultiplicity=0;
    T_coin_hEL_LO_ROC2_tdcTime=0;
    T_coin_hEL_LO_ROC2_tdcTimeRaw=0;
    T_coin_hEL_REAL_ROC1_tdcMultiplicity=0;
    T_coin_hEL_REAL_ROC1_tdcTime=0;
    T_coin_hEL_REAL_ROC1_tdcTimeRaw=0;
    T_coin_hEL_REAL_ROC2_tdcMultiplicity=0;
    T_coin_hEL_REAL_ROC2_tdcTime=0;
    T_coin_hEL_REAL_ROC2_tdcTimeRaw=0;
    T_coin_hFADC_TREF_ROC1_adcMultiplicity=0;
    T_coin_hFADC_TREF_ROC1_adcPed=0;
    T_coin_hFADC_TREF_ROC1_adcPedRaw=0;
    T_coin_hFADC_TREF_ROC1_adcPulseAmp=0;
    T_coin_hFADC_TREF_ROC1_adcPulseAmpRaw=0;
    T_coin_hFADC_TREF_ROC1_adcPulseInt=0;
    T_coin_hFADC_TREF_ROC1_adcPulseIntRaw=0;
    T_coin_hFADC_TREF_ROC1_adcPulseTime=0;
    T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw=0;
    T_coin_hHODO_RF_tdcMultiplicity=0;
    T_coin_hHODO_RF_tdcTime=0;
    T_coin_hHODO_RF_tdcTimeRaw=0;
    T_coin_hPRE100_ROC1_tdcMultiplicity=0;
    T_coin_hPRE100_ROC1_tdcTime=0;
    T_coin_hPRE100_ROC1_tdcTimeRaw=0;
    T_coin_hPRE100_ROC2_tdcMultiplicity=0;
    T_coin_hPRE100_ROC2_tdcTime=0;
    T_coin_hPRE100_ROC2_tdcTimeRaw=0;
    T_coin_hPRE150_ROC1_tdcMultiplicity=0;
    T_coin_hPRE150_ROC1_tdcTime=0;
    T_coin_hPRE150_ROC1_tdcTimeRaw=0;
    T_coin_hPRE150_ROC2_tdcMultiplicity=0;
    T_coin_hPRE150_ROC2_tdcTime=0;
    T_coin_hPRE150_ROC2_tdcTimeRaw=0;
    T_coin_hPRE200_ROC1_tdcMultiplicity=0;
    T_coin_hPRE200_ROC1_tdcTime=0;
    T_coin_hPRE200_ROC1_tdcTimeRaw=0;
    T_coin_hPRE200_ROC2_tdcMultiplicity=0;
    T_coin_hPRE200_ROC2_tdcTime=0;
    T_coin_hPRE200_ROC2_tdcTimeRaw=0;
    T_coin_hPRE40_ROC1_tdcMultiplicity=0;
    T_coin_hPRE40_ROC1_tdcTime=0;
    T_coin_hPRE40_ROC1_tdcTimeRaw=0;
    T_coin_hPRE40_ROC2_tdcMultiplicity=0;
    T_coin_hPRE40_ROC2_tdcTime=0;
    T_coin_hPRE40_ROC2_tdcTimeRaw=0;
    T_coin_hPRHI_tdcMultiplicity=0;
    T_coin_hPRHI_tdcTime=0;
    T_coin_hPRHI_tdcTimeRaw=0;
    T_coin_hPRLO_tdcMultiplicity=0;
    T_coin_hPRLO_tdcTime=0;
    T_coin_hPRLO_tdcTimeRaw=0;
    T_coin_hPSHWR_adcMultiplicity=0;
    T_coin_hPSHWR_adcPed=0;
    T_coin_hPSHWR_adcPedRaw=0;
    T_coin_hPSHWR_adcPulseAmp=0;
    T_coin_hPSHWR_adcPulseAmpRaw=0;
    T_coin_hPSHWR_adcPulseInt=0;
    T_coin_hPSHWR_adcPulseIntRaw=0;
    T_coin_hPSHWR_adcPulseTime=0;
    T_coin_hPSHWR_adcPulseTimeRaw=0;
    T_coin_hRF_tdcMultiplicity=0;
    T_coin_hRF_tdcTime=0;
    T_coin_hRF_tdcTimeRaw=0;
    T_coin_hSHWR_adcMultiplicity=0;
    T_coin_hSHWR_adcPed=0;
    T_coin_hSHWR_adcPedRaw=0;
    T_coin_hSHWR_adcPulseAmp=0;
    T_coin_hSHWR_adcPulseAmpRaw=0;
    T_coin_hSHWR_adcPulseInt=0;
    T_coin_hSHWR_adcPulseIntRaw=0;
    T_coin_hSHWR_adcPulseTime=0;
    T_coin_hSHWR_adcPulseTimeRaw=0;
    T_coin_hSHWR_tdcMultiplicity=0;
    T_coin_hSHWR_tdcTime=0;
    T_coin_hSHWR_tdcTimeRaw=0;
    T_coin_hSTOF_ROC1_tdcMultiplicity=0;
    T_coin_hSTOF_ROC1_tdcTime=0;
    T_coin_hSTOF_ROC1_tdcTimeRaw=0;
    T_coin_hSTOF_ROC2_tdcMultiplicity=0;
    T_coin_hSTOF_ROC2_tdcTime=0;
    T_coin_hSTOF_ROC2_tdcTimeRaw=0;
    T_coin_hT1_tdcMultiplicity=0;
    T_coin_hT1_tdcTime=0;
    T_coin_hT1_tdcTimeRaw=0;
    T_coin_hT2_tdcMultiplicity=0;
    T_coin_hT2_tdcTime=0;
    T_coin_hT2_tdcTimeRaw=0;
    T_coin_hTRIG1_ROC1_tdcMultiplicity=0;
    T_coin_hTRIG1_ROC1_tdcTime=0;
    T_coin_hTRIG1_ROC1_tdcTimeRaw=0;
    T_coin_hTRIG1_ROC2_tdcMultiplicity=0;
    T_coin_hTRIG1_ROC2_tdcTime=0;
    T_coin_hTRIG1_ROC2_tdcTimeRaw=0;
    T_coin_hTRIG2_ROC1_tdcMultiplicity=0;
    T_coin_hTRIG2_ROC1_tdcTime=0;
    T_coin_hTRIG2_ROC1_tdcTimeRaw=0;
    T_coin_hTRIG2_ROC2_tdcMultiplicity=0;
    T_coin_hTRIG2_ROC2_tdcTime=0;
    T_coin_hTRIG2_ROC2_tdcTimeRaw=0;
    T_coin_hTRIG3_ROC1_tdcMultiplicity=0;
    T_coin_hTRIG3_ROC1_tdcTime=0;
    T_coin_hTRIG3_ROC1_tdcTimeRaw=0;
    T_coin_hTRIG3_ROC2_tdcMultiplicity=0;
    T_coin_hTRIG3_ROC2_tdcTime=0;
    T_coin_hTRIG3_ROC2_tdcTimeRaw=0;
    T_coin_hTRIG4_ROC1_tdcMultiplicity=0;
    T_coin_hTRIG4_ROC1_tdcTime=0;
    T_coin_hTRIG4_ROC1_tdcTimeRaw=0;
    T_coin_hTRIG4_ROC2_tdcMultiplicity=0;
    T_coin_hTRIG4_ROC2_tdcTime=0;
    T_coin_hTRIG4_ROC2_tdcTimeRaw=0;
    T_coin_hTRIG5_ROC1_tdcMultiplicity=0;
    T_coin_hTRIG5_ROC1_tdcTime=0;
    T_coin_hTRIG5_ROC1_tdcTimeRaw=0;
    T_coin_hTRIG5_ROC2_tdcMultiplicity=0;
    T_coin_hTRIG5_ROC2_tdcTime=0;
    T_coin_hTRIG5_ROC2_tdcTimeRaw=0;
    T_coin_hTRIG6_ROC1_tdcMultiplicity=0;
    T_coin_hTRIG6_ROC1_tdcTime=0;
    T_coin_hTRIG6_ROC1_tdcTimeRaw=0;
    T_coin_hTRIG6_ROC2_tdcMultiplicity=0;
    T_coin_hTRIG6_ROC2_tdcTime=0;
    T_coin_hTRIG6_ROC2_tdcTimeRaw=0;
    T_coin_p1T_tdcMultiplicity=0;
    T_coin_p1T_tdcTime=0;
    T_coin_p1T_tdcTimeRaw=0;
    T_coin_p1X_tdcMultiplicity=0;
    T_coin_p1X_tdcTime=0;
    T_coin_p1X_tdcTimeRaw=0;
    T_coin_p1Y_tdcMultiplicity=0;
    T_coin_p1Y_tdcTime=0;
    T_coin_p1Y_tdcTimeRaw=0;
    T_coin_p2T_tdcMultiplicity=0;
    T_coin_p2T_tdcTime=0;
    T_coin_p2T_tdcTimeRaw=0;
    T_coin_p2X_tdcMultiplicity=0;
    T_coin_p2X_tdcTime=0;
    T_coin_p2X_tdcTimeRaw=0;
    T_coin_p2Y_tdcMultiplicity=0;
    T_coin_p2Y_tdcTime=0;
    T_coin_p2Y_tdcTimeRaw=0;
    T_coin_pAER_adcMultiplicity=0;
    T_coin_pAER_adcPed=0;
    T_coin_pAER_adcPedRaw=0;
    T_coin_pAER_adcPulseAmp=0;
    T_coin_pAER_adcPulseAmpRaw=0;
    T_coin_pAER_adcPulseInt=0;
    T_coin_pAER_adcPulseIntRaw=0;
    T_coin_pAER_adcPulseTime=0;
    T_coin_pAER_adcPulseTimeRaw=0;
    T_coin_pAER_tdcMultiplicity=0;
    T_coin_pAER_tdcTime=0;
    T_coin_pAER_tdcTimeRaw=0;
    T_coin_pDCREF10_tdcMultiplicity=0;
    T_coin_pDCREF10_tdcTime=0;
    T_coin_pDCREF10_tdcTimeRaw=0;
    T_coin_pDCREF1_tdcMultiplicity=0;
    T_coin_pDCREF1_tdcTime=0;
    T_coin_pDCREF1_tdcTimeRaw=0;
    T_coin_pDCREF2_tdcMultiplicity=0;
    T_coin_pDCREF2_tdcTime=0;
    T_coin_pDCREF2_tdcTimeRaw=0;
    T_coin_pDCREF3_tdcMultiplicity=0;
    T_coin_pDCREF3_tdcTime=0;
    T_coin_pDCREF3_tdcTimeRaw=0;
    T_coin_pDCREF4_tdcMultiplicity=0;
    T_coin_pDCREF4_tdcTime=0;
    T_coin_pDCREF4_tdcTimeRaw=0;
    T_coin_pDCREF5_tdcMultiplicity=0;
    T_coin_pDCREF5_tdcTime=0;
    T_coin_pDCREF5_tdcTimeRaw=0;
    T_coin_pDCREF6_tdcMultiplicity=0;
    T_coin_pDCREF6_tdcTime=0;
    T_coin_pDCREF6_tdcTimeRaw=0;
    T_coin_pDCREF7_tdcMultiplicity=0;
    T_coin_pDCREF7_tdcTime=0;
    T_coin_pDCREF7_tdcTimeRaw=0;
    T_coin_pDCREF8_tdcMultiplicity=0;
    T_coin_pDCREF8_tdcTime=0;
    T_coin_pDCREF8_tdcTimeRaw=0;
    T_coin_pDCREF9_tdcMultiplicity=0;
    T_coin_pDCREF9_tdcTime=0;
    T_coin_pDCREF9_tdcTimeRaw=0;
    T_coin_pEDTM_tdcMultiplicity=0;
    T_coin_pEDTM_tdcTime=0;
    T_coin_pEDTM_tdcTimeRaw=0;
    T_coin_pEL_CLEAN_ROC1_tdcMultiplicity=0;
    T_coin_pEL_CLEAN_ROC1_tdcTime=0;
    T_coin_pEL_CLEAN_ROC1_tdcTimeRaw=0;
    T_coin_pEL_CLEAN_ROC2_tdcMultiplicity=0;
    T_coin_pEL_CLEAN_ROC2_tdcTime=0;
    T_coin_pEL_CLEAN_ROC2_tdcTimeRaw=0;
    T_coin_pEL_HI_ROC1_tdcMultiplicity=0;
    T_coin_pEL_HI_ROC1_tdcTime=0;
    T_coin_pEL_HI_ROC1_tdcTimeRaw=0;
    T_coin_pEL_HI_ROC2_tdcMultiplicity=0;
    T_coin_pEL_HI_ROC2_tdcTime=0;
    T_coin_pEL_HI_ROC2_tdcTimeRaw=0;
    T_coin_pEL_LO_LO_ROC1_tdcMultiplicity=0;
    T_coin_pEL_LO_LO_ROC1_tdcTime=0;
    T_coin_pEL_LO_LO_ROC1_tdcTimeRaw=0;
    T_coin_pEL_LO_LO_ROC2_tdcMultiplicity=0;
    T_coin_pEL_LO_LO_ROC2_tdcTime=0;
    T_coin_pEL_LO_LO_ROC2_tdcTimeRaw=0;
    T_coin_pEL_LO_ROC1_tdcMultiplicity=0;
    T_coin_pEL_LO_ROC1_tdcTime=0;
    T_coin_pEL_LO_ROC1_tdcTimeRaw=0;
    T_coin_pEL_LO_ROC2_tdcMultiplicity=0;
    T_coin_pEL_LO_ROC2_tdcTime=0;
    T_coin_pEL_LO_ROC2_tdcTimeRaw=0;
    T_coin_pEL_REAL_ROC1_tdcMultiplicity=0;
    T_coin_pEL_REAL_ROC1_tdcTime=0;
    T_coin_pEL_REAL_ROC1_tdcTimeRaw=0;
    T_coin_pEL_REAL_ROC2_tdcMultiplicity=0;
    T_coin_pEL_REAL_ROC2_tdcTime=0;
    T_coin_pEL_REAL_ROC2_tdcTimeRaw=0;
    T_coin_pFADC_TREF_ROC2_adcMultiplicity=0;
    T_coin_pFADC_TREF_ROC2_adcPed=0;
    T_coin_pFADC_TREF_ROC2_adcPedRaw=0;
    T_coin_pFADC_TREF_ROC2_adcPulseAmp=0;
    T_coin_pFADC_TREF_ROC2_adcPulseAmpRaw=0;
    T_coin_pFADC_TREF_ROC2_adcPulseInt=0;
    T_coin_pFADC_TREF_ROC2_adcPulseIntRaw=0;
    T_coin_pFADC_TREF_ROC2_adcPulseTime=0;
    T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw=0;
    T_coin_pHEL_MPS_adcMultiplicity=0;
    T_coin_pHEL_MPS_adcPed=0;
    T_coin_pHEL_MPS_adcPedRaw=0;
    T_coin_pHEL_MPS_adcPulseAmp=0;
    T_coin_pHEL_MPS_adcPulseAmpRaw=0;
    T_coin_pHEL_MPS_adcPulseInt=0;
    T_coin_pHEL_MPS_adcPulseIntRaw=0;
    T_coin_pHEL_MPS_adcPulseTime=0;
    T_coin_pHEL_MPS_adcPulseTimeRaw=0;
    T_coin_pHEL_NEG_adcMultiplicity=0;
    T_coin_pHEL_NEG_adcPed=0;
    T_coin_pHEL_NEG_adcPedRaw=0;
    T_coin_pHEL_NEG_adcPulseAmp=0;
    T_coin_pHEL_NEG_adcPulseAmpRaw=0;
    T_coin_pHEL_NEG_adcPulseInt=0;
    T_coin_pHEL_NEG_adcPulseIntRaw=0;
    T_coin_pHEL_NEG_adcPulseTime=0;
    T_coin_pHEL_NEG_adcPulseTimeRaw=0;
    T_coin_pHEL_POS_adcMultiplicity=0;
    T_coin_pHEL_POS_adcPed=0;
    T_coin_pHEL_POS_adcPedRaw=0;
    T_coin_pHEL_POS_adcPulseAmp=0;
    T_coin_pHEL_POS_adcPulseAmpRaw=0;
    T_coin_pHEL_POS_adcPulseInt=0;
    T_coin_pHEL_POS_adcPulseIntRaw=0;
    T_coin_pHEL_POS_adcPulseTime=0;
    T_coin_pHEL_POS_adcPulseTimeRaw=0;
    T_coin_pHGCER_MOD_adcMultiplicity=0;
    T_coin_pHGCER_MOD_adcPed=0;
    T_coin_pHGCER_MOD_adcPedRaw=0;
    T_coin_pHGCER_MOD_adcPulseAmp=0;
    T_coin_pHGCER_MOD_adcPulseAmpRaw=0;
    T_coin_pHGCER_MOD_adcPulseInt=0;
    T_coin_pHGCER_MOD_adcPulseIntRaw=0;
    T_coin_pHGCER_MOD_adcPulseTime=0;
    T_coin_pHGCER_MOD_adcPulseTimeRaw=0;
    T_coin_pHGCER_adcMultiplicity=0;
    T_coin_pHGCER_adcPed=0;
    T_coin_pHGCER_adcPedRaw=0;
    T_coin_pHGCER_adcPulseAmp=0;
    T_coin_pHGCER_adcPulseAmpRaw=0;
    T_coin_pHGCER_adcPulseInt=0;
    T_coin_pHGCER_adcPulseIntRaw=0;
    T_coin_pHGCER_adcPulseTime=0;
    T_coin_pHGCER_adcPulseTimeRaw=0;
    T_coin_pHGCER_tdcMultiplicity=0;
    T_coin_pHGCER_tdcTime=0;
    T_coin_pHGCER_tdcTimeRaw=0;
    T_coin_pHODO_RF_tdcMultiplicity=0;
    T_coin_pHODO_RF_tdcTime=0;
    T_coin_pHODO_RF_tdcTimeRaw=0;
    T_coin_pNGCER_MOD_adcMultiplicity=0;
    T_coin_pNGCER_MOD_adcPed=0;
    T_coin_pNGCER_MOD_adcPedRaw=0;
    T_coin_pNGCER_MOD_adcPulseAmp=0;
    T_coin_pNGCER_MOD_adcPulseAmpRaw=0;
    T_coin_pNGCER_MOD_adcPulseInt=0;
    T_coin_pNGCER_MOD_adcPulseIntRaw=0;
    T_coin_pNGCER_MOD_adcPulseTime=0;
    T_coin_pNGCER_MOD_adcPulseTimeRaw=0;
    T_coin_pNGCER_adcMultiplicity=0;
    T_coin_pNGCER_adcPed=0;
    T_coin_pNGCER_adcPedRaw=0;
    T_coin_pNGCER_adcPulseAmp=0;
    T_coin_pNGCER_adcPulseAmpRaw=0;
    T_coin_pNGCER_adcPulseInt=0;
    T_coin_pNGCER_adcPulseIntRaw=0;
    T_coin_pNGCER_adcPulseTime=0;
    T_coin_pNGCER_adcPulseTimeRaw=0;
    T_coin_pNGCER_tdcMultiplicity=0;
    T_coin_pNGCER_tdcTime=0;
    T_coin_pNGCER_tdcTimeRaw=0;
    T_coin_pPRE100_ROC1_tdcMultiplicity=0;
    T_coin_pPRE100_ROC1_tdcTime=0;
    T_coin_pPRE100_ROC1_tdcTimeRaw=0;
    T_coin_pPRE100_ROC2_tdcMultiplicity=0;
    T_coin_pPRE100_ROC2_tdcTime=0;
    T_coin_pPRE100_ROC2_tdcTimeRaw=0;
    T_coin_pPRE150_ROC1_tdcMultiplicity=0;
    T_coin_pPRE150_ROC1_tdcTime=0;
    T_coin_pPRE150_ROC1_tdcTimeRaw=0;
    T_coin_pPRE150_ROC2_tdcMultiplicity=0;
    T_coin_pPRE150_ROC2_tdcTime=0;
    T_coin_pPRE150_ROC2_tdcTimeRaw=0;
    T_coin_pPRE200_ROC1_tdcMultiplicity=0;
    T_coin_pPRE200_ROC1_tdcTime=0;
    T_coin_pPRE200_ROC1_tdcTimeRaw=0;
    T_coin_pPRE200_ROC2_tdcMultiplicity=0;
    T_coin_pPRE200_ROC2_tdcTime=0;
    T_coin_pPRE200_ROC2_tdcTimeRaw=0;
    T_coin_pPRE40_ROC1_tdcMultiplicity=0;
    T_coin_pPRE40_ROC1_tdcTime=0;
    T_coin_pPRE40_ROC1_tdcTimeRaw=0;
    T_coin_pPRE40_ROC2_tdcMultiplicity=0;
    T_coin_pPRE40_ROC2_tdcTime=0;
    T_coin_pPRE40_ROC2_tdcTimeRaw=0;
    T_coin_pPRHI_tdcMultiplicity=0;
    T_coin_pPRHI_tdcTime=0;
    T_coin_pPRHI_tdcTimeRaw=0;
    T_coin_pPRLO_tdcMultiplicity=0;
    T_coin_pPRLO_tdcTime=0;
    T_coin_pPRLO_tdcTimeRaw=0;
    T_coin_pPSHWR_adcMultiplicity=0;
    T_coin_pPSHWR_adcPed=0;
    T_coin_pPSHWR_adcPedRaw=0;
    T_coin_pPSHWR_adcPulseAmp=0;
    T_coin_pPSHWR_adcPulseAmpRaw=0;
    T_coin_pPSHWR_adcPulseInt=0;
    T_coin_pPSHWR_adcPulseIntRaw=0;
    T_coin_pPSHWR_adcPulseTime=0;
    T_coin_pPSHWR_adcPulseTimeRaw=0;
    T_coin_pRF_tdcMultiplicity=0;
    T_coin_pRF_tdcTime=0;
    T_coin_pRF_tdcTimeRaw=0;
    T_coin_pSTOF_ROC1_tdcMultiplicity=0;
    T_coin_pSTOF_ROC1_tdcTime=0;
    T_coin_pSTOF_ROC1_tdcTimeRaw=0;
    T_coin_pSTOF_ROC2_tdcMultiplicity=0;
    T_coin_pSTOF_ROC2_tdcTime=0;
    T_coin_pSTOF_ROC2_tdcTimeRaw=0;
    T_coin_pT1_tdcMultiplicity=0;
    T_coin_pT1_tdcTime=0;
    T_coin_pT1_tdcTimeRaw=0;
    T_coin_pT2_tdcMultiplicity=0;
    T_coin_pT2_tdcTime=0;
    T_coin_pT2_tdcTimeRaw=0;
    T_coin_pT3_tdcMultiplicity=0;
    T_coin_pT3_tdcTime=0;
    T_coin_pT3_tdcTimeRaw=0;
    T_coin_pTRIG1_ROC1_tdcMultiplicity=0;
    T_coin_pTRIG1_ROC1_tdcTime=0;
    T_coin_pTRIG1_ROC1_tdcTimeRaw=0;
    T_coin_pTRIG1_ROC2_tdcMultiplicity=0;
    T_coin_pTRIG1_ROC2_tdcTime=0;
    T_coin_pTRIG1_ROC2_tdcTimeRaw=0;
    T_coin_pTRIG2_ROC1_tdcMultiplicity=0;
    T_coin_pTRIG2_ROC1_tdcTime=0;
    T_coin_pTRIG2_ROC1_tdcTimeRaw=0;
    T_coin_pTRIG2_ROC2_tdcMultiplicity=0;
    T_coin_pTRIG2_ROC2_tdcTime=0;
    T_coin_pTRIG2_ROC2_tdcTimeRaw=0;
    T_coin_pTRIG3_ROC1_tdcMultiplicity=0;
    T_coin_pTRIG3_ROC1_tdcTime=0;
    T_coin_pTRIG3_ROC1_tdcTimeRaw=0;
    T_coin_pTRIG3_ROC2_tdcMultiplicity=0;
    T_coin_pTRIG3_ROC2_tdcTime=0;
    T_coin_pTRIG3_ROC2_tdcTimeRaw=0;
    T_coin_pTRIG4_ROC1_tdcMultiplicity=0;
    T_coin_pTRIG4_ROC1_tdcTime=0;
    T_coin_pTRIG4_ROC1_tdcTimeRaw=0;
    T_coin_pTRIG4_ROC2_tdcMultiplicity=0;
    T_coin_pTRIG4_ROC2_tdcTime=0;
    T_coin_pTRIG4_ROC2_tdcTimeRaw=0;
    T_coin_pTRIG5_ROC1_tdcMultiplicity=0;
    T_coin_pTRIG5_ROC1_tdcTime=0;
    T_coin_pTRIG5_ROC1_tdcTimeRaw=0;
    T_coin_pTRIG5_ROC2_tdcMultiplicity=0;
    T_coin_pTRIG5_ROC2_tdcTime=0;
    T_coin_pTRIG5_ROC2_tdcTimeRaw=0;
    T_coin_pTRIG6_ROC1_tdcMultiplicity=0;
    T_coin_pTRIG6_ROC1_tdcTime=0;
    T_coin_pTRIG6_ROC1_tdcTimeRaw=0;
    T_coin_pTRIG6_ROC2_tdcMultiplicity=0;
    T_coin_pTRIG6_ROC2_tdcTime=0;
    T_coin_pTRIG6_ROC2_tdcTimeRaw=0;

    H_kin_secondary_Erecoil=0;
    H_kin_secondary_MandelS=0;
    H_kin_secondary_MandelT=0;
    H_kin_secondary_MandelU=0;
    H_kin_secondary_Mrecoil=0;
    H_kin_secondary_Prec_x=0;
    H_kin_secondary_Prec_y=0;
    H_kin_secondary_Prec_z=0;
    H_kin_secondary_emiss=0;
    H_kin_secondary_emiss_nuc=0;
    H_kin_secondary_ph_bq=0;
    H_kin_secondary_ph_xq=0;
    H_kin_secondary_phb_cm=0;
    H_kin_secondary_phx_cm=0;
    H_kin_secondary_pmiss=0;
    H_kin_secondary_pmiss_x=0;
    H_kin_secondary_pmiss_y=0;
    H_kin_secondary_pmiss_z=0;
    H_kin_secondary_px_cm=0;
    H_kin_secondary_t_tot_cm=0;
    H_kin_secondary_tb=0;
    H_kin_secondary_tb_cm=0;
    H_kin_secondary_th_bq=0;
    H_kin_secondary_th_xq=0;
    H_kin_secondary_thb_cm=0;
    H_kin_secondary_thx_cm=0;
    H_kin_secondary_tx=0;
    H_kin_secondary_tx_cm=0;
    H_kin_secondary_xangle=0;

    P_kin_primary_Q2=0;
    P_kin_primary_W=0;
    P_kin_primary_W2=0;
    P_kin_primary_epsilon=0;
    P_kin_primary_nu=0;
    P_kin_primary_omega=0;  
    P_kin_primary_ph_q=0;
    P_kin_primary_q3m=0;
    P_kin_primary_q_x=0;
    P_kin_primary_q_y=0;
    P_kin_primary_q_z=0;
    P_kin_primary_scat_ang_deg=0;
    P_kin_primary_scat_ang_rad=0;
    P_kin_primary_th_q=0;
    P_kin_primary_x_bj=0;
    P_react_ok=0;
    P_react_x=0;
    P_react_y=0;

    return 1;

}

////////////////////////////////////////////////////////////////////:

int ReadHallCData::InitBranch(string process){
	if (!iChain) return 0;

   iChain->SetBranchAddress("H.gtr.xptar", &H_gtr_xptar, &b_H_gtr_xptar);
   iChain->SetBranchAddress("H.gtr.yptar", &H_gtr_yptar, &b_H_gtr_yptar);
   iChain->SetBranchAddress("P.gtr.xptar", &P_gtr_xptar, &b_P_gtr_xptar);
   iChain->SetBranchAddress("P.gtr.yptar", &P_gtr_yptar, &b_P_gtr_yptar);
   iChain->SetBranchAddress("mmiss", &mmiss, &b_mmiss);
   iChain->SetBranchAddress("CTime.epCoinTime_ROC1", &CTime_epCoinTime_ROC1, &b_CTime_epCoinTime_ROC1);
   iChain->SetBranchAddress("CTime.epCoinTime_ROC2", &CTime_epCoinTime_ROC2, &b_CTime_epCoinTime_ROC2);
   iChain->SetBranchAddress("CTime.epCoinTime_TRIG1", &CTime_epCoinTime_TRIG1, &b_CTime_epCoinTime_TRIG1);
   iChain->SetBranchAddress("CTime.epCoinTime_TRIG4", &CTime_epCoinTime_TRIG4, &b_CTime_epCoinTime_TRIG4);
   iChain->SetBranchAddress("H.cal.eprtrack", &H_cal_eprtrack, &b_H_cal_eprtrack);
   iChain->SetBranchAddress("H.cal.eprtracknorm", &H_cal_eprtracknorm, &b_H_cal_eprtracknorm);
   iChain->SetBranchAddress("H.cal.etot", &H_cal_etot, &b_H_cal_etot);
   iChain->SetBranchAddress("H.cal.etotnorm", &H_cal_etotnorm, &b_H_cal_etotnorm);
   iChain->SetBranchAddress("H.cal.etottracknorm", &H_cal_etottracknorm, &b_H_cal_etottracknorm);
   iChain->SetBranchAddress("H.cal.etrack", &H_cal_etrack, &b_H_cal_etrack);
   iChain->SetBranchAddress("H.cal.etracknorm", &H_cal_etracknorm, &b_H_cal_etracknorm);
   iChain->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum, &b_H_cer_npeSum);
   iChain->SetBranchAddress("H.dc.ntrack", &H_dc_ntrack, &b_H_dc_ntrack);
   iChain->SetBranchAddress("H.dc.x_fp", &H_dc_x_fp, &b_H_dc_x_fp);
   iChain->SetBranchAddress("H.dc.xp_fp", &H_dc_xp_fp, &b_H_dc_xp_fp);
   iChain->SetBranchAddress("H.dc.y_fp", &H_dc_y_fp, &b_H_dc_y_fp);
   iChain->SetBranchAddress("H.dc.yp_fp", &H_dc_yp_fp, &b_H_dc_yp_fp);
   iChain->SetBranchAddress("H.extcor.xsieve", &H_extcor_xsieve, &b_H_extcor_xsieve);
   iChain->SetBranchAddress("H.extcor.ysieve", &H_extcor_ysieve, &b_H_extcor_ysieve);
   iChain->SetBranchAddress("H.gtr.beta", &H_gtr_beta, &b_H_gtr_beta);
   iChain->SetBranchAddress("H.gtr.dp", &H_gtr_dp, &b_H_gtr_dp);
   iChain->SetBranchAddress("H.gtr.index", &H_gtr_index, &b_H_gtr_index);
   iChain->SetBranchAddress("H.gtr.ok", &H_gtr_ok, &b_H_gtr_ok);
   iChain->SetBranchAddress("H.gtr.p", &H_gtr_p, &b_H_gtr_p);
   iChain->SetBranchAddress("H.gtr.ph", &H_gtr_ph, &b_H_gtr_ph);
   iChain->SetBranchAddress("H.gtr.px", &H_gtr_px, &b_H_gtr_px);
   iChain->SetBranchAddress("H.gtr.py", &H_gtr_py, &b_H_gtr_py);
   iChain->SetBranchAddress("H.gtr.pz", &H_gtr_pz, &b_H_gtr_pz);
   iChain->SetBranchAddress("H.gtr.th", &H_gtr_th, &b_H_gtr_th);
   iChain->SetBranchAddress("H.gtr.x", &H_gtr_x, &b_H_gtr_x);
   iChain->SetBranchAddress("H.gtr.y", &H_gtr_y, &b_H_gtr_y);
   iChain->SetBranchAddress("H.hod.beta", &H_hod_beta, &b_H_hod_beta);
   iChain->SetBranchAddress("H.hod.betachisqnotrack", &H_hod_betachisqnotrack, &b_H_hod_betachisqnotrack);
   iChain->SetBranchAddress("H.hod.betanotrack", &H_hod_betanotrack, &b_H_hod_betanotrack);
   iChain->SetBranchAddress("H.hod.fpHitsTime", &H_hod_fpHitsTime, &b_H_hod_fpHitsTime);
   iChain->SetBranchAddress("H.hod.starttime", &H_hod_starttime, &b_H_hod_starttime);
   iChain->SetBranchAddress("P.aero.npeSum", &P_aero_npeSum, &b_P_aero_npeSum);
   iChain->SetBranchAddress("P.cal.eprtrack", &P_cal_eprtrack, &b_P_cal_eprtrack);
   iChain->SetBranchAddress("P.cal.eprtracknorm", &P_cal_eprtracknorm, &b_P_cal_eprtracknorm);
   iChain->SetBranchAddress("P.cal.etot", &P_cal_etot, &b_P_cal_etot);
   iChain->SetBranchAddress("P.cal.etotnorm", &P_cal_etotnorm, &b_P_cal_etotnorm);
   iChain->SetBranchAddress("P.cal.etottracknorm", &P_cal_etottracknorm, &b_P_cal_etottracknorm);
   iChain->SetBranchAddress("P.cal.etracknorm", &P_cal_etracknorm, &b_P_cal_etracknorm);
   iChain->SetBranchAddress("P.cal.pr.eplane", &P_cal_pr_eplane, &b_P_cal_pr_eplane);
   iChain->SetBranchAddress("P.dc.ntrack", &P_dc_ntrack, &b_P_dc_ntrack);
   iChain->SetBranchAddress("P.dc.x_fp", &P_dc_x_fp, &b_P_dc_x_fp);
   iChain->SetBranchAddress("P.dc.xp_fp", &P_dc_xp_fp, &b_P_dc_xp_fp);
   iChain->SetBranchAddress("P.dc.y_fp", &P_dc_y_fp, &b_P_dc_y_fp);
   iChain->SetBranchAddress("P.dc.yp_fp", &P_dc_yp_fp, &b_P_dc_yp_fp);
   iChain->SetBranchAddress("P.extcor.xsieve", &P_extcor_xsieve, &b_P_extcor_xsieve);
   iChain->SetBranchAddress("P.extcor.ysieve", &P_extcor_ysieve, &b_P_extcor_ysieve);
   iChain->SetBranchAddress("P.gtr.beta", &P_gtr_beta, &b_P_gtr_beta);
   iChain->SetBranchAddress("P.gtr.dp", &P_gtr_dp, &b_P_gtr_dp);
   iChain->SetBranchAddress("P.gtr.index", &P_gtr_index, &b_P_gtr_index);
   iChain->SetBranchAddress("P.gtr.ok", &P_gtr_ok, &b_P_gtr_ok);
   iChain->SetBranchAddress("P.gtr.p", &P_gtr_p, &b_P_gtr_p);
   iChain->SetBranchAddress("P.gtr.ph", &P_gtr_ph, &b_P_gtr_ph);
   iChain->SetBranchAddress("P.gtr.px", &P_gtr_px, &b_P_gtr_px);
   iChain->SetBranchAddress("P.gtr.py", &P_gtr_py, &b_P_gtr_py);
   iChain->SetBranchAddress("P.gtr.pz", &P_gtr_pz, &b_P_gtr_pz);
   iChain->SetBranchAddress("P.gtr.th", &P_gtr_th, &b_P_gtr_th);
   iChain->SetBranchAddress("P.gtr.x", &P_gtr_x, &b_P_gtr_x);
   iChain->SetBranchAddress("P.gtr.y", &P_gtr_y, &b_P_gtr_y);
   iChain->SetBranchAddress("P.hgcer.npeSum", &P_hgcer_npeSum, &b_P_hgcer_npeSum);
   iChain->SetBranchAddress("P.hod.beta", &P_hod_beta, &b_P_hod_beta);
   iChain->SetBranchAddress("P.hod.betachisqnotrack", &P_hod_betachisqnotrack, &b_P_hod_betachisqnotrack);
   iChain->SetBranchAddress("P.hod.betanotrack", &P_hod_betanotrack, &b_P_hod_betanotrack);
   iChain->SetBranchAddress("P.hod.fpHitsTime", &P_hod_fpHitsTime, &b_P_hod_fpHitsTime);
   iChain->SetBranchAddress("P.hod.goodscinhit", &P_hod_goodscinhit, &b_P_hod_goodscinhit);
   iChain->SetBranchAddress("P.hod.starttime", &P_hod_starttime, &b_P_hod_starttime);
   iChain->SetBranchAddress("P.ngcer.npeSum", &P_ngcer_npeSum, &b_P_ngcer_npeSum);
   iChain->SetBranchAddress("P.react.z", &P_react_z, &b_P_react_z);
   iChain->SetBranchAddress("T.coin.h1T_tdcMultiplicity", &T_coin_h1T_tdcMultiplicity, &b_T_coin_h1T_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.h1T_tdcTime", &T_coin_h1T_tdcTime, &b_T_coin_h1T_tdcTime);
   iChain->SetBranchAddress("T.coin.h1T_tdcTimeRaw", &T_coin_h1T_tdcTimeRaw, &b_T_coin_h1T_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.h1X_tdcMultiplicity", &T_coin_h1X_tdcMultiplicity, &b_T_coin_h1X_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.h1X_tdcTime", &T_coin_h1X_tdcTime, &b_T_coin_h1X_tdcTime);
   iChain->SetBranchAddress("T.coin.h1X_tdcTimeRaw", &T_coin_h1X_tdcTimeRaw, &b_T_coin_h1X_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.h1Y_tdcMultiplicity", &T_coin_h1Y_tdcMultiplicity, &b_T_coin_h1Y_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.h1Y_tdcTime", &T_coin_h1Y_tdcTime, &b_T_coin_h1Y_tdcTime);
   iChain->SetBranchAddress("T.coin.h1Y_tdcTimeRaw", &T_coin_h1Y_tdcTimeRaw, &b_T_coin_h1Y_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.h2T_tdcMultiplicity", &T_coin_h2T_tdcMultiplicity, &b_T_coin_h2T_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.h2T_tdcTime", &T_coin_h2T_tdcTime, &b_T_coin_h2T_tdcTime);
   iChain->SetBranchAddress("T.coin.h2T_tdcTimeRaw", &T_coin_h2T_tdcTimeRaw, &b_T_coin_h2T_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.h2X_tdcMultiplicity", &T_coin_h2X_tdcMultiplicity, &b_T_coin_h2X_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.h2X_tdcTime", &T_coin_h2X_tdcTime, &b_T_coin_h2X_tdcTime);
   iChain->SetBranchAddress("T.coin.h2X_tdcTimeRaw", &T_coin_h2X_tdcTimeRaw, &b_T_coin_h2X_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.h2Y_tdcMultiplicity", &T_coin_h2Y_tdcMultiplicity, &b_T_coin_h2Y_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.h2Y_tdcTime", &T_coin_h2Y_tdcTime, &b_T_coin_h2Y_tdcTime);
   iChain->SetBranchAddress("T.coin.h2Y_tdcTimeRaw", &T_coin_h2Y_tdcTimeRaw, &b_T_coin_h2Y_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hAER_adcMultiplicity", &T_coin_hAER_adcMultiplicity, &b_T_coin_hAER_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hAER_adcPed", &T_coin_hAER_adcPed, &b_T_coin_hAER_adcPed);
   iChain->SetBranchAddress("T.coin.hAER_adcPedRaw", &T_coin_hAER_adcPedRaw, &b_T_coin_hAER_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hAER_adcPulseAmp", &T_coin_hAER_adcPulseAmp, &b_T_coin_hAER_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hAER_adcPulseAmpRaw", &T_coin_hAER_adcPulseAmpRaw, &b_T_coin_hAER_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hAER_adcPulseInt", &T_coin_hAER_adcPulseInt, &b_T_coin_hAER_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hAER_adcPulseIntRaw", &T_coin_hAER_adcPulseIntRaw, &b_T_coin_hAER_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hAER_adcPulseTime", &T_coin_hAER_adcPulseTime, &b_T_coin_hAER_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hAER_adcPulseTimeRaw", &T_coin_hAER_adcPulseTimeRaw, &b_T_coin_hAER_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hASUM_adcMultiplicity", &T_coin_hASUM_adcMultiplicity, &b_T_coin_hASUM_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hASUM_adcPed", &T_coin_hASUM_adcPed, &b_T_coin_hASUM_adcPed);
   iChain->SetBranchAddress("T.coin.hASUM_adcPedRaw", &T_coin_hASUM_adcPedRaw, &b_T_coin_hASUM_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hASUM_adcPulseAmp", &T_coin_hASUM_adcPulseAmp, &b_T_coin_hASUM_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hASUM_adcPulseAmpRaw", &T_coin_hASUM_adcPulseAmpRaw, &b_T_coin_hASUM_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hASUM_adcPulseInt", &T_coin_hASUM_adcPulseInt, &b_T_coin_hASUM_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hASUM_adcPulseIntRaw", &T_coin_hASUM_adcPulseIntRaw, &b_T_coin_hASUM_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hASUM_adcPulseTime", &T_coin_hASUM_adcPulseTime, &b_T_coin_hASUM_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hASUM_adcPulseTimeRaw", &T_coin_hASUM_adcPulseTimeRaw, &b_T_coin_hASUM_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hASUM_tdcMultiplicity", &T_coin_hASUM_tdcMultiplicity, &b_T_coin_hASUM_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hASUM_tdcTime", &T_coin_hASUM_tdcTime, &b_T_coin_hASUM_tdcTime);
   iChain->SetBranchAddress("T.coin.hASUM_tdcTimeRaw", &T_coin_hASUM_tdcTimeRaw, &b_T_coin_hASUM_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hBSUM_adcMultiplicity", &T_coin_hBSUM_adcMultiplicity, &b_T_coin_hBSUM_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hBSUM_adcPed", &T_coin_hBSUM_adcPed, &b_T_coin_hBSUM_adcPed);
   iChain->SetBranchAddress("T.coin.hBSUM_adcPedRaw", &T_coin_hBSUM_adcPedRaw, &b_T_coin_hBSUM_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hBSUM_adcPulseAmp", &T_coin_hBSUM_adcPulseAmp, &b_T_coin_hBSUM_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hBSUM_adcPulseAmpRaw", &T_coin_hBSUM_adcPulseAmpRaw, &b_T_coin_hBSUM_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hBSUM_adcPulseInt", &T_coin_hBSUM_adcPulseInt, &b_T_coin_hBSUM_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hBSUM_adcPulseIntRaw", &T_coin_hBSUM_adcPulseIntRaw, &b_T_coin_hBSUM_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hBSUM_adcPulseTime", &T_coin_hBSUM_adcPulseTime, &b_T_coin_hBSUM_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hBSUM_adcPulseTimeRaw", &T_coin_hBSUM_adcPulseTimeRaw, &b_T_coin_hBSUM_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hBSUM_tdcMultiplicity", &T_coin_hBSUM_tdcMultiplicity, &b_T_coin_hBSUM_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hBSUM_tdcTime", &T_coin_hBSUM_tdcTime, &b_T_coin_hBSUM_tdcTime);
   iChain->SetBranchAddress("T.coin.hBSUM_tdcTimeRaw", &T_coin_hBSUM_tdcTimeRaw, &b_T_coin_hBSUM_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hCER_adcMultiplicity", &T_coin_hCER_adcMultiplicity, &b_T_coin_hCER_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hCER_adcPed", &T_coin_hCER_adcPed, &b_T_coin_hCER_adcPed);
   iChain->SetBranchAddress("T.coin.hCER_adcPedRaw", &T_coin_hCER_adcPedRaw, &b_T_coin_hCER_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hCER_adcPulseAmp", &T_coin_hCER_adcPulseAmp, &b_T_coin_hCER_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hCER_adcPulseAmpRaw", &T_coin_hCER_adcPulseAmpRaw, &b_T_coin_hCER_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hCER_adcPulseInt", &T_coin_hCER_adcPulseInt, &b_T_coin_hCER_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hCER_adcPulseIntRaw", &T_coin_hCER_adcPulseIntRaw, &b_T_coin_hCER_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hCER_adcPulseTime", &T_coin_hCER_adcPulseTime, &b_T_coin_hCER_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hCER_adcPulseTimeRaw", &T_coin_hCER_adcPulseTimeRaw, &b_T_coin_hCER_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hCER_tdcMultiplicity", &T_coin_hCER_tdcMultiplicity, &b_T_coin_hCER_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hCER_tdcTime", &T_coin_hCER_tdcTime, &b_T_coin_hCER_tdcTime);
   iChain->SetBranchAddress("T.coin.hCER_tdcTimeRaw", &T_coin_hCER_tdcTimeRaw, &b_T_coin_hCER_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hCSUM_adcMultiplicity", &T_coin_hCSUM_adcMultiplicity, &b_T_coin_hCSUM_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hCSUM_adcPed", &T_coin_hCSUM_adcPed, &b_T_coin_hCSUM_adcPed);
   iChain->SetBranchAddress("T.coin.hCSUM_adcPedRaw", &T_coin_hCSUM_adcPedRaw, &b_T_coin_hCSUM_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hCSUM_adcPulseAmp", &T_coin_hCSUM_adcPulseAmp, &b_T_coin_hCSUM_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hCSUM_adcPulseAmpRaw", &T_coin_hCSUM_adcPulseAmpRaw, &b_T_coin_hCSUM_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hCSUM_adcPulseInt", &T_coin_hCSUM_adcPulseInt, &b_T_coin_hCSUM_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hCSUM_adcPulseIntRaw", &T_coin_hCSUM_adcPulseIntRaw, &b_T_coin_hCSUM_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hCSUM_adcPulseTime", &T_coin_hCSUM_adcPulseTime, &b_T_coin_hCSUM_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hCSUM_adcPulseTimeRaw", &T_coin_hCSUM_adcPulseTimeRaw, &b_T_coin_hCSUM_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hCSUM_tdcMultiplicity", &T_coin_hCSUM_tdcMultiplicity, &b_T_coin_hCSUM_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hCSUM_tdcTime", &T_coin_hCSUM_tdcTime, &b_T_coin_hCSUM_tdcTime);
   iChain->SetBranchAddress("T.coin.hCSUM_tdcTimeRaw", &T_coin_hCSUM_tdcTimeRaw, &b_T_coin_hCSUM_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hDCREF1_tdcMultiplicity", &T_coin_hDCREF1_tdcMultiplicity, &b_T_coin_hDCREF1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hDCREF1_tdcTime", &T_coin_hDCREF1_tdcTime, &b_T_coin_hDCREF1_tdcTime);
   iChain->SetBranchAddress("T.coin.hDCREF1_tdcTimeRaw", &T_coin_hDCREF1_tdcTimeRaw, &b_T_coin_hDCREF1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hDCREF2_tdcMultiplicity", &T_coin_hDCREF2_tdcMultiplicity, &b_T_coin_hDCREF2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hDCREF2_tdcTime", &T_coin_hDCREF2_tdcTime, &b_T_coin_hDCREF2_tdcTime);
   iChain->SetBranchAddress("T.coin.hDCREF2_tdcTimeRaw", &T_coin_hDCREF2_tdcTimeRaw, &b_T_coin_hDCREF2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hDCREF3_tdcMultiplicity", &T_coin_hDCREF3_tdcMultiplicity, &b_T_coin_hDCREF3_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hDCREF3_tdcTime", &T_coin_hDCREF3_tdcTime, &b_T_coin_hDCREF3_tdcTime);
   iChain->SetBranchAddress("T.coin.hDCREF3_tdcTimeRaw", &T_coin_hDCREF3_tdcTimeRaw, &b_T_coin_hDCREF3_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hDCREF4_tdcMultiplicity", &T_coin_hDCREF4_tdcMultiplicity, &b_T_coin_hDCREF4_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hDCREF4_tdcTime", &T_coin_hDCREF4_tdcTime, &b_T_coin_hDCREF4_tdcTime);
   iChain->SetBranchAddress("T.coin.hDCREF4_tdcTimeRaw", &T_coin_hDCREF4_tdcTimeRaw, &b_T_coin_hDCREF4_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hDCREF5_tdcMultiplicity", &T_coin_hDCREF5_tdcMultiplicity, &b_T_coin_hDCREF5_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hDCREF5_tdcTime", &T_coin_hDCREF5_tdcTime, &b_T_coin_hDCREF5_tdcTime);
   iChain->SetBranchAddress("T.coin.hDCREF5_tdcTimeRaw", &T_coin_hDCREF5_tdcTimeRaw, &b_T_coin_hDCREF5_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hDSUM_adcMultiplicity", &T_coin_hDSUM_adcMultiplicity, &b_T_coin_hDSUM_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hDSUM_adcPed", &T_coin_hDSUM_adcPed, &b_T_coin_hDSUM_adcPed);
   iChain->SetBranchAddress("T.coin.hDSUM_adcPedRaw", &T_coin_hDSUM_adcPedRaw, &b_T_coin_hDSUM_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hDSUM_adcPulseAmp", &T_coin_hDSUM_adcPulseAmp, &b_T_coin_hDSUM_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hDSUM_adcPulseAmpRaw", &T_coin_hDSUM_adcPulseAmpRaw, &b_T_coin_hDSUM_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hDSUM_adcPulseInt", &T_coin_hDSUM_adcPulseInt, &b_T_coin_hDSUM_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hDSUM_adcPulseIntRaw", &T_coin_hDSUM_adcPulseIntRaw, &b_T_coin_hDSUM_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hDSUM_adcPulseTime", &T_coin_hDSUM_adcPulseTime, &b_T_coin_hDSUM_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hDSUM_adcPulseTimeRaw", &T_coin_hDSUM_adcPulseTimeRaw, &b_T_coin_hDSUM_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hDSUM_tdcMultiplicity", &T_coin_hDSUM_tdcMultiplicity, &b_T_coin_hDSUM_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hDSUM_tdcTime", &T_coin_hDSUM_tdcTime, &b_T_coin_hDSUM_tdcTime);
   iChain->SetBranchAddress("T.coin.hDSUM_tdcTimeRaw", &T_coin_hDSUM_tdcTimeRaw, &b_T_coin_hDSUM_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEDTM_tdcMultiplicity", &T_coin_hEDTM_tdcMultiplicity, &b_T_coin_hEDTM_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEDTM_tdcTime", &T_coin_hEDTM_tdcTime, &b_T_coin_hEDTM_tdcTime);
   iChain->SetBranchAddress("T.coin.hEDTM_tdcTimeRaw", &T_coin_hEDTM_tdcTimeRaw, &b_T_coin_hEDTM_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_CLEAN_ROC1_tdcMultiplicity", &T_coin_hEL_CLEAN_ROC1_tdcMultiplicity, &b_T_coin_hEL_CLEAN_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_CLEAN_ROC1_tdcTime", &T_coin_hEL_CLEAN_ROC1_tdcTime, &b_T_coin_hEL_CLEAN_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_CLEAN_ROC1_tdcTimeRaw", &T_coin_hEL_CLEAN_ROC1_tdcTimeRaw, &b_T_coin_hEL_CLEAN_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_CLEAN_ROC2_tdcMultiplicity", &T_coin_hEL_CLEAN_ROC2_tdcMultiplicity, &b_T_coin_hEL_CLEAN_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_CLEAN_ROC2_tdcTime", &T_coin_hEL_CLEAN_ROC2_tdcTime, &b_T_coin_hEL_CLEAN_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_CLEAN_ROC2_tdcTimeRaw", &T_coin_hEL_CLEAN_ROC2_tdcTimeRaw, &b_T_coin_hEL_CLEAN_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_HI_ROC1_tdcMultiplicity", &T_coin_hEL_HI_ROC1_tdcMultiplicity, &b_T_coin_hEL_HI_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_HI_ROC1_tdcTime", &T_coin_hEL_HI_ROC1_tdcTime, &b_T_coin_hEL_HI_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_HI_ROC1_tdcTimeRaw", &T_coin_hEL_HI_ROC1_tdcTimeRaw, &b_T_coin_hEL_HI_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_HI_ROC2_tdcMultiplicity", &T_coin_hEL_HI_ROC2_tdcMultiplicity, &b_T_coin_hEL_HI_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_HI_ROC2_tdcTime", &T_coin_hEL_HI_ROC2_tdcTime, &b_T_coin_hEL_HI_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_HI_ROC2_tdcTimeRaw", &T_coin_hEL_HI_ROC2_tdcTimeRaw, &b_T_coin_hEL_HI_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_LO_LO_ROC1_tdcMultiplicity", &T_coin_hEL_LO_LO_ROC1_tdcMultiplicity, &b_T_coin_hEL_LO_LO_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_LO_LO_ROC1_tdcTime", &T_coin_hEL_LO_LO_ROC1_tdcTime, &b_T_coin_hEL_LO_LO_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_LO_LO_ROC1_tdcTimeRaw", &T_coin_hEL_LO_LO_ROC1_tdcTimeRaw, &b_T_coin_hEL_LO_LO_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_LO_LO_ROC2_tdcMultiplicity", &T_coin_hEL_LO_LO_ROC2_tdcMultiplicity, &b_T_coin_hEL_LO_LO_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_LO_LO_ROC2_tdcTime", &T_coin_hEL_LO_LO_ROC2_tdcTime, &b_T_coin_hEL_LO_LO_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_LO_LO_ROC2_tdcTimeRaw", &T_coin_hEL_LO_LO_ROC2_tdcTimeRaw, &b_T_coin_hEL_LO_LO_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_LO_ROC1_tdcMultiplicity", &T_coin_hEL_LO_ROC1_tdcMultiplicity, &b_T_coin_hEL_LO_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_LO_ROC1_tdcTime", &T_coin_hEL_LO_ROC1_tdcTime, &b_T_coin_hEL_LO_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_LO_ROC1_tdcTimeRaw", &T_coin_hEL_LO_ROC1_tdcTimeRaw, &b_T_coin_hEL_LO_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_LO_ROC2_tdcMultiplicity", &T_coin_hEL_LO_ROC2_tdcMultiplicity, &b_T_coin_hEL_LO_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_LO_ROC2_tdcTime", &T_coin_hEL_LO_ROC2_tdcTime, &b_T_coin_hEL_LO_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_LO_ROC2_tdcTimeRaw", &T_coin_hEL_LO_ROC2_tdcTimeRaw, &b_T_coin_hEL_LO_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_REAL_ROC1_tdcMultiplicity", &T_coin_hEL_REAL_ROC1_tdcMultiplicity, &b_T_coin_hEL_REAL_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_REAL_ROC1_tdcTime", &T_coin_hEL_REAL_ROC1_tdcTime, &b_T_coin_hEL_REAL_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_REAL_ROC1_tdcTimeRaw", &T_coin_hEL_REAL_ROC1_tdcTimeRaw, &b_T_coin_hEL_REAL_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hEL_REAL_ROC2_tdcMultiplicity", &T_coin_hEL_REAL_ROC2_tdcMultiplicity, &b_T_coin_hEL_REAL_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hEL_REAL_ROC2_tdcTime", &T_coin_hEL_REAL_ROC2_tdcTime, &b_T_coin_hEL_REAL_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hEL_REAL_ROC2_tdcTimeRaw", &T_coin_hEL_REAL_ROC2_tdcTimeRaw, &b_T_coin_hEL_REAL_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcMultiplicity", &T_coin_hFADC_TREF_ROC1_adcMultiplicity, &b_T_coin_hFADC_TREF_ROC1_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcPed", &T_coin_hFADC_TREF_ROC1_adcPed, &b_T_coin_hFADC_TREF_ROC1_adcPed);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcPedRaw", &T_coin_hFADC_TREF_ROC1_adcPedRaw, &b_T_coin_hFADC_TREF_ROC1_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcPulseAmp", &T_coin_hFADC_TREF_ROC1_adcPulseAmp, &b_T_coin_hFADC_TREF_ROC1_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcPulseAmpRaw", &T_coin_hFADC_TREF_ROC1_adcPulseAmpRaw, &b_T_coin_hFADC_TREF_ROC1_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcPulseInt", &T_coin_hFADC_TREF_ROC1_adcPulseInt, &b_T_coin_hFADC_TREF_ROC1_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcPulseIntRaw", &T_coin_hFADC_TREF_ROC1_adcPulseIntRaw, &b_T_coin_hFADC_TREF_ROC1_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcPulseTime", &T_coin_hFADC_TREF_ROC1_adcPulseTime, &b_T_coin_hFADC_TREF_ROC1_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hFADC_TREF_ROC1_adcPulseTimeRaw", &T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw, &b_T_coin_hFADC_TREF_ROC1_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hHODO_RF_tdcMultiplicity", &T_coin_hHODO_RF_tdcMultiplicity, &b_T_coin_hHODO_RF_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hHODO_RF_tdcTime", &T_coin_hHODO_RF_tdcTime, &b_T_coin_hHODO_RF_tdcTime);
   iChain->SetBranchAddress("T.coin.hHODO_RF_tdcTimeRaw", &T_coin_hHODO_RF_tdcTimeRaw, &b_T_coin_hHODO_RF_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRE100_ROC1_tdcMultiplicity", &T_coin_hPRE100_ROC1_tdcMultiplicity, &b_T_coin_hPRE100_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRE100_ROC1_tdcTime", &T_coin_hPRE100_ROC1_tdcTime, &b_T_coin_hPRE100_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRE100_ROC1_tdcTimeRaw", &T_coin_hPRE100_ROC1_tdcTimeRaw, &b_T_coin_hPRE100_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRE100_ROC2_tdcMultiplicity", &T_coin_hPRE100_ROC2_tdcMultiplicity, &b_T_coin_hPRE100_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRE100_ROC2_tdcTime", &T_coin_hPRE100_ROC2_tdcTime, &b_T_coin_hPRE100_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRE100_ROC2_tdcTimeRaw", &T_coin_hPRE100_ROC2_tdcTimeRaw, &b_T_coin_hPRE100_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRE150_ROC1_tdcMultiplicity", &T_coin_hPRE150_ROC1_tdcMultiplicity, &b_T_coin_hPRE150_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRE150_ROC1_tdcTime", &T_coin_hPRE150_ROC1_tdcTime, &b_T_coin_hPRE150_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRE150_ROC1_tdcTimeRaw", &T_coin_hPRE150_ROC1_tdcTimeRaw, &b_T_coin_hPRE150_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRE150_ROC2_tdcMultiplicity", &T_coin_hPRE150_ROC2_tdcMultiplicity, &b_T_coin_hPRE150_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRE150_ROC2_tdcTime", &T_coin_hPRE150_ROC2_tdcTime, &b_T_coin_hPRE150_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRE150_ROC2_tdcTimeRaw", &T_coin_hPRE150_ROC2_tdcTimeRaw, &b_T_coin_hPRE150_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRE200_ROC1_tdcMultiplicity", &T_coin_hPRE200_ROC1_tdcMultiplicity, &b_T_coin_hPRE200_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRE200_ROC1_tdcTime", &T_coin_hPRE200_ROC1_tdcTime, &b_T_coin_hPRE200_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRE200_ROC1_tdcTimeRaw", &T_coin_hPRE200_ROC1_tdcTimeRaw, &b_T_coin_hPRE200_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRE200_ROC2_tdcMultiplicity", &T_coin_hPRE200_ROC2_tdcMultiplicity, &b_T_coin_hPRE200_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRE200_ROC2_tdcTime", &T_coin_hPRE200_ROC2_tdcTime, &b_T_coin_hPRE200_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRE200_ROC2_tdcTimeRaw", &T_coin_hPRE200_ROC2_tdcTimeRaw, &b_T_coin_hPRE200_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRE40_ROC1_tdcMultiplicity", &T_coin_hPRE40_ROC1_tdcMultiplicity, &b_T_coin_hPRE40_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRE40_ROC1_tdcTime", &T_coin_hPRE40_ROC1_tdcTime, &b_T_coin_hPRE40_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRE40_ROC1_tdcTimeRaw", &T_coin_hPRE40_ROC1_tdcTimeRaw, &b_T_coin_hPRE40_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRE40_ROC2_tdcMultiplicity", &T_coin_hPRE40_ROC2_tdcMultiplicity, &b_T_coin_hPRE40_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRE40_ROC2_tdcTime", &T_coin_hPRE40_ROC2_tdcTime, &b_T_coin_hPRE40_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRE40_ROC2_tdcTimeRaw", &T_coin_hPRE40_ROC2_tdcTimeRaw, &b_T_coin_hPRE40_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRHI_tdcMultiplicity", &T_coin_hPRHI_tdcMultiplicity, &b_T_coin_hPRHI_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRHI_tdcTime", &T_coin_hPRHI_tdcTime, &b_T_coin_hPRHI_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRHI_tdcTimeRaw", &T_coin_hPRHI_tdcTimeRaw, &b_T_coin_hPRHI_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPRLO_tdcMultiplicity", &T_coin_hPRLO_tdcMultiplicity, &b_T_coin_hPRLO_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPRLO_tdcTime", &T_coin_hPRLO_tdcTime, &b_T_coin_hPRLO_tdcTime);
   iChain->SetBranchAddress("T.coin.hPRLO_tdcTimeRaw", &T_coin_hPRLO_tdcTimeRaw, &b_T_coin_hPRLO_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcMultiplicity", &T_coin_hPSHWR_adcMultiplicity, &b_T_coin_hPSHWR_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcPed", &T_coin_hPSHWR_adcPed, &b_T_coin_hPSHWR_adcPed);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcPedRaw", &T_coin_hPSHWR_adcPedRaw, &b_T_coin_hPSHWR_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcPulseAmp", &T_coin_hPSHWR_adcPulseAmp, &b_T_coin_hPSHWR_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcPulseAmpRaw", &T_coin_hPSHWR_adcPulseAmpRaw, &b_T_coin_hPSHWR_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcPulseInt", &T_coin_hPSHWR_adcPulseInt, &b_T_coin_hPSHWR_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcPulseIntRaw", &T_coin_hPSHWR_adcPulseIntRaw, &b_T_coin_hPSHWR_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcPulseTime", &T_coin_hPSHWR_adcPulseTime, &b_T_coin_hPSHWR_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hPSHWR_adcPulseTimeRaw", &T_coin_hPSHWR_adcPulseTimeRaw, &b_T_coin_hPSHWR_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hRF_tdcMultiplicity", &T_coin_hRF_tdcMultiplicity, &b_T_coin_hRF_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hRF_tdcTime", &T_coin_hRF_tdcTime, &b_T_coin_hRF_tdcTime);
   iChain->SetBranchAddress("T.coin.hRF_tdcTimeRaw", &T_coin_hRF_tdcTimeRaw, &b_T_coin_hRF_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hSHWR_adcMultiplicity", &T_coin_hSHWR_adcMultiplicity, &b_T_coin_hSHWR_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.hSHWR_adcPed", &T_coin_hSHWR_adcPed, &b_T_coin_hSHWR_adcPed);
   iChain->SetBranchAddress("T.coin.hSHWR_adcPedRaw", &T_coin_hSHWR_adcPedRaw, &b_T_coin_hSHWR_adcPedRaw);
   iChain->SetBranchAddress("T.coin.hSHWR_adcPulseAmp", &T_coin_hSHWR_adcPulseAmp, &b_T_coin_hSHWR_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.hSHWR_adcPulseAmpRaw", &T_coin_hSHWR_adcPulseAmpRaw, &b_T_coin_hSHWR_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.hSHWR_adcPulseInt", &T_coin_hSHWR_adcPulseInt, &b_T_coin_hSHWR_adcPulseInt);
   iChain->SetBranchAddress("T.coin.hSHWR_adcPulseIntRaw", &T_coin_hSHWR_adcPulseIntRaw, &b_T_coin_hSHWR_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.hSHWR_adcPulseTime", &T_coin_hSHWR_adcPulseTime, &b_T_coin_hSHWR_adcPulseTime);
   iChain->SetBranchAddress("T.coin.hSHWR_adcPulseTimeRaw", &T_coin_hSHWR_adcPulseTimeRaw, &b_T_coin_hSHWR_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.hSHWR_tdcMultiplicity", &T_coin_hSHWR_tdcMultiplicity, &b_T_coin_hSHWR_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hSHWR_tdcTime", &T_coin_hSHWR_tdcTime, &b_T_coin_hSHWR_tdcTime);
   iChain->SetBranchAddress("T.coin.hSHWR_tdcTimeRaw", &T_coin_hSHWR_tdcTimeRaw, &b_T_coin_hSHWR_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hSTOF_ROC1_tdcMultiplicity", &T_coin_hSTOF_ROC1_tdcMultiplicity, &b_T_coin_hSTOF_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hSTOF_ROC1_tdcTime", &T_coin_hSTOF_ROC1_tdcTime, &b_T_coin_hSTOF_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hSTOF_ROC1_tdcTimeRaw", &T_coin_hSTOF_ROC1_tdcTimeRaw, &b_T_coin_hSTOF_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hSTOF_ROC2_tdcMultiplicity", &T_coin_hSTOF_ROC2_tdcMultiplicity, &b_T_coin_hSTOF_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hSTOF_ROC2_tdcTime", &T_coin_hSTOF_ROC2_tdcTime, &b_T_coin_hSTOF_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hSTOF_ROC2_tdcTimeRaw", &T_coin_hSTOF_ROC2_tdcTimeRaw, &b_T_coin_hSTOF_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hT1_tdcMultiplicity", &T_coin_hT1_tdcMultiplicity, &b_T_coin_hT1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hT1_tdcTime", &T_coin_hT1_tdcTime, &b_T_coin_hT1_tdcTime);
   iChain->SetBranchAddress("T.coin.hT1_tdcTimeRaw", &T_coin_hT1_tdcTimeRaw, &b_T_coin_hT1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hT2_tdcMultiplicity", &T_coin_hT2_tdcMultiplicity, &b_T_coin_hT2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hT2_tdcTime", &T_coin_hT2_tdcTime, &b_T_coin_hT2_tdcTime);
   iChain->SetBranchAddress("T.coin.hT2_tdcTimeRaw", &T_coin_hT2_tdcTimeRaw, &b_T_coin_hT2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG1_ROC1_tdcMultiplicity", &T_coin_hTRIG1_ROC1_tdcMultiplicity, &b_T_coin_hTRIG1_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG1_ROC1_tdcTime", &T_coin_hTRIG1_ROC1_tdcTime, &b_T_coin_hTRIG1_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG1_ROC1_tdcTimeRaw", &T_coin_hTRIG1_ROC1_tdcTimeRaw, &b_T_coin_hTRIG1_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG1_ROC2_tdcMultiplicity", &T_coin_hTRIG1_ROC2_tdcMultiplicity, &b_T_coin_hTRIG1_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG1_ROC2_tdcTime", &T_coin_hTRIG1_ROC2_tdcTime, &b_T_coin_hTRIG1_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG1_ROC2_tdcTimeRaw", &T_coin_hTRIG1_ROC2_tdcTimeRaw, &b_T_coin_hTRIG1_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG2_ROC1_tdcMultiplicity", &T_coin_hTRIG2_ROC1_tdcMultiplicity, &b_T_coin_hTRIG2_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG2_ROC1_tdcTime", &T_coin_hTRIG2_ROC1_tdcTime, &b_T_coin_hTRIG2_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG2_ROC1_tdcTimeRaw", &T_coin_hTRIG2_ROC1_tdcTimeRaw, &b_T_coin_hTRIG2_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG2_ROC2_tdcMultiplicity", &T_coin_hTRIG2_ROC2_tdcMultiplicity, &b_T_coin_hTRIG2_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG2_ROC2_tdcTime", &T_coin_hTRIG2_ROC2_tdcTime, &b_T_coin_hTRIG2_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG2_ROC2_tdcTimeRaw", &T_coin_hTRIG2_ROC2_tdcTimeRaw, &b_T_coin_hTRIG2_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG3_ROC1_tdcMultiplicity", &T_coin_hTRIG3_ROC1_tdcMultiplicity, &b_T_coin_hTRIG3_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG3_ROC1_tdcTime", &T_coin_hTRIG3_ROC1_tdcTime, &b_T_coin_hTRIG3_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG3_ROC1_tdcTimeRaw", &T_coin_hTRIG3_ROC1_tdcTimeRaw, &b_T_coin_hTRIG3_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG3_ROC2_tdcMultiplicity", &T_coin_hTRIG3_ROC2_tdcMultiplicity, &b_T_coin_hTRIG3_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG3_ROC2_tdcTime", &T_coin_hTRIG3_ROC2_tdcTime, &b_T_coin_hTRIG3_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG3_ROC2_tdcTimeRaw", &T_coin_hTRIG3_ROC2_tdcTimeRaw, &b_T_coin_hTRIG3_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG4_ROC1_tdcMultiplicity", &T_coin_hTRIG4_ROC1_tdcMultiplicity, &b_T_coin_hTRIG4_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG4_ROC1_tdcTime", &T_coin_hTRIG4_ROC1_tdcTime, &b_T_coin_hTRIG4_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG4_ROC1_tdcTimeRaw", &T_coin_hTRIG4_ROC1_tdcTimeRaw, &b_T_coin_hTRIG4_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG4_ROC2_tdcMultiplicity", &T_coin_hTRIG4_ROC2_tdcMultiplicity, &b_T_coin_hTRIG4_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG4_ROC2_tdcTime", &T_coin_hTRIG4_ROC2_tdcTime, &b_T_coin_hTRIG4_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG4_ROC2_tdcTimeRaw", &T_coin_hTRIG4_ROC2_tdcTimeRaw, &b_T_coin_hTRIG4_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG5_ROC1_tdcMultiplicity", &T_coin_hTRIG5_ROC1_tdcMultiplicity, &b_T_coin_hTRIG5_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG5_ROC1_tdcTime", &T_coin_hTRIG5_ROC1_tdcTime, &b_T_coin_hTRIG5_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG5_ROC1_tdcTimeRaw", &T_coin_hTRIG5_ROC1_tdcTimeRaw, &b_T_coin_hTRIG5_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG5_ROC2_tdcMultiplicity", &T_coin_hTRIG5_ROC2_tdcMultiplicity, &b_T_coin_hTRIG5_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG5_ROC2_tdcTime", &T_coin_hTRIG5_ROC2_tdcTime, &b_T_coin_hTRIG5_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG5_ROC2_tdcTimeRaw", &T_coin_hTRIG5_ROC2_tdcTimeRaw, &b_T_coin_hTRIG5_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG6_ROC1_tdcMultiplicity", &T_coin_hTRIG6_ROC1_tdcMultiplicity, &b_T_coin_hTRIG6_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG6_ROC1_tdcTime", &T_coin_hTRIG6_ROC1_tdcTime, &b_T_coin_hTRIG6_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG6_ROC1_tdcTimeRaw", &T_coin_hTRIG6_ROC1_tdcTimeRaw, &b_T_coin_hTRIG6_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.hTRIG6_ROC2_tdcMultiplicity", &T_coin_hTRIG6_ROC2_tdcMultiplicity, &b_T_coin_hTRIG6_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.hTRIG6_ROC2_tdcTime", &T_coin_hTRIG6_ROC2_tdcTime, &b_T_coin_hTRIG6_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.hTRIG6_ROC2_tdcTimeRaw", &T_coin_hTRIG6_ROC2_tdcTimeRaw, &b_T_coin_hTRIG6_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.p1T_tdcMultiplicity", &T_coin_p1T_tdcMultiplicity, &b_T_coin_p1T_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.p1T_tdcTime", &T_coin_p1T_tdcTime, &b_T_coin_p1T_tdcTime);
   iChain->SetBranchAddress("T.coin.p1T_tdcTimeRaw", &T_coin_p1T_tdcTimeRaw, &b_T_coin_p1T_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.p1X_tdcMultiplicity", &T_coin_p1X_tdcMultiplicity, &b_T_coin_p1X_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.p1X_tdcTime", &T_coin_p1X_tdcTime, &b_T_coin_p1X_tdcTime);
   iChain->SetBranchAddress("T.coin.p1X_tdcTimeRaw", &T_coin_p1X_tdcTimeRaw, &b_T_coin_p1X_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.p1Y_tdcMultiplicity", &T_coin_p1Y_tdcMultiplicity, &b_T_coin_p1Y_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.p1Y_tdcTime", &T_coin_p1Y_tdcTime, &b_T_coin_p1Y_tdcTime);
   iChain->SetBranchAddress("T.coin.p1Y_tdcTimeRaw", &T_coin_p1Y_tdcTimeRaw, &b_T_coin_p1Y_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.p2T_tdcMultiplicity", &T_coin_p2T_tdcMultiplicity, &b_T_coin_p2T_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.p2T_tdcTime", &T_coin_p2T_tdcTime, &b_T_coin_p2T_tdcTime);
   iChain->SetBranchAddress("T.coin.p2T_tdcTimeRaw", &T_coin_p2T_tdcTimeRaw, &b_T_coin_p2T_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.p2X_tdcMultiplicity", &T_coin_p2X_tdcMultiplicity, &b_T_coin_p2X_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.p2X_tdcTime", &T_coin_p2X_tdcTime, &b_T_coin_p2X_tdcTime);
   iChain->SetBranchAddress("T.coin.p2X_tdcTimeRaw", &T_coin_p2X_tdcTimeRaw, &b_T_coin_p2X_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.p2Y_tdcMultiplicity", &T_coin_p2Y_tdcMultiplicity, &b_T_coin_p2Y_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.p2Y_tdcTime", &T_coin_p2Y_tdcTime, &b_T_coin_p2Y_tdcTime);
   iChain->SetBranchAddress("T.coin.p2Y_tdcTimeRaw", &T_coin_p2Y_tdcTimeRaw, &b_T_coin_p2Y_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pAER_adcMultiplicity", &T_coin_pAER_adcMultiplicity, &b_T_coin_pAER_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pAER_adcPed", &T_coin_pAER_adcPed, &b_T_coin_pAER_adcPed);
   iChain->SetBranchAddress("T.coin.pAER_adcPedRaw", &T_coin_pAER_adcPedRaw, &b_T_coin_pAER_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pAER_adcPulseAmp", &T_coin_pAER_adcPulseAmp, &b_T_coin_pAER_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pAER_adcPulseAmpRaw", &T_coin_pAER_adcPulseAmpRaw, &b_T_coin_pAER_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pAER_adcPulseInt", &T_coin_pAER_adcPulseInt, &b_T_coin_pAER_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pAER_adcPulseIntRaw", &T_coin_pAER_adcPulseIntRaw, &b_T_coin_pAER_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pAER_adcPulseTime", &T_coin_pAER_adcPulseTime, &b_T_coin_pAER_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pAER_adcPulseTimeRaw", &T_coin_pAER_adcPulseTimeRaw, &b_T_coin_pAER_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pAER_tdcMultiplicity", &T_coin_pAER_tdcMultiplicity, &b_T_coin_pAER_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pAER_tdcTime", &T_coin_pAER_tdcTime, &b_T_coin_pAER_tdcTime);
   iChain->SetBranchAddress("T.coin.pAER_tdcTimeRaw", &T_coin_pAER_tdcTimeRaw, &b_T_coin_pAER_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF10_tdcMultiplicity", &T_coin_pDCREF10_tdcMultiplicity, &b_T_coin_pDCREF10_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF10_tdcTime", &T_coin_pDCREF10_tdcTime, &b_T_coin_pDCREF10_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF10_tdcTimeRaw", &T_coin_pDCREF10_tdcTimeRaw, &b_T_coin_pDCREF10_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF1_tdcMultiplicity", &T_coin_pDCREF1_tdcMultiplicity, &b_T_coin_pDCREF1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF1_tdcTime", &T_coin_pDCREF1_tdcTime, &b_T_coin_pDCREF1_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF1_tdcTimeRaw", &T_coin_pDCREF1_tdcTimeRaw, &b_T_coin_pDCREF1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF2_tdcMultiplicity", &T_coin_pDCREF2_tdcMultiplicity, &b_T_coin_pDCREF2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF2_tdcTime", &T_coin_pDCREF2_tdcTime, &b_T_coin_pDCREF2_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF2_tdcTimeRaw", &T_coin_pDCREF2_tdcTimeRaw, &b_T_coin_pDCREF2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF3_tdcMultiplicity", &T_coin_pDCREF3_tdcMultiplicity, &b_T_coin_pDCREF3_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF3_tdcTime", &T_coin_pDCREF3_tdcTime, &b_T_coin_pDCREF3_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF3_tdcTimeRaw", &T_coin_pDCREF3_tdcTimeRaw, &b_T_coin_pDCREF3_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF4_tdcMultiplicity", &T_coin_pDCREF4_tdcMultiplicity, &b_T_coin_pDCREF4_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF4_tdcTime", &T_coin_pDCREF4_tdcTime, &b_T_coin_pDCREF4_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF4_tdcTimeRaw", &T_coin_pDCREF4_tdcTimeRaw, &b_T_coin_pDCREF4_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF5_tdcMultiplicity", &T_coin_pDCREF5_tdcMultiplicity, &b_T_coin_pDCREF5_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF5_tdcTime", &T_coin_pDCREF5_tdcTime, &b_T_coin_pDCREF5_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF5_tdcTimeRaw", &T_coin_pDCREF5_tdcTimeRaw, &b_T_coin_pDCREF5_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF6_tdcMultiplicity", &T_coin_pDCREF6_tdcMultiplicity, &b_T_coin_pDCREF6_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF6_tdcTime", &T_coin_pDCREF6_tdcTime, &b_T_coin_pDCREF6_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF6_tdcTimeRaw", &T_coin_pDCREF6_tdcTimeRaw, &b_T_coin_pDCREF6_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF7_tdcMultiplicity", &T_coin_pDCREF7_tdcMultiplicity, &b_T_coin_pDCREF7_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF7_tdcTime", &T_coin_pDCREF7_tdcTime, &b_T_coin_pDCREF7_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF7_tdcTimeRaw", &T_coin_pDCREF7_tdcTimeRaw, &b_T_coin_pDCREF7_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF8_tdcMultiplicity", &T_coin_pDCREF8_tdcMultiplicity, &b_T_coin_pDCREF8_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF8_tdcTime", &T_coin_pDCREF8_tdcTime, &b_T_coin_pDCREF8_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF8_tdcTimeRaw", &T_coin_pDCREF8_tdcTimeRaw, &b_T_coin_pDCREF8_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pDCREF9_tdcMultiplicity", &T_coin_pDCREF9_tdcMultiplicity, &b_T_coin_pDCREF9_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pDCREF9_tdcTime", &T_coin_pDCREF9_tdcTime, &b_T_coin_pDCREF9_tdcTime);
   iChain->SetBranchAddress("T.coin.pDCREF9_tdcTimeRaw", &T_coin_pDCREF9_tdcTimeRaw, &b_T_coin_pDCREF9_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEDTM_tdcMultiplicity", &T_coin_pEDTM_tdcMultiplicity, &b_T_coin_pEDTM_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEDTM_tdcTime", &T_coin_pEDTM_tdcTime, &b_T_coin_pEDTM_tdcTime);
   iChain->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw", &T_coin_pEDTM_tdcTimeRaw, &b_T_coin_pEDTM_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_CLEAN_ROC1_tdcMultiplicity", &T_coin_pEL_CLEAN_ROC1_tdcMultiplicity, &b_T_coin_pEL_CLEAN_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_CLEAN_ROC1_tdcTime", &T_coin_pEL_CLEAN_ROC1_tdcTime, &b_T_coin_pEL_CLEAN_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_CLEAN_ROC1_tdcTimeRaw", &T_coin_pEL_CLEAN_ROC1_tdcTimeRaw, &b_T_coin_pEL_CLEAN_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_CLEAN_ROC2_tdcMultiplicity", &T_coin_pEL_CLEAN_ROC2_tdcMultiplicity, &b_T_coin_pEL_CLEAN_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_CLEAN_ROC2_tdcTime", &T_coin_pEL_CLEAN_ROC2_tdcTime, &b_T_coin_pEL_CLEAN_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_CLEAN_ROC2_tdcTimeRaw", &T_coin_pEL_CLEAN_ROC2_tdcTimeRaw, &b_T_coin_pEL_CLEAN_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_HI_ROC1_tdcMultiplicity", &T_coin_pEL_HI_ROC1_tdcMultiplicity, &b_T_coin_pEL_HI_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_HI_ROC1_tdcTime", &T_coin_pEL_HI_ROC1_tdcTime, &b_T_coin_pEL_HI_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_HI_ROC1_tdcTimeRaw", &T_coin_pEL_HI_ROC1_tdcTimeRaw, &b_T_coin_pEL_HI_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_HI_ROC2_tdcMultiplicity", &T_coin_pEL_HI_ROC2_tdcMultiplicity, &b_T_coin_pEL_HI_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_HI_ROC2_tdcTime", &T_coin_pEL_HI_ROC2_tdcTime, &b_T_coin_pEL_HI_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_HI_ROC2_tdcTimeRaw", &T_coin_pEL_HI_ROC2_tdcTimeRaw, &b_T_coin_pEL_HI_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_LO_LO_ROC1_tdcMultiplicity", &T_coin_pEL_LO_LO_ROC1_tdcMultiplicity, &b_T_coin_pEL_LO_LO_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_LO_LO_ROC1_tdcTime", &T_coin_pEL_LO_LO_ROC1_tdcTime, &b_T_coin_pEL_LO_LO_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_LO_LO_ROC1_tdcTimeRaw", &T_coin_pEL_LO_LO_ROC1_tdcTimeRaw, &b_T_coin_pEL_LO_LO_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_LO_LO_ROC2_tdcMultiplicity", &T_coin_pEL_LO_LO_ROC2_tdcMultiplicity, &b_T_coin_pEL_LO_LO_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_LO_LO_ROC2_tdcTime", &T_coin_pEL_LO_LO_ROC2_tdcTime, &b_T_coin_pEL_LO_LO_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_LO_LO_ROC2_tdcTimeRaw", &T_coin_pEL_LO_LO_ROC2_tdcTimeRaw, &b_T_coin_pEL_LO_LO_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_LO_ROC1_tdcMultiplicity", &T_coin_pEL_LO_ROC1_tdcMultiplicity, &b_T_coin_pEL_LO_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_LO_ROC1_tdcTime", &T_coin_pEL_LO_ROC1_tdcTime, &b_T_coin_pEL_LO_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_LO_ROC1_tdcTimeRaw", &T_coin_pEL_LO_ROC1_tdcTimeRaw, &b_T_coin_pEL_LO_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_LO_ROC2_tdcMultiplicity", &T_coin_pEL_LO_ROC2_tdcMultiplicity, &b_T_coin_pEL_LO_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_LO_ROC2_tdcTime", &T_coin_pEL_LO_ROC2_tdcTime, &b_T_coin_pEL_LO_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_LO_ROC2_tdcTimeRaw", &T_coin_pEL_LO_ROC2_tdcTimeRaw, &b_T_coin_pEL_LO_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_REAL_ROC1_tdcMultiplicity", &T_coin_pEL_REAL_ROC1_tdcMultiplicity, &b_T_coin_pEL_REAL_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_REAL_ROC1_tdcTime", &T_coin_pEL_REAL_ROC1_tdcTime, &b_T_coin_pEL_REAL_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_REAL_ROC1_tdcTimeRaw", &T_coin_pEL_REAL_ROC1_tdcTimeRaw, &b_T_coin_pEL_REAL_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pEL_REAL_ROC2_tdcMultiplicity", &T_coin_pEL_REAL_ROC2_tdcMultiplicity, &b_T_coin_pEL_REAL_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pEL_REAL_ROC2_tdcTime", &T_coin_pEL_REAL_ROC2_tdcTime, &b_T_coin_pEL_REAL_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pEL_REAL_ROC2_tdcTimeRaw", &T_coin_pEL_REAL_ROC2_tdcTimeRaw, &b_T_coin_pEL_REAL_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcMultiplicity", &T_coin_pFADC_TREF_ROC2_adcMultiplicity, &b_T_coin_pFADC_TREF_ROC2_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcPed", &T_coin_pFADC_TREF_ROC2_adcPed, &b_T_coin_pFADC_TREF_ROC2_adcPed);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcPedRaw", &T_coin_pFADC_TREF_ROC2_adcPedRaw, &b_T_coin_pFADC_TREF_ROC2_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcPulseAmp", &T_coin_pFADC_TREF_ROC2_adcPulseAmp, &b_T_coin_pFADC_TREF_ROC2_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcPulseAmpRaw", &T_coin_pFADC_TREF_ROC2_adcPulseAmpRaw, &b_T_coin_pFADC_TREF_ROC2_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcPulseInt", &T_coin_pFADC_TREF_ROC2_adcPulseInt, &b_T_coin_pFADC_TREF_ROC2_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcPulseIntRaw", &T_coin_pFADC_TREF_ROC2_adcPulseIntRaw, &b_T_coin_pFADC_TREF_ROC2_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcPulseTime", &T_coin_pFADC_TREF_ROC2_adcPulseTime, &b_T_coin_pFADC_TREF_ROC2_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pFADC_TREF_ROC2_adcPulseTimeRaw", &T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw, &b_T_coin_pFADC_TREF_ROC2_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcMultiplicity", &T_coin_pHEL_MPS_adcMultiplicity, &b_T_coin_pHEL_MPS_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcPed", &T_coin_pHEL_MPS_adcPed, &b_T_coin_pHEL_MPS_adcPed);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcPedRaw", &T_coin_pHEL_MPS_adcPedRaw, &b_T_coin_pHEL_MPS_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcPulseAmp", &T_coin_pHEL_MPS_adcPulseAmp, &b_T_coin_pHEL_MPS_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcPulseAmpRaw", &T_coin_pHEL_MPS_adcPulseAmpRaw, &b_T_coin_pHEL_MPS_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcPulseInt", &T_coin_pHEL_MPS_adcPulseInt, &b_T_coin_pHEL_MPS_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcPulseIntRaw", &T_coin_pHEL_MPS_adcPulseIntRaw, &b_T_coin_pHEL_MPS_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcPulseTime", &T_coin_pHEL_MPS_adcPulseTime, &b_T_coin_pHEL_MPS_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pHEL_MPS_adcPulseTimeRaw", &T_coin_pHEL_MPS_adcPulseTimeRaw, &b_T_coin_pHEL_MPS_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcMultiplicity", &T_coin_pHEL_NEG_adcMultiplicity, &b_T_coin_pHEL_NEG_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcPed", &T_coin_pHEL_NEG_adcPed, &b_T_coin_pHEL_NEG_adcPed);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcPedRaw", &T_coin_pHEL_NEG_adcPedRaw, &b_T_coin_pHEL_NEG_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcPulseAmp", &T_coin_pHEL_NEG_adcPulseAmp, &b_T_coin_pHEL_NEG_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcPulseAmpRaw", &T_coin_pHEL_NEG_adcPulseAmpRaw, &b_T_coin_pHEL_NEG_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcPulseInt", &T_coin_pHEL_NEG_adcPulseInt, &b_T_coin_pHEL_NEG_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcPulseIntRaw", &T_coin_pHEL_NEG_adcPulseIntRaw, &b_T_coin_pHEL_NEG_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcPulseTime", &T_coin_pHEL_NEG_adcPulseTime, &b_T_coin_pHEL_NEG_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pHEL_NEG_adcPulseTimeRaw", &T_coin_pHEL_NEG_adcPulseTimeRaw, &b_T_coin_pHEL_NEG_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcMultiplicity", &T_coin_pHEL_POS_adcMultiplicity, &b_T_coin_pHEL_POS_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcPed", &T_coin_pHEL_POS_adcPed, &b_T_coin_pHEL_POS_adcPed);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcPedRaw", &T_coin_pHEL_POS_adcPedRaw, &b_T_coin_pHEL_POS_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcPulseAmp", &T_coin_pHEL_POS_adcPulseAmp, &b_T_coin_pHEL_POS_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcPulseAmpRaw", &T_coin_pHEL_POS_adcPulseAmpRaw, &b_T_coin_pHEL_POS_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcPulseInt", &T_coin_pHEL_POS_adcPulseInt, &b_T_coin_pHEL_POS_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcPulseIntRaw", &T_coin_pHEL_POS_adcPulseIntRaw, &b_T_coin_pHEL_POS_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcPulseTime", &T_coin_pHEL_POS_adcPulseTime, &b_T_coin_pHEL_POS_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pHEL_POS_adcPulseTimeRaw", &T_coin_pHEL_POS_adcPulseTimeRaw, &b_T_coin_pHEL_POS_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcMultiplicity", &T_coin_pHGCER_MOD_adcMultiplicity, &b_T_coin_pHGCER_MOD_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcPed", &T_coin_pHGCER_MOD_adcPed, &b_T_coin_pHGCER_MOD_adcPed);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcPedRaw", &T_coin_pHGCER_MOD_adcPedRaw, &b_T_coin_pHGCER_MOD_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcPulseAmp", &T_coin_pHGCER_MOD_adcPulseAmp, &b_T_coin_pHGCER_MOD_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcPulseAmpRaw", &T_coin_pHGCER_MOD_adcPulseAmpRaw, &b_T_coin_pHGCER_MOD_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcPulseInt", &T_coin_pHGCER_MOD_adcPulseInt, &b_T_coin_pHGCER_MOD_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcPulseIntRaw", &T_coin_pHGCER_MOD_adcPulseIntRaw, &b_T_coin_pHGCER_MOD_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcPulseTime", &T_coin_pHGCER_MOD_adcPulseTime, &b_T_coin_pHGCER_MOD_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pHGCER_MOD_adcPulseTimeRaw", &T_coin_pHGCER_MOD_adcPulseTimeRaw, &b_T_coin_pHGCER_MOD_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_adcMultiplicity", &T_coin_pHGCER_adcMultiplicity, &b_T_coin_pHGCER_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pHGCER_adcPed", &T_coin_pHGCER_adcPed, &b_T_coin_pHGCER_adcPed);
   iChain->SetBranchAddress("T.coin.pHGCER_adcPedRaw", &T_coin_pHGCER_adcPedRaw, &b_T_coin_pHGCER_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_adcPulseAmp", &T_coin_pHGCER_adcPulseAmp, &b_T_coin_pHGCER_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pHGCER_adcPulseAmpRaw", &T_coin_pHGCER_adcPulseAmpRaw, &b_T_coin_pHGCER_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_adcPulseInt", &T_coin_pHGCER_adcPulseInt, &b_T_coin_pHGCER_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pHGCER_adcPulseIntRaw", &T_coin_pHGCER_adcPulseIntRaw, &b_T_coin_pHGCER_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_adcPulseTime", &T_coin_pHGCER_adcPulseTime, &b_T_coin_pHGCER_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pHGCER_adcPulseTimeRaw", &T_coin_pHGCER_adcPulseTimeRaw, &b_T_coin_pHGCER_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pHGCER_tdcMultiplicity", &T_coin_pHGCER_tdcMultiplicity, &b_T_coin_pHGCER_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pHGCER_tdcTime", &T_coin_pHGCER_tdcTime, &b_T_coin_pHGCER_tdcTime);
   iChain->SetBranchAddress("T.coin.pHGCER_tdcTimeRaw", &T_coin_pHGCER_tdcTimeRaw, &b_T_coin_pHGCER_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pHODO_RF_tdcMultiplicity", &T_coin_pHODO_RF_tdcMultiplicity, &b_T_coin_pHODO_RF_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pHODO_RF_tdcTime", &T_coin_pHODO_RF_tdcTime, &b_T_coin_pHODO_RF_tdcTime);
   iChain->SetBranchAddress("T.coin.pHODO_RF_tdcTimeRaw", &T_coin_pHODO_RF_tdcTimeRaw, &b_T_coin_pHODO_RF_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcMultiplicity", &T_coin_pNGCER_MOD_adcMultiplicity, &b_T_coin_pNGCER_MOD_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcPed", &T_coin_pNGCER_MOD_adcPed, &b_T_coin_pNGCER_MOD_adcPed);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcPedRaw", &T_coin_pNGCER_MOD_adcPedRaw, &b_T_coin_pNGCER_MOD_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcPulseAmp", &T_coin_pNGCER_MOD_adcPulseAmp, &b_T_coin_pNGCER_MOD_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcPulseAmpRaw", &T_coin_pNGCER_MOD_adcPulseAmpRaw, &b_T_coin_pNGCER_MOD_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcPulseInt", &T_coin_pNGCER_MOD_adcPulseInt, &b_T_coin_pNGCER_MOD_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcPulseIntRaw", &T_coin_pNGCER_MOD_adcPulseIntRaw, &b_T_coin_pNGCER_MOD_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcPulseTime", &T_coin_pNGCER_MOD_adcPulseTime, &b_T_coin_pNGCER_MOD_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pNGCER_MOD_adcPulseTimeRaw", &T_coin_pNGCER_MOD_adcPulseTimeRaw, &b_T_coin_pNGCER_MOD_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_adcMultiplicity", &T_coin_pNGCER_adcMultiplicity, &b_T_coin_pNGCER_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pNGCER_adcPed", &T_coin_pNGCER_adcPed, &b_T_coin_pNGCER_adcPed);
   iChain->SetBranchAddress("T.coin.pNGCER_adcPedRaw", &T_coin_pNGCER_adcPedRaw, &b_T_coin_pNGCER_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_adcPulseAmp", &T_coin_pNGCER_adcPulseAmp, &b_T_coin_pNGCER_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pNGCER_adcPulseAmpRaw", &T_coin_pNGCER_adcPulseAmpRaw, &b_T_coin_pNGCER_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_adcPulseInt", &T_coin_pNGCER_adcPulseInt, &b_T_coin_pNGCER_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pNGCER_adcPulseIntRaw", &T_coin_pNGCER_adcPulseIntRaw, &b_T_coin_pNGCER_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_adcPulseTime", &T_coin_pNGCER_adcPulseTime, &b_T_coin_pNGCER_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pNGCER_adcPulseTimeRaw", &T_coin_pNGCER_adcPulseTimeRaw, &b_T_coin_pNGCER_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pNGCER_tdcMultiplicity", &T_coin_pNGCER_tdcMultiplicity, &b_T_coin_pNGCER_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pNGCER_tdcTime", &T_coin_pNGCER_tdcTime, &b_T_coin_pNGCER_tdcTime);
   iChain->SetBranchAddress("T.coin.pNGCER_tdcTimeRaw", &T_coin_pNGCER_tdcTimeRaw, &b_T_coin_pNGCER_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRE100_ROC1_tdcMultiplicity", &T_coin_pPRE100_ROC1_tdcMultiplicity, &b_T_coin_pPRE100_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRE100_ROC1_tdcTime", &T_coin_pPRE100_ROC1_tdcTime, &b_T_coin_pPRE100_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRE100_ROC1_tdcTimeRaw", &T_coin_pPRE100_ROC1_tdcTimeRaw, &b_T_coin_pPRE100_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRE100_ROC2_tdcMultiplicity", &T_coin_pPRE100_ROC2_tdcMultiplicity, &b_T_coin_pPRE100_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRE100_ROC2_tdcTime", &T_coin_pPRE100_ROC2_tdcTime, &b_T_coin_pPRE100_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRE100_ROC2_tdcTimeRaw", &T_coin_pPRE100_ROC2_tdcTimeRaw, &b_T_coin_pPRE100_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRE150_ROC1_tdcMultiplicity", &T_coin_pPRE150_ROC1_tdcMultiplicity, &b_T_coin_pPRE150_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRE150_ROC1_tdcTime", &T_coin_pPRE150_ROC1_tdcTime, &b_T_coin_pPRE150_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRE150_ROC1_tdcTimeRaw", &T_coin_pPRE150_ROC1_tdcTimeRaw, &b_T_coin_pPRE150_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRE150_ROC2_tdcMultiplicity", &T_coin_pPRE150_ROC2_tdcMultiplicity, &b_T_coin_pPRE150_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRE150_ROC2_tdcTime", &T_coin_pPRE150_ROC2_tdcTime, &b_T_coin_pPRE150_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRE150_ROC2_tdcTimeRaw", &T_coin_pPRE150_ROC2_tdcTimeRaw, &b_T_coin_pPRE150_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRE200_ROC1_tdcMultiplicity", &T_coin_pPRE200_ROC1_tdcMultiplicity, &b_T_coin_pPRE200_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRE200_ROC1_tdcTime", &T_coin_pPRE200_ROC1_tdcTime, &b_T_coin_pPRE200_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRE200_ROC1_tdcTimeRaw", &T_coin_pPRE200_ROC1_tdcTimeRaw, &b_T_coin_pPRE200_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRE200_ROC2_tdcMultiplicity", &T_coin_pPRE200_ROC2_tdcMultiplicity, &b_T_coin_pPRE200_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRE200_ROC2_tdcTime", &T_coin_pPRE200_ROC2_tdcTime, &b_T_coin_pPRE200_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRE200_ROC2_tdcTimeRaw", &T_coin_pPRE200_ROC2_tdcTimeRaw, &b_T_coin_pPRE200_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRE40_ROC1_tdcMultiplicity", &T_coin_pPRE40_ROC1_tdcMultiplicity, &b_T_coin_pPRE40_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRE40_ROC1_tdcTime", &T_coin_pPRE40_ROC1_tdcTime, &b_T_coin_pPRE40_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRE40_ROC1_tdcTimeRaw", &T_coin_pPRE40_ROC1_tdcTimeRaw, &b_T_coin_pPRE40_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRE40_ROC2_tdcMultiplicity", &T_coin_pPRE40_ROC2_tdcMultiplicity, &b_T_coin_pPRE40_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRE40_ROC2_tdcTime", &T_coin_pPRE40_ROC2_tdcTime, &b_T_coin_pPRE40_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRE40_ROC2_tdcTimeRaw", &T_coin_pPRE40_ROC2_tdcTimeRaw, &b_T_coin_pPRE40_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRHI_tdcMultiplicity", &T_coin_pPRHI_tdcMultiplicity, &b_T_coin_pPRHI_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRHI_tdcTime", &T_coin_pPRHI_tdcTime, &b_T_coin_pPRHI_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRHI_tdcTimeRaw", &T_coin_pPRHI_tdcTimeRaw, &b_T_coin_pPRHI_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPRLO_tdcMultiplicity", &T_coin_pPRLO_tdcMultiplicity, &b_T_coin_pPRLO_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPRLO_tdcTime", &T_coin_pPRLO_tdcTime, &b_T_coin_pPRLO_tdcTime);
   iChain->SetBranchAddress("T.coin.pPRLO_tdcTimeRaw", &T_coin_pPRLO_tdcTimeRaw, &b_T_coin_pPRLO_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcMultiplicity", &T_coin_pPSHWR_adcMultiplicity, &b_T_coin_pPSHWR_adcMultiplicity);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcPed", &T_coin_pPSHWR_adcPed, &b_T_coin_pPSHWR_adcPed);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcPedRaw", &T_coin_pPSHWR_adcPedRaw, &b_T_coin_pPSHWR_adcPedRaw);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcPulseAmp", &T_coin_pPSHWR_adcPulseAmp, &b_T_coin_pPSHWR_adcPulseAmp);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcPulseAmpRaw", &T_coin_pPSHWR_adcPulseAmpRaw, &b_T_coin_pPSHWR_adcPulseAmpRaw);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcPulseInt", &T_coin_pPSHWR_adcPulseInt, &b_T_coin_pPSHWR_adcPulseInt);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcPulseIntRaw", &T_coin_pPSHWR_adcPulseIntRaw, &b_T_coin_pPSHWR_adcPulseIntRaw);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcPulseTime", &T_coin_pPSHWR_adcPulseTime, &b_T_coin_pPSHWR_adcPulseTime);
   iChain->SetBranchAddress("T.coin.pPSHWR_adcPulseTimeRaw", &T_coin_pPSHWR_adcPulseTimeRaw, &b_T_coin_pPSHWR_adcPulseTimeRaw);
   iChain->SetBranchAddress("T.coin.pRF_tdcMultiplicity", &T_coin_pRF_tdcMultiplicity, &b_T_coin_pRF_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pRF_tdcTime", &T_coin_pRF_tdcTime, &b_T_coin_pRF_tdcTime);
   iChain->SetBranchAddress("T.coin.pRF_tdcTimeRaw", &T_coin_pRF_tdcTimeRaw, &b_T_coin_pRF_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pSTOF_ROC1_tdcMultiplicity", &T_coin_pSTOF_ROC1_tdcMultiplicity, &b_T_coin_pSTOF_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pSTOF_ROC1_tdcTime", &T_coin_pSTOF_ROC1_tdcTime, &b_T_coin_pSTOF_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pSTOF_ROC1_tdcTimeRaw", &T_coin_pSTOF_ROC1_tdcTimeRaw, &b_T_coin_pSTOF_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pSTOF_ROC2_tdcMultiplicity", &T_coin_pSTOF_ROC2_tdcMultiplicity, &b_T_coin_pSTOF_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pSTOF_ROC2_tdcTime", &T_coin_pSTOF_ROC2_tdcTime, &b_T_coin_pSTOF_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pSTOF_ROC2_tdcTimeRaw", &T_coin_pSTOF_ROC2_tdcTimeRaw, &b_T_coin_pSTOF_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pT1_tdcMultiplicity", &T_coin_pT1_tdcMultiplicity, &b_T_coin_pT1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pT1_tdcTime", &T_coin_pT1_tdcTime, &b_T_coin_pT1_tdcTime);
   iChain->SetBranchAddress("T.coin.pT1_tdcTimeRaw", &T_coin_pT1_tdcTimeRaw, &b_T_coin_pT1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pT2_tdcMultiplicity", &T_coin_pT2_tdcMultiplicity, &b_T_coin_pT2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pT2_tdcTime", &T_coin_pT2_tdcTime, &b_T_coin_pT2_tdcTime);
   iChain->SetBranchAddress("T.coin.pT2_tdcTimeRaw", &T_coin_pT2_tdcTimeRaw, &b_T_coin_pT2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pT3_tdcMultiplicity", &T_coin_pT3_tdcMultiplicity, &b_T_coin_pT3_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pT3_tdcTime", &T_coin_pT3_tdcTime, &b_T_coin_pT3_tdcTime);
   iChain->SetBranchAddress("T.coin.pT3_tdcTimeRaw", &T_coin_pT3_tdcTimeRaw, &b_T_coin_pT3_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcMultiplicity", &T_coin_pTRIG1_ROC1_tdcMultiplicity, &b_T_coin_pTRIG1_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTime", &T_coin_pTRIG1_ROC1_tdcTime, &b_T_coin_pTRIG1_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTimeRaw", &T_coin_pTRIG1_ROC1_tdcTimeRaw, &b_T_coin_pTRIG1_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcMultiplicity", &T_coin_pTRIG1_ROC2_tdcMultiplicity, &b_T_coin_pTRIG1_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTime", &T_coin_pTRIG1_ROC2_tdcTime, &b_T_coin_pTRIG1_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw", &T_coin_pTRIG1_ROC2_tdcTimeRaw, &b_T_coin_pTRIG1_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG2_ROC1_tdcMultiplicity", &T_coin_pTRIG2_ROC1_tdcMultiplicity, &b_T_coin_pTRIG2_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG2_ROC1_tdcTime", &T_coin_pTRIG2_ROC1_tdcTime, &b_T_coin_pTRIG2_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG2_ROC1_tdcTimeRaw", &T_coin_pTRIG2_ROC1_tdcTimeRaw, &b_T_coin_pTRIG2_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG2_ROC2_tdcMultiplicity", &T_coin_pTRIG2_ROC2_tdcMultiplicity, &b_T_coin_pTRIG2_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG2_ROC2_tdcTime", &T_coin_pTRIG2_ROC2_tdcTime, &b_T_coin_pTRIG2_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG2_ROC2_tdcTimeRaw", &T_coin_pTRIG2_ROC2_tdcTimeRaw, &b_T_coin_pTRIG2_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG3_ROC1_tdcMultiplicity", &T_coin_pTRIG3_ROC1_tdcMultiplicity, &b_T_coin_pTRIG3_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG3_ROC1_tdcTime", &T_coin_pTRIG3_ROC1_tdcTime, &b_T_coin_pTRIG3_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG3_ROC1_tdcTimeRaw", &T_coin_pTRIG3_ROC1_tdcTimeRaw, &b_T_coin_pTRIG3_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG3_ROC2_tdcMultiplicity", &T_coin_pTRIG3_ROC2_tdcMultiplicity, &b_T_coin_pTRIG3_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG3_ROC2_tdcTime", &T_coin_pTRIG3_ROC2_tdcTime, &b_T_coin_pTRIG3_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG3_ROC2_tdcTimeRaw", &T_coin_pTRIG3_ROC2_tdcTimeRaw, &b_T_coin_pTRIG3_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcMultiplicity", &T_coin_pTRIG4_ROC1_tdcMultiplicity, &b_T_coin_pTRIG4_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTime", &T_coin_pTRIG4_ROC1_tdcTime, &b_T_coin_pTRIG4_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTimeRaw", &T_coin_pTRIG4_ROC1_tdcTimeRaw, &b_T_coin_pTRIG4_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcMultiplicity", &T_coin_pTRIG4_ROC2_tdcMultiplicity, &b_T_coin_pTRIG4_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTime", &T_coin_pTRIG4_ROC2_tdcTime, &b_T_coin_pTRIG4_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw", &T_coin_pTRIG4_ROC2_tdcTimeRaw, &b_T_coin_pTRIG4_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG5_ROC1_tdcMultiplicity", &T_coin_pTRIG5_ROC1_tdcMultiplicity, &b_T_coin_pTRIG5_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG5_ROC1_tdcTime", &T_coin_pTRIG5_ROC1_tdcTime, &b_T_coin_pTRIG5_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG5_ROC1_tdcTimeRaw", &T_coin_pTRIG5_ROC1_tdcTimeRaw, &b_T_coin_pTRIG5_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcMultiplicity", &T_coin_pTRIG5_ROC2_tdcMultiplicity, &b_T_coin_pTRIG5_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcTime", &T_coin_pTRIG5_ROC2_tdcTime, &b_T_coin_pTRIG5_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcTimeRaw", &T_coin_pTRIG5_ROC2_tdcTimeRaw, &b_T_coin_pTRIG5_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcMultiplicity", &T_coin_pTRIG6_ROC1_tdcMultiplicity, &b_T_coin_pTRIG6_ROC1_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTime", &T_coin_pTRIG6_ROC1_tdcTime, &b_T_coin_pTRIG6_ROC1_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTimeRaw", &T_coin_pTRIG6_ROC1_tdcTimeRaw, &b_T_coin_pTRIG6_ROC1_tdcTimeRaw);
   iChain->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcMultiplicity", &T_coin_pTRIG6_ROC2_tdcMultiplicity, &b_T_coin_pTRIG6_ROC2_tdcMultiplicity);
   iChain->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTime", &T_coin_pTRIG6_ROC2_tdcTime, &b_T_coin_pTRIG6_ROC2_tdcTime);
   iChain->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTimeRaw", &T_coin_pTRIG6_ROC2_tdcTimeRaw, &b_T_coin_pTRIG6_ROC2_tdcTimeRaw);
   iChain->SetBranchAddress("g.evtyp", &g_evtyp, &b_g_evtyp);
   iChain->SetBranchAddress("IBC3H00CRCUR4", &IBC3H00CRCUR4, &b_IBC3H00CRCUR4);
   iChain->SetBranchAddress("IGL1I00OD16_16", &IGL1I00OD16_16, &b_IGL1I00OD16_16);
   iChain->SetBranchAddress("HELFREQ", &HELFREQ, &b_HELFREQ);
   iChain->SetBranchAddress("HELDELAYs", &HELDELAYs, &b_HELDELAYs);
   iChain->SetBranchAddress("hac_bcm_average", &hac_bcm_average, &b_hac_bcm_average);
   iChain->SetBranchAddress("ibcm1", &ibcm1, &b_ibcm1);
   iChain->SetBranchAddress("ibcm2", &ibcm2, &b_ibcm2);
   iChain->SetBranchAddress("iunser", &iunser, &b_iunser);
   iChain->SetBranchAddress("itov3out", &itov3out, &b_itov3out);
   iChain->SetBranchAddress("itov4out", &itov4out, &b_itov4out);
   iChain->SetBranchAddress("IPM3H07A.XPOS", &IPM3H07A_XPOS, &b_IPM3H07A_XPOS);
   iChain->SetBranchAddress("IPM3H07A.YPOS", &IPM3H07A_YPOS, &b_IPM3H07A_YPOS);
   iChain->SetBranchAddress("IPM3H07B.XPOS", &IPM3H07B_XPOS, &b_IPM3H07B_XPOS);
   iChain->SetBranchAddress("IPM3H07B.YPOS", &IPM3H07B_YPOS, &b_IPM3H07B_YPOS);
   iChain->SetBranchAddress("IPM3H07C.XPOS", &IPM3H07C_XPOS, &b_IPM3H07C_XPOS);
   iChain->SetBranchAddress("IPM3H07C.YPOS", &IPM3H07C_YPOS, &b_IPM3H07C_YPOS);
   iChain->SetBranchAddress("IPM3H07A.XRAW", &IPM3H07A_XRAW, &b_IPM3H07A_XRAW);
   iChain->SetBranchAddress("IPM3H07A.YRAW", &IPM3H07A_YRAW, &b_IPM3H07A_YRAW);
   iChain->SetBranchAddress("IPM3H07B.XRAW", &IPM3H07B_XRAW, &b_IPM3H07B_XRAW);
   iChain->SetBranchAddress("IPM3H07B.YRAW", &IPM3H07B_YRAW, &b_IPM3H07B_YRAW);
   iChain->SetBranchAddress("IPM3H07C.XRAW", &IPM3H07C_XRAW, &b_IPM3H07C_XRAW);
   iChain->SetBranchAddress("IPM3H07C.YRAW", &IPM3H07C_YRAW, &b_IPM3H07C_YRAW);
   iChain->SetBranchAddress("IPM3H07A.XSOF", &IPM3H07A_XSOF, &b_IPM3H07A_XSOF);
   iChain->SetBranchAddress("IPM3H07A.YSOF", &IPM3H07A_YSOF, &b_IPM3H07A_YSOF);
   iChain->SetBranchAddress("IPM3H07B.XSOF", &IPM3H07B_XSOF, &b_IPM3H07B_XSOF);
   iChain->SetBranchAddress("IPM3H07B.YSOF", &IPM3H07B_YSOF, &b_IPM3H07B_YSOF);
   iChain->SetBranchAddress("IPM3H07C.XSOF", &IPM3H07C_XSOF, &b_IPM3H07C_XSOF);
   iChain->SetBranchAddress("IPM3H07C.YSOF", &IPM3H07C_YSOF, &b_IPM3H07C_YSOF);
   iChain->SetBranchAddress("fEvtHdr.fEvtTime", &fEvtHdr_fEvtTime, &b_Event_Branch_fEvtHdr_fEvtTime);
   iChain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtHdr_fEvtNum, &b_Event_Branch_fEvtHdr_fEvtNum);
   iChain->SetBranchAddress("fEvtHdr.fEvtType", &fEvtHdr_fEvtType, &b_Event_Branch_fEvtHdr_fEvtType);
   iChain->SetBranchAddress("fEvtHdr.fEvtLen", &fEvtHdr_fEvtLen, &b_Event_Branch_fEvtHdr_fEvtLen);
   iChain->SetBranchAddress("fEvtHdr.fHelicity", &fEvtHdr_fHelicity, &b_Event_Branch_fEvtHdr_fHelicity);
   iChain->SetBranchAddress("fEvtHdr.fTargetPol", &fEvtHdr_fTargetPol, &b_Event_Branch_fEvtHdr_fTargetPol);
   iChain->SetBranchAddress("fEvtHdr.fRun", &fEvtHdr_fRun, &b_Event_Branch_fEvtHdr_fRun);
  

  if (process.compare("elasticLT")==0 || process.compare("pi0LT")==0 || process.compare("vcsLT")==0 
      || process.compare("pi+LT")==0 || process.compare("K+LT")==0){
   iChain->SetBranchAddress("P.kin.secondary.Erecoil", &P_kin_secondary_Erecoil, &b_P_kin_secondary_Erecoil);
   iChain->SetBranchAddress("P.kin.secondary.MandelS", &P_kin_secondary_MandelS, &b_P_kin_secondary_MandelS);
   iChain->SetBranchAddress("P.kin.secondary.MandelT", &P_kin_secondary_MandelT, &b_P_kin_secondary_MandelT);
   iChain->SetBranchAddress("P.kin.secondary.MandelU", &P_kin_secondary_MandelU, &b_P_kin_secondary_MandelU);
   iChain->SetBranchAddress("P.kin.secondary.Mrecoil", &P_kin_secondary_Mrecoil, &b_P_kin_secondary_Mrecoil);
   iChain->SetBranchAddress("P.kin.secondary.Prec_x", &P_kin_secondary_Prec_x, &b_P_kin_secondary_Prec_x);
   iChain->SetBranchAddress("P.kin.secondary.Prec_y", &P_kin_secondary_Prec_y, &b_P_kin_secondary_Prec_y);
   iChain->SetBranchAddress("P.kin.secondary.Prec_z", &P_kin_secondary_Prec_z, &b_P_kin_secondary_Prec_z);
   iChain->SetBranchAddress("P.kin.secondary.emiss", &P_kin_secondary_emiss, &b_P_kin_secondary_emiss);
   iChain->SetBranchAddress("P.kin.secondary.emiss_nuc", &P_kin_secondary_emiss_nuc, &b_P_kin_secondary_emiss_nuc);
   iChain->SetBranchAddress("P.kin.secondary.ph_bq", &P_kin_secondary_ph_bq, &b_P_kin_secondary_ph_bq);
   iChain->SetBranchAddress("P.kin.secondary.ph_xq", &P_kin_secondary_ph_xq, &b_P_kin_secondary_ph_xq);
   iChain->SetBranchAddress("P.kin.secondary.phb_cm", &P_kin_secondary_phb_cm, &b_P_kin_secondary_phb_cm);
   iChain->SetBranchAddress("P.kin.secondary.phx_cm", &P_kin_secondary_phx_cm, &b_P_kin_secondary_phx_cm);
   iChain->SetBranchAddress("P.kin.secondary.pmiss", &P_kin_secondary_pmiss, &b_P_kin_secondary_pmiss);
   iChain->SetBranchAddress("P.kin.secondary.pmiss_x", &P_kin_secondary_pmiss_x, &b_P_kin_secondary_pmiss_x);
   iChain->SetBranchAddress("P.kin.secondary.pmiss_y", &P_kin_secondary_pmiss_y, &b_P_kin_secondary_pmiss_y);
   iChain->SetBranchAddress("P.kin.secondary.pmiss_z", &P_kin_secondary_pmiss_z, &b_P_kin_secondary_pmiss_z);
   iChain->SetBranchAddress("P.kin.secondary.px_cm", &P_kin_secondary_px_cm, &b_P_kin_secondary_px_cm);
   iChain->SetBranchAddress("P.kin.secondary.t_tot_cm", &P_kin_secondary_t_tot_cm, &b_P_kin_secondary_t_tot_cm);
   iChain->SetBranchAddress("P.kin.secondary.tb", &P_kin_secondary_tb, &b_P_kin_secondary_tb);
   iChain->SetBranchAddress("P.kin.secondary.tb_cm", &P_kin_secondary_tb_cm, &b_P_kin_secondary_tb_cm);
   iChain->SetBranchAddress("P.kin.secondary.th_bq", &P_kin_secondary_th_bq, &b_P_kin_secondary_th_bq);
   iChain->SetBranchAddress("P.kin.secondary.th_xq", &P_kin_secondary_th_xq, &b_P_kin_secondary_th_xq);
   iChain->SetBranchAddress("P.kin.secondary.thb_cm", &P_kin_secondary_thb_cm, &b_P_kin_secondary_thb_cm);
   iChain->SetBranchAddress("P.kin.secondary.thx_cm", &P_kin_secondary_thx_cm, &b_P_kin_secondary_thx_cm);
   iChain->SetBranchAddress("P.kin.secondary.tx", &P_kin_secondary_tx, &b_P_kin_secondary_tx);
   iChain->SetBranchAddress("P.kin.secondary.tx_cm", &P_kin_secondary_tx_cm, &b_P_kin_secondary_tx_cm);
   iChain->SetBranchAddress("P.kin.secondary.xangle", &P_kin_secondary_xangle, &b_P_kin_secondary_xangle);
   iChain->SetBranchAddress("H.kin.primary.Q2", &H_kin_primary_Q2, &b_H_kin_primary_Q2);
   iChain->SetBranchAddress("H.kin.primary.W", &H_kin_primary_W, &b_H_kin_primary_W);
   iChain->SetBranchAddress("H.kin.primary.W2", &H_kin_primary_W2, &b_H_kin_primary_W2);
   iChain->SetBranchAddress("H.kin.primary.epsilon", &H_kin_primary_epsilon, &b_H_kin_primary_epsilon);
   iChain->SetBranchAddress("H.kin.primary.nu", &H_kin_primary_nu, &b_H_kin_primary_nu);
   iChain->SetBranchAddress("H.kin.primary.omega", &H_kin_primary_omega, &b_H_kin_primary_omega);
   iChain->SetBranchAddress("H.kin.primary.ph_q", &H_kin_primary_ph_q, &b_H_kin_primary_ph_q);
   iChain->SetBranchAddress("H.kin.primary.q3m", &H_kin_primary_q3m, &b_H_kin_primary_q3m);
   iChain->SetBranchAddress("H.kin.primary.q_x", &H_kin_primary_q_x, &b_H_kin_primary_q_x);
   iChain->SetBranchAddress("H.kin.primary.q_y", &H_kin_primary_q_y, &b_H_kin_primary_q_y);
   iChain->SetBranchAddress("H.kin.primary.q_z", &H_kin_primary_q_z, &b_H_kin_primary_q_z);
   iChain->SetBranchAddress("H.kin.primary.scat_ang_deg", &H_kin_primary_scat_ang_deg, &b_H_kin_primary_scat_ang_deg);
   iChain->SetBranchAddress("H.kin.primary.scat_ang_rad", &H_kin_primary_scat_ang_rad, &b_H_kin_primary_scat_ang_rad);
   iChain->SetBranchAddress("H.kin.primary.th_q", &H_kin_primary_th_q, &b_H_kin_primary_th_q);
   iChain->SetBranchAddress("H.kin.primary.x_bj", &H_kin_primary_x_bj, &b_H_kin_primary_x_bj);
   iChain->SetBranchAddress("H.rb.raster.fr_xbpm_tar", &H_rb_raster_fr_xbpm_tar, &b_H_rb_raster_fr_xbpm_tar);
   iChain->SetBranchAddress("H.rb.raster.fr_ybpm_tar", &H_rb_raster_fr_ybpm_tar, &b_H_rb_raster_fr_ybpm_tar);
   iChain->SetBranchAddress("H.react.x", &H_react_x, &b_H_react_x);
   iChain->SetBranchAddress("H.react.y", &H_react_y, &b_H_react_y);
   iChain->SetBranchAddress("H.react.z", &H_react_z, &b_H_react_z);
  }
   // if process !!!
  if (process.compare("elastic")==0 || process.compare("pi0")==0 || process.compare("vcs")==0 
      || process.compare("pi+")==0 || process.compare("K+")==0){
   iChain->SetBranchAddress("H.kin.secondary.Erecoil", &H_kin_secondary_Erecoil, &b_H_kin_secondary_Erecoil);
   iChain->SetBranchAddress("H.kin.secondary.MandelS", &H_kin_secondary_MandelS, &b_H_kin_secondary_MandelS);
   iChain->SetBranchAddress("H.kin.secondary.MandelT", &H_kin_secondary_MandelT, &b_H_kin_secondary_MandelT);
   iChain->SetBranchAddress("H.kin.secondary.MandelU", &H_kin_secondary_MandelU, &b_H_kin_secondary_MandelU);
   iChain->SetBranchAddress("H.kin.secondary.Mrecoil", &H_kin_secondary_Mrecoil, &b_H_kin_secondary_Mrecoil);
   iChain->SetBranchAddress("H.kin.secondary.Prec_x", &H_kin_secondary_Prec_x, &b_H_kin_secondary_Prec_x);
   iChain->SetBranchAddress("H.kin.secondary.Prec_y", &H_kin_secondary_Prec_y, &b_H_kin_secondary_Prec_y);
   iChain->SetBranchAddress("H.kin.secondary.Prec_z", &H_kin_secondary_Prec_z, &b_H_kin_secondary_Prec_z);
   iChain->SetBranchAddress("H.kin.secondary.emiss", &H_kin_secondary_emiss, &b_H_kin_secondary_emiss);
   iChain->SetBranchAddress("H.kin.secondary.emiss_nuc", &H_kin_secondary_emiss_nuc, &b_H_kin_secondary_emiss_nuc);
   iChain->SetBranchAddress("H.kin.secondary.ph_bq", &H_kin_secondary_ph_bq, &b_H_kin_secondary_ph_bq);
   iChain->SetBranchAddress("H.kin.secondary.ph_xq", &H_kin_secondary_ph_xq, &b_H_kin_secondary_ph_xq);
   iChain->SetBranchAddress("H.kin.secondary.phb_cm", &H_kin_secondary_phb_cm, &b_H_kin_secondary_phb_cm);
   iChain->SetBranchAddress("H.kin.secondary.phx_cm", &H_kin_secondary_phx_cm, &b_H_kin_secondary_phx_cm);
   iChain->SetBranchAddress("H.kin.secondary.pmiss", &H_kin_secondary_pmiss, &b_H_kin_secondary_pmiss);
   iChain->SetBranchAddress("H.kin.secondary.pmiss_x", &H_kin_secondary_pmiss_x, &b_H_kin_secondary_pmiss_x);
   iChain->SetBranchAddress("H.kin.secondary.pmiss_y", &H_kin_secondary_pmiss_y, &b_H_kin_secondary_pmiss_y);
   iChain->SetBranchAddress("H.kin.secondary.pmiss_z", &H_kin_secondary_pmiss_z, &b_H_kin_secondary_pmiss_z);
   iChain->SetBranchAddress("H.kin.secondary.px_cm", &H_kin_secondary_px_cm, &b_H_kin_secondary_px_cm);
   iChain->SetBranchAddress("H.kin.secondary.t_tot_cm", &H_kin_secondary_t_tot_cm, &b_H_kin_secondary_t_tot_cm);
   iChain->SetBranchAddress("H.kin.secondary.tb", &H_kin_secondary_tb, &b_H_kin_secondary_tb);
   iChain->SetBranchAddress("H.kin.secondary.tb_cm", &H_kin_secondary_tb_cm, &b_H_kin_secondary_tb_cm);
   iChain->SetBranchAddress("H.kin.secondary.th_bq", &H_kin_secondary_th_bq, &b_H_kin_secondary_th_bq);
   iChain->SetBranchAddress("H.kin.secondary.th_xq", &H_kin_secondary_th_xq, &b_H_kin_secondary_th_xq);
   iChain->SetBranchAddress("H.kin.secondary.thb_cm", &H_kin_secondary_thb_cm, &b_H_kin_secondary_thb_cm);
   iChain->SetBranchAddress("H.kin.secondary.thx_cm", &H_kin_secondary_thx_cm, &b_H_kin_secondary_thx_cm);
   iChain->SetBranchAddress("H.kin.secondary.tx", &H_kin_secondary_tx, &b_H_kin_secondary_tx);
   iChain->SetBranchAddress("H.kin.secondary.tx_cm", &H_kin_secondary_tx_cm, &b_H_kin_secondary_tx_cm);
   iChain->SetBranchAddress("H.kin.secondary.xangle", &H_kin_secondary_xangle, &b_H_kin_secondary_xangle);
   iChain->SetBranchAddress("P.kin.primary.Q2", &P_kin_primary_Q2, &b_P_kin_primary_Q2);
   iChain->SetBranchAddress("P.kin.primary.W", &P_kin_primary_W, &b_P_kin_primary_W);
   iChain->SetBranchAddress("P.kin.primary.W2", &P_kin_primary_W2, &b_P_kin_primary_W2);
   iChain->SetBranchAddress("P.kin.primary.epsilon", &P_kin_primary_epsilon, &b_P_kin_primary_epsilon);
   iChain->SetBranchAddress("P.kin.primary.nu", &P_kin_primary_nu, &b_P_kin_primary_nu);
   iChain->SetBranchAddress("P.kin.primary.omega", &P_kin_primary_omega, &b_P_kin_primary_omega);
   iChain->SetBranchAddress("P.kin.primary.ph_q", &P_kin_primary_ph_q, &b_P_kin_primary_ph_q);
   iChain->SetBranchAddress("P.kin.primary.q3m", &P_kin_primary_q3m, &b_P_kin_primary_q3m);
   iChain->SetBranchAddress("P.kin.primary.q_x", &P_kin_primary_q_x, &b_P_kin_primary_q_x);
   iChain->SetBranchAddress("P.kin.primary.q_y", &P_kin_primary_q_y, &b_P_kin_primary_q_y);
   iChain->SetBranchAddress("P.kin.primary.q_z", &P_kin_primary_q_z, &b_P_kin_primary_q_z);
   iChain->SetBranchAddress("P.kin.primary.scat_ang_deg", &P_kin_primary_scat_ang_deg, &b_P_kin_primary_scat_ang_deg);
   iChain->SetBranchAddress("P.kin.primary.scat_ang_rad", &P_kin_primary_scat_ang_rad, &b_P_kin_primary_scat_ang_rad);
   iChain->SetBranchAddress("P.kin.primary.th_q", &P_kin_primary_th_q, &b_P_kin_primary_th_q);
   iChain->SetBranchAddress("P.kin.primary.x_bj", &P_kin_primary_x_bj, &b_P_kin_primary_x_bj);
  }
   
   Notify();
  return 1;
}




