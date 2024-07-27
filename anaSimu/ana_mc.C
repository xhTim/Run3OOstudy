/*
Main script for running the analysis on MC.

	1. Apply the cuts on the MC sample
	2. Create the histograms for various kinematic variables

*/

#include "inFiles/headers/ParticleTree.C"
#include "../common/functions.C"
#include "../common/EvtSel.C"
#include "../common/matching.C"

//* Ka GEN level
TH3D *hPtEtaPhi_KaPos_gen	= new TH3D("hPtEtaPhi_KaPos_gen", "hPtEtaPhi_KaPos_gen; p_{T} [GeV/c]; #eta; #phi", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig, mHistPhiBins, mHistPhiLow, mHistPhiHig);
TH3D *hPtEtaPhi_KaPos_acc	= new TH3D("hPtEtaPhi_KaPos_acc", "hPtEtaPhi_KaPos_acc; p_{T} [GeV/c]; #eta; #phi", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig, mHistPhiBins, mHistPhiLow, mHistPhiHig);
TH3D *hPtEtaPhi_KaNeg_gen	= new TH3D("hPtEtaPhi_KaNeg_gen", "hPtEtaPhi_KaNeg_gen; p_{T} [GeV/c]; #eta; #phi", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig, mHistPhiBins, mHistPhiLow, mHistPhiHig);
TH3D *hPtEtaPhi_KaNeg_acc	= new TH3D("hPtEtaPhi_KaNeg_acc", "hPtEtaPhi_KaNeg_acc; p_{T} [GeV/c]; #eta; #phi", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig, mHistPhiBins, mHistPhiLow, mHistPhiHig);
TH3D *hPtEtaPhi_Ka_gen		= new TH3D("hPtEtaPhi_Ka_gen", "hPtEtaPhi_Ka_gen; p_{T} [GeV/c]; #eta; #phi", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig, mHistPhiBins, mHistPhiLow, mHistPhiHig);
TH3D *hPtEtaPhi_Ka_acc		= new TH3D("hPtEtaPhi_Ka_acc", "hPtEtaPhi_Ka_acc; p_{T} [GeV/c]; #eta; #phi", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig, mHistPhiBins, mHistPhiLow, mHistPhiHig);

//* DiKa GEN level
TH3D *hPtRapPhi_DiKa_gen	= new TH3D("hPtRapPhi_DiKa_gen", "hPtRapPhi_DiKa_gen; p_{T} [GeV/c]; y; #phi", mHistPtBins, mHistPtLow, mHistPtHig, mHistRapBins, mHistRapLow, mHistRapHig, mHistPhiBins, mHistPhiLow, mHistPhiHig);
TH3D *hPtRapMass_DiKa_gen	= new TH3D("hPtRapMass_DiKa_gen", "hPtRapMass_DiKa_gen; p_{T} [GeV/c]; y; M_{KK} [GeV/c^{2}]", mHistPtBins, mHistPtLow, mHistPtHig, mHistRapBins, mHistRapLow, mHistRapHig, mHistMassBins, mHistMassLow, mHistMassHig);
//* DiKa GEN level within acceptance
TH3D *hPtRapPhi_DiKa_acc 		= new TH3D("hPtRapPhi_DiKa_acc", "hPtRapPhi_DiKa_acc; p_{T} [GeV/c]; y; #phi", mHistPtBins, mHistPtLow, mHistPtHig, mHistRapBins, mHistRapLow, mHistRapHig, mHistPhiBins, mHistPhiLow, mHistPhiHig);
TH3D *hPtRapMass_DiKa_acc 		= new TH3D("hPtRapMass_DiKa_acc", "hPtRapMass_DiKa_acc; p_{T} [GeV/c]; y; M_{KK} [GeV/c^{2}]", mHistPtBins, mHistPtLow, mHistPtHig, mHistRapBins, mHistRapLow, mHistRapHig, mHistMassBins, mHistMassLow, mHistMassHig);

TH2D *hPtResDeltaR_GenReco = new TH2D("hPtResDeltaR_GenReco", "hPtResDeltaR_GenReco; (p_{T}^{reco} - p_{T}^{gen}) / p_{T}^{gen}; #DeltaR", 200, -1.0, 1.0, 200, 0, 1.0);
TH2D *hPtResDeltaR_RecoGen = new TH2D("hPtResDeltaR_RecoGen", "hPtResDeltaR_RecoGen; (p_{T}^{gen} - p_{T}^{reco}) / p_{T}^{reco}; #DeltaR", 200, -1.0, 1.0, 200, 0, 1.0);
TH3D *hAbsPtResDeltaRMass_pos = new TH3D("hAbsPtResDeltaRMass_pos", "hPtResDeltaRMass_pos; |p_{T}^{reco} - p_{T}^{gen}| / p_{T}^{gen}; #DeltaR; M_{KK} [GeV/c^{2}]", 200, 0, 1.0, 200, 0, 1.0, mHistMassBins, mHistMassLow, mHistMassHig);
TH3D *hAbsPtResDeltaRMass_neg = new TH3D("hAbsPtResDeltaRMass_neg", "hPtResDeltaRMass_neg; |p_{T}^{reco} - p_{T}^{gen}| / p_{T}^{gen}; #DeltaR; M_{KK} [GeV/c^{2}]", 200, 0, 1.0, 200, 0, 1.0, mHistMassBins, mHistMassLow, mHistMassHig);

TH2D *hPtResDeltaR_RecoGen1 = new TH2D("hPtResDeltaR_RecoGen1", "hPtResDeltaR_RecoGen1; (p_{T}^{gen} - p_{T}^{reco}) / p_{T}^{reco}; #DeltaR", 200, -1.0, 1.0, 200, 0, 1.0);
TH2D *hPtResDeltaR_RecoGen2 = new TH2D("hPtResDeltaR_RecoGen2", "hPtResDeltaR_RecoGen2; (p_{T}^{gen} - p_{T}^{reco}) / p_{T}^{reco}; #DeltaR", 200, -1.0, 1.0, 200, 0, 1.0);
TH2D *hPtResDeltaR_RecoGen3 = new TH2D("hPtResDeltaR_RecoGen3", "hPtResDeltaR_RecoGen3; (p_{T}^{gen} - p_{T}^{reco}) / p_{T}^{reco}; #DeltaR", 200, -1.0, 1.0, 200, 0, 1.0);
TH2D *hPtResDeltaR_RecoGen4 = new TH2D("hPtResDeltaR_RecoGen4", "hPtResDeltaR_RecoGen4; (p_{T}^{gen} - p_{T}^{reco}) / p_{T}^{reco}; #DeltaR", 200, -1.0, 1.0, 200, 0, 1.0);
TH2D *hPtResDeltaR_RecoGen5 = new TH2D("hPtResDeltaR_RecoGen5", "hPtResDeltaR_RecoGen5; (p_{T}^{gen} - p_{T}^{reco}) / p_{T}^{reco}; #DeltaR", 200, -1.0, 1.0, 200, 0, 1.0);

TH2D *hPtEta1 = new TH2D("hPtEta1", "hPtEta1; p_{T} [GeV/c]; #eta", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig);
TH2D *hPtEta2 = new TH2D("hPtEta2", "hPtEta2; p_{T} [GeV/c]; #eta", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig);
TH2D *hPtEta3 = new TH2D("hPtEta3", "hPtEta3; p_{T} [GeV/c]; #eta", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig);
TH2D *hPtEta4 = new TH2D("hPtEta4", "hPtEta4; p_{T} [GeV/c]; #eta", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig);
TH2D *hPtEta5 = new TH2D("hPtEta5", "hPtEta5; p_{T} [GeV/c]; #eta", mHistPtBins, mHistPtLow, mHistPtHig, mHistEtaBins, mHistEtaLow, mHistEtaHig);

TH3D *hAbsPtResDeltaRMass_trk_pos = new TH3D("hAbsPtResDeltaRMass_trk_pos", "hPtResDeltaRMass_trk_pos; |p_{T}^{reco} - p_{T}^{gen}| / p_{T}^{gen}; #DeltaR; M_{KK} [GeV/c^{2}]", 200, 0, 1.0, 200, 0, 1.0, mHistMassBins, mHistMassLow, mHistMassHig);
TH3D *hAbsPtResDeltaRMass_trk_neg = new TH3D("hAbsPtResDeltaRMass_trk_neg", "hPtResDeltaRMass_trk_neg; |p_{T}^{reco} - p_{T}^{gen}| / p_{T}^{gen}; #DeltaR; M_{KK} [GeV/c^{2}]", 200, 0, 1.0, 200, 0, 1.0, mHistMassBins, mHistMassLow, mHistMassHig);

TH3D *hMvsPtvsRap = new TH3D("hMvsPtvsRap", "hMvsPtvsRap; Rapidity; p_{T} [GeV/c]; M_{MuMu} [GeV/c^{2}]", mHistRapBins, mHistRapLow, mHistRapHig, mHistPtBins, mHistPtLow, mHistPtHig, mHistMassBins, mHistMassLow, mHistMassHig);

void ana_mc(TString FileName, int HFVeto_option = mHFVetoOpt)
{
	//= Read the input file ==================================================================

	// TString FileName = mInFileName_mc;
	// TString FileName = "diPi_ana_mc";

	TString inFileName = "inFiles/diMu_ana_mc_" + FileName + ".root";
	TString outFileName = "outFiles/ana_mc_" + FileName + ".root";

	TFile *inFile = new TFile(inFileName, "READ");
	TTree *inTree = (TTree*)inFile->Get("diMuAna/ParticleTree");
	ParticleTree *tree = new ParticleTree(inTree);

	//= Loop over the tree ====================================================================
	int nentries = inTree->GetEntries();
	int max_entries = 6000000;
	cout << "Total Entries: " << inTree->GetEntries() << endl;
	for (int ientry = 0; ientry < nentries; ientry++)
	{
		if (ientry > max_entries) break;

		update_progress(ientry, nentries, 10);
		tree->GetEntry(ientry);

		bool pass_acc = true;

		int ngen_pdgId = tree->gen_pdgId->size();
		//= GEN level ============================================================================
		for (int igen = 0; igen < ngen_pdgId; igen++)
		{
			int pdgId = tree->gen_pdgId->at(igen);
			int status = tree->gen_status->at(igen);
			int charge = tree->gen_charge->at(igen);
			double pT = tree->gen_pT->at(igen);
			double eta = tree->gen_eta->at(igen);
			double phi = tree->gen_phi->at(igen);
			double y = tree->gen_y->at(igen);
			double mass = tree->gen_mass->at(igen);

			if (isPhi(pdgId))
			{
				hPtRapPhi_DiKa_gen->Fill(pT, y, phi);
				hPtRapMass_DiKa_gen->Fill(pT, y, mass);

				//* find the daughter kaons
				auto dau_Idx = tree->gen_dauIdx->at(igen);
				if (dau_Idx.size() != 2)
				{
					cout << "Warning: the number of GEN daughters is not 2" << endl;
					continue;
				}

				//*
				pass_gen_match_reco_cand(tree, igen, hPtResDeltaR_GenReco);


				double pT_pos = tree->gen_charge->at(dau_Idx[0]) > 0 ? tree->gen_pT->at(dau_Idx[0]) : tree->gen_pT->at(dau_Idx[1]);
				double pT_neg = tree->gen_charge->at(dau_Idx[0]) < 0 ? tree->gen_pT->at(dau_Idx[0]) : tree->gen_pT->at(dau_Idx[1]);
				double eta_pos = tree->gen_charge->at(dau_Idx[0]) > 0 ? tree->gen_eta->at(dau_Idx[0]) : tree->gen_eta->at(dau_Idx[1]);
				double eta_neg = tree->gen_charge->at(dau_Idx[0]) < 0 ? tree->gen_eta->at(dau_Idx[0]) : tree->gen_eta->at(dau_Idx[1]);
				double phi_pos = tree->gen_charge->at(dau_Idx[0]) > 0 ? tree->gen_phi->at(dau_Idx[0]) : tree->gen_phi->at(dau_Idx[1]);
				double phi_neg = tree->gen_charge->at(dau_Idx[0]) < 0 ? tree->gen_phi->at(dau_Idx[0]) : tree->gen_phi->at(dau_Idx[1]);

				hPtEtaPhi_KaPos_gen->Fill(pT_pos, eta_pos, phi_pos);
				hPtEtaPhi_KaNeg_gen->Fill(pT_neg, eta_neg, phi_neg);

				bool pass_acc_pos = pass_acc_gen(pT_pos, eta_pos);
				bool pass_acc_neg = pass_acc_gen(pT_neg, eta_neg);
				if (pass_acc_pos) hPtEtaPhi_KaPos_acc->Fill(pT_pos, eta_pos, phi_pos);
				if (pass_acc_neg) hPtEtaPhi_KaNeg_acc->Fill(pT_neg, eta_neg, phi_neg);
				pass_acc = pass_acc_pos && pass_acc_neg;
				if (pass_acc)
				{
					hPtRapPhi_DiKa_acc->Fill(pT, y, phi);
					hPtRapMass_DiKa_acc->Fill(pT, y, mass);
				}
			}
		}

		//= Event Selection ======================================================================
		bool passEvtSel = pass_EvtSel_standard(tree, nullptr, mPVFilterIdx, HFVeto_option);
		if (!passEvtSel) continue;

		//= RECO level Phi ======================================================================
		const int ncand = tree->cand_pdgId->size();
		for (int icand = 0; icand < ncand; icand++)
		{
			//! only consider phi!!!!!!!!!!!!!!!!!
			int pdgId = tree->cand_pdgId->at(icand);
			//if (!isPhi(pdgId)) continue;
			if (!check_dau_cand(tree, icand)) continue;
			cout<<"if any? "<<icand<<endl;
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			reco_match_gen_fillMass_cand(tree, icand, hAbsPtResDeltaRMass_pos, hAbsPtResDeltaRMass_neg);

			//if (! pass_acc_gen_cand(tree, icand)) continue;
			// if (! pass_TrkQuality_cand(tree, icand)) continue;

			//* fill the candidate information
			//double mass = tree->cand_mass->at(icand);
			//double y = tree->cand_y->at(icand);
			//double pT = tree->cand_pT->at(icand);
			//hPtRapMass->Fill(pT, y, mass);

			pass_reco_match_gen_cand(tree, icand, hPtResDeltaR_RecoGen);
			reco_match_gen_fillMass_cand(tree, icand, hAbsPtResDeltaRMass_trk_pos, hAbsPtResDeltaRMass_trk_neg);

			auto daughters = tree->cand_dauIdx->at(icand);
			int dauIdx_pos = (tree->cand_charge->at(daughters[0]) > 0) ? daughters[0] : daughters[1];
			int dauIdx_neg = (tree->cand_charge->at(daughters[0]) < 0) ? daughters[0] : daughters[1];
			double pT_pos = tree->cand_pT->at(dauIdx_pos);
			double pT_neg = tree->cand_pT->at(dauIdx_neg);
			double eta_pos = tree->cand_eta->at(dauIdx_pos);
			double eta_neg = tree->cand_eta->at(dauIdx_neg);
			double phi_pos = tree->cand_phi->at(dauIdx_pos);
			double phi_neg = tree->cand_phi->at(dauIdx_neg);
			
			TLorentzVector posFourMom, negFourMom, pairFourMom;
			posFourMom.SetPtEtaPhiM(pT_pos, eta_pos, phi_pos, MUON_MASS);
			negFourMom.SetPtEtaPhiM(pT_neg, eta_neg, phi_neg, MUON_MASS);
			pairFourMom = posFourMom + negFourMom;

			Double_t pt   = pairFourMom.Pt();
			Double_t eta  = pairFourMom.Eta();
			Double_t phi  = pairFourMom.Phi();
			Double_t mass = pairFourMom.M();
			Double_t y    = pairFourMom.Rapidity();
			hMvsPtvsRap->Fill(y, pt, mass);

			if (abs(eta_pos) < 1.2)
			{
				reco_match_gen(tree, dauIdx_pos, hPtResDeltaR_RecoGen1);
			}
			if (1.2 < abs(eta_pos) && abs(eta_pos) < 1.8)
			{
				reco_match_gen(tree, dauIdx_pos, hPtResDeltaR_RecoGen2);
			}
			if (1.8 < abs(eta_pos) && abs(eta_pos) < 2.1)
			{
				reco_match_gen(tree, dauIdx_pos, hPtResDeltaR_RecoGen3);
			}
			if (2.1 < abs(eta_pos) && abs(eta_pos) < 2.4)
			{
				reco_match_gen(tree, dauIdx_pos, hPtResDeltaR_RecoGen4);
			}

			if (abs(eta_neg) < 1.2)
			{
				reco_match_gen(tree, dauIdx_neg, hPtResDeltaR_RecoGen1);
			}
			if (1.2 < abs(eta_neg) && abs(eta_neg) < 1.8)
			{
				reco_match_gen(tree, dauIdx_neg, hPtResDeltaR_RecoGen2);
			}
			if (1.8 < abs(eta_neg) && abs(eta_neg) < 2.1)
			{
				reco_match_gen(tree, dauIdx_neg, hPtResDeltaR_RecoGen3);
			}
			if (2.1 < abs(eta_neg) && abs(eta_neg) < 2.4)
			{
				reco_match_gen(tree, dauIdx_neg, hPtResDeltaR_RecoGen4);
			}


			if ( 0.5 < abs(y) < 1.0)
			{
				// file pt eta of the kaons
				hPtEta1->Fill(pT_pos, eta_pos);
				hPtEta1->Fill(pT_neg, eta_neg);
			}
			if ( 0.3 < abs(y) < 0.5)
			{
				// file pt eta of the kaons
				hPtEta2->Fill(pT_pos, eta_pos);
				hPtEta2->Fill(pT_neg, eta_neg);
			}
			if ( 0.2 < abs(y) < 0.3)
			{
				// file pt eta of the kaons
				hPtEta3->Fill(pT_pos, eta_pos);
				hPtEta3->Fill(pT_neg, eta_neg);
			}
			if ( 0.1 < abs(y) < 0.2)
			{
				// file pt eta of the kaons
				hPtEta4->Fill(pT_pos, eta_pos);
				hPtEta4->Fill(pT_neg, eta_neg);
			}
			if ( 0.0 < abs(y) < 0.1)
			{
				// file pt eta of the kaons
				hPtEta5->Fill(pT_pos, eta_pos);
				hPtEta5->Fill(pT_neg, eta_neg);
			}
		}
	}

	//* Combine the Kaon+ and Kaon- to get the Kaon
	hPtEtaPhi_Ka_gen->Add(hPtEtaPhi_KaPos_gen, hPtEtaPhi_KaNeg_gen);
	hPtEtaPhi_Ka_acc->Add(hPtEtaPhi_KaPos_acc, hPtEtaPhi_KaNeg_acc);

	//= Write the output file ==================================================================
	TFile *outFile = new TFile(outFileName, "RECREATE");
	outFile->cd();

	// write_sequence( hPtEtaPhi_KaPos_gen );
	// write_sequence( hPtEtaPhi_KaNeg_gen );
	// write_sequence( hPtEtaPhi_KaPos_acc );
	// write_sequence( hPtEtaPhi_KaNeg_acc );
	// write_sequence( hPtEtaPhi_Ka_gen, {"xy"}, 1);
	// write_sequence( hPtEtaPhi_Ka_acc, {"xy"}, 1);
	// write_sequence( hPtRapPhi_DiKa_gen, {"xy"}, 1);
	// write_sequence( hPtRapMass_DiKa_gen, {"xy"}, 1);
	// write_sequence( hPtRapPhi_DiKa_acc, {"xy"}, 1);
	// write_sequence( hPtRapMass_DiKa_acc, {"xy"}, 1);
	// write_sequence( hPtResDeltaR_GenReco);
	// write_sequence( hPtResDeltaR_RecoGen);
	// write_sequence( hAbsPtResDeltaRMass_pos, {"xy"}, 1);
	// write_sequence( hAbsPtResDeltaRMass_neg, {"xy"}, 1);
	// write_sequence( hAbsPtResDeltaRMass_trk_pos, {"xy"}, 1);
	// write_sequence( hAbsPtResDeltaRMass_trk_neg, {"xy"}, 1);
	write_sequence( hMvsPtvsRap, {}, 1);
	write_sequence( hPtResDeltaR_RecoGen1);
	write_sequence( hPtResDeltaR_RecoGen2);
	write_sequence( hPtResDeltaR_RecoGen3);
	write_sequence( hPtResDeltaR_RecoGen4);
	write_sequence( hPtResDeltaR_RecoGen5);
	write_sequence( hPtEta1, 1);
	write_sequence( hPtEta2, 1);
	write_sequence( hPtEta3, 1);
	write_sequence( hPtEta4, 1);
	write_sequence( hPtEta5, 1);

	outFile->Close();	
}
