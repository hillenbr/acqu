#include "TA2MesonPhysics.h"
int GLOBAL_COUNT = 0;
ClassImp(TA2MesonPhysics)

//-----------------------------------------------------------------------------

TA2MesonPhysics::TA2MesonPhysics(const char* Name, TA2Analysis* Analysis) : TA2BasePhysics(Name, Analysis)
{

  fVeto			= NULL; // TAPS Vetos
  fPID			= NULL; // PID

  PID_Hits                = new Int_t[24];
  Veto_Hits               = new Int_t[437];

  coveredCrystals.push_back(23);     // BaF2 which could be behind the Pizza Detector
  coveredCrystals.push_back(29);
  coveredCrystals.push_back(35);
  coveredCrystals.push_back(36);
  coveredCrystals.push_back(41);
  coveredCrystals.push_back(42);
  coveredCrystals.push_back(43);
  coveredCrystals.push_back(44);
  coveredCrystals.push_back(50);
  coveredCrystals.push_back(51);
  coveredCrystals.push_back(52);
  coveredCrystals.push_back(53);
  coveredCrystals.push_back(61);
  coveredCrystals.push_back(62);
  coveredCrystals.push_back(63);
  coveredCrystals.push_back(86);
  coveredCrystals.push_back(87);
  coveredCrystals.push_back(88);
  coveredCrystals.push_back(89);
  coveredCrystals.push_back(90);
  coveredCrystals.push_back(91);
  coveredCrystals.push_back(92);
  coveredCrystals.push_back(93);
  coveredCrystals.push_back(94);
  coveredCrystals.push_back(95);
  coveredCrystals.push_back(97);
  coveredCrystals.push_back(98);
  coveredCrystals.push_back(99);
  coveredCrystals.push_back(100);
  coveredCrystals.push_back(103);
  coveredCrystals.push_back(104);
  coveredCrystals.push_back(110);
  coveredCrystals.push_back(158);
  coveredCrystals.push_back(161);
  coveredCrystals.push_back(162);
  fcoveredCrystals.push_back(36);    // BaF2 which should be completely (at least 80-90%) behind the Pizza Detector
  fcoveredCrystals.push_back(43);
  fcoveredCrystals.push_back(44);
  fcoveredCrystals.push_back(90);
  fcoveredCrystals.push_back(92);
  fcoveredCrystals.push_back(93);
  fcoveredCrystals.push_back(97);
  dCrystals.push_back(310);          // BaF2 which are definitely not behind the Pizza detector
  dCrystals.push_back(311);
  dCrystals.push_back(312);
  dCrystals.push_back(313);
  dCrystals.push_back(314);
  dCrystals.push_back(315);
  dCrystals.push_back(316);
  dCrystals.push_back(317);
  dCrystals.push_back(318);
  dCrystals.push_back(319);
  dCrystals.push_back(320);
  dCrystals.push_back(400);
  dCrystals.push_back(401);
  dCrystals.push_back(402);
  dCrystals.push_back(403);
  dCrystals.push_back(404);
  dCrystals.push_back(405);
  dCrystals.push_back(406);
  dCrystals.push_back(407);
  dCrystals.push_back(408);
}

//-----------------------------------------------------------------------------

TA2MesonPhysics::~TA2MesonPhysics()
{

}

//---------------------------------------------------------------------------

void TA2MesonPhysics::SetConfig(Char_t* line, Int_t key)
{
  TA2BasePhysics::SetConfig(line, key);
}

//---------------------------------------------------------------------------

void TA2MesonPhysics::LoadVariable()
{
  TA2BasePhysics::LoadVariable();
}

//---------------------------------------------------------------------------

void TA2MesonPhysics::PostInit()
{

    // PID
    fPID = (TA2PlasticPID*)((TA2Analysis*)fParent)->GetGrandChild("PID");
    if ( !fPID) PrintError( "", "<No NaI class found>", EErrFatal);

    // Veto
    fVeto = (TA2PlasticPID*)((TA2Analysis*)fParent)->GetGrandChild("Veto");
    if ( !fVeto) PrintError( "", "<No NaI class found>", EErrFatal);

  //Call default PostInit()
  TA2BasePhysics::PostInit();

  hist[0]                       = new TH1I("ADC_PMTlinks", "ADC_PMTlinks", 4100, 0, 4100);
  hist[1]                       = new TH1I("ADC_PMTrechts", "ADC_PMTrechts", 4100, 0, 4100);
  hist[2]                       = new TH1I("ADC_beidePMT", "ADC_beidePMT", 8200, 0, 8200);
  hist[3]                       = new TH1I("TDC_links", "TDC_links", 15000, -10000, 5000);
  hist[4]                       = new TH1I("TDC_rechts", "TDC_rechts", 15000, -10000, 5000);
  hist[5]                       = new TH1I("TDC_beide", "TDC_beide", 15000, -10000, 5000);
  hist[6]                       = new TH1I("ADC_Proton_in_TAPS", "ADC_Proton_in_TAPS", 10000, 0, 10000);
  ADC_moeglicheKristalle	= new TH1I("ADC_BaF2area", "ADC_BaF2area", 10000, 0, 10000);
  ADC_PizzaHitBaF2area          = new TH1I("ADC_PizzaHitBaF2area", "ADC_PizzaHitBaF2area", 10000, 0, 10000);
  VetoADC_BaF2area	        = new TH1I("VetoADC_BaF2area", "VetoADC_BaF2area", 100, 0, 10);
  VetoADC_BaF2full	        = new TH1I("VetoADC_BaF2full", "VetoADC_BaF2full", 100, 0, 10);
  ADC_CutCrystals               = new TH1I("ADC_BaF2area_Cut", "ADC_BaF2area_Cut", 10000, 0, 10000);
  adcarea                       = new TH1I("ADC_BaF2area_Proton", "ADC_BaF2area_Proton", 10000, 0, 10000);
  adcarea1                      = new TH1I("ADC_BaF2area_Photon", "ADC_BaF2area_Photon", 10000, 0, 10000);
  adcarea2                      = new TH1I("ADC_BaF2area_PiPLus", "ADC_BaF2area_PiPLus", 10000, 0, 10000);
  adcfull                       = new TH1I("ADC_BaF2full_Proton", "ADC_BaF2full_Proton", 10000, 0, 10000);
  adcdiff                       = new TH1I("ADC_BaF2diff_Proton", "ADC_BaF2diff_Proton", 10000, 0, 10000);
  EnergyVeto                    = new TH1D("VetoE", "VetoE", 1000, 0, 10);
  dEvE_OR                       = new TH2D("dEvE_Veto", "dEvE_Veto", 700, 0, 700, 100, 0, 10);
  dEvE_Pizza                    = new TH2D("dEvE_Pizza", "dEvE_Pizza", 700, 0, 700, 410, 0, 8200);
  dEvE_Pizza1                   = new TH2D("dEvE_PMTlinks", "dEvE_PMTlinks", 700, 0, 700, 210, 0, 4200);
  dEvE_Pizza2                   = new TH2D("dEvE_PMTrechts", "dEvE_PMTrechts", 700, 0, 700, 210, 0, 4200);
  dEvE_Veto29                   = new TH2D("dEvE_Veto29", "dEvE_Veto29", 700, 0, 700, 300, 0, 1500);
  VetoCut_dEvE_Pizza            = new TH2D("VetoCut_dEvE_Pizza", "VetoCut_dEvE_Pizza", 700, 0, 700, 410, 0, 8200);
  VetoCut_dEvE_Veto             = new TH2D("VetoCut_dEvE_Veto", "VetoCut_dEvE_Veto", 700, 0, 700, 100, 0, 10);
  PizzaCut_dEvE_Pizza           = new TH2D("PizzaCut_dEvE_Pizza", "PizzaCut_dEvE_Pizza", 700, 0, 700, 410, 0, 8200);
  PizzaCut_dEvE_Veto            = new TH2D("PizzaCut_dEvE_Veto", "PizzaCut_dEvE_Veto", 700, 0, 700, 100, 0, 10);
  dEvE_BaF2area                 = new TH2D("dEvE_BaF2area_Proton", "dEvE_BaF2area_Proton", 900, 0, 900, 410, 0, 8200);
  dEvE_BaF2diff                 = new TH2D("dEvE_BaF2diff_Proton", "dEvE_BaF2diff_Proton", 900, 0, 900, 410, 0, 8200);
  dEvE_BaF2full                 = new TH2D("dEvE_BaF2full_Proton", "dEvE_BaF2full_Proton", 900, 0, 900, 410, 0, 8200);
  dEvE_adcarea1                 = new TH2D("dEvE_BaF2area_Photon", "dEvE_BaF2area_Photon", 900, 0, 900, 410, 0, 8200);
  dEvE_adcarea2                 = new TH2D("dEvE_BaF2area_PiPlus", "dEvE_BaF2area_PiPlus", 900, 0, 900, 410, 0, 8200);
  precuthist                    = new TH1D("ADC_PMT_BaF2_Cluster", "ADC_PMT_BaF2_Cluster", 8200, 0, 8200);
  cuthist1                      = new TH1D("ADC_PMT_BaF2_Cluster_PizzaCut", "ADC_PMT_BaF2_Cluster_PizzaCut", 8200, 0, 8200);
  cuthist                       = new TH1D("ADC_PMT_BaF2_Cluster_VetoCut", "ADC_PMT_BaF2_Cluster_VetoCut", 8200, 0, 8200);
  ADCvE_PizzaHitBaF2area        = new TH2D("ADCvE_BaF2area", "ADCvE_BaF2area", 800, 0, 800, 410, 0, 8200);
  ADCvE_ProtonPizzaCutBaF2area  = new TH2D("ADCvE_BaF2area_PizzaCut", "ADCvE_BaF2area_PizzaCut", 800, 0, 800, 410, 0, 8200);
  VetovE_ProtPizzaCutBaF2area   = new TH2D("VetoADCvE_BaF2Area", "VetoADCvE_BaF2Area", 800, 0, 800, 100, 0, 10);
  VetovE_ProtPizzaCutBaF2full   = new TH2D("VetoADCvE_BaF2Full", "VetoADCvE_BaF2Full", 800, 0, 800, 100, 0, 10);
  // Count Histograms
  count                         = new TH1I("Counts","Counts",10,0,10);
  Treffer                       = new TH1I("VetoHits_when_PizzaHit","VetoHits_when_PizzaHit",438,0,438);
  TrefferProt                   = new TH1I("VetoHits_when_PizzaProton_and_BaF2_Cluster","VetoHits_when_PizzaProton_and_BaF2_Cluster",438,0,438);
  TrefferProt1                  = new TH1I("VetoHits_when_PizzaHit_and_BaF2_Cluster","VetoHits_when_PizzaHit_and_BaF2_Cluster",438,0,438);
  TrefferProtVeto               = new TH1I("VetoHits_when_VetoProton","VetoHits_when_VetoProton",438,0,438);

  //Call master default PostInit()
  TA2Physics::PostInit();
}

//-----------------------------------------------------------------------------

void TA2MesonPhysics::Reconstruct()
{
 VarInit();

 //Perform basic physics tasks
 TA2BasePhysics::Reconstruct();

 // Get Detector Hits
 nVeto_Hits      = fVeto->GetNhits();
 nPID_Hits       = fPID->GetNhits();
 Bool_t pizzahit = false;
 for(int i=0; i<nPID_Hits; i++) {
     PID_Hits[i] = fPID->GetHits(i);
     if(PID_Hits[i] == 22 || PID_Hits[i] == 23) pizzahit=true;
 }
 //skip other stuff if no event is detected in my pizza slice

 //if(!pizzahit) return;

 // Introducing: "The Energies" and "The Clusters"

 Double_t* vetoE     = TAPS->GetVeto()->GetEnergyOR();
 Double_t* baFlE     = BaF2->GetEnergyOR();
 Double_t* ClE       = BaF2->GetClEnergyOR();
 UInt_t* testclust   = BaF2->GetClustHit();
 Bool_t pizzahitprot = false;
 HitCluster_t* tmpclust;

 // ADC Spectra from Pizza Detector

        hist[0]->Fill((gAR->GetADC())[124]);
        hist[1]->Fill((gAR->GetADC())[125]);
        hist[2]->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);

     // Protons in TAPS, nProtonMarc from Basephysics
     if(nProtonMarc > 0)
        hist[6]->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);


 // TDC Spectra from Pizza Detector

        for(int i=0; i<gAR->GetMulti(2002)->GetNstore(); i++)
                hist[3]->Fill(gAR->GetMulti(2002)->GetHit(i));
        for(int i=0; i<gAR->GetMulti(2003)->GetNstore(); i++)
                hist[4]->Fill(gAR->GetMulti(2003)->GetHit(i));
        count->Fill(1);

        if(gAR->GetMulti(2002)->GetNstore() == 1)
        {
                if(gAR->GetMulti(2003)->GetNstore() == 1)
                {
                        count->Fill(8);
                        hist[5]->Fill(gAR->GetMulti(2002)->GetHit(0) - gAR->GetMulti(2003)->GetHit(0));
                }
        }

 // Loop over all Cluster in TAPS

 for(int i=0; i<BaF2->GetNCluster(); i++)
        {
            tmpclust = BaF2->GetCluster(i);
            dEvE_OR->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
            EnergyVeto->Fill(vetoE[tmpclust->GetIndex()]);
            dEvE_Pizza->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
            dEvE_Pizza1->Fill(ClE[i],(gAR->GetADC())[124]);
            dEvE_Pizza2->Fill(ClE[i],(gAR->GetADC())[125]);
            dEvE_Veto29->Fill(ClE[i],(gAR->GetADC())[26075]);
            precuthist->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);

 // Loop BaF2 Cluster - Only Cluster in BaF2, which are possibly behind the Pizza Detector

            if(std::find(coveredCrystals.begin(), coveredCrystals.end(), testclust[i]) != coveredCrystals.end()){
            VetovE_ProtPizzaCutBaF2area->Fill(ClE[i],TAPS->GetVeto()->GetEnergy(testclust[i]));
            VetoADC_BaF2area->Fill(TAPS->GetVeto()->GetEnergy(testclust[i]));
            ADC_moeglicheKristalle->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
            }

 // Loop BaF2 Cluster - Only Cluster in BaF2, which are definitely behind the Pizza Detector

            if(std::find(fcoveredCrystals.begin(), fcoveredCrystals.end(), testclust[i]) != fcoveredCrystals.end()){
            VetovE_ProtPizzaCutBaF2full->Fill(ClE[i],TAPS->GetVeto()->GetEnergy(testclust[i]));
            VetoADC_BaF2full->Fill(TAPS->GetVeto()->GetEnergy(testclust[i]));
            }

 // Loop BaF2 Cluster - Hit in Pizza Detector

            if(pizzahit){
                for(int j=0; j<nVeto_Hits; j++)
                    TrefferProt1->Fill(fVeto->GetHits(j));

 // Loop BaF2 Cluster - Hit Pizza Detector - Only Cluster in BaF2, which are possibly behind the Pizza Detector

                if(std::find(coveredCrystals.begin(), coveredCrystals.end(), testclust[i]) != coveredCrystals.end()){
                    ADCvE_PizzaHitBaF2area->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    ADC_PizzaHitBaF2area->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                }
            }

 // Loop BaF2 Cluster - Simple Cuts on dEvE Veto-TAPS and Pizza-TAPS

            // Cut on Veto Energy
            if((2.9*(TMath::Exp(-0.0105*ClE[i]))+2.1) < vetoE[tmpclust->GetIndex()])
            {
                VetoCut_dEvE_Veto->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
                VetoCut_dEvE_Pizza->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                cuthist->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
            }
            else
            {
            }
            // Cut on Pizza ADC
            if((2200*(TMath::Exp(-0.008524*ClE[i]))+2300) < (gAR->GetADC())[124] + (gAR->GetADC())[125])
            {
                PizzaCut_dEvE_Pizza->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                PizzaCut_dEvE_Veto->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
                cuthist1->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                // Hit in Pizza Detector
                if(pizzahit){
                    for(int j=0; j<nVeto_Hits; j++)
                        TrefferProt->Fill(fVeto->GetHits(j));
                    // Cluster in BaF2, which is possibly covered by the Pizza Detector
                    if(std::find(coveredCrystals.begin(), coveredCrystals.end(), testclust[i]) != coveredCrystals.end()){
                        ADCvE_ProtonPizzaCutBaF2area->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                        ADC_CutCrystals->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    }
                }
            }
            else
            {
            }
        }

 // End Loop
 //

 // Counts Veto Hits
 for(int i=0; i<nVeto_Hits; i++) {
               // Proton Hits in Veto
               if(nProtonMarc > 0){
                   TrefferProtVeto->Fill(fVeto->GetHits(i));
               }
           }
 // Counts Veto Hits, when Pizza Detector was hit
 if(pizzahit){
    for(int i=0; i<nVeto_Hits; i++) {
        Treffer->Fill(fVeto->GetHits(i));
    }
}

 // Counts Protons in BaF2 regions ; Protons from BasePhysics
 for(int i = 0; i < nProton; i++){
            if(std::find(dCrystals.begin(), dCrystals.end(), Proton[i].GetCentralIndex()) != dCrystals.end())
                count->Fill(2);
            if(std::find(coveredCrystals.begin(), coveredCrystals.end(), Proton[i].GetCentralIndex()) != coveredCrystals.end())
                count->Fill(4);
            if(std::find(fcoveredCrystals.begin(), fcoveredCrystals.end(), Proton[i].GetCentralIndex()) != fcoveredCrystals.end())
                count->Fill(6);
        }

 // ADC Spectra Pizza Detector for different Particles from BasePhysics with Loop over different BaF2 regions

        bool alreadyFilled = false;
        // Proton in TAPS
        if(nProtonMarc > 0)
            for(int i = 0; i < nProtonMarc; i++)
                // Only BaF2, which could possibly be behind the Pizza Detector
                if(std::find(coveredCrystals.begin(), coveredCrystals.end(), ProtonMarc[i].GetCentralIndex()) != coveredCrystals.end()){
                     if(!alreadyFilled){
                        adcarea->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                        dEvE_BaF2area->Fill(ProtonMarc[i].GetP4().E()-ProtonMarc[i].GetP4().M(),(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                        alreadyFilled = true;
                    }
                     count->Fill(5);
                }

        bool alreadyFilled1 = false;
        // Proton in TAPS
        if(nProtonMarc1 > 0)
            for(int i = 0; i < nProtonMarc1; i++)
                // Only BaF2, which are definitely not behind the Pizza Detector
                if(std::find(dCrystals.begin(), dCrystals.end(), ProtonMarc1[i].GetCentralIndex()) != dCrystals.end()){
                    if(!alreadyFilled1){
                    adcdiff->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    dEvE_BaF2diff->Fill(ProtonMarc1[i].GetP4().E()-ProtonMarc1[i].GetP4().M(),(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    alreadyFilled1 = true;
                }
                 count->Fill(3);
                }

        bool alreadyFilled2 = false;
        // Proton in TAPS
        if(nProtonMarc2 > 0)
            for(int i = 0; i < nProtonMarc2; i++)
                // Only BaF2, which are definitely (at least 80-90%)  behind the Pizza Detector
                if(std::find(fcoveredCrystals.begin(), fcoveredCrystals.end(), ProtonMarc2[i].GetCentralIndex()) != fcoveredCrystals.end()){
                    if(!alreadyFilled2){
                    adcfull->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    dEvE_BaF2full->Fill(ProtonMarc2[i].GetP4().E()-ProtonMarc2[i].GetP4().M(),(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    alreadyFilled2 = true;
                }
                 count->Fill(9);
                }


        bool alreadyFilled3 = false;
        // Photon in TAPS
        if(nPhotonMarc > 0)
            for(int i = 0; i < nPhotonMarc; i++)
                // Only BaF2, which could possibly be behind the Pizza Detector
                if(std::find(coveredCrystals.begin(), coveredCrystals.end(), PhotonMarc[i].GetCentralIndex()) != coveredCrystals.end()){
                     if(!alreadyFilled3){
                        adcarea1->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                        dEvE_adcarea1->Fill(PhotonMarc[i].GetP4().E()-PhotonMarc[i].GetP4().M(),(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                        alreadyFilled3 = true;
                    }
                }

        bool alreadyFilled4 = false;
        // PiPLus in TAPS
        if(nPiPlusMarc > 0)
            for(int i = 0; i < nPiPlusMarc; i++)
                // Only BaF2, which could possibly be behind the Pizza Detector
                if(std::find(coveredCrystals.begin(), coveredCrystals.end(), PiPlusMarc[i].GetCentralIndex()) != coveredCrystals.end()){
                     if(!alreadyFilled4){
                        adcarea2->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                        dEvE_adcarea2->Fill(PiPlusMarc[i].GetP4().E()-PiPlusMarc[i].GetP4().M(),(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                        alreadyFilled4 = true;
                    }
                }

  //Clean up and set EBufferEnd end markers for arrays
  TermArrays();
}

//-----------------------------------------------------------------------------

void TA2MesonPhysics::VarInit()
{
  //Initally, set single value histograms to EBufferEnd (this value will not be plotted)
 
}

//-----------------------------------------------------------------------------

void TA2MesonPhysics::TermArrays()
{
  //Set EBufferEnd end markers on multi value histograms
  
}

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

void TA2MesonPhysics::ParseMisc(char* line)
{

  //Get keyword
  /*if(sscanf(line, "%s", sWord)!=1) return;

  if(!strcmp("Threshold", sWord))
  {
    sscanf(line, "%*s %lf", &Threshold);
    printf("Reaction threshold is %f\n", Threshold);
    return;
  }

  if(!strcmp("InvariantMass", sWord))
  {
    sscanf(line, "%*s %lf %lf", &MinvLo, &MinvHi);
    printf("Invariant mass cut window from %f to %f\n", MinvLo, MinvHi);
    return;
  }*/


  TA2BasePhysics::ParseMisc(line);
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
