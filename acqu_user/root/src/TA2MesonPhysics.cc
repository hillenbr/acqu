#include "TA2MesonPhysics.h"

ClassImp(TA2MesonPhysics)

//-----------------------------------------------------------------------------

TA2MesonPhysics::TA2MesonPhysics(const char* Name, TA2Analysis* Analysis) : TA2BasePhysics(Name, Analysis)
{

  fVeto			= NULL; // TAPS Vetos
  fPID			= NULL; // PID

  PID_Hits                = new Int_t[24];
  Veto_Hits               = new Int_t[437];

  hist[0]=0;
  hist[1]=0;
  hist[2]=0;
  hist[6]=0;
  //hist[3]=0;
  //hist[4]=0;
  //hist[5]=0;
  coveredCrystals.push_back(23);
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
  coveredCrystals.push_back(103);
  coveredCrystals.push_back(104);
  coveredCrystals.push_back(158);
  coveredCrystals.push_back(161);
  coveredCrystals.push_back(162);
  fcoveredCrystals.push_back(29);
  fcoveredCrystals.push_back(36);
  fcoveredCrystals.push_back(42);
  fcoveredCrystals.push_back(43);
  fcoveredCrystals.push_back(90);
  fcoveredCrystals.push_back(92);
  fcoveredCrystals.push_back(93);
  fcoveredCrystals.push_back(97);
  fcoveredCrystal.push_back(28);
  fcoveredCrystal.push_back(35);
  fcoveredCrystal.push_back(41);
  fcoveredCrystal.push_back(42);
  fcoveredCrystal.push_back(89);
  fcoveredCrystal.push_back(91);
  fcoveredCrystal.push_back(92);
  fcoveredCrystal.push_back(96);
  dCrystals.push_back(310);
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

    // NaI
    fPID = (TA2PlasticPID*)((TA2Analysis*)fParent)->GetGrandChild("PID");
    if ( !fPID) PrintError( "", "<No NaI class found>", EErrFatal);

    // NaI
    fVeto = (TA2PlasticPID*)((TA2Analysis*)fParent)->GetGrandChild("Veto");
    if ( !fVeto) PrintError( "", "<No NaI class found>", EErrFatal);

  //Call default PostInit()
  TA2BasePhysics::PostInit();

  hist[0]	= new TH1I("ADC_PMTlinks", "ADC_PMTlinks", 5000, 0, 5000);
  hist[1]	= new TH1I("ADC_PMTrechts", "ADC_PMTrechts", 5000, 0, 5000);
  hist[2]	= new TH1I("ADC_beidePMT", "ADC_beidePMT", 10000, 0, 10000);
  hist[6]	= new TH1I("ADC_Proton_in_TAPS", "ADC_Proton_in_TAPS", 10000, 0, 10000);
  adcarea	= new TH1I("ADC_moeglicheBaF2", "ADC_moeglicheBaF2", 10000, 0, 10000);
  adcfull	= new TH1I("ADC_definitiveBaF2", "ADC_definitiveBaF2", 10000, 0, 10000);
  adcfull1	= new TH1I("ADC_definitiveBaF2_minus1", "ADC_definitiveBaF2_minus1", 10000, 0, 10000);
  adcdiff	= new TH1I("ADC_unmoeglicheBaF2", "ADC_unmoeglicheBaF2", 10000, 0, 10000);
  hist[3]	= new TH1I("TDC_links", "TDC_links", 150000, -145000, 5000);
  hist[4]	= new TH1I("TDC_rechts", "TDC_rechts", 150000, -145000, 5000);
  hist[5]	= new TH1I("TDC_beide", "TDC_beide", 15000, -10000, 5000);
  //hist[6]	= new TH1I("CombinedTDCAdd", "CombinedTDCAdd", 15000, -10000, 5000);
  /*TAPS_adcarea	= new TH1I("TAPS_adcarea", "TAPS_adcarea", 10000, 0, 1000000000);
  TAPS_adcfull	= new TH1I("TAPS_adcfull", "TAPS_adcfull", 10000, 0, 1000000000);
  TAPS_adcfull1	= new TH1I("TAPS_adcfull1", "TAPS_adcfull1", 10000, 0, 1000000000);
  TAPS_adcdiff	= new TH1I("TAPS_adcdiff", "TAPS_adcdiff", 10000, 0, 1000000000);*/

        count	           = new TH1I("Zaehler","Zaehler",10,0,10);
        Treffer	           = new TH1I("Counts","Counts",438,0,438);
        TrefferProt	   = new TH1I("CountsProt","CountsProt",438,0,438);
        TrefferVeto        = new TH1I("CountsVeto","CountsVeto",438,0,438);
        TrefferProtVeto    = new TH1I("CountsProtVeto","CountsProtVeto",438,0,438);
        //TrefferPCProtVeto  = new TH1I("CountsPCProtVeto","CountsPCProtVeto",438,0,438);

  char    str[64];
    /*for(int i=0; i<BaF2->GetNelem(); i++)
    {
      sprintf(str,"dEvE_%d", i);
        dEvE[i]          = new TH2D(str, str, 1000, 0, 1000, 100, 0, 10);
    }*/
    EnergyVeto           = new TH1D("VetoE", "VetoE", 1000, 0, 10);
    dEvE_OR              = new TH2D("dEvE_Veto", "dEvE_Veto", 1000, 0, 1000, 100, 0, 10);
    //dEvE_OR1             = new TH2D("dEvE_OR1", "dEvE_OR1", 1000, 0, 1000, 8200, 0, 10);
    dEvE_Pizza           = new TH2D("dEvE_Pizza", "dEvE_Pizza", 1000, 0, 1000, 8200, 0, 8200);
    dEvE_Pizza1          = new TH2D("dEvE_PMTlinks", "dEvE_PMTlinks", 1000, 0, 1000, 4200, 0, 4200);
    dEvE_Pizza2          = new TH2D("dEvE_PMTrechts", "dEvE_PMTrechts", 1000, 0, 1000, 4200, 0, 4200);
    dEvE_Veto28          = new TH2D("dEvE_Veto28", "dEvE_Veto28", 1000, 0, 1000, 15000, 0, 1500);
    dEvE_Veto29          = new TH2D("dEvE_Veto29", "dEvE_Veto29", 1000, 0, 1000, 15000, 0, 1500);
    dEvE_Veto30          = new TH2D("dEvE_Veto30", "dEvE_Veto30", 1000, 0, 1000, 15000, 0, 1500);
    VetoCut_dEvE_Pizza   = new TH2D("VetoCut_dEvE_Pizza", "VetoCut_dEvE_Pizza", 1000, 0, 1000, 8200, 0, 8200);
    VetoCut_dEvE_Veto    = new TH2D("VetoCut_dEvE_Veto", "VetoCut_dEvE_Veto", 1000, 0, 1000, 100, 0, 10);
    PizzaCut_dEvE_Pizza  = new TH2D("PizzaCut_dEvE_Pizza", "PizzaCut_dEvE_Pizza", 1000, 0, 1000, 8200, 0, 8200);
    PizzaCut_dEvE_Veto   = new TH2D("PizzaCut_dEvE_Veto", "PizzaCut_dEvE_Veto", 1000, 0, 1000, 100, 0, 10);
    precuthist           = new TH1D("Precut", "Precut", 8200, 0, 8200);
    cuthist1             = new TH1D("PizzaCut", "PizzaCut", 8200, 0, 8200);
    cuthist              = new TH1D("VetoCut", "VetoCut", 8200, 0, 8200);
  //Call master default PostInit()
  TA2Physics::PostInit();
}

//-----------------------------------------------------------------------------

void TA2MesonPhysics::Reconstruct()
{
  //Perform basic physics tasks
 TA2BasePhysics::Reconstruct();


        Double_t* vetoE = TAPS->GetVeto()->GetEnergyOR();
        Double_t* baFlE = BaF2->GetEnergyOR();
        //Double_t* baFlE = BaF2->GetEnergyAll()

        Double_t* ClE = BaF2->GetClEnergyOR();
        Bool_t pizzahitprot = false;

        for(int i=0; i<BaF2->GetNCluster(); i++)
        {
            HitCluster_t* tmpclust = BaF2->GetCluster(i);

            //dEvE[i]->Fill(baFlE[i],vetoE[i]);
            //printf("dEvE_OR->Fill(%f,%f)\n", BaF2->GetCluster(i)->GetEnergy(), vetoE[(BaF2->GetCluster(i)->GetHits())[0]]);
            //unsigned int tmp = (BaF2->GetCluster(i)->GetHits())[0];
            //UInt_t* temp = BaF2->GetCluster(i)->GetHits();
            dEvE_OR->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
            EnergyVeto->Fill(vetoE[tmpclust->GetIndex()]);
            //dEvE_OR1->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
            dEvE_Pizza->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
            dEvE_Pizza1->Fill(ClE[i],(gAR->GetADC())[124]);
            dEvE_Pizza2->Fill(ClE[i],(gAR->GetADC())[125]);
            dEvE_Veto29->Fill(ClE[i],(gAR->GetADC())[26075]);
            dEvE_Veto28->Fill(ClE[i],(gAR->GetADC())[26074]);
            dEvE_Veto30->Fill(ClE[i],(gAR->GetADC())[26076]);
            precuthist->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);

            /*if(4.2-(0.01*ClE[i]) < vetoE[tmpclust->GetIndex()])
            {
                VetoCut_dEvE_Veto->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
                VetoCut_dEvE_Pizza->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                cuthist->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
            }
            else
            {
            }*/
            if((2.9*(TMath::Exp(-0.0105*ClE[i]))+2.1) < vetoE[tmpclust->GetIndex()])
            {
                VetoCut_dEvE_Veto->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
                VetoCut_dEvE_Pizza->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                cuthist->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
            }
            else
            {
            }
            if((2200*(TMath::Exp(-0.008523741*ClE[i]))+2300) < (gAR->GetADC())[124] + (gAR->GetADC())[125])
            {
                PizzaCut_dEvE_Pizza->Fill(ClE[i],(gAR->GetADC())[124] + (gAR->GetADC())[125]);
                PizzaCut_dEvE_Veto->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
                cuthist1->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                pizzahitprot = true;
            }
            else
            {
            }
        }


        nVeto_Hits = fVeto->GetNhits();

           for(int i=0; i<nVeto_Hits; i++) {
               TrefferVeto->Fill(fVeto->GetHits(i));
               if(nProtonMarc > 0){
                   TrefferProtVeto->Fill(fVeto->GetHits(i));
               }
           }


 // Get Detector Hits
 nPID_Hits = fPID->GetNhits();
 Bool_t pizzahit = false;
 for(int i=0; i<nPID_Hits; i++) {
     PID_Hits[i] = fPID->GetHits(i);
     if(PID_Hits[i] == 22 || PID_Hits[i] == 23) pizzahit=true;
 }

 nVeto_Hits = fVeto->GetNhits();

 /*for(int i=0; i<nVeto_Hits; i++) {
     if(pizzahitprot == true){
         TrefferPCProtVeto->Fill(fVeto->GetHits(i));
     }
 }*/

 if(pizzahit == true){
    for(int i=0; i<nVeto_Hits; i++) {
        Treffer->Fill(fVeto->GetHits(i));
        if(pizzahitprot == true){
            TrefferProt->Fill(fVeto->GetHits(i));
        }
    }
}
 //printf("PID: %d\t Veto: %d\n", nPID_Hits, nVeto_Hits);

	VarInit();

	//printf("maxADC: %d\t adc124: %d\n", gAR->GetMaxADC(), (gAR->GetADC())[124]);
	hist[0]->Fill((gAR->GetADC())[124]);
	hist[1]->Fill((gAR->GetADC())[125]);
        hist[2]->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
	//cout << "NStore: " << gAR->GetMulti(2002)->GetNstore() << endl;
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
                        //hist[6]->Fill(gAR->GetMulti(2002)->GetHit(0) + gAR->GetMulti(2003)->GetHit(0));
		}
	}

        for(int i = 0; i < nProton; i++){
            if(std::find(dCrystals.begin(), dCrystals.end(), Proton[i].GetCentralIndex()) != dCrystals.end())
                count->Fill(2);
            if(std::find(coveredCrystals.begin(), coveredCrystals.end(), Proton[i].GetCentralIndex()) != coveredCrystals.end())
                count->Fill(4);
            if(std::find(fcoveredCrystals.begin(), fcoveredCrystals.end(), Proton[i].GetCentralIndex()) != fcoveredCrystals.end())
                count->Fill(6);
        }

        if(nProtonMarc > 0)
            hist[6]->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);


        UInt_t TAPS_ADC_sum;
        if(nProtonMarc > 0){
           TAPS_ADC_sum = 0;
           for (std::vector<Int_t>::iterator it = coveredCrystals.begin() ; it != coveredCrystals.end(); ++it)
                TAPS_ADC_sum += (BaF2->GetADC())[*it];
          // TAPS_adcarea->Fill(TAPS_ADC_sum);
       }

        bool alreadyFilled = false;
        if(nProtonMarc > 0)
            for(int i = 0; i < nProtonMarc; i++)
                if(std::find(coveredCrystals.begin(), coveredCrystals.end(), ProtonMarc[i].GetCentralIndex()) != coveredCrystals.end()){
                     if(!alreadyFilled){
                        adcarea->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                        alreadyFilled = true;
                    }
                     count->Fill(5);
                }

        if(nProtonMarc > 0){
           TAPS_ADC_sum = 0;
           for (std::vector<Int_t>::iterator it = fcoveredCrystals.begin() ; it != fcoveredCrystals.end(); ++it)
                TAPS_ADC_sum += (BaF2->GetADC())[*it];
          // TAPS_adcfull->Fill(TAPS_ADC_sum);
       }

        alreadyFilled = false;
        if(nProtonMarc > 0)
            for(int i = 0; i < nProtonMarc; i++)
                if(std::find(fcoveredCrystals.begin(), fcoveredCrystals.end(), ProtonMarc[i].GetCentralIndex()) != fcoveredCrystals.end()){
                    if(!alreadyFilled){
                    adcfull->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    alreadyFilled = true;
                }
                 count->Fill(7);
                }

        if(nProtonMarc > 0){
           TAPS_ADC_sum = 0;
           for (std::vector<Int_t>::iterator it = dCrystals.begin(); it != dCrystals.end(); ++it)
                TAPS_ADC_sum += (BaF2->GetADC())[*it];
          // TAPS_adcdiff->Fill(TAPS_ADC_sum);
       }

        alreadyFilled = false;
        if(nProtonMarc > 0)
            for(int i = 0; i < nProtonMarc; i++)
                if(std::find(dCrystals.begin(), dCrystals.end(), ProtonMarc[i].GetCentralIndex()) != dCrystals.end()){
                    if(!alreadyFilled){
                    adcdiff->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    alreadyFilled = true;
                }
                 count->Fill(3);
                }

        if(nProtonMarc > 0){
           TAPS_ADC_sum = 0;
           for (std::vector<Int_t>::iterator it = fcoveredCrystal.begin() ; it != fcoveredCrystal.end(); ++it)
                TAPS_ADC_sum += (BaF2->GetADC())[*it];
          // TAPS_adcfull1->Fill(TAPS_ADC_sum);
       }

        alreadyFilled = false;
        if(nProtonMarc > 0)
            for(int i = 0; i < nProtonMarc; i++)
                if(std::find(fcoveredCrystal.begin(), fcoveredCrystal.end(), ProtonMarc[i].GetCentralIndex()) != fcoveredCrystal.end()){
                    if(!alreadyFilled){
                    adcfull1->Fill((gAR->GetADC())[124] + (gAR->GetADC())[125]);
                    alreadyFilled = true;
                }
                 count->Fill(9);
                }
/*        for(int i = 0; i < nProtonMarc; i++){
            ProtonMarc[i].GetCentralIndex()
        }*/
	

  
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
