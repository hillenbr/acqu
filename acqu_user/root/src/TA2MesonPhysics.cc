#include "TA2MesonPhysics.h"

ClassImp(TA2MesonPhysics)

//-----------------------------------------------------------------------------

TA2MesonPhysics::TA2MesonPhysics(const char* Name, TA2Analysis* Analysis) : TA2BasePhysics(Name, Analysis)
{
  hist[0]=0;
  hist[1]=0;
  hist[2]=0;
  //hist[3]=0;
  //hist[4]=0;
  //hist[5]=0;
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
  //Call default PostInit()
  TA2BasePhysics::PostInit();

  hist[0]	= new TH1D("MarcADC124", "MarcADC124", 10000, 0, 10000);
  hist[1]	= new TH1D("MarcADC125", "MarcADC125", 10000, 0, 10000);
  hist[2]	= new TH1D("CombinedADC", "CombinedADC", 10000, 0, 10000);
  hist[3]	= new TH1D("MarcTDC2002", "MarcTDC2002", 15000, -10000, 5000);
  hist[4]	= new TH1D("MarcTDC2003", "MarcTDC2003", 15000, -10000, 5000);
  hist[5]	= new TH1D("CombinedTDC", "CombinedTDC", 15000, -10000, 5000);
  hist[6]	= new TH1D("CombinedTDCAdd", "CombinedTDCAdd", 15000, -10000, 5000);

	count	= new TH1I("Countdzjudtj","Countdzjudtj",10,0,10);

  char    str[64];
    for(int i=0; i<BaF2->GetNelem(); i++)
    {
      sprintf(str,"dEvE_%d", i);
        dEvE[i]          = new TH2D(str, str, 1000, 0, 1000, 100, 0, 10);
    }
    dEvE_OR          = new TH2D("dEvE_OR", "dEvE_OR", 1000, 0, 1000, 100, 0, 10);
    Cut_dEvE_OR      = new TH2D("Cut_dEvE_OR", "Cut_dEvE_OR", 1000, 0, 1000, 100, 0, 10);
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

        for(int i=0; i<BaF2->GetNCluster(); i++)
        {
            HitCluster_t* tmpclust = BaF2->GetCluster(i);

            //dEvE[i]->Fill(baFlE[i],vetoE[i]);
            //printf("dEvE_OR->Fill(%f,%f)\n", BaF2->GetCluster(i)->GetEnergy(), vetoE[(BaF2->GetCluster(i)->GetHits())[0]]);
            //unsigned int tmp = (BaF2->GetCluster(i)->GetHits())[0];
            //UInt_t* temp = BaF2->GetCluster(i)->GetHits();
            dEvE_OR->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);

            if(4-(0.01*ClE[i]) < vetoE[tmpclust->GetIndex()])
            {
                //Cut_dEvE_OR->Fill(ClE[i],vetoE[tmpclust->GetIndex()]);
                hist[6]->Fill(gAR->GetMulti(2002)->GetHit(0) + gAR->GetMulti(2003)->GetHit(0));
            }
            else
            {
            }
        }


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
			count->Fill(2);
			hist[5]->Fill(gAR->GetMulti(2002)->GetHit(0) - gAR->GetMulti(2003)->GetHit(0));
			hist[6]->Fill(gAR->GetMulti(2002)->GetHit(0) + gAR->GetMulti(2003)->GetHit(0));
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
