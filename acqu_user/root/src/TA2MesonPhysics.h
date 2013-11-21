#ifndef __TA2MesonPhysics_h__
#define __TA2MesonPhysics_h__

#include "TA2BasePhysics.h"
#include "TA2KFParticle.h"
#include "TA2CBKinematicFitter.h"

#include <vector>
#include <algorithm>

#define MASS_PIZERO    134.9766
#define MASS_ETA       547.7500

class TA2Apparatus;

class TA2MesonPhysics : public TA2BasePhysics
{
private:
        //TA2TAPS_Veto*	veto;

  protected:
        TH1*	hist[7];
        TH1*	cuthist;
        TH1*	cuthist1;
        TH1*	adcarea;    // TAPS Crystals + ADC maybe covered
        TH1*	adcarea1;    // TAPS Crystals + ADC maybe covered
        TH1*	adcarea2;    // TAPS Crystals + ADC maybe covered
        TH1*	adcfull;    // TAPS Crystals + ADC fully covered
        TH1*	adcdiff;    // TAPS Crystals + ADC far away from module
        TH1*	precuthist;
        TH1*    ADC_PizzaHitBaF2area;
        TH1*    VetoADC_BaF2area;
        TH1*    VetoADC_BaF2full;
        TH1*	ADC_moeglicheKristalle;
        TH1*	Test;
        TH1*	count;
        TH1*    ADC_CutCrystals;
        TH1*	Treffer;
        TH1*	TrefferProt;
        TH1*	TrefferProt1;
        TH1*	TrefferProtVeto;
        TH1*	EnergyVeto;
        TH2D*   dEvE_OR;
        TH2D*   dEvE_Pizza;
        TH2D*   dEvE_Pizza1;
        TH2D*   dEvE_Pizza2;
        TH2D*   dEvE_Veto28;
        TH2D*   dEvE_Veto29;
        TH2D*   dEvE_Veto30;
        TH2D*   VetoCut_dEvE_Pizza;
        TH2D*   VetoCut_dEvE_Veto;
        TH2D*   PizzaCut_dEvE_Pizza;
        TH2D*   PizzaCut_dEvE_Veto;
        TH2D*   dEvE_BaF2area;
        TH2D*   dEvE_BaF2diff;
        TH2D*   dEvE_BaF2full;
        TH2D*   dEvE_adcarea1;
        TH2D*   dEvE_adcarea2;
        TH2D*   ADCvE_PizzaHitBaF2area;
        TH2D*   ADCvE_ProtonPizzaCutBaF2area;
        TH2D*   VetovE_ProtPizzaCutBaF2area;
        TH2D*   VetovE_ProtPizzaCutBaF2full;
        Int_t   nPID_Hits;
        Int_t   nVeto_Hits;
        Int_t*  PID_Hits;
        Int_t*  Veto_Hits;
        std::vector<Int_t> coveredCrystals;
        std::vector<Int_t> fcoveredCrystals;
        std::vector<Int_t> dCrystals;

    // Begin by initialising Detectors
    TA2PlasticPID	*fPID; 	        // PID
    TA2PlasticPID	*fVeto; 	// TAPS Vetos

    //General functions
    void VarInit();                    //Clear all used variables
    void TermArrays();                 //Terminate arrays with EBufferEnd markers
  
  public:
    TA2MesonPhysics(const char*, TA2Analysis*);
    virtual ~TA2MesonPhysics();
    virtual void LoadVariable();            //Creates histograms
    virtual void SetConfig(Char_t*, Int_t); //Parses general information from configuration file
    virtual void ParseMisc(char* line);     //Parses additional information from configuration file
    virtual void Reconstruct();             //Event reconstruction
    virtual void PostInit();                //Initialisation etc.

  ClassDef(TA2MesonPhysics,1)
};

//-----------------------------------------------------------------------------

#endif

