#ifndef __TA2MesonPhysics_h__
#define __TA2MesonPhysics_h__

#include "TA2BasePhysics.h"
#include "TA2KFParticle.h"
#include "TA2CBKinematicFitter.h"

#define MASS_PIZERO    134.9766
#define MASS_ETA       547.7500

class TA2Apparatus;

class TA2MesonPhysics : public TA2BasePhysics
{
private:
        //TA2TAPS_Veto*	veto;

  protected:
        TH1D*	hist[7];
	TH1I*	count;
        TH2D*   dEvE[300];
        TH2D*   dEvE_OR;
        TH2D*   Cut_dEvE_OR;

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

