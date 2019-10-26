//SimpleReactions.cpp
#include "ReactionScheme.h"
/***************************************************************************
							SIMPLE DECAY
***************************************************************************/
CSimpleDecay::CSimpleDecay()
             :CReactionScheme(false){
	DecayRates=NULL;
}
//---------------------------------------------------------------------------
CSimpleDecay::CSimpleDecay(double *decayrates)
             :CReactionScheme(false){
	DecayRates=new double [nspecies];
	for (int s=0; s<nspecies; s++){DecayRates[s]=decayrates[s];}
}
//---------------------------------------------------------------------------
CSimpleDecay::~CSimpleDecay(){
	if (globaldebug){cout <<"DESTROYING SIMPLE DECAY SCHEME "<<endl;}
	delete [] DecayRates;
}
/***************************************************************************/
void CSimpleDecay::React(Ironclad1DArray  Caq,
												 Ironclad1DArray  Caq_end,
												 Writeable1DArray Caq_new,
												 Ironclad1DArray  Cs,
												 Ironclad1DArray  Cs_end,
												 Writeable1DArray Cs_new,	
												 const int        soil,
												 const double    &n,  
												 const double    &tstep){
	int s;

	if (Caq==NULL){RxNExitGracefully("CSimpleDecay::React: provided NULL concentration array",RUNTIME_ERR);}

	for (s=0; s<nspecies; s++){
		  Caq_new[s]=Caq_end[s]-0.5*(Caq_end[s]+Caq[s])*DecayRates[s]*tstep;

		if (Cs!=NULL){
			 Cs_new[s]= Cs_end[s]-0.5*( Cs_end[s]+ Cs[s])*DecayRates[s]*tstep;
		}
	}

}
/***************************************************************************
							PARENT DAUGHTER DECAY
***************************************************************************/