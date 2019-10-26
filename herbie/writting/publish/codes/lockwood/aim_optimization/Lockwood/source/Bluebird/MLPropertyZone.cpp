#include "MLPropZone.h"
//************************************************************************
//                           CmlPropZone Constructors
//************************************************************************
CmlPropZone::CmlPropZone(const proptype typ,const double *val){
	ptype  =typ;
	for (int l=0; l<MAX_MLAYERS;l++){
		values[l]=val[l];
	}
	pAllPropZones[TotalPropZones]=this;
	TotalPropZones++;
}
//------------------------------------------------------------------------
CmlPropZone::~CmlPropZone(){
	if (globaldebug){cout<<"  DESTROYING MULTILAYER PROPERTY ZONE"<<endl;}
}
//************************************************************************
//									ACCESSORS
//************************************************************************
Ironclad1DArray CmlPropZone::GetValues()          const{return values;}
//------------------------------------------------------------------------
proptype CmlPropZone::GetType ()                const{return ptype;}
//------------------------------------------------------------------------
double   CmlPropZone::GetArea ()                const{return area;}
//------------------------------------------------------------------------
Ironclad1DArray   CmlPropZone::GetValues(const cmplex &z) const{
	static double novalue[MAX_MLAYERS];

	for (int l=0; l<MAX_MLAYERS;l++){
		novalue[l]=NO_VALUE;
	}
	if (IsInside(z)){return values;   }
	else            {return novalue;}
}
//************************************************************************
bool CmlPropZone::IsInside  (const cmplex &z                ) const {cout <<"PZVirBug1"<<endl;return false; }
//------------------------------------------------------------------------
bool CmlPropZone::IsInSquare(const cmplex &zc,const double w) const {cout <<"PZVirBug2"<<endl;return false; }
/*************************************************************************
							STATIC MEMBER FUNCTIONS
*************************************************************************/
int         CmlPropZone::TotalPropZones=0;
//------------------------------------------------------------------------
CmlPropZone  *CmlPropZone::pAllPropZones[];
//------------------------------------------------------------------------
void CmlPropZone::DestroyAllMLPropZones(){
	if (globaldebug){cout <<"DESTROYING ALL PROPERTY ZONES"<<endl;}
	for (int i=0; i<TotalPropZones; i++){delete pAllPropZones[i];}
}

/************************************************************************
					Nested Sift
------------------------------------------------------------------------
Searches through an array (size=nzones) of pointers to abstract property zones (pZones)
returns the value of the zone in which z is located 
returns back_val if z is outside of all zones
MUST OPTIMIZE -this would be optimized if zonearray was sorted by size
************************************************************************/
Ironclad1DArray CmlPropZone::NestedSift(const cmplex	 &z, 
																				CmlPropZone   **pZones, 
																				const int				nzones, 
																				Ironclad1DArray back_val){
	if (nzones==0){return back_val;}   //quick out

	Unchangeable1DArray vlocal;
	static double smallArea;
	static double thisval[MAX_MLAYERS];
	static double tmpArea;	
	int l;
	for (l=0;l<MAX_MLAYERS;l++){
		thisval[l]=NO_VALUE;
	}

	smallArea=ALMOST_INF;
	for (int i=0; i<nzones; i++){      //sift thru zones
    vlocal=pZones[i]->GetValues(z);   //get value due to zone at z
		if (vlocal[0]!=NO_VALUE){                
      tmpArea=pZones[i]->GetArea();  //get area of zone
			if (tmpArea<smallArea){        //if zone is smaller than smallest zone 
				for (l=0;l<MAX_MLAYERS;l++){
					thisval[l]=vlocal[l];            //then this is true value (accounts for nesting)
				}
        smallArea=tmpArea;                 
			} 
		}
	}
	if (thisval[0]==NO_VALUE){return back_val;}       //if no value is found, background is used
	else										 {return thisval; }
}
