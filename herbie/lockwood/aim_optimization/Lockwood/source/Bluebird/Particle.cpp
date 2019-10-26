#include "Particle.h"
#include "FluxGrid.h"
/**********************************************************************
  Constructors
**********************************************************************/
CParticle::CParticle(){}
//--------------------------------------------------------------------
CParticle::CParticle(int							   ID, 
										 const CAquiferABC  *pAquifer,
										 trackdir            direct){

  partID=ID; 
	pAq=pAquifer;
	if (direct==TRACK_FORWARD){dir= 1;} 
	else                      {dir=-1;}
}
//--------------------------------------------------------------------
CParticle::~CParticle(){}
/*********************************************************************
   STATIC MEMBERS
*********************************************************************/
double CParticle::sm_change=     0.01;		
double CParticle::lg_change=     0.02;		
double CParticle::bad_change=    0.03;		
double CParticle::sm_adjust=     1.20;		
double CParticle::lg_adjust=     0.80;			
double CParticle::bad_adjust=    0.50;			
//---------------------------------------------------------------------
void   CParticle::SetTrackPrecision(const double sm_ch,const double lg_ch,const double bd_ch,
																		const double sm_ad,const double lg_ad,const double bd_ad){
	sm_change =sm_ch; 
	sm_adjust =sm_ad;
	lg_change =lg_ch; 
	lg_adjust =lg_ad;
	bad_change=bd_ch; 
	bad_adjust=bd_ad;
}
//---------------------------------------------------------------------
void CParticle::SetPrecision(int Precision){
  switch(Precision){
		case(0):{CParticle::SetTrackPrecision(0.10 ,0.15,0.30 ,1.3,0.9,0.6);break;}
		case(1):{CParticle::SetTrackPrecision(0.07 ,0.10,0.14 ,1.3,0.9,0.6);break;}
		case(2):{CParticle::SetTrackPrecision(0.05 ,0.07,0.10 ,1.2,0.8,0.5);break;}
		case(3):{CParticle::SetTrackPrecision(0.03 ,0.04,0.06 ,1.2,0.8,0.5);break;}
		case(4):{CParticle::SetTrackPrecision(0.02 ,0.03,0.045,1.2,0.8,0.5);break;}
		case(5):{CParticle::SetTrackPrecision(0.01 ,0.02,0.03 ,1.2,0.8,0.5);break;}
		case(9):{CParticle::SetTrackPrecision(0.005,0.01,0.017,1.2,0.8,0.5);break;}
	}
}
/*********************************************************************
   Manipulators 
*********************************************************************/
void CParticle::SetDirection(trackdir	direct){
	if (direct==TRACK_FORWARD){dir=1; } 
	else                      {dir=-1;}
}
//--------------------------------------------------------------------
void CParticle::SetAquifer(const CAquiferABC  *pAquifer){
	ExitGracefullyIf((pAquifer==NULL),"CParticle::SetAquifer: NULL Aquifer",RUNTIME_ERR);
	pAq=pAquifer;
}
/*********************************************************************
   Accessors 
*********************************************************************/
int  CParticle::GetID() const                  {return partID;}
/*********************************************************************
   TRACK
----------------------------------------------------------------------
	tracks an abstract particle through a three dimensional domain
*********************************************************************/
void CParticle::Track(const double    &timeperiod, 
											const advtype    atype, 
											      double     tstep, //if atype=ADAPTIVE_TIME_STEP, tstep=min tstep
											      double     sstep, //if atype=ADAPTIVE_TIME_STEP, sstep=min step
											const bool       eff_vel, 
											const disp_info &disp, 
											const bool       track3D){

	//MUST OPTIMIZE
  bool            step_at_min(true),done(false);
  int             counter(0),NumSteps(0);
  static pt3D     pt;
  static vector   v1,v2,v3,v4, K1,K2,K3,K4,vstart,v1old;
  static double   vx,vy,vt,vz,change,starttime;
	static double   min_time_step, min_step,t,lastt;

	if      (IsCaptured()){ return;} 

	pt        =GetLocation();
	t         =GetCurrentT();
	lastt     =GetLastT();
  starttime=t;

	//set duration of first time step
	if (track3D)
	{
		if (eff_vel){vstart=pAq->GetEffVelocity3D(pt,t,disp);}
		else        {vstart=pAq->   GetVelocity3D(pt,t     );}
	}
	else        
	{
		if (eff_vel){vstart=pAq->GetEffVelocity2D(pt,t,disp);}
		else        {vstart=pAq->   GetVelocity2D(pt,t     );}
	}

  vt=abs(vstart); 

	/*if      (atype==ADAPTIVE_TIME_STEP) {cout <<"adaptive tstep:  "<<tstep<< " sstep: "<<sstep<<endl;}
	else if (atype==CONSTANT_SPACE_STEP){cout <<"const space step:"<<tstep<< " sstep: "<<sstep<<endl;}
	else                                {cout <<"const time step: "<<tstep<< " sstep: "<<sstep<<endl;}*/

	if (atype==ADAPTIVE_TIME_STEP)
	{
		min_time_step=max(G_MIN_TIME_STEP,tstep);
		min_step     =max(G_MIN_STEP,     sstep);

		if   (vt==0.0){tstep= min_time_step/DELTA0;}   
		else          {tstep= min_step     /vt    ;}
		tstep=max(max(tstep,(t-lastt)/2.0),G_MIN_TIME_STEP); //corrects for particles which already have good time steps
	}

  v1   =0.0;
	done =false;
  //track particle until time is up or it is captured or has exceeded MAX_STEPS
  while (!done){

   	pt=GetLocation();
		t =GetCurrentT();
   
		//Runge-Kutta integration of particle location---------------------------------
		v1old=v1;		
		if (track3D){ v1=double(dir)*pAq->GetVelocity3D(pt,t); }
		else        { v1=double(dir)*pAq->GetVelocity2D(pt,t); }

		if (atype==CONSTANT_SPACE_STEP){
			if(abs(v1)!=0.0){tstep=sstep/abs(v1);}
			else            {tstep=sstep/REALSMALL;}
		}
		
		//insures particle moves no further than requested
		if (tstep>(starttime+timeperiod-t)){tstep=(starttime+timeperiod-t);} 

		if (track3D)
		{
			if (eff_vel)
			{
																																	K1=tstep*v1;
				v2=double(dir)*pAq->GetEffVelocity3D(pt+(0.5*K1),t,disp); K2=tstep*v2;
				v3=double(dir)*pAq->GetEffVelocity3D(pt+(0.5*K2),t,disp); K3=tstep*v3;
				v4=double(dir)*pAq->GetEffVelocity3D(pt+ K3,     t,disp); K4=tstep*v4;
			}
			else
			{
																																  K1=tstep*v1;
				v2=double(dir)*pAq->GetVelocity3D(pt+(0.5*K1),t);			    K2=tstep*v2;
				v3=double(dir)*pAq->GetVelocity3D(pt+(0.5*K2),t);				  K3=tstep*v3;
				v4=double(dir)*pAq->GetVelocity3D(pt+ K3,     t);				  K4=tstep*v4;		
			}
		}
		else
		{
			if (eff_vel)
			{
																																	K1=tstep*v1;
				v2=double(dir)*pAq->GetEffVelocity2D(pt+(0.5*K1),t,disp); K2=tstep*v2;
				v3=double(dir)*pAq->GetEffVelocity2D(pt+(0.5*K2),t,disp); K3=tstep*v3;
				v4=double(dir)*pAq->GetEffVelocity2D(pt+ K3,     t,disp); K4=tstep*v4;
			}
			else
			{
																																  K1=tstep*v1;
				v2=double(dir)*pAq->GetVelocity2D(pt+(0.5*K1),t);			    K2=tstep*v2;
				v3=double(dir)*pAq->GetVelocity2D(pt+(0.5*K2),t);				  K3=tstep*v3;
				v4=double(dir)*pAq->GetVelocity2D(pt+ K3,     t);				  K4=tstep*v4;		
			}		
		} /* end if (track3D) */
   	
		vx=(v1.x+(2.0*v2.x)+(2.0*v3.x)+v4.x)/6.0;
    vy=(v1.y+(2.0*v2.y)+(2.0*v3.y)+v4.y)/6.0;
    if (track3D){vz=(v1.z+(2.0*v2.z)+(2.0*v3.z)+v4.z)/6.0; }else{vz=0.0;}
    vt=sqrt(vx*vx+vy*vy+vz*vz);

		Advect((1.0/6.0*(K1+(2.0*K2)+(2.0*K3)+K4)),tstep);

		if (vt==0.0){Capture(timeperiod);/*cout <<"captured0*"<<endl;*/} 

		//check for sinks,final time, or max steps--------------------------------------
		if ((((v1old.x*v1.x)<0.0) && (v1old.x*v1old.x>0.5*abs(v1old))) ||
		    (((v1old.y*v1.y)<0.0) && (v1old.y*v1old.y>0.5*abs(v1old))) || 
				(((v1old.z*v1.z)<0.0) && (v1old.z*v1old.z>0.5*abs(v1old)))){
			Capture(timeperiod); //capture @ drastic changes in velocity
		}

		if      (HasBeenCaptured())          { Capture(timeperiod); /*cout <<"captured"<<endl;*/  done=true;}//for now->should be overall duration
		else if (NumSteps>=MAX_STEPS-1)			 { cout <<"exceeded max_steps!!"<<endl;               done=true;}
		else if (t>=starttime+timeperiod)		 { /*cout <<"tracked."            <<endl; */          done=true;}//must set final time
		//else if      (ProgramAborted())      {                                                    done=true;}//really slow
   	
    //adjust time step, if desired -------------------------------------------------
		if (!done)
		{
			if ((atype==ADAPTIVE_TIME_STEP) && (vt!=0.0))//if non-stagnant point
			{  
				//cout <<endl<<"tstep:"<<tstep;
        //calculate dimensionless change
				change=0.0;
				change+=(fabs(v1.x-v2.x)+fabs(v1.x-v3.x)+fabs(v1.x-v4.x))/vt/6.0;
				change+=(fabs(v1.y-v2.y)+fabs(v1.y-v3.y)+fabs(v1.y-v4.y))/vt/6.0;
				change+=(fabs(v1.z-v2.z)+fabs(v1.z-v3.z)+fabs(v1.z-v4.z))/vt/6.0;
				
				//if change is way too big, redo with much better time step
				if ((change>bad_change) && (!step_at_min)) 
				{  
					tstep      =max((tstep * bad_adjust), (min_step/vt));
					step_at_min=((fabs(tstep-min_step/vt)<REALSMALL) || (tstep<min_time_step));
					BackTrack();                                                     /*cout <<"B";*/
				} 
				else
				{
					//modify time step without redo
					if (fabs(tstep-min_step/vt)< REALSMALL)
					{
						step_at_min=true;
					}
					else
					{   

		 				//if step is too small, make bigger  but keep prior calc of position
						//if step is too big,   make smaller but keep prior calc of position
  					if      (change<sm_change){tstep=max(min((tstep*sm_adjust),(MAX_STEP/vt)),min_time_step);/*cout <<"s";*/}
						else if (change>lg_change){tstep=max(max((tstep*lg_adjust),(min_step/vt)),min_time_step);/*cout <<"b";*/}

  					//if time step creates minimum allowable movement, set step_at_min to true
						step_at_min=(fabs(tstep-min_step/vt)<REALSMALL);
					}
				}
			}//end if (atype==ADAPTIVE_TIME_STEP)...

			NumSteps++;

		}//end if !done...
  }//end while !done...
}
//*********************************************************************
void CParticle::TrackPollock(const double    &timeperiod,
														 const int        intervals,
														 CFluxGrid       *pFluxGrid,
											       const bool       track3D,
														 const bool       forward){

	double tstep=timeperiod/(double)(intervals);

	pt3D here;
	for (double t=0; t<timeperiod-0.5*tstep; t+=tstep){
		here=this->GetLocation();
		this->Advect(pFluxGrid->PollockTrack(here,tstep,forward)-here,tstep);
	}

}