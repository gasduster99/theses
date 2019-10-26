//Observations.cpp
#include "AnalysisLocation.h"

/**********************************************************************
													HEAD OBSERVATION POINTS
**********************************************************************/
CHeadObsArray::CHeadObsArray():CAnalysisLocation(){
	pAq=NULL;
	pHeadObservations=NULL;
	nHeadObs=0;
}
//---------------------------------------------------------------------
CHeadObsArray::CHeadObsArray(const CAquiferABC *pAquifer):CAnalysisLocation(){
	pAq=pAquifer;
	pHeadObservations=NULL;
	nHeadObs=0;
}
//---------------------------------------------------------------------
CHeadObsArray::~CHeadObsArray(){
	if (pHeadObservations!=NULL){
		for (int i=0; i<nHeadObs; i++){delete pHeadObservations[i];}
	}
	delete [] pHeadObservations;
}
//*********************************************************************
void CHeadObsArray::AddObservation(const cmplex &z, const double &h, char *Name){
	HeadObs *obsptr=new HeadObs;
	obsptr->head   =0.0;
	obsptr->obshead=h;
	obsptr->pt     =c2Dto3D(z);
	if (Name!=NULL){
    obsptr->name= new char[strlen(Name)+1];
    obsptr->name=strcpy(obsptr->name,Name);
	}
	if (!DynArrayAppend((void**&)(pHeadObservations),(void*)(obsptr),nHeadObs)){
	  ExitGracefully("CHeadObsArray::AddObservation: adding NULL Head Observation",BAD_DATA);};
}
//*********************************************************************
void CHeadObsArray::CalculateStatistics(const double &t){
	for (int i=0; i<nHeadObs; i++){
		pHeadObservations[i]->head=pAq->GetHead(pHeadObservations[i]->pt,t);
	}
}
//*********************************************************************
void CHeadObsArray::WriteOutput   (const double &t) const{
	ofstream HEADS;
	double RMSerror(0.0);
	double err;
	int    i;
	if (nHeadObs>0){
		cout <<"    Head Observations("<<nHeadObs<<")..."<<endl;
		HEADS.open("headerr.csv");
		for (i=0; i<nHeadObs; i++){
			err=pHeadObservations[i]->head-pHeadObservations[i]->obshead;
			HEADS << pHeadObservations[i]->pt.x     <<",";
			HEADS << pHeadObservations[i]->pt.y     <<",";
			HEADS << pHeadObservations[i]->head     <<",";
			HEADS << pHeadObservations[i]->obshead  <<",";
			HEADS << err                            <<endl;
			RMSerror+=err*err;
		}	
		HEADS<<endl<<"SUM_SQ_ERROR: "<<RMSerror<<endl;	
		RMSerror=pow(RMSerror,0.5);
		HEADS<<endl<<"RMS_ERROR: "<<RMSerror<<endl;	

		ofstream OBSOUT;
		OBSOUT.open("obs_errors.csv",ios::app);

		//"HEAD_OBS",[index],[name],[x],[y],[z],[modeled head],[observed head],[error] 
		for (i=0; i<nHeadObs; i++){
			err=pHeadObservations[i]->head-pHeadObservations[i]->obshead;
			OBSOUT << "HEAD_OBS,"<< i                <<",";
			OBSOUT << pHeadObservations[i]->name     <<",";
			OBSOUT << pHeadObservations[i]->pt.x     <<",";
			OBSOUT << pHeadObservations[i]->pt.y     <<",";
			OBSOUT << pHeadObservations[i]->pt.z     <<",";
			OBSOUT << pHeadObservations[i]->head     <<",";
			OBSOUT << pHeadObservations[i]->obshead  <<",";
			OBSOUT << err                        <<endl;
		}	
		OBSOUT.close();
	}
}









/**********************************************************************
													BASE FLOW OBSERVATION POINTS
**********************************************************************/
CBFObsArray::CBFObsArray():CAnalysisLocation(){
	pAq=NULL;
	pObsArray=NULL;
	nObs=0;
}
//---------------------------------------------------------------------
CBFObsArray::CBFObsArray(const CAquiferABC *pAquifer):CAnalysisLocation(){
	pAq=pAquifer;
	pObsArray=NULL;
	nObs=0;
}
//---------------------------------------------------------------------
CBFObsArray::~CBFObsArray(){
	if (pObsArray!=NULL){
		for (int i=0; i<nObs; i++){delete pObsArray[i];}
	}
	delete [] pObsArray;
}
//*********************************************************************
void CBFObsArray::AddObservation(const cmplex &z, const double &flow, char *Name){
	BFObs *obsptr=new BFObs;
	obsptr->flow   =0.0;
	obsptr->obsflow=flow;
	obsptr->pt     =c2Dto3D(z);
	if (Name!=NULL){
    obsptr->name= new char[strlen(Name)+1];
    obsptr->name=strcpy(obsptr->name,Name);
	}
	if (!DynArrayAppend((void**&)(pObsArray),(void*)(obsptr),nObs)){
	  ExitGracefully("CBFObsArray::AddObservation: adding NULL BF Observation",BAD_DATA);};
}
//*********************************************************************
void CBFObsArray::CalculateStatistics(const double &t){
	for (int i=0; i<nObs; i++){
		pObsArray[i]->flow=pAq->GetBaseflow(pObsArray[i]->pt,t);
	}
}
//*********************************************************************
void CBFObsArray::WriteOutput   (const double &t) const{
	ofstream BFS;
	double RMSerror(0.0);
	double err;
	int    i;
	if (nObs>0){
		cout <<"    Baseflow Observations ("<<nObs<<")..."<<endl;
		for (i=0; i<nObs; i++){
			err=pObsArray[i]->flow-pObsArray[i]->obsflow;
			RMSerror+=err*err;
		}	
		//BFS<<endl<<"SUM_SQ_ERROR: "<<RMSerror<<endl;	
		//RMSerror=pow(RMSerror,0.5);
		//BFS<<endl<<"RMS_ERROR: "<<RMSerror<<endl;	

		ofstream OBSOUT;
		OBSOUT.open("obs_errors.csv",ios::app);

		//"BF_OBS",[index],[name],[x],[y],[z],[modeled base flow],[observed BF],[error] 
		for (i=0; i<nObs; i++){
			err=pObsArray[i]->flow-pObsArray[i]->obsflow;
			OBSOUT << "BF_OBS,"<< i          <<",";
			OBSOUT << pObsArray[i]->name     <<",";
			OBSOUT << pObsArray[i]->pt.x     <<",";
			OBSOUT << pObsArray[i]->pt.y     <<",";
			OBSOUT << pObsArray[i]->pt.z     <<",";
			OBSOUT << pObsArray[i]->flow     <<",";
			OBSOUT << pObsArray[i]->obsflow  <<",";
			OBSOUT << err                    <<endl;
		}	
		OBSOUT.close();
		}
}









/**********************************************************************
													GRADIENT OBSERVATION POINTS
**********************************************************************/
CGradObsArray::CGradObsArray():CAnalysisLocation(){
	pAq=NULL;
	pObsArray=NULL;
	nObs=0;
}
//---------------------------------------------------------------------
CGradObsArray::CGradObsArray(const CAquiferABC *pAquifer):CAnalysisLocation(){
	pAq=pAquifer;
	pObsArray=NULL;
	nObs=0;
}
//---------------------------------------------------------------------
CGradObsArray::~CGradObsArray(){
	if (pObsArray!=NULL){
		for (int i=0; i<nObs; i++){delete pObsArray[i];}
	}
	delete [] pObsArray;
}
//*********************************************************************
void CGradObsArray::AddObservation(const cmplex &z, const double &grad, const double &angle, char *Name){
	//angle in radians
	GradObs *obsptr=new GradObs;
	obsptr->gradient   =0.0;
	obsptr->obsgradient=grad;
	obsptr->angle      =0.0;
	obsptr->obsangle   =angle;
	obsptr->pt         =c2Dto3D(z);
	if (Name!=NULL){
    obsptr->name= new char[strlen(Name)+1];
    obsptr->name=strcpy(obsptr->name,Name);
	}
	if (!DynArrayAppend((void**&)(pObsArray),(void*)(obsptr),nObs)){
	  ExitGracefully("CGradObsArray::AddObservation: adding NULL Grad Observation",BAD_DATA);};
}
//*********************************************************************
void CGradObsArray::CalculateStatistics(const double &t){
	cmplex v;
	double n,K;
	for (int i=0; i<nObs; i++){
		v=c3Dto2D(pAq->GetVelocity2D(pObsArray[i]->pt,t));
		n=pAq->GetPoro(pObsArray[i]->pt);
		K=pAq->GetCond(pObsArray[i]->pt);
		pObsArray[i]->angle=adjarg(v);
		if (K>0){
			pObsArray[i]->gradient=-abs(v*n/K); //actual gradient (steepest descent)
			v=cmplex(v.real()*cos(pObsArray[i]->obsangle)+v.imag()*sin(pObsArray[i]->obsangle), //rotate coordinate system
               v.imag()*cos(pObsArray[i]->obsangle)-v.real()*sin(pObsArray[i]->obsangle));
			pObsArray[i]->gradcomponent=-v.real()*n/K; //gradient in specified direction
		}
		else    {pObsArray[i]->gradient=pObsArray[i]->gradcomponent=0;}
		
	}
}
//*********************************************************************
void CGradObsArray::WriteOutput   (const double &t) const{
	ofstream GradS;
	double RMSerror(0.0),RMSerror2(0.0),RMSerror3(0.0);
	double err,err2,err3;
	int    i;
  
	if (nObs>0){
		cout <<"    Gradient Observations("<<nObs<<")..."<<endl;
		for (i=0; i<nObs; i++){
			err =pObsArray[i]->gradient     -pObsArray[i]->obsgradient;
      err2=pObsArray[i]->gradcomponent-pObsArray[i]->obsgradient;
			err3=pObsArray[i]->angle        -pObsArray[i]->obsangle;
			err3=min(fabs(err3),2*PI-fabs(err3));

			RMSerror+=err*err;
			RMSerror2+=err2*err2;
			RMSerror3+=err3*err3;
		}	
		//GradS<<endl<<"SUM_SQ_ERROR: "<<RMSerror<<endl;	
		//RMSerror=pow(RMSerror,0.5);
		//GradS<<endl<<"RMS_ERROR: "<<RMSerror<<endl;	

		ofstream OBSOUT;
		OBSOUT.open("obs_errors.csv",ios::app);

		//"Grad_OBS",[index],[name],[x],[y],[z],
		//           [modeled gradient (steepest descent)],[modeled gradient (in observed direction)],
		//           [observed Grad],[error (in magnitude)], [error (in this direction)], 
		//           [modeled angle], [observed angle], [error in angle] 
		for (i=0; i<nObs; i++){
			err =pObsArray[i]->gradient     -pObsArray[i]->obsgradient;
      err2=pObsArray[i]->gradcomponent-pObsArray[i]->obsgradient;
			err3=pObsArray[i]->angle        -pObsArray[i]->obsangle;
			err3=min(fabs(err3),2*PI-fabs(err3));

			OBSOUT << "GRAD_OBS,"<< i            <<",";
			OBSOUT << pObsArray[i]->name         <<",";
			OBSOUT << pObsArray[i]->pt.x         <<",";
			OBSOUT << pObsArray[i]->pt.y         <<",";
			OBSOUT << pObsArray[i]->pt.z         <<",";
			OBSOUT << pObsArray[i]->gradient     <<",";//modeled gradient 
			OBSOUT << pObsArray[i]->gradcomponent<<",";//modeled gradient in specified direction
			OBSOUT << pObsArray[i]->obsgradient  <<",";//observed gradient
			OBSOUT << err                        <<",";
			OBSOUT << err2                       <<",";
			OBSOUT << pObsArray[i]->angle*180.0/PI    <<",";//modeled angle
			OBSOUT << pObsArray[i]->obsangle*180.0/PI <<",";//observed angle
			OBSOUT << err3*180.0/PI                   <<endl;
		}	
		OBSOUT.close();
	}
}









/**********************************************************************
													BASE FLOW OBSERVATION POINTS
**********************************************************************/
CLakeObsArray::CLakeObsArray():CAnalysisLocation(){
	pAq=NULL;
	pObsArray=NULL;
	nObs=0;
}
//---------------------------------------------------------------------
CLakeObsArray::CLakeObsArray(const CAquiferABC *pAquifer):CAnalysisLocation(){
	pAq=pAquifer;
	pObsArray=NULL;
	nObs=0;
}
//---------------------------------------------------------------------
CLakeObsArray::~CLakeObsArray(){
	if (pObsArray!=NULL){
		for (int i=0; i<nObs; i++){delete pObsArray[i];}
	}
	delete [] pObsArray;
}
//*********************************************************************
void CLakeObsArray::AddObservation(const cmplex &z, const double &flux, char *Name){
	LakeObs *obsptr=new LakeObs;
	obsptr->flux   =0.0;
	obsptr->obsflux=flux;
	obsptr->pt     =c2Dto3D(z);
	if (Name!=NULL){
    obsptr->name= new char[strlen(Name)+1];
    obsptr->name=strcpy(obsptr->name,Name);
	}
	if (!DynArrayAppend((void**&)(pObsArray),(void*)(obsptr),nObs)){
	  ExitGracefully("CLakeObsArray::AddObservation: adding NULL BF Observation",BAD_DATA);};
}
//*********************************************************************
void CLakeObsArray::CalculateStatistics(const double &t){
	for (int i=0; i<nObs; i++){
		pObsArray[i]->flux=0.0;//pAq->GetLakeFlux(pObsArray[i]->pt,t);
	}
}
//*********************************************************************
void CLakeObsArray::WriteOutput   (const double &t) const{
	ofstream BFS;
	double RMSerror(0.0);
	double err;
	int    i;
	if (nObs>0){
		cout <<"    Lake Flux Observations ("<<nObs<<")..."<<endl;
		for (i=0; i<nObs; i++){
			err=pObsArray[i]->flux-pObsArray[i]->obsflux;
			RMSerror+=err*err;
		}	
		//BFS<<endl<<"SUM_SQ_ERROR: "<<RMSerror<<endl;	
		//RMSerror=pow(RMSerror,0.5);
		//BFS<<endl<<"RMS_ERROR: "<<RMSerror<<endl;	

		ofstream OBSOUT;
		OBSOUT.open("obs_errors.csv",ios::app);

		//"LAKEFLUX_OBS",[index],[name],[x],[y],[z],[modeled lake flux],[observed lake flux],[error] 
		for (i=0; i<nObs; i++){
			err=pObsArray[i]->flux-pObsArray[i]->obsflux;
			OBSOUT << "LAKEFLUX_OBS,"<< i          <<",";
			OBSOUT << pObsArray[i]->name     <<",";
			OBSOUT << pObsArray[i]->pt.x     <<",";
			OBSOUT << pObsArray[i]->pt.y     <<",";
			OBSOUT << pObsArray[i]->pt.z     <<",";
			OBSOUT << pObsArray[i]->flux     <<",";
			OBSOUT << pObsArray[i]->obsflux  <<",";
			OBSOUT << err                    <<endl;
		}	
		OBSOUT.close();
		}
}