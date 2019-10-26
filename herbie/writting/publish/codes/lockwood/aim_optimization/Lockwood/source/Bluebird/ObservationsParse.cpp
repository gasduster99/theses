//ObservationsParse.cpp
#include "Bluebird.h"
class CObservations;

/******************************************************************************
		OBSERVATIONS FILE PARSING ROUTINE
-------------------------------------------------------------------------------
IMPORTANT: order of commands:
					-
******************************************************************************/
bool ObservationsParse(char             *filename,
									     CAquifer        *&pAq){

  //parsing variables----------------------------------------------------
  char  *s[MAXINPUTITEMS];
	char  *thisname=new char  [255];

  int    Len,code,l(0),i;
	bool   eof (false);
	bool   done(false);
	bool   noisy(parserdebug); //FOR PARSER DEBUG, turn parserdebug to true (in "MasterInclude.h")

	//creation variables-------------------------------------------------
	CHeadObsArray *pHeadObsArray=new CHeadObsArray(pAq); 
	CBFObsArray   *pBFObsArray  =new CBFObsArray(pAq); 
	CGradObsArray *pGradObsArray=new CGradObsArray(pAq); 
  CLakeObsArray *pLakeObsArray=new CLakeObsArray(pAq);

 	char             blank[1]; blank[0]='\0';
	char             space[1]; space[0]=' ';

  ifstream OBS(filename);  
  if (OBS.fail()){cout << "Cannot find observation file "<<filename <<endl; return false;}

	ExitGracefullyIf(pAq==NULL,"ObservationsParse: NULL Aquifer",RUNTIME_ERR);

  //-----------------------------------------------------------------
  cout << "Parsing Observations File " << filename <<"..."<<endl;
	cout << "------------------------------------------------------------"<<endl;

	//--sift through file-----------------------------------------------
	while ((!TokenizeLine(OBS,s,Len)) && (!eof)){

		l++;  
		if (noisy){ cout << "reading line " << l << ": ";}

		if      (Len==0)                                                            {code=-1;}
		else if ((!strcmp(s[0],"HeadObs"               )) )                         {code=101;}
		else if ((!strcmp(s[0],"BaseflowObs"           )) )                         {code=102;}
		else if ((!strcmp(s[0],"GradientObs"           )) )                         {code=103;}
		else if ((!strcmp(s[0],"SWFluxObs"             )) )                         {code=104;}

		//disabled or commented out
    else if ((!strcmp(s[0],"#"										 )) || 
			       (!strcmp(s[0],"rem"									 )) || 
						 (!strcmp(s[0],"*"										 )) || 
						 (!strcmp(s[0],"&"										 )) )													{code=-1;}
		//end file
		else if ((!strcmp(s[0],"end"									 )) || 
						 (!strcmp(s[0],"End"									 )) || 
			       (!strcmp(s[0],"EndInput"							 )) )							            {code=-1;eof=true;}
		//unrecognized
    else																																			  {code=-3;}


		//disable 

    switch(code){
		case(101)://-----------------------------------------------------------------------------------
			{ //Head Observation (single-layer)
		    //"HeadObs" {double x-location} {double y-location} {double observed head [L]} [string name]
        if (noisy){cout<<"Head Observation"<<endl;}
				if (Len>=4){
					thisname=strcpy(thisname,blank);
					for(i=4; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
					pHeadObsArray->AddObservation(s_to_c(s[1],s[2]),s_to_d(s[3]),thisname);
				}
		    break;
		  }
		case(102)://-----------------------------------------------------------------------------------
			{ //Base Flow Observation (single-layer)
		    //"BaseFlowObs" {double x-location} {double y-location} {double observed cumulative baseflow [L^3/T]} [string name]
			  //&
        if (noisy){cout<<"Base Flow Observation"<<endl;}
				if (Len>=4){
					thisname=strcpy(thisname,blank);
					for(i=4; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
					pBFObsArray->AddObservation(s_to_c(s[1],s[2]),s_to_d(s[3]),thisname);
				}
		    break;
		  }
		case(103)://-----------------------------------------------------------------------------------
			{ //Gradient Observation (single-layer)
		    //"GradientObs"{double x-location} {double y-location} {double gradient(dimless)} {double orientation (degrees)}  [string name]
        if (noisy){cout<<"Gradient Observation"<<endl;}
				if (Len>=5){
					thisname=strcpy(thisname,blank);
					for(i=5; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
					pGradObsArray->AddObservation(s_to_c(s[1],s[2]),s_to_d(s[3]),s_to_d(s[4])*PI/180,thisname);
				}
		    break;
		  }
		case(104)://-----------------------------------------------------------------------------------
			{//SWFlux Observation
		   //{double x-location} {double y-location} {char gaining/losing}
			 //&
       if (noisy){cout<<"Surface water flux Observation"<<endl;}
		   break;
		  }
		case(105)://-----------------------------------------------------------------------------------
			{ //Lake Flux Observation (single-layer)
		    //"LakeFluxObs"{double x-location} {double y-location} {double orientation (degrees)} {double gradient(dimless)} [string name]
        if (noisy){cout<<"Lake Flux Observation"<<endl;}
				if (Len>=4){
					thisname=strcpy(thisname,blank);
					for(i=4; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
					pLakeObsArray->AddObservation(s_to_c(s[1],s[2]),s_to_d(s[3]),thisname);
				}
		    break;
		  }
		case(201)://-----------------------------------------------------------------------------------
			{//Travel Time Observation
		   //{double particle x-location} {double particle y-location} {double time-to-capture}
			 //&
       if (noisy){cout<<"Travel Time Observation"<<endl;}
		   break;
		  }
	  case(-1)://------------------------------------------------------------------------------------
		  {//comment out
				if (noisy) {cout <<"##"<< endl;}
				break;   
		  }
	  case(-2)://------------------------------------------------------------------------------------
		  {//disabled
				if (noisy) {cout <<"disabled"<< endl;} 
				break;  
		  }
		case(-3)://------------------------------------------------------------------------------------
		  {//unrecognize
				cout <<"observation input file: command unrecognized at line "<< l << endl;
				break;   
		  }
		default: {
			  if (noisy) {cout << "##"<<endl;} 
				break;
			}//{cout<< "BAD OBSERVATIONS FILE INPUT, line "<< l+1 <<endl;return false;}

		}/*end switch*/
	}/* end while Tokenize & !eof*/



  //*******************************************************************************************
	//                       DONE PARSING- CREATE & BUILD SYSTEM
  //*******************************************************************************************

	OBS.close();

	cout<<"...done parsing observations input file"<<endl<<endl;

	//delete [] thisname;

	return true;
}
/*---------------------------------------------------------------------
  OldHead File Parse
	For backward compatibility
	should be only called if observations.dat doesnt exist
---------------------------------------------------------------------*/
bool OldHeadFileParse(CAquifer        *&pAq){
  int    Len,i;
	char  *thisname=new char  [255];

	char   blank[1]; blank[0]='\0';
	char   space[1]; space[0]=' ';

	ExitGracefullyIf(pAq==NULL,"Parsing Head.dat: NULL Aquifer",RUNTIME_ERR);

	CHeadObsArray *pObsArray=new CHeadObsArray(pAq); 
	ifstream HEAD;
	HEAD.open("head.dat");
	char  *hs[MAXINPUTITEMS];			 
	if (HEAD.fail()){ExitGracefully("Cannot find file head.dat",BAD_DATA);return false;}
	while (!TokenizeLine(HEAD,hs,Len)){
		if (Len>3){
			thisname=strcpy(thisname,blank);
			for(i=3; i<Len; i++){thisname=strcat(thisname,hs[i]);thisname=strcat(thisname,space);}
			pObsArray->AddObservation(s_to_c(hs[0],hs[1]),s_to_d(hs[2]),thisname);
		}
	}
	HEAD.close();
	return true;
}