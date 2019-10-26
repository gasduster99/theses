#include "Particle.h"
#include "FluxGrid.h"
//*********************************************************************
//                     CONSTRUCTORS
//*********************************************************************
CPathline::CPathline(){}
//---------------------------------------------------------------------
CPathline::CPathline(int                ID, 
										 const CAquiferABC *pAquifer,
										 pt3D					      startpt, 
										 trackdir           direct):
					 CParticle(ID,pAquifer,direct){

  numsteps=0;

	firstnode=new pathnode;
  firstnode->pt=startpt;
	firstnode->t=0.0;
	firstnode->s=0.0;
	firstnode->next=PATHEND;
	firstnode->last=PATHSTART;
	lastnode=firstnode;
  captured=false;

	buffercount=0;
}
//---------------------------------------------------------------------
CPathline::~CPathline(){
 if (globaldebug){cout <<"   DESTROYING PATHLINE "<<endl;}
 while (firstnode->next!=PATHEND){
	 firstnode=firstnode->next;
	 delete firstnode->last;
 }
 delete firstnode;
 for (int i=0; i<buffercount;i++){delete buffarray[i];}
}
/*********************************************************************
                   STATIC MEMBERS
*********************************************************************/
pathnode *CPathline::PATHSTART =NULL;
pathnode *CPathline::PATHEND   =NULL;
int       CPathline::Resolution=100;
double    CPathline::Duration  =100;
CFluxGrid *CPathline::pPollockGrid=NULL;
void      CPathline::Destroy(){delete pPollockGrid;}
//---------------------------------------------------------------------
void   CPathline::SetDuration  (double time){CPathline::Duration  =time;}
void   CPathline::SetResolution(int    res ){CPathline::Resolution=res ;}
//---------------------------------------------------------------------
double CPathline::GetDuration  ()           {return CPathline::Duration;}
//---------------------------------------------------------------------
void   CPathline::SetPollockGrid   (CFluxGrid *pFluxGrid){pPollockGrid=pFluxGrid;}
/*********************************************************************
                   ACCESSOR FUNCTIONS
*********************************************************************/
int CPathline::GetNumOfSteps() const {return numsteps;}
//---------------------------------------------------------------------
pt3D CPathline::GetLocation(const double &T) const{
	pathnode *node;
	node=firstnode;
  while (node->next!=PATHEND){
		if      (node->t<T)            {node=node->next;}  //not there yet
		else if (node->last==PATHSTART){break;}            //start of the path
		else                           {                   //interpolate
			return (node->last->pt)+(T-node->last->t)/((node->t)-(node->last->t))*((node->pt)-(node->last->pt));
		}    																    	
	}
  return node->pt;                                     //start/end of the path
}
//*********************************************************************
double CPathline::GetCurrentT() const {return lastnode->t;}
//---------------------------------------------------------------------
pt3D   CPathline::GetLocation() const {return lastnode->pt;}
//---------------------------------------------------------------------
double CPathline::GetLastT   () const {
	if ((lastnode->last)==PATHSTART){return lastnode->t;}
	else                            {return lastnode->last->t;}
}

//*********************************************************************
bool CPathline::IsCaptured() const {return captured;}
//---------------------------------------------------------------------
bool CPathline::HasBeenCaptured(){
	//geometric manipulation to identify sporadic pathline (near sink)
	pathnode *l=lastnode;
  //cout <<"HasBeenCaptured "<<numsteps<<endl;
  if (captured){return true;}

	//TMP DEBUG - capture if exiting window==========================================
	/*window cap_win;
	cap_win.e= 20.000001;
	cap_win.w=-20.000001;
	cap_win.n=  0.000001;
	cap_win.s=- 1.000001;

	if ((l->pt.x<cap_win.w) ||
		  (l->pt.x>cap_win.e) ||
		  (l->pt.y<cap_win.s) ||
			(l->pt.y>cap_win.n)){return true;}//TMP DEBUG*/
	//TMP DEBUG - capture if exiting window==========================================

	if (((l->last            !=PATHSTART) &&
    	 (l->last->last      !=PATHSTART) &&
		   (l->last->last->last!=PATHSTART) &&

		 ((abs((l->pt)-(l->last->last->pt))/
			(abs((l->pt)-(l->last->pt)     )+abs((l->last->pt)-(l->last->last->pt)  )))<=BAD_TURN) &&

		 ((abs((l->last->pt)-(l->last->last->last->pt))/
		  (abs((l->last->pt)-(l->last->last->pt      ))+abs((l->last->last->pt)-(l->last->last->last->pt) )))<=BAD_TURN)) ||

	   ((l->last      !=PATHSTART) &&
			(l->last->last!=PATHSTART) &&
			(abs(l->pt-l->last->pt)>10.0*abs(l->last->pt-l->last->last->pt)))){

			 return true;	
	}
	return false;
}
//*********************************************************************
void CPathline::Capture(const double &endtime){
	
  if (!captured){
    double capture_time;
		pt3D   capture_loc;
    
		if (lastnode->last==PATHSTART){
			capture_time=lastnode->t;
			capture_loc =lastnode->pt;
		}
		else if (lastnode->last->last==PATHSTART){
			capture_time=lastnode->last->t;
			capture_loc =lastnode->last->pt;
		}
		else{
			capture_time=lastnode->last->last->t;
			capture_loc =lastnode->last->last->pt;
		}

		//TMP DEBUG - capture if exiting window==========================================
		/*window cap_win;
		cap_win.e= 20.000001;
		cap_win.w=-20.000001;
		cap_win.n=  0.000001;
		cap_win.s=- 1.000001;
		double x,y;
		pathnode *out=lastnode;
		pathnode *in =lastnode->last;
    if (out->pt.x<cap_win.w){
		  y           =LinInterp(out->pt.y,in->pt.y,out->pt.x,in->pt.x,cap_win.w);
			x=cap_win.w;
      capture_time=LinInterp(out->t   ,in->t   ,out->pt.x,in->pt.x,cap_win.w);
			capture_loc=pt3D(x,y,0);
		}
    if (out->pt.x> cap_win.e){
		  y           =LinInterp(out->pt.y,in->pt.y,out->pt.x,in->pt.x,cap_win.e);
			x=cap_win.e;
      capture_time=LinInterp(out->t   ,in->t   ,out->pt.x,in->pt.x,cap_win.e);
			capture_loc=pt3D(x,y,0);
		}
		else if (out->pt.y>cap_win.n){
			x           =LinInterp(out->pt.x,in->pt.x,out->pt.y,in->pt.y,cap_win.n);
			y=cap_win.n;
      capture_time=LinInterp(out->t   ,in->t   ,out->pt.y,in->pt.y,cap_win.n);
			capture_loc=pt3D(x,y,0);
		}
		else if (out->pt.y<cap_win.s){
			x           =LinInterp(out->pt.x,in->pt.x,out->pt.y,in->pt.y,cap_win.s);
			y=cap_win.s;
      capture_time=LinInterp(out->t   ,in->t   ,out->pt.y,in->pt.y,cap_win.s);
			capture_loc=pt3D(x,y,0);
		}*/
		//TMP DEBUG - capture if exiting window==========================================

		lastnode=lastnode->last;           //move back one node (final two are junk from tracking)
	 
		if (lastnode->next!=PATHEND){
			delete lastnode->next;             //delete last node, replace with end of pathline
			lastnode->next=PATHEND;
		}
		numsteps-=1;

		lastnode->t      =endtime;				 //modify time and location accordingly
		lastnode->pt     =capture_loc;
		captured=true;	
	}
		
}

//*********************************************************************
void CPathline::Advect(const vector &movement, const double &tstep){

	if (buffercount==0){
    //create new array of pathnodes
		for (int i=0;i<BUFFERSIZE;i++){buffarray[i]=new pathnode;}
		buffercount=BUFFERSIZE;
	}
	
	lastnode->next=buffarray[buffercount-1];      //create new node
	buffercount--;

	//lastnode->next=new pathnode;                //create new node
  numsteps++;

  //initialize node
	lastnode->next->pt=lastnode->pt+movement;     //move in space
	lastnode->next->t =lastnode->t +tstep;        //move in time
	lastnode->next->s =lastnode->s +abs(movement); 

	lastnode->next->next=PATHEND;   
  lastnode->next->last=lastnode;                //change lastnode

	lastnode=lastnode->next;	                    //reset value for last node as this new node

  if ((numsteps>100) && (numsteps%100==0)){cout <<".";}
	//WriteEllipsisToScreen(numsteps,100,1);
}
//---------------------------------------------------------------------
void CPathline::BackTrack(){
	if (lastnode->last!=PATHSTART){
		numsteps--;
		lastnode=lastnode->last; //back up one
		delete lastnode->next;   //delete last node
		lastnode->next=PATHEND;  
	}
}
//*********************************************************************
void CPathline::ReversePath(){
	//switches direction of path from t=T to 0 to t=0 to T
	pathnode *node;
	pathnode *tmpstart;
	pathnode *tmp;
	double lasttime, lastdist;
	
	node=lastnode;
	tmpstart=lastnode;
	lasttime=lastnode->t;
	lastdist=lastnode->s;
	while (node->next !=PATHSTART){
    tmp=node->last;
		node->last=node->next;
		node->next=tmp;
		node->t=lasttime-node->t;
		node->s=lastdist-node->s;
	}
	lastnode=firstnode;
	firstnode=lastnode;
	dir=-dir;

}
//*********************************************************************
void CPathline::CleanPath(){
  //removes unneccesary nodes in particle path, counts nodes, calculates cumulative distance s
	pathnode *a,*node,*c;  //->a->node->c->

	int numnodes;
		
	if      (firstnode->next      ==PATHEND){numnodes=1;return;}
	else if (firstnode->next->next==PATHEND){numnodes=2;return;}

	node=firstnode->next;  //start at second node
	numnodes=2;

	while (node->next!=PATHEND){
		a=node->last;
		c=node->next;
		//if node may be represented by linearization and distance between neighbors not too small, remove node
		if (((abs( (node->t-a->t)/(c->t-a->t)*(c->pt-a->pt)+a->pt-node->pt)/
			   abs(c->pt-a->pt))<LIN_TOLERANCE) &&
				 (abs(a->pt-c->pt)<abs(firstnode->pt-lastnode->pt)/CPathline::Resolution )){
			a->next=c;       //redirect pointers forward and backward of node to be deleted
			c->last=a;
			delete node;     //delete bad node
			node=c;          //goto next node
		}
		else{
			node->s=(node->last->s)+abs((node->pt)-(node->last->pt)); //calculate cumulative distance
			node=node->next;          //goto next node
			numnodes++;               //count number of steps
		}
	}

  //cout <<endl<< "Path cleaned. " <<numsteps<<" before. " << numnodes <<"after."<<endl;
	numsteps=numnodes;
}
//*********************************************************************
void CPathline::WriteOutput() const{
	int j;
	pathnode *node;
	ofstream TRACK;
  TRACK.open("tracks.bna",ios::app);
  TRACK.precision(16);
	//create starting point graphic
	TRACK << "\" Point " << partID << " \" 2" <<endl;
	TRACK	<< (double)(firstnode->pt.x)<< " ";
	TRACK	<< (double)(firstnode->pt.y)<< endl;
	TRACK	<< "0.1 0.1"         << endl;

	//Create pathline graphic
	node=firstnode;
	if (numsteps>2){
		TRACK << "\" track"<< partID << " \"  " << -numsteps <<endl;
		for (j=0; j<numsteps; j++){
			TRACK << node->pt.x <<" ";
			TRACK << node->pt.y <<endl;
			node=node->next;
		}  
	}
  TRACK.close();

	//stores more detailed information on particle tracks
	ofstream TRACK2;
  TRACK2.open("tracks.dat",ios::app);
  TRACK2.precision(16);
	node=firstnode;
	TRACK2<<numsteps<<endl;
	j=1;
	while (node!=lastnode){
		TRACK2<<node->pt.x<<" "<<node->pt.y<<" "<<node->pt.z<<" "<<node->t<<" "<<node->s<<" "<<j<< endl;
		node=node->next;
		j++;
	}
  TRACK2.close();

	ofstream TRACK3;
  TRACK3.open("mattstracks.csv",ios::app);
  TRACK3.precision(16);
	node=firstnode;
	TRACK3<<numsteps<<endl;
	j=1;
	while (node!=lastnode){
		TRACK3<<this->GetID()<<","<<node->t<<","<<node->pt.x<<","<<node->pt.y<<","<<node->s<<" "<<j<< endl;
		node=node->next;
		j++;
	}
  TRACK3.close();
}
//*********************************************************************
void CPathline::WriteOutput(const double &t) const{}
//*********************************************************************
void CPathline::TrackAllPathlines(CPathline   **Particles,
																	const int     numpart, 
																	const double  timeperiod, 
																	ofstream     &PROGRESS){
	clock_t  t1,t2;  
 	disp_info junk;  junk.al=junk.at=junk.D =0.0;

  cout <<endl<< "Tracking "<<numpart<<" particles..." <<endl;
  cout << "------------------------------------------------------------"<<endl;

	//TMP DEBUG - endtimes plotted
	//ofstream ENDTIME;
	//ENDTIME.open("endtimes.csv");
	
	PROGRESS <<"track"<<endl;
  for (int i=0; i<numpart; i++){
		if (ProgramAborted()){
			CPathline::WriteAllPathlines(Particles,numpart);
			break;
		}
	  else{
			cout<<endl<<" -Particle "<< Particles[i]->GetID()<< "...  ";
			//cout<<Particles[i]->GetLocation(0).x<<endl; 
			//cout<<Particles[i]->GetLocation(0).y<<endl; 
	    t1=clock();

			if (CPathline::pPollockGrid==NULL)
			{
				Particles[i]->Track(timeperiod,ADAPTIVE_TIME_STEP,G_MIN_TIME_STEP,G_MIN_STEP,false,junk,false);  //highest possible precision
			  //Particles[i]->Track(timeperiod,CONSTANT_SPACE_STEP,G_MIN_TIME_STEP,0.01,false,junk,false);
			}
			else 
			{
				Particles[i]->TrackPollock(timeperiod,100,pPollockGrid,false,true);//Particles[i]->TrackedForward()
			}
			t2=clock();
			cout <<float(t2-t1)/CLK_TCK << " seconds elapsed. ";
			PROGRESS << i+1 <<" " << numpart <<endl;
		
			//ENDTIME<<Particles[i]->GetLocation(0.0).x<<",";//start location (x) 
			//ENDTIME<<Particles[i]->GetLocation(0.0).y<<",";//start location (y)
			//ENDTIME<<Particles[i]->GetLastT() <<endl;      //tracking duration (T)

		}
	}
	PROGRESS<<"done"<<endl;
  //ENDTIME.close();
  cout << endl<<"             ...Done. "<< endl;
}
//*********************************************************************
void CPathline::WriteAllPathlines(CPathline **Particles,const int numpart){
  ofstream TRACK;
  TRACK.open("tracks.bna");
	TRACK.close();
	ofstream TRACK2;
  TRACK2.open("tracks.dat");
	TRACK2.close();
	ofstream TRACK3;
  TRACK3.open("mattstracks.csv");
	TRACK3<<"ID,t,x,y,s,step"<< endl;
  TRACK3.close();
  for (int i=0; i<numpart; i++){
		Particles[i]->CleanPath();
    Particles[i]->WriteOutput();
	}
}
