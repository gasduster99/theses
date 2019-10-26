#include "Grid.h"

//global members- visible to all gridding routines (would like to make this less messy)
//-------------------------------------------------------------------
const CSingleLayerABC *pPlotLayer;  
const CAquicludeABC   *pPlotAquiclude;
const CAnalyticElem   *pPlotElem;

struct GridProgWrite{
	ofstream             *PROG;
	int                  percentdone;
	int                  ps;
	int                  pe;
	bool                 writing;
};
GridProgWrite PlotProg;
//-------------------------------------------------------------------
void SetPlotLayer(const CSingleLayerABC *lay){
	pPlotLayer=lay;
}
/*******************************************************************
       COMPLEX-VALUED GRIDDING FUNCTIONS
*******************************************************************/
cmplex grdHeadStream		(const cmplex &z, const double &t){
	cmplex om;
	double head,B,K;
	pPlotLayer->GetHeadAndPotential(z,om,head,t);
	B=pPlotLayer->GetBase(z);
	K=pPlotLayer->GetCond(z);
	om=cmplex(max(head,B),om.imag());
	if ((K==0.0) || (om.real()<=B)){om=cmplex(SURFERMAGIC,om.imag());}
	return om;
}
//-------------------------------------------------------------------
cmplex grdQxQy					(const cmplex &z, const double &t){
	return conj(pPlotLayer->GetW(z,t));     //returns Qx+iQy      
}
//-------------------------------------------------------------------
cmplex grdGx				  	(const cmplex &z, const double &t){
	return conj(pPlotLayer->GetGx(z,t));    //returns dQx/dx+idQy/dx      
}
//-------------------------------------------------------------------
cmplex grdGxNum				  (const cmplex &z, const double &t){
	double step=0.000001;
	cmplex W1,W2; 
	W1=conj(pPlotLayer->GetW(z,t));
	W2=conj(pPlotLayer->GetW(z+step,t));
  //returns dQx/dx+idQy/dx 
	return (W2-W1)/step;
	//return cmplex ((W2.real()-W1.real())/step,(W2.imag()-W1.imag())/step);
}
//-------------------------------------------------------------------
cmplex grdGy				  	(const cmplex &z, const double &t){
	//returns dQx/dy+idQy/dy   
	return conj(IM*pPlotLayer->GetGx(z,t))-(IM*-pPlotLayer->GetLeakage(z,t,FROMTOPANDBOTTOM));//+pPlotLayer->GetCurl(z,t);
} 
//-------------------------------------------------------------------
cmplex grdGyNum				  (const cmplex &z, const double &t){
	double step=0.000001;
	cmplex W1,W2; 
	W1=conj(pPlotLayer->GetW(z,t));        //Qx+iQy
	W2=conj(pPlotLayer->GetW(z+IM*step,t));//Qx+iQy
  //returns dQx/dy+idQy/dy  
	return (W2-W1)/step;
	//return cmplex ((W2.real()-W1.real())/step,(W2.imag()-W1.imag())/step);
}
//-------------------------------------------------------------------
cmplex grdPotStream			(const cmplex &z, const double &t){return pPlotLayer->GetDischargePotential(z,t);}
//-------------------------------------------------------------------
cmplex grdElemPotStream (const cmplex &z, const double &t){return pPlotElem ->GetDischargePotential(z,t);}
//-------------------------------------------------------------------
cmplex grdvxvy					(const cmplex &z, const double &t){return pPlotLayer->GetVelocity2D(z,t);        }
//-------------------------------------------------------------------
cmplex grdEffVel			  (const cmplex &z, const double &t){
	disp_info disp;
	disp.al=1.0; disp.at=1.0; disp.D=0.0; 
	cmplex v=pPlotLayer->GetVelocity2D(z,t);
	return abs(v);//pPlotLayer->GetEffVelocity2D(z,t,disp)-v)*100/abs(v);           
}

//-------------------------------------------------------------------
cmplex grdDxxDxDyyDy			      (const cmplex &z, const double &t){
	//Numerical evaluation of Dxx/dx and Dyy/dy
	double step=0.000001;
	double Dxx1,Dyy1,Dxx2,Dyy2,vx,vy,aL,aT,vnrm;
	cmplex v1,v2;
	disp_info disp;
	aL=10.0; aT=7.0; disp.D=0.0; 

  v1=pPlotLayer->GetVelocity2D(z,t);
	vx=v1.real();
	vy=v1.imag();
	vnrm=abs(v1);
  Dxx1=					aL*vx*vx/vnrm+aT*vy*vy/vnrm+disp.D;
  Dyy1=					aT*vx*vx/vnrm+aL*vy*vy/vnrm+disp.D;

  v2=pPlotLayer->GetVelocity2D(z+step,t);
	vx=v2.real();
	vy=v2.imag();
	vnrm=abs(v2);
  Dxx2=					aL*vx*vx/vnrm+aT*vy*vy/vnrm+disp.D;
  
	v2=pPlotLayer->GetVelocity2D(z+IM*step,t);
	vx=v2.real();
	vy=v2.imag();
	vnrm=abs(v2);
  Dyy2=					aT*vx*vx/vnrm+aL*vy*vy/vnrm+disp.D;

	return cmplex((Dxx2-Dxx1)/step,(Dyy2-Dyy1)/step);
	
}
//-------------------------------------------------------------------
cmplex grdDxydxDxydy			      (const cmplex &z, const double &t){
	//Numerical evaluation of Dxy/dx and Dxy/dy
	double step=0.000001;
	double Dxy1,Dxy2,Dxy3,vx,vy,aL,aT,vnrm;
	cmplex v1,v2;
	disp_info disp;
	aL=10.0; aT=7.0; disp.D=0.0; 

  v1=pPlotLayer->GetVelocity2D(z,t);
	vx=v1.real();
	vy=v1.imag();
	vnrm=abs(v1);
  Dxy1=					(aL-aT)*vx*vy/vnrm+disp.D;

  v2=pPlotLayer->GetVelocity2D(z+step,t);
	vx=v2.real();
	vy=v2.imag();
	vnrm=abs(v2);
  Dxy2=					(aL-aT)*vx*vy/vnrm+disp.D;
  
	v2=pPlotLayer->GetVelocity2D(z+IM*step,t);
	vx=v2.real();
	vy=v2.imag();
	vnrm=abs(v2);
  Dxy3=					(aL-aT)*vx*vy/vnrm+disp.D;

	return cmplex((Dxy2-Dxy1)/step,(Dxy3-Dxy1)/step);
	
}
//-------------------------------------------------------------------
cmplex grdVmagDxDy				  (const cmplex &z, const double &t){
	double step=0.000001;
	cmplex v1,v2,v3; 
	v1=pPlotLayer->GetVelocity2D(z,t) ;
	v2=pPlotLayer->GetVelocity2D(z+step,t);
	v3=pPlotLayer->GetVelocity2D(z+IM*step,t);
  //returns d|v|/dx + i d|v|/dy
	return cmplex ((abs(v2)-abs(v1))/step,(abs(v3)-abs(v1))/step);
}
//-------------------------------------------------------------------
cmplex grdvderx				  (const cmplex &z, const double &t){
	double step=0.000001;
	cmplex v1,v2; 
	v1=pPlotLayer->GetVelocity2D(z,t) ;
	v2=pPlotLayer->GetVelocity2D(z+step,t);
  //returns dvx/dx + i dvy/dx
	return (v2-v1)/step;
	//return cmplex ((v2.real()-v1.real())/step,(v2.imag()-v1.imag())/step);
}
//-------------------------------------------------------------------
cmplex grdvdery				  (const cmplex &z, const double &t){
	double step=0.000001;
	cmplex v1,v2; 
	v1=pPlotLayer->GetVelocity2D(z,t) ;
	v2=pPlotLayer->GetVelocity2D(z+IM*step,t);
  //returns dvx/dy + i dvy/dy
	return (v2-v1)/step;
	//return cmplex ((v2.real()-v1.real())/step,(v2.imag()-v1.imag())/step);

}
/*******************************************************************
       REAL-VALUED GRIDDING FUNCTIONS
*******************************************************************/
double grdLeak					(const cmplex &z, const double &t){
	double leak;
	leak=pPlotLayer->GetLeakage(z,t,FROMTOPANDBOTTOM);
	if (leak==0.0){leak=SURFERMAGIC;}
	return leak;
}
//-------------------------------------------------------------------
double grdSaltwater     (const cmplex &z, const double &t){
	double B(pPlotLayer->GetBase(z));
	return GetSaltwaterElev(pPlotLayer->GetHead(z,t)-B,pPlotLayer->GetSeaLevel()-B,pPlotLayer->GetSaltwaterSG())+B;	
}
//-------------------------------------------------------------------
double grdPot				  	(const cmplex &z, const double &t){return pPlotLayer->GetDischargePotential(z,t).real();}
//-------------------------------------------------------------------
double grdCond					(const cmplex &z, const double &t){return pPlotLayer->GetCond(z);}
//-------------------------------------------------------------------
double grdBase					(const cmplex &z, const double &t){return pPlotLayer->GetBase(z);}
//-------------------------------------------------------------------
double grdThick			    (const cmplex &z, const double &t){return pPlotLayer->GetThick(z);}
//-------------------------------------------------------------------
double grdSatThick			(const cmplex &z, const double &t){return min(grdThick(z,t),pPlotLayer->GetHead(z,t));}
//-------------------------------------------------------------------
double grdLayerTop			(const cmplex &z, const double &t){return pPlotLayer->GetThick(z)+pPlotLayer->GetBase(z);}
//-------------------------------------------------------------------
double grdCurl   			  (const cmplex &z, const double &t){return pPlotLayer->GetCurl(z,0.0);}
//-------------------------------------------------------------------
double grdConduct			  (const cmplex &z, const double &t){return pPlotAquiclude->GetConductance(z);}

//-------------------------------------------------------------------
cmplex grdEffVelNum     (const cmplex &z, const double &t){
	double step=0.000001;
	double aL,aT,vnrm,vx,vy;
	disp_info disp;
	aL=10.0; aT=7.0; disp.D=0.0; 

	cmplex vel=grdvxvy(z,t); 
	vx=vel.real();
	vy=vel.imag();
	vnrm=abs(vel);

  double Dxx=     aL*vx*vx/vnrm+aT*vy*vy/vnrm+disp.D;
  double Dyy=     aT*vx*vx/vnrm+aL*vy*vy/vnrm+disp.D;
  double Dxy=(aL-aT)*vx*vy/vnrm              +disp.D;

	double h=grdSatThick(z,t);
	double dhdx=(grdSatThick(z+step,t)-h)/step;
	double dhdy=(grdSatThick(z+IM*step,t)-h)/step;

	cmplex tmp=grdDxxDxDyyDy(z,t);
	double DxxDx=tmp.real();
	double DyyDy=tmp.imag();
         tmp=grdDxydxDxydy(z,t);
	double DxyDx=tmp.real();
	double DxyDy=tmp.imag();
	return vel-
				 cmplex(DxxDx,     DyyDy)-
				 cmplex(DxyDy,     DxyDx)-
				 cmplex(Dxx/h*dhdx,Dyy/h*dhdy)-
				 cmplex(Dxy/h*dhdy,Dxy/h*dhdx); //Porosity not yet included
}
//*******************************************************************
void   GridDriver     (const CAquifer *Aq, 
											 GridCommands    G, 
											 double         &t, 
											 ofstream       &PROGRESS){

	int    NLayers(Aq->GetNumLayers());
	double numgrids(0);
	int    ngridded(0);

	if (pPlotAquiclude==NULL){G.conduct=false;}

  if (G.base)     {numgrids++; }
	if (G.cond)     {numgrids++; }
	if (G.conduct)  {numgrids++; }
	if (G.head)     {numgrids++; }
	if (G.leak)     {numgrids++; }
	if (G.saltH20)  {numgrids++; }
	if (G.pot)      {numgrids++; }
	if (G.stream)   {numgrids++; }
	if (G.satthick) {numgrids++; }
	if (G.thick)    {numgrids++; }
	if (G.QxQy)     {numgrids+=2;}
	if (G.GxGy)     {numgrids+=4;} 
	if (G.GxGynum)  {numgrids+=4;} 
	if (G.vxvy)     {numgrids+=2;}
	if (G.effvel)   {numgrids+=2;}
	if (G.elem)     {numgrids+=2;}
	if (G.DxxDx)    {numgrids+=2;}
	if (G.DxyDx)    {numgrids+=2;}
	if (G.VmagDer)  {numgrids+=2;}
	if (G.vder)     {numgrids+=2;}
	if (G.effvelnum){numgrids+=2;}

	numgrids*=NLayers;
	
	PlotProg.writing=true;
	PlotProg.PROG   =&PROGRESS; //Make this a pointer to an ofStream? (Clearwater compiler bug)

	char filename1[30];
	char filename2[30];

	PROGRESS << "grid"<<endl;
  ngridded=0;

	CAnalyticElem::SetCosmetic(true); 

	//CAnalyticElem::SetBranchCutLocus(true,cmplex(860,627.5)); 

	for (int L=0; L<NLayers; L++){
		pPlotLayer    =Aq->GetLayer(L);
		pPlotAquiclude=Aq->GetAquiclude(L);
		//-------Grid Head and Stream Function--------------------
		if ((G.head) && (G.stream)){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			if (NLayers==1){
				cGrid("Head.grd", "Stream.grd",G.W, G.res, grdHeadStream,t);
			}
			else {
				sprintf(filename1,"Head_%d.grd"  ,L+1);   
				sprintf(filename2,"Stream_%d.grd",L+1);
				cGrid  (filename1,filename2,G.W,G.res,	grdHeadStream,t);	 
			}
			ngridded+=2;
		}
		//-------Grid Qx and Qy ----------------------------------
		if (G.QxQy){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			if (NLayers==1){
				cGrid("Qx.grd","Qy.grd",G.W,G.res, grdQxQy,t);
			}
			else{
				sprintf(filename1,"Qx_%d.grd",L+1);   
				sprintf(filename2,"Qy_%d.grd",L+1);
				cGrid  (filename1,filename2,G.W,G.res,grdQxQy,t);	            
			}
			ngridded+=2;
		}
		//-------Grid Gxx, Gxy, Gyy and Gyx ----------------------------------
		if (G.GxGy){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"Qxdx_%d.grd",L+1);   
			sprintf(filename2,"Qydx_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdGx,t);	 
			ngridded+=2;

		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"Qxdy_%d.grd",L+1);   
			sprintf(filename2,"Qydy_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdGy,t);	 
			ngridded+=2;
		}
		//-------Grid Gxx, Gxy, Gyy and Gyx (NUMERICAL)---------------------------
		if (G.GxGy){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"Qxdx(num)_%d.grd",L+1);   
			sprintf(filename2,"Qydx(num)_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdGxNum,t);	 
			ngridded+=2;

		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"Qxdy(num)_%d.grd",L+1);   
			sprintf(filename2,"Qydy(num)_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdGyNum,t);	 
			ngridded+=2;
		}
		//-------Grid Leakage ------------------------------------
 		if (G.leak){   
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);
			if (NLayers==1){
				rGrid("Leakage.grd",G.W, G.res, grdLeak,t);
			}
			else{  
				sprintf(filename1,"Leakage_%d.grd",L+1);
				rGrid  (filename1,G.W, G.res, grdLeak,t);                        
			}
			ngridded++;
		}
		//-------Grid Saltwater elevation ------------------------
 		if (G.saltH20){   
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);                   
			sprintf(filename1,"Saltwater_%d.grd",L+1);
			rGrid  (filename1,G.W, G.res, grdSaltwater,t);
			ngridded++;
		}
		//-------Grid Potential ----------------------------------
		if (G.pot){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);  
			sprintf(filename1,"Potential_%d.grd",L+1);
			rGrid  (filename1,G.W, G.res, grdPot,t);
			ngridded++;
		}		
		//-------Grid Conductivity -------------------------------
		if (G.cond){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0); 
			sprintf(filename1,"Conductivity_%d.grd",L+1);
			rGrid  (filename1,G.W, G.res, grdCond,t);                
			ngridded++;	
		}	
		//-------Grid Base -------------------------------
		if (G.base){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);
			sprintf(filename1,"Base_%d.grd",L+1);  
			rGrid  (filename1,G.W, G.res, grdBase,t);                         
			ngridded++;
		}		
		//-------Grid Thickness -------------------------------
		if (G.thick){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);
			sprintf(filename1,"Thickness_%d.grd",L+1);  
			rGrid  (filename1,G.W, G.res, grdThick,t);                  
			ngridded++;
		}		
		//-------Grid Thickness -------------------------------
		if (G.thick){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);
			sprintf(filename1,"SaturatedThickness_%d.grd",L+1);  
			rGrid  (filename1,G.W, G.res, grdSatThick,t);                  
			ngridded++;
		}		
		//-------Grid Layer Top -------------------------------
		if (G.top){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);
			sprintf(filename1,"Top_%d.grd",L+1); 
			rGrid  (filename1,G.W, G.res, grdLayerTop,t);                 
			ngridded++;
		}		
		//-------Grid Conductance -------------------------------
		if (G.conduct){    
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);
			sprintf(filename1,"Conductance_%d.grd",L+1); 
			rGrid  (filename1,G.W, G.res, grdConduct,t);                   
			ngridded++;
		}		
		//-------Grid velocity  ----------------------------------
		if (G.vxvy){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"vx_%d.grd",L+1);   
			sprintf(filename2,"vy_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdvxvy,t);	 
			ngridded+=2;
		}
		//-------Grid effective velocity  ------------------------
		if (G.effvel){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"effective_vx_%d.grd",L+1);   
			sprintf(filename2,"effective_vy_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdEffVel,t);	 
			ngridded+=2;
		}
		//-------Grid Effective Velocity (Numerical)  ------------------
		if (G.effvelnum){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"EffVelx(num)_%d.grd",L+1);   
			sprintf(filename2,"EffVely(num)_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdEffVelNum,t);	 
			ngridded+=2;
		}
		//-------Grid Curl -------------------------------
		if (G.curl){    
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+1)/(double)(numgrids)*100.0);
			sprintf(filename1,"Curl_%d.grd",L+1); 
			rGrid  (filename1,G.W, G.res, grdCurl,t);                   
			ngridded++;
		}	
		//-------Grid dDxx/Dx ,dDyy/dy  ------------------
		if (G.DxxDx){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"DxxDx_%d.grd",L+1);   
			sprintf(filename2,"DyyDy_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdDxxDxDyyDy,t);	 
			ngridded+=2;
		}
		//-------Grid dDxy/Dx ,dDxy/dy  ------------------
		if (G.DxyDx){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"DxyDx(num)_%d.grd",L+1);   
			sprintf(filename2,"DxyDy(num)_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdDxydxDxydy,t);	 
			ngridded+=2;
		}
		//-------Grid d|v|/Dx ,d|v|/dy  ------------------
		if (G.VmagDer){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"dvmag-dx_%d.grd",L+1);   
			sprintf(filename2,"dvmag-dy_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdVmagDxDy,t);	 
			ngridded+=2;
		}
		//-------Grid vx/dx, vy/dx, vx/dy and vy/dy (NUMERICAL)---------------------------
		if (G.vder){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"vxdx(num)_%d.grd",L+1);   
			sprintf(filename2,"vydx(num)_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdvderx,t);	 
			ngridded+=2;

		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			sprintf(filename1,"vxdy(num)_%d.grd",L+1);   
			sprintf(filename2,"vydy(num)_%d.grd",L+1);
			cGrid  (filename1,filename2,G.W,G.res,grdvdery,t);	 
			ngridded+=2;
		}

	}
  //-------Grid Element -------------------------------
	if (G.elem){
		pPlotElem=G.pGrdElement;
		if (pPlotElem!=NULL){
		  PlotProg.ps=(int)((double)(ngridded)/(double)(numgrids)*100.0);
			PlotProg.pe=(int)((double)(ngridded+2)/(double)(numgrids)*100.0);
			cGrid("ElemPot.grd", "ElemStream.grd",G.W,G.res,grdElemPotStream,t);
			ngridded+=2;
		}
	}
  CAnalyticElem::SetCosmetic(false); 
	//---------------------
	PROGRESS<<100<<endl;
  PROGRESS<<"done"<<endl;
	PlotProg.writing=false;
}
//*******************************************************************
void WriteProgress(const double pctdone){
	if (PlotProg.writing){
		if ((int)((PlotProg.pe-PlotProg.ps)*pctdone+PlotProg.ps) > PlotProg.percentdone){
			PlotProg.percentdone=(int)((PlotProg.pe-PlotProg.ps)*pctdone+PlotProg.ps); 
			*PlotProg.PROG << PlotProg.percentdone<<endl;
		}
	}
}
/*******************************************************************
  cGrid (COMPLEXGRID)
	Plots time (t) and space (z)-dependent complex function (funct(z,t)) 
	to two *.grd files specified by filename1 and filename2
	The plotting window(W)  and resolution(resolution) are specified
*******************************************************************/
bool cGrid(       char		*filename1, 
					        char		*filename2, 
					 const  window	 W, 
					 const  int			 resolution, 
					        cmplex (*funct)(const cmplex &z, const double &t),
					 const  double	&time){

  double          xincrem((W.e-W.w)/resolution);          //x-increment between plot points
	double          yincrem((W.n-W.s)/resolution);          //y-increment between plot points
  double          minV1( ALMOST_INF),minV2( ALMOST_INF);  //minimum values (Real and imag)
	double          maxV1(-ALMOST_INF),maxV2(-ALMOST_INF);  //maximum values (Real and imag)
	double          xgrid,ygrid;                            //current x & y
  int             v1,v2,progcount;                        //counters for writing to file and screen
  ofstream        V1GRD,V2GRD;                            //output streams 
	bool            stopped(false);                         //true if stopped
	int             stoprow(resolution);                    //row at which gridding stopped
  cmplex        **plotpts=NULL;                           //array of plot points
  int             i,j;

	V1GRD.open(filename1);							
	V2GRD.open(filename2);
	if ((xincrem<=0) || (yincrem<=0)){ExitGracefully("cGrid::Improper window specified",RUNTIME_ERR);}
	if (resolution<=0)               {ExitGracefully("cGrid::Improper resolution specified",RUNTIME_ERR);}

	//allot dynamic memory, initialize to SurferMagic
  plotpts=new cmplex *[resolution+2];
  if (plotpts  ==NULL){ExitGracefully("cGrid::plotpoints",OUT_OF_MEMORY);}
  for(i=0; i<=resolution+1; i++){
    plotpts    [i]=new cmplex[resolution+2];
    if (plotpts[i]==NULL){ExitGracefully("cGrid::plotpoints",OUT_OF_MEMORY);}
		for(j=0; j<=resolution+1; j++){
			plotpts  [i][j]=SURFERMAGIC;
		}
	}

  cout << endl << "Gridding files " <<filename1 << " and "<< filename2;
	cout << " on "<< resolution << " by "<< resolution <<" window... "<< endl;
	cout << "------------------------------------------------------------"<<endl;

  v1=v2=progcount=0; 
 
	/*ofstream TMP;
	TMP.open("mapgrid.csv");
	TMP<<"x,y,phi,psi,vx,vy"<<endl;*/
	//generate grid
	//---------------------------------------------------------------------------------------------
	for (xgrid=W.w; xgrid<=W.e+(xincrem/4.0); xgrid+=xincrem){

		if (ProgramAborted()){stopped=true; xgrid=W.e+xincrem;}

		if (!stopped){
			for (ygrid=W.s; ygrid<=W.n+(yincrem/4.0); ygrid+=yincrem){

				//Get values at points
				plotpts[v1][v2]=funct(cmplex(xgrid,ygrid),time);

				//TMP DEBUG-----------------------------
				/*cmplex z=cmplex(xgrid,ygrid);
        cmplex om=pPlotLayer->GetDischargePotential(z,time);
				cmplex v=pPlotLayer->GetVelocity2D(z,time);
				om=exp(-z)/z;
				TMP << xgrid<<","<<ygrid<<","<<om.real()<<","<<om.imag()<<","<<v.real()<<","<<v.imag()<<endl;*/
				//TMP DEBUG-----------------------------

	      if (((plotpts[v1][v2].real()> 1e8) || 
					   (plotpts[v1][v2].real()<-1e8)) && 
						 (plotpts[v1][v2].real()!=SURFERMAGIC)){plotpts[v1][v2]-=   plotpts[v1][v2].real();}
	      if (((plotpts[v1][v2].imag()> 1e8) || 
					   (plotpts[v1][v2].imag()<-1e8)) && 
						 (plotpts[v1][v2].imag()!=SURFERMAGIC)){plotpts[v1][v2]-=IM*plotpts[v1][v2].imag();}

				//calc min & max
				if (plotpts[v1][v2].real()!=SURFERMAGIC){
					upperswap(maxV1,plotpts[v1][v2].real()); lowerswap(minV1,plotpts[v1][v2].real());}
				if (plotpts[v1][v2].imag()!=SURFERMAGIC){
					upperswap(maxV2,plotpts[v1][v2].imag()); lowerswap(minV2,plotpts[v1][v2].imag());}

				progcount++;
				WriteProgress((double)(progcount)/(double)(resolution+1)/(double)(resolution+1));

				v2++;if ((resolution>=16)&&(v2%(resolution/16)==0)&&((v1%(resolution/8))==0)) {cout<<"#";}
		  }
			v1++; v2=0; if ((resolution>=8)&&((v1%(resolution/8))==0)){cout<<endl;}
		}
		else {
			stoprow=v1;
		}
  }
	//TMP.close();

  //create ASCII files
	//---------------------------------------------------------------------------------------------
  V1GRD<< "DSAA"												 <<endl;	V2GRD<< "DSAA"												 <<endl;
  V1GRD<< resolution+1 <<" "<< stoprow+1 <<endl;	V2GRD<< resolution+1 <<" "<< stoprow+1 <<endl;
  V1GRD<< W.w					 <<" "<< W.e			 <<endl;	V2GRD<< W.w					 <<" "<< W.e			 <<endl;
  V1GRD<< W.s				 	 <<" "<< W.n			 <<endl;	V2GRD<< W.s					 <<" "<< W.n			 <<endl;
  V1GRD<< minV1				 <<" "<< maxV1		 <<endl;	V2GRD<< minV2				 <<" "<< maxV2		 <<endl; 

  //fill with array
  for (v2=0; v2<stoprow+1; v2++){
  for (v1=0; v1<resolution+1; v1++){
	V1GRD << plotpts[v1][v2].real() << " ";					V2GRD << plotpts[v1][v2].imag() << " ";
	} 
	V1GRD<<endl;														        V2GRD<<endl; 
  }
  //close files
  V1GRD.close();																	V2GRD.close();

  for(i=0; i<=resolution+1; i++){delete [] plotpts  [i];}
  delete [] plotpts;  

  cout <<"             ...Done. "<< endl;
	//ExitGracefully("TMP DEBUG-cGrid",BAD_DATA);
	return stopped;
}
/*******************************************************************
  rGrid (REALGRID)
	Plots time (t) and space (z)-dependent real-valued function (funct(z,t)) 
	to one *.grd file specified by filename
	The plotting window(W)  and resolution(resolution) are specified
*******************************************************************/
bool rGrid(       char		*filename, 
					 const  window	 W, 
					 const  int			 resolution, 
					        double (*funct)(const cmplex &z, const double &t),
					 const  double	&time){

  double          xincrem((W.e-W.w)/resolution);					//x-increment between plot points
	double          yincrem((W.n-W.s)/resolution);					//y-increment between plot points
  double          minV( ALMOST_INF);											//minimum value
	double          maxV(-ALMOST_INF);											//maximum value
	double          xgrid,ygrid;														//current x & y
  int             v1,v2,progcount;												//counters for writing to file and screen
  ofstream        GRD;																		//output stream
	bool            stopped(false);													//true if stopped
	int             stoprow(resolution);                    //row at which gridding stopped
  double        **plotpts=NULL;														//array of plot points                        
  int             i,j;

	GRD.open(filename);							

	if ((xincrem<=0) || (yincrem<=0)){ExitGracefully("rGrid::Improper window specified",RUNTIME_ERR);}
	if (resolution<=0)               {ExitGracefully("rGrid::Improper resolution specified",RUNTIME_ERR);}
	//allot dynamic memory, initialize to SurferMagic
  plotpts=new double *[resolution+2];
  if (plotpts  ==NULL){ExitGracefully("rGrid::Out of memory",OUT_OF_MEMORY);}
  for(i=0; i<=resolution+1; i++){
    plotpts    [i]=new double[resolution+2];
    if (plotpts[i]==NULL){ExitGracefully("rGrid::Out of memory",OUT_OF_MEMORY);}
		for(j=0; j<=resolution+1; j++){plotpts[i][j]=SURFERMAGIC;}
	}

  cout << endl << "Gridding file " <<filename;
	cout << " on "<< resolution << " by "<< resolution <<" window... "<< endl;
	cout << "------------------------------------------------------------"<<endl;

  v1=v2=progcount=0; 
 
	//generate grid
	//---------------------------------------------------------------------------------------------
	for (xgrid=W.w; xgrid<=W.e+(xincrem/4.0); xgrid+=xincrem){


		if (ProgramAborted())  {stopped=true; xgrid=W.e+xincrem; stoprow=v1;}

		if (!stopped){
			for (ygrid=W.s; ygrid<=W.n+(yincrem/4.0); ygrid+=yincrem){

				//Get values at points
				plotpts[v1][v2]=funct(cmplex(xgrid,ygrid),time);
	      if (((plotpts[v1][v2]> 1e8) || 
					   (plotpts[v1][v2]<-1e8)) &&
						(plotpts[v1][v2]!=SURFERMAGIC)){plotpts[v1][v2]=0.0;}

				if (plotpts[v1][v2]!=SURFERMAGIC){
					upperswap(maxV,plotpts[v1][v2]); 
					lowerswap(minV,plotpts[v1][v2]);
				}

				progcount++;
				WriteProgress((double)(progcount)/(double)(resolution+1)/(double)(resolution+1));

				v2++;if ((resolution>=16)&&(v2%(resolution/16)==0)&&((v1%(resolution/8))==0)) {cout<<"#";}
		  }
			v1++; v2=0; if ((resolution>=8)&&((v1%(resolution/8))==0)){cout<<endl;}
		}
  }

  //create ASCII file
	//---------------------------------------------------------------------------------------------
  GRD<< "DSAA"													<<endl; 
  GRD<< resolution+1 <<" "<< stoprow+1  <<endl; 
  GRD<< W.w					 <<" "<< W.e				<<endl; 
  GRD<< W.s					 <<" "<< W.n				<<endl; 
  GRD<< minV				 <<" "<< maxV				<<endl; 

  //fill with array
  for (v2=0; v2<stoprow+1; v2++){
		for (v1=0; v1<resolution+1; v1++){
			GRD << plotpts[v1][v2] << " ";	
		} 
		GRD<<endl;														     
  }
  //close file
  GRD.close();																	

  for(i=0; i<=resolution+1; i++){delete [] plotpts  [i];}
  delete [] plotpts;  

  cout <<"             ...Done. "<< endl;
	return stopped;
}
/*******************************************************************
  vGrid (VECTORGRID)
	Plots time (t) and space (z)-dependent numval-valued function (funct(z,t,v[],nv)) 
	to one *.grd file specified by filename
	The plotting window(W)  and resolution(resolution) are specified
*******************************************************************/
bool vGrid(       char		**filenames, 
					 const  window	 W, 
					 const  int			 resolution, 
					 const  int      numval,
					        void (*funct)(const cmplex &z, const double &t, double *v, const int nv),
					 const  double	&time){

  double          xincrem((W.e-W.w)/resolution);					//x-increment between plot points
	double          yincrem((W.n-W.s)/resolution);					//y-increment between plot points
  double         *minV;																		//minimum values
	double         *maxV;																		//maximum values
	double          xgrid,ygrid;														//current x & y
  int             v1,v2,progcount;												//counters for writing to file and screen
  ofstream       *GRD;																		//array of output streams
	bool            stopped(false);													//true if stopped
	int             stoprow(resolution);                    //row at which gridding stopped
  double       ***plotpts=NULL;														//array of plot points  
	double         *vals;                                   //temporary array of values at plot points
  int             i,j,n;

	GRD =new ofstream [numval];
	vals=new double   [numval];
	minV=new double   [numval];
	maxV=new double   [numval];
	for (n=0;n<numval;n++){
		GRD[n].open(filenames[n]);
		minV[n]= ALMOST_INF;
		maxV[n]=-ALMOST_INF;
	}

	if ((xincrem<=0) || (yincrem<=0)){ExitGracefully("vGrid::Improper window specified",RUNTIME_ERR);}
	if (resolution<=0)               {ExitGracefully("vGrid::Improper resolution specified",RUNTIME_ERR);}
	//allot dynamic memory, initialize to SurferMagic
  plotpts=new double **[resolution+2];
  if (plotpts  ==NULL){ExitGracefully("vGrid::Out of memory",OUT_OF_MEMORY);}
  for(i=0; i<=resolution+1; i++){
    plotpts    [i]=new double *[resolution+2];
    if (plotpts[i]==NULL){ExitGracefully("vGrid::Out of memory",OUT_OF_MEMORY);}
		for(j=0; j<=resolution+1; j++){
			plotpts    [i][j]=new double [numval];
			if (plotpts[i][j]==NULL){ExitGracefully("vGrid::Out of memory",OUT_OF_MEMORY);}
			for (n=0;n<numval; n++){
				plotpts[i][j][n]=SURFERMAGIC;
			}
		}
	}

  cout << endl << "Gridding multiple files (" <<filenames[0];
	cout << " through " <<filenames[numval-1]<<") on "<< resolution << " by "<< resolution <<" window... "<< endl;
	cout << "------------------------------------------------------------"<<endl;

  v1=v2=progcount=0; 
 
	//generate grid
	//---------------------------------------------------------------------------------------------
	for (xgrid=W.w; xgrid<=W.e+(xincrem/4.0); xgrid+=xincrem){


		if (ProgramAborted())  {stopped=true; xgrid=W.e+xincrem; stoprow=v1;}

		if (!stopped){
			for (ygrid=W.s; ygrid<=W.n+(yincrem/4.0); ygrid+=yincrem){

				//Get values at points
        funct(cmplex(xgrid,ygrid),time,vals,numval);
				for (n=0;n<numval;n++){
					plotpts[v1][v2][n]=vals[n];
					if (((plotpts[v1][v2][n]> 1e8) || 
							 (plotpts[v1][v2][n]<-1e8)) &&
							(plotpts[v1][v2][n]!=SURFERMAGIC)){plotpts[v1][v2][n]=0.0;}

					if (plotpts[v1][v2][n]!=SURFERMAGIC){
						upperswap(maxV[n],plotpts[v1][v2][n]); 
						lowerswap(minV[n],plotpts[v1][v2][n]);
					}
				}

				progcount++;
				WriteProgress((double)(progcount)/(double)(resolution+1)/(double)(resolution+1));

				v2++;if ((resolution>=16)&&(v2%(resolution/16)==0)&&((v1%(resolution/8))==0)) {cout<<"#";}
		  }
			v1++; v2=0; if ((resolution>=8)&&((v1%(resolution/8))==0)){cout<<endl;}
		}
  }

  //create ASCII file
	//---------------------------------------------------------------------------------------------
  for (n=0;n<numval;n++){
		GRD[n]<< "DSAA"													<<endl; 
		GRD[n]<< resolution+1 <<" "<< stoprow+1 <<endl; 
		GRD[n]<< W.w					<<" "<< W.e				<<endl; 
		GRD[n]<< W.s					<<" "<< W.n				<<endl; 
		GRD[n]<< minV[n]			<<" "<< maxV[n]	  <<endl; 

		//fill with array
		for (v2=0; v2<stoprow+1; v2++){
			for (v1=0; v1<resolution+1; v1++){
				GRD[n] << plotpts[v1][v2][n] << " ";	
			} 
			GRD[n]<<endl;														     
		}
		//close file
		GRD[n].close();																	
	}
  for(i=0; i<=resolution+1; i++){
		for(j=0; j<=resolution+1; j++){
			delete [] plotpts  [i][j];}
		delete [] plotpts  [i];}
  delete [] plotpts;  
  delete [] GRD;
	delete [] vals;
	
  cout <<"             ...Done. "<< endl;
	return stopped;
}


