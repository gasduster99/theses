#include "Flownode.h"

//************************************************************************
//                           Constructors
//************************************************************************
CFlowNode::CFlowNode(){}
//------------------------------------------------------------------------
CFlowNode::CFlowNode(cmplex z){
  znode=z;
  Qin=0;
  Qout=0;
	NumInElems=0;
	isonlake=false;
  LakeElem=NULL;
	for (int i=0; i<MAX_FLOW_INPUTS; i++){
	  InElems[i]=NULL;
		Inflow [i]=NULL;
		//InNodes[i]=NULL;
	}
	OutElem=NULL;
  //OutNode=NULL;
}
void CFlowNode::Destroy(){
	if (globaldebug){cout<<"DESTROYING ALL FLOW NODES"<<endl;}
	for(int i=0; i<nflownodes; i++){
		delete NodeArray[i];
	}
}
//************************************************************************
//          STATIC MEMBER INITIALIZATION
//************************************************************************
CFlowNode *CFlowNode::NodeArray[]={NULL};
int        CFlowNode::nflownodes=0;

//************************************************************************
//         Accessor/modifier functions
//************************************************************************
cmplex            CFlowNode::GetZ         () const     {return znode;}
CAnalyticElem    *CFlowNode::GetOutElem   () const     {return OutElem;}
CAnalyticElem    *CFlowNode::GetInElem    (int i) const{return InElems[i];}
int               CFlowNode::GetNumInElems() const     {return NumInElems;}
double            CFlowNode::GetHeadWater () const     {return Qin;}
void              CFlowNode::SetHeadWater (double Q)   {Qin=Q;Qout+=Qin;}
//************************************************************************
//         BuildFlowNetwork
//************************************************************************	
void CFlowNode::BuildFlowNetwork(){
	for (int i=0; i<nflownodes; i++){
	  NodeArray[i]->Adjust();
	}
}
//------------------------------------------------------------------------
void CFlowNode::Adjust(){
	int junk;
  if ((OutElem==NULL) && (NumInElems==0) && (Qin==0)){
    //unused lake node, should delete
 	}
	else if ((OutElem==NULL) && (isonlake)){
		SpecifyOutflowElem(LakeElem);
	}
	else if (isonlake) {
	  AddInflowElem(LakeElem,junk);
	}
}
//************************************************************************
//         AddInflowElem, OutFlowElem
//************************************************************************	
void CFlowNode::SpecifyOutflowElem(CAnalyticElem *elem){OutElem =elem;}
//------------------------------------------------------------------------
void CFlowNode::SpecifyLakeElem		(CAnalyticElem *elem){LakeElem=elem;isonlake=true;}
//------------------------------------------------------------------------
void CFlowNode::AddInflowElem			(CAnalyticElem *elem,int &ID){
 InElems[NumInElems]=elem;
 Inflow[NumInElems]=0.0;
 ID=NumInElems; 
 NumInElems++;
 //cout << "******* "<<ID<<endl;
}
//************************************************************************
//        Add Flow Node /Add Headwater
//************************************************************************
CFlowNode *CFlowNode::AddFlowNode(CAnalyticElem *elem, cmplex z,flowdir dir,double Q,int &ID){
	bool nodeexists(false);
	int  index(0);

	for (int i=0; i<nflownodes; i++){
		if (abs(NodeArray[i]->GetZ()-z)<REALSMALL){
			nodeexists=true; index=i;
		} //node exists
	}
  if (!nodeexists){                                  //if node doesn't exist, create new node
		if (nflownodes>=MAX_FLOWNODES){ExitGracefully("Add Flow Nodes: Too many river/lake nodes",TOO_MANY);}
	  NodeArray[nflownodes]=new CFlowNode(z);
		index=nflownodes;
		nflownodes++;
	}
	if (elem!=NULL){
		if      (dir==UPSTREAM  ){                       //water body is upstream
			NodeArray[index]->AddInflowElem     (elem,ID);
		}
		else if (dir==DOWNSTREAM){                       //water body is downstream
			NodeArray[index]->SpecifyOutflowElem(elem);
		}
		else if (dir==FLAT      ){                       //water body dir is unknown (lake)
			NodeArray[index]->SpecifyLakeElem   (elem);
		}
	}
  else if (elem==NULL){                              //water body is user-defined inflow (headwater)
		NodeArray[index]->SetHeadWater(Q);
	}

	return NodeArray[index];
}
//************************************************************************
//           Public Member Functions
//************************************************************************
double CFlowNode::GetOutflow(double t){return Qout;}
//------------------------------------------------------------------------
void CFlowNode::UpdateFluxes(double cumflow, int ID,  double t){ 
  if ((ID<0) || (ID>=nflownodes)){ExitGracefully("CFlowNode::UpdateFluxes: bad node ID",RUNTIME_ERR);}
	Inflow[ID]=cumflow;
	Qout=Qin;
	for (int i=0; i<NumInElems; i++){Qout+=Inflow[i];}
}
