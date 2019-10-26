#ifndef FLOWNODE_H
#define FLOWNODE_H

#include "BluebirdLibrary.h"
#include "AnalyticElem.h"

enum flowdir {UPSTREAM, DOWNSTREAM, FLAT};
const int    MAX_FLOWNODES=  50000;    // maximum number of flow nodes in system
const int    MAX_FLOW_INPUTS= 4;        // maximum number of inputs to a flow node

/****************************************************
  Class CFlowNode
  Surface Water Juncture Data Abstraction
  used as connector between surface water features and AEM analogues
****************************************************/
/*
   * in[0]   *in[1]


          *   this node


          *   downstream
*/
class CFlowNode {
 private:/*----------------------------------------------------------*/
	//Static Variables:
	static CFlowNode *NodeArray[MAX_FLOWNODES]; //array of all flow nodes in system
	static int        nflownodes;                //number of flow nodes

	//Member Variables:
	CAnalyticElem    *InElems[MAX_FLOW_INPUTS];  //pointers to elements which provide water
	//CFlowNode        *InNodes[MAX_FLOW_INPUTS];  //pointer to upstream nodes
	double            Inflow [MAX_FLOW_INPUTS];  //current input flow 
	int               NumInElems;						     //number of input elements

	CAnalyticElem    *OutElem;								   //pointer to element which recieves water
	//CFlowNode        *OutNode;                   //pointer to downstream node
	double            Qout;                      //Cumulative base flow to node
	double            Qbranch;                   //Cumulative base flow downstream of node

  cmplex            znode;									   //location of flow junction

	double            Qin;										   //given headwater (for non-element input)

	bool              isonlake;
	CAnalyticElem    *LakeElem;
 
	//Private Accessor Functions
	CAnalyticElem    *GetOutElem        () const;
	CAnalyticElem    *GetInElem         (int i) const;
  int               GetNumInElems     () const;
	double            GetHeadWater      () const;
  cmplex            GetZ              () const;

	//Private Manipulator Functions
	void              SpecifyOutflowElem(CAnalyticElem *elem);
  void              SpecifyLakeElem   (CAnalyticElem *elem); 
  void              AddInflowElem     (CAnalyticElem *elem,int &ID);
  void              SetHeadWater      (double Q);
  
	//Private Member Functions
	void              Adjust            ();

 public:/*----------------------------------------------------------*/
  //Constructors
	CFlowNode();
	CFlowNode(cmplex z);
 
	//Static Member Functions
	static CFlowNode *AddFlowNode (CAnalyticElem *elem, cmplex z, flowdir dir,double Q,int &ID);
//	static CFlowNode *AddFlowLink
	static void       BuildFlowNetwork();
	static void       Destroy();

	//Member functions
	double            GetBranchFlow(double t);
	double            GetOutflow  (double t);
	void              UpdateFluxes(double cumflow, int ID,  double t);
};

/*class CFlowNode{
	
};
class CFlowLink{
	CAnalyticElem **pContibutingElems;
  CFlowNode     **pUpstream;
	CFlowNode     * pDownstream;
}


class CFlowNetwork {
 private:
  
};*/

#endif