#ifndef STAGE_H
#define STAGE_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"
#include "Flownode.h"

/************************************************************************
  Class CStage
  Analytic Stage-Specified Leaky Bottom Linesink Element Data Abstraction
  parents: CAnalyticElem, CStringElem
************************************************************************/

class CStage : public CStringElem {
 protected:/*----------------------------------------------------------*/
  double           *depth;        //pointer to array of depths at endpoints
  double          **depthctrl;    //array of depths at all control points
  double           *head;         //pointer to array of Head at endpoints
  double          **headctrl;     //array of heads at all control points

  double						conductivity; //conductivity of leaky bottom
	double						thickness;    //thickness of leaky bottom
	double						width;        //width of stream
  double            resistance;   //resistance of leaky bottom (thickness/cond)
				
	CFlowNode       **flownodes;    //array of pointers to surface flow nodes
	int              *outflow_ID;   //array of surface flow IDs
	flowdir           direction;    //direction of surface flow

  bool            **connected;    //2D array [i][m] - true if hydraulically connected at this ctrl point [i][m]

	bool							fresh;        //true if not yet solved

	//Member Functions
  void            WriteXSection(double t) const;

 public:/*-------------------------------------------------------------*/
  //Constructor
  CStage();
  CStage(char									 *Name,        //name of element
		     const CSingleLayerABC *pLay,        //layer of element
				 const cmplex					 *points,      //vertices of element polyline
				 const int							NumOfLines,  //number of polyline segments
				 const double					 *heads,       //specified heads at element vertices
				 const double					 *depths,      //specified depth at element vertices
		     const double						bot_cond,    //conductivity of leaky layer 
				 const double						bot_thick,   //thickness of leaky layer
				 const double						bot_width,   //width of stream
				 const int							prec);       //precision of element
 ~CStage();
  
  static CStage *Parse         (ifstream &input, int &l, CSingleLayerABC *pLay, char * Name);
	
	//Accessor Functions:
  double          GetRiverHead (double dist) const;
  virtual void    GetDepth     (int i, const double inflow, const double t);            //blank for stage

	//Manipulator functions (inherited from CStringElem)
	void            SetSegCoeff  (const int i, double *coeff);

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself  (double &change, double &objective,const double t);
	double          GetMaxError  (const double &t) const;
	void            WriteOutput  (const double &t) const;

};

/************************************************************************
 *  Class CVariableStage
 *  Analytic Stage-Specified Leaky Bottom Linesink Element
 *  attached to 1-D mannings flow model
 *  ****IN DEVELOPMENT****
 *  parents: CAnalyticElem, CStringElem, CStage
 ***********************************************************************/
class CVariableStage : public CStage {

 private:/*----------------------------------------------------------*/
  //double          **flowctrl;    //pointer to array of Flows in river. 
  double           *pRunoff;      
	double           *pSlope;
	double           *pBase;
	double          **basectrl;
  double            roughness;

 public:/*-----------------------------------------------------------*/
  //Constructors
  CVariableStage();
  CVariableStage(char			 *Name, 
								 CSingleLayerABC *pLay, 
								 cmplex		 *points,
								 int				NumOfLines,
								 double		 *base_elevs, 
								 double		 *zeros, 
		             double		 *runoff,
								 double			bot_cond, 
								 double			bot_thick, 
								 double			bot_width, 
								 double			mannings_n, 
								 int				prec);

  static CVariableStage *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Member Functions (Unique):
  void   GetDepth(int i, const double inflow, const double t);            //changes the head, depth in the river

  //Member Functions (Inherited from CAnalyticElem via CStringElem via stageelem):
};
#endif
