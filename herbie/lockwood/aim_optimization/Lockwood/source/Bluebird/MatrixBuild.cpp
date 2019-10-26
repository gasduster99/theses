//MatrixBuilder.cpp
//Builds and solves Explicit Matrix for AEM solution
#ifndef MATRIX_BUILD
#define MATRIX_BUILD

#include "AnalyticElem.h"
#include "MatrixInclude.h"

/***********************************************************************
				BuildMatrix
************************************************************************
Builds explicit matrix for solution of Analytic element coefficients
-----------------------------------------------------------------------*/

void BuildMatrix(CAnalyticElem **pElems,        //array of pointers to elements/element segments 
								 int            *segindices,    //indices of element segments
								 int             numelems,      //number of elements
								 CAnalyticElem  *OtherElems,    //aggregate pointer representing all outer elements
								 const double    t){            //times

	int             row,col;                //counters
	int             g,i,j,m,n,s;            //counters
	double        **AA,*BB,*SS;             //Master matrix, RHS, and solution vector
	int						 *DOF,totalDOF;           //degrees of freedom (indexed by elem ind), total degrees of freedom
	int            *FirstCol;								//first column which has local influence of this element
	CAnalyticElem **GivenElems;             //array of "given" elements (elements that dont need to be solved for
	int            *GivenSegs, numgiven;    //array of given segment indices, number of given segments
	int            *SegsByCol,*ElemsByCol;  //indexes of elements/segments indexed by row/column in matrix
	double          alpha;
	MatrixInfo      MI;                     
	double          uPhi [MAXCONTROLPTS];		//Unit potential influences at ctrl point of current(col) element/elem seg
  cmplex          uQ   [MAXCONTROLPTS];		//Unit discharge influences at ctrl point of current(col) element/elem seg


	DOF       =new int            [numelems];
	FirstCol  =new int            [numelems];
	GivenElems=new CAnalyticElem *[numelems]; // obviously, a maximum estimate
	GivenSegs =new int            [numelems];
	if (GivenSegs==NULL){ExitGracefully("BuildMatrix: Out of memory",OUT_OF_MEMORY);}

	totalDOF=0;
	numgiven=0;

	//for each element/element segment, get # of degrees of freedom, identify given elems/elem segs
	for (i=0; i<numelems; i++){
		//cout << "Processing "<< pElems[i]->GetName() << " (segment "<<segindices[i]<<")"<< endl;

		FirstCol[i]=totalDOF;

		if (segindices[i]==NOT_A_STRING){DOF[i]=pElems[i]->GetDegreesOfFreedom   ();             }
		else                            {DOF[i]=pElems[i]->GetSegDegreesOfFreedom(segindices[i]);}

		totalDOF+=DOF[i];
		if (pElems[i]->IsGiven()){
			GivenElems[numgiven]=pElems[i];
			GivenSegs [numgiven]=segindices[i];
			numgiven++;
		}
	}

	//build referential index of system by column (element indices and segment indices)
	//create memory for solver matrix

	ElemsByCol= new int     [totalDOF];   
	SegsByCol = new int     [totalDOF];
	AA        = new double *[totalDOF];
	BB        = new double  [totalDOF];
	SS        = new double  [totalDOF];	
	if (SS==NULL){ExitGracefully("BuildMatrix: Out of memory(2)",OUT_OF_MEMORY);}
	col=0;
	for (i=0; i<numelems; i++){
		for (j=0;j<DOF[i];j++){
			ElemsByCol[col]=i;
			SegsByCol [col]=segindices[i];
			col++;
		}
	}
	for (i=0;i<totalDOF;i++){
		AA[i]=new double [totalDOF];
		for (j=0; j<totalDOF; j++){
			AA[i][j]=0.0;
		}
		BB[i]=0.0;
		SS[i]=0.0;
	}

	//Go through, element by element (group of rows by group of rows), building matrix
	row=col=0;
	for (i=0; i<numelems; i++){ //for each element/element segment

		//get local unit influence coefficients, control points, coefficents, etc. 
		if (segindices[i]==NOT_A_STRING){pElems[i]->GetMatrixBuildInfo   (               MI);}
		else                            {pElems[i]->GetSegMatrixBuildInfo(segindices[i], MI);}				

		for (s=0; s<DOF[i]; s++){
			
			//cout << "building matrix...row "<< row <<endl;
			
			//go thru column by column
			for (col=0; col<totalDOF; col++){		
				
				//Fill AA Matrix entries-----------------------------------------------------

				//Matrix Coefficients from this element
				if ((col>=FirstCol[i]) && (col<FirstCol[i]+DOF[i])){
			  	n=col-FirstCol[i];
				  AA[row][col]=0.0;
					for (m=0; m<MI.nctrl; m++){
						AA[row][col]+=MI.unit[s][m]*MI.unit[n][m]; 
						//AA[row][col]=s;
					}
				}

				//Matrix Coefficients from other elements
				else { 
					n=col-FirstCol[ElemsByCol[col]];
					if (SegsByCol[col]==NOT_A_STRING){
						pElems[ElemsByCol[col]]->GetUnitInfluences   (               n,MI.zctrl,MI.nctrl,uPhi,uQ,t);
					}
					else                             {
						pElems[ElemsByCol[col]]->GetSegUnitInfluences(SegsByCol[col],n,MI.zctrl,MI.nctrl,uPhi,uQ,t);
					}
          AA[row][col]=0.0;
					for (m=0; m<MI.nctrl; m++){
						alpha=MI.phiCoeff*uPhi[m]+MI.QxCoeff*uQ[m].real()+MI.QyCoeff*uQ[m].imag();
						AA[row][col]-=alpha*MI.unit[s][m];
						//AA[row][col]=0.0;
					}
				}
			}//end column (col)

			//Fill BB Vector (RHS) entries-----------------------------------------------------
			BB[row]=0.0;

			//get RHS from outer elements (non-explicit elements)
			if (OtherElems!=NULL){
				for (m=0; m<MI.nctrl; m++){
					if (MI.phiCoeff!=0){
						uPhi[m]=OtherElems->GetDischargePotential(MI.zctrl[m],t).real();
					}
					if ((MI.QxCoeff!=0) || (MI.QyCoeff!=0)){
						uQ  [m]=OtherElems->GetW                 (MI.zctrl[m],t);
					}
					alpha=MI.phiCoeff*uPhi[m]+MI.QxCoeff*uQ[m].real()+MI.QyCoeff*uQ[m].imag();
					BB[row]+=alpha*MI.unit[s][m];
				}
			}

			//get RHS from given elements
			for (g=0;g<numgiven;g++){
				
				for (m=0; m<MI.nctrl; m++){
					if (MI.phiCoeff!=0){
						uPhi[m]=GivenElems[g]->GetDischargePotential(MI.zctrl[m],t).real();
					}
					if ((MI.QxCoeff!=0) || (MI.QyCoeff!=0)){
						uQ  [m]=GivenElems[g]->GetW                 (MI.zctrl[m],t);
					}
					alpha=MI.phiCoeff*uPhi[m]+MI.QxCoeff*uQ[m].real()+MI.QyCoeff*uQ[m].imag();
					BB[row]+=alpha*MI.unit[s][m];
				}
			} //end given (g)

			g=4;
			for (m=0; m<MI.nctrl; m++){
				//constants from element 
				BB[row]+=MI.elemrhs[m]*MI.unit[s][m];
			}
			row++;
		} //end DOF (s)
	}	//end element (i) 



	//Print Matrix
	/*ofstream TED;
	TED.open("explicitmatrix.txt");
  cout <<endl<< "A matrix"<<endl;
  for(s=0; s<totalDOF; s++){ 
		cout<<" | ";  
		for(n=0; n<totalDOF; n++){ 
			if(fabs(AA[s][n])>REALSMALL){cout.width(12);cout.precision(4); cout<<AA[s][n]<<" ,";} 
      else{cout.width(12); cout<<0.0<<" ,";}
		}
		cout<<" | "<<BB[s]<<endl;
	}  
   cout << "b vector"<<endl; for(n=0; n<totalDOF; n++){cout <<BB[n]<<" ";} 
  cout <<endl;
	TED.close();*/


	//Solve Matrix
	if (!SVD(AA,BB,SS,totalDOF)){
		ExitGracefully("BuildMatrix: SVD routine failed",SINGMAT);}

	/*int numiter(300);
	PCG(AA,BB,SS,totalDOF, REALSMALL, numiter, PRECONDITIONED_CG);*/

  cout << "sol vector: "<<endl; for(n=0; n<totalDOF; n++){cout <<SS[n]<<" ";} cout <<endl;

	//go through and assign values to coefficents 
	double tmp[MAXDOF];
	for (i=0; i<numelems; i++){
		for (s=0; s<DOF[i]; s++){
			tmp[s]=SS[FirstCol[i]+s];
		}		
		if (segindices[i]==NOT_A_STRING){
		  pElems[i]->SetCoeff(tmp);
		}
		else {
			pElems[i]->SetSegCoeff(segindices[i],tmp);
		}
	}

	delete [] SegsByCol;
	delete [] ElemsByCol;
	for (i=0;i<totalDOF;i++){delete [] AA[i];} delete [] AA;
	delete [] BB;
	delete [] DOF;
	delete [] FirstCol;
	delete [] GivenElems;
	delete [] GivenSegs;
}

#endif
//Functions I must now build (this is the real beast)

/*
int  GetDegreesOfFreedom   ();
	//easy
void GetMatrixBuildInfo    (int &nctrl,cmplex *zctrl,double **unitmat, double *elemrhs, double &phiCoeff,double &Qxcoeff,double &QyCoeff);
	//neccessary for every BC
	//nctrl easy, zctrl easy, phiCoeff, QxCoeff, QyCoeff easy
	//unit hard, but already done for GenSolve

void GetUnitInfluences     (const int n,const cmplex *zctrl,const int NumCtrl,double *uPhi, complex *uQ);
	//save coeff
	//set pJumpCoeff[n] to 1, others to zero
	//set FFcoeff
	//for each zctrl,
		//uPhi[m]=GetDischargePot(zctrl,t).real()
		//uQ[m]=GetW(zctrl,t);
	//set JumpCoeff Back
	//Set FFCoeff Back

void SetCoeff              (double *coeff);
	//easy

//strings only (and an AE virtual)
int  GetSegDegreesOfFreedom   (const int i);
void GetMatrixBuildInfo       (const int i, int &nctrl,cmplex *zctrl,double **unitmat,double *elemrhs, double &phiCoeff,double &Qxcoeff,double &QyCoeff);
void GetSegUnitInfluences     (const int i,const int n,const cmplex *zctrl,const int NumCtrl,double *uPhi, complex *uQ);
void SetSegCoeff              (const int i,double *coeff);
*/
