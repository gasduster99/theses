//ParentDaughter.h

#ifndef PARENTDAUGHER_H
#define PARENTDAUGHER_H
/**************************************************************
	Class CParentDaughter
	Data Abstraction for Parent Daghter sequential decay
**************************************************************/
class CParentDaughter:public CReactionScheme{
 private:/*--------------------------------------------------*/
	int     *Indices;    //array of indices [NS] of decay species
	int      NS;         //number of species

	double  *K;          //Array of decay coefficients 
	double  *y;          //Array of yield coefficients (species i to child (i+1))
	
	int      soiltype;   //global index of applicable soil

	double  *C;          //Caq [NS]
	double  *S;          //Csorbed [NS]
	double  *a;          //transformed concentrations [NS]
	double  *Cneg;       //negative concentrations
	double  *Sneg; 

	double               Validate(const double &pb, const double &n) const;

 public:/*---------------------------------------------------*/
	//Constructors
	CParentDaughter();
  CParentDaughter::CParentDaughter(int     NumCations, 
																   int    *indices, 
																	 double *decay_coeff,
																	 double *yield_fractions,
																	 int     soiltype_index);
 ~CParentDaughter();

  static CParentDaughter *CParentDaughter::Parse(ifstream &input, int &l);
  
	void Initialize();
	void React(Ironclad1DArray  Caq,
						 Ironclad1DArray  Caq_end,
						 Writeable1DArray Caq_new,
		         Ironclad1DArray  Cs,
						 Ironclad1DArray  Cs_end,
						 Writeable1DArray Cs_new,	
						 const int        soil,
						 const double    &n,  
						 const double    &tstep);
};
void TestParentDaughter();

#endif