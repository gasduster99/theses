#ifndef ISOTHERM_H
#define ISOTHERM_H

#include "RxNLibraryInclude.h"

/****************************************************
 *  Class CIsotherm
 *  Isotherm Data Abstraction
 ***************************************************/
class CIsotherm{
 private:/*----------------------------------------*/

 public: /*----------------------------------------*/
	CIsotherm(){}
	virtual ~CIsotherm(){}
	virtual bool   IsEquilibrium    () const {return true;}
	virtual double TranslateToSorbed(const double &C) const {return 0.0;}
	virtual double GetRetardation   (const double &pb, const double &n, const double &C) const {return 1.0;}
};

/****************************************************
 *  Class CLinearIsotherm
 ***************************************************/
class CLinearIsotherm:public CIsotherm{
 private:/*----------------------------------------*/
	
	 double Kd; //Linear isotherm coeff. (L/kg)

 public:/*----------------------------------------*/
	CLinearIsotherm(const double K){
		Kd=K;
	}
	double TranslateToSorbed(const double &C) const {
		return Kd*C; 
	}
	double GetRetardation   (const double &pb, const double &n, const double &C) const {
		return 1.0+pb/n*Kd;
	}
};
/****************************************************
 *  Class CFreundlichIsotherm
 ***************************************************/
class CFreundlichIsotherm:public CIsotherm{
 private:/*----------------------------------------*/

	double K;   //Freundlich isotherm coeff. (L/kg)
	double a;		//Freundlich exponent

 public:/*----------------------------------------*/

	CFreundlichIsotherm(const double Kf, const double exponent){
		K=Kf;
		a=exponent;
	}
	double TranslateToSorbed(const double &C) const {
		return K*pow(C,a);
	}
	double GetRetardation   (const double &pb, const double &n, const double &C) const {
		return 1.0+pb/n*a*K*pow(C,a-1.0);
	}
};
/****************************************************
 *  Class CLangmuirIsotherm
 ***************************************************/
class CLangmuirIsotherm:public CIsotherm{
 private:/*----------------------------------------*/

	double K;   //Langmuir isotherm coeff. (L/kg)
	double a;   //
	double S;   //Aqueous Solubility (mg/L)

 public:/*----------------------------------------*/

	CLangmuirIsotherm(const double Kf, const double solubility, const double exponent){
		K=Kf;
		a=exponent;
		S=solubility;
	}
	double TranslateToSorbed(const double &C) const {
		return (C*K*S)/(1.0+K*C);
	}
	double GetRetardation   (const double &pb, const double &n, const double &C) const {
		return 1.0+pb/n*(K*S)/pow(1.0+K*C,2);
	}
};
/****************************************************
 *  Class CPDMIsotherm
 ***************************************************/
class CPDMIsotherm:public CIsotherm{
 private:/*----------------------------------------*/

	double qmax;          //maximum adsorption capacity (cm^3/g-1)
	double d;             //adsorption isotherm exponent

 public:/*----------------------------------------*/

	CPDMIsotherm(const double qmax_, const double exponent){
		qmax=qmax_;
		d=exponent;
	}
	double TranslateToSorbed(const double &C) const {
		return 0.0;
		//return qmax*exp(-pow((R*T*log(I.S/C))/(I.b1*I.b2*I.E0,I.a));
	}
	double GetRetardation   (const double &pb, const double &n, const double &C) const {
		return 1.0;//TMP DEBUG
	}
};

#endif 