//TimeSeries.h
#ifndef TIMESERIES_H
#define TIMESERIES_H

#include "MasterInclude.h"

const int MAX_TIMESERIES=100;
const int MAX_DISCRETE_TIMES=1000;
/****************************************************
 *  Class CTimeSeries
 *  Time Series Data Abstraction
 ***************************************************/
class CTimeSeries{
 protected:/*----------------------------------------------------------*/
	static CTimeSeries  *pAllTimeSeries[MAX_TIMESERIES];
	static int           TotalNumTimeSeries;

	double starttime;
	double endtime;
	double novalue;    //value returned when time series is off (zero by default)

 public:/*-------------------------------------------------------------*/

  CTimeSeries();
	CTimeSeries(double startt, double endt, double no_val);
  virtual ~CTimeSeries();

	virtual double   GetValue        (const double &t) const;

	static void DestroyAllTimeSeries();
};
/****************************************************
 *  Class CHeavisideTimeSeries
 *  Heaviside Time Series Data Abstraction
 ***************************************************/
class CHeavisideTimeSeries: public CTimeSeries {
 private:/*----------------------------------------------------------*/

	double  val;  //pulse value

 public:/*-----------------------------------------------------------*/

	CHeavisideTimeSeries();
	CHeavisideTimeSeries(const double    val,
											 const double    startt,
											 const double    endt,
										 	 const double    novalue);
  ~CHeavisideTimeSeries();

	static CHeavisideTimeSeries *Parse(ifstream &input, int &l);

	double GetValue(const double &t) const;

};
/****************************************************
 *  Class CDiscreteTimeSeries
 *  Piecewise linear Time Series Data Abstraction
 ***************************************************/
class CDiscreteTimeSeries: public CTimeSeries {
 private:/*----------------------------------------------------------*/

	double *ts;  //time values
	double *vs;  //series values
	int     N;   //number of values

 public:/*-----------------------------------------------------------*/
	//Constructors
	CDiscreteTimeSeries();
	CDiscreteTimeSeries(const double   *times,
										  const double   *vals,
											const int       numvals,
											const double    novalue);
  ~CDiscreteTimeSeries();

	//Static Member Functions
	static CDiscreteTimeSeries *Parse(ifstream &input, int &l);

	//Accessor Functions
	double   GetValue(const double &t) const;

};
/****************************************************
 *  Class CPolynomialTimeSeries
 *  Polynomial Time Series Data Abstraction
 ***************************************************/
class CPolynomialTimeSeries: public CTimeSeries {
 private:/*----------------------------------------------------------*/
	//Member Variables
	double *a;  //coefficients
	int     N;  //number of terms

 public:/*-----------------------------------------------------------*/
	//Constructors
	CPolynomialTimeSeries();
	CPolynomialTimeSeries(const double   *coeff,
												const int       numterms,
												const double    startt,
												const double    endt,
												const double    novalue);
  ~CPolynomialTimeSeries();

	//Static Member Functions
	static CPolynomialTimeSeries *Parse(ifstream &input, int &l);

	//Accessor Functions
	double   GetValue      (const double &t) const;

};

#endif