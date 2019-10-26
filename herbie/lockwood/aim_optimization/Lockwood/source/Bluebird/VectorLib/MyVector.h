//My Vector.h
#ifndef MYVECTOR_H
#define MYVECTOR_H

#include <math.h>
#include <complex>

/************************************************************************
  Class CVector
  vector number Data Abstraction
************************************************************************/

class CVector { 
  private:/*----------------------------------------------------------*/
	public:/*----------------------------------------------------------*/
		double x;
    double y;
		double z;

		CVector(){}
		CVector(const CVector &in)                                 {x=in.x; y=in.y; z=in.z;}
		//CVector(const complex <double> &z)                       {x=z.real();y=z.imag();z=0.0;}
		CVector(const double &d1,const double &d2,const double &d3){x=d1;y=d2;z=d3;}
    CVector(const double &d1)                                  {x=d1;y=0.0;z=0.0;}
	 ~CVector(){}
    
		inline CVector& CVector::operator=(const CVector &v){x=v.x;y=v.y;z=v.z;return *this;}
		inline CVector& CVector::operator=(const double &d) {x=d;  y=0.0;z=0.0;return *this;}
};

    /*inline*/ double abs(const CVector &v);
		
    //addition-----------------------------------------------------
		/*inline*/ CVector operator+(const CVector &v1, const CVector &v2);
    /*inline*/ void operator+=(CVector &v1, const CVector &v2);

    //subtraction-----------------------------------------------------
    /*inline*/ CVector operator-(const CVector &v1, const CVector &v2);
    /*inline*/ void operator-=(CVector &v1,const CVector &v2);

		//multiplication-----------------------------------------------------
    /*inline CVector operator*(const CVector &v1, const CVector &v2){
		  return CVector(v1.x * v2.x - v1.y * v2.y,
				             v1.x * v2.y + v1.y * v2.x);
		}*/
    /*inline*/ CVector operator*(const double &d1, const CVector &v2);
    /*inline*/ CVector operator*(const CVector &v1, const double &d2);
    /*inlinevoid operator*=(CVector &v1, const CVector &v2){
		  v1= v1*v2;
		}*/
    /*inline*/ void operator*=(CVector &v1, const double &d2);

#endif