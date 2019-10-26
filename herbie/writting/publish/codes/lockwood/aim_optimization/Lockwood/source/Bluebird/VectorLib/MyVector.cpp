#include "MyVector.h"

    /*inline*/ double abs(const CVector &v){   
		  return pow(v.x*v.x + v.y*v.y+v.z*v.z,0.5);
		}
		
    //addition-----------------------------------------------------
		/*inline*/ CVector operator+(const CVector &v1, const CVector &v2){
		  return CVector(v1.x + v2.x,v1.y + v2.y,v1.z+v2.z);
		}
    /*inline*/ void operator+=(CVector &v1, const CVector &v2){
		  v1.x=v1.x+v2.x;
		  v1.y=v1.y+v2.y;
			v1.z=v1.z+v2.z;
		}

    //subtraction-----------------------------------------------------
    /*inline*/ CVector operator-(const CVector &v1, const CVector &v2){
		  return CVector(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z);
		}
    /*inline*/ void operator-=(CVector &v1,const CVector &v2){
		  v1.x=v1.x-v2.x;
		  v1.y=v1.y-v2.y;
		  v1.z=v1.z-v2.z;
		}

		//multiplication-----------------------------------------------------
    /*inline CVector operator*(const CVector &v1, const CVector &v2){
		  return CVector(v1.x * v2.x - v1.y * v2.y,
				             v1.x * v2.y + v1.y * v2.x);
		}*/
    /*inline*/ CVector operator*(const double &d1, const CVector &v2){
		  return CVector(d1 * v2.x, d1 * v2.y, d1 * v2.z);
		}
    /*inline*/ CVector operator*(const CVector &v1, const double &d2){
		  return CVector(d2 * v1.x, d2 * v1.y, d2*v1.z);
		}
    /*inline void operator*=(CVector &v1, const CVector &v2){
		  v1= v1*v2;
		}*/
    /*inline*/ void operator*=(CVector &v1, const double &d2){
		  v1.x=v1.x*d2;
			v1.y=v1.y*d2;
			v1.z=v1.z*d2;
		}