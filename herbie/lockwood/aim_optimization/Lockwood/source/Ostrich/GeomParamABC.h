/******************************************************************************
File      : GeomParamABC.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a 'geometry' parameter. Geometry parameters are variables in the model 
which are composed on one or more spatial vertices. The ABC for the geometry 
parameters encapsulates the interface used by other Ostrich modules, allowing 
various specific geometry parameter relationships (line3, poly2, poly3, etc.) to be 
implemented as needed with minimal code change (just need to add the specific 
geometry parameter class and some additional input file parsing). The purpose of a
geometry parameter is to facilitate changes to geometric model properties as part of
a calibration or optimization exercise. Ostrich will ensure that all geometry parameters
are topologically correct in that:
   1) vertices will be automatically inserted if elements overlap
   2) if a given ordering of vertices is not valid, the polygon vertices will be 
   randomly reordered until a valid polygon is found

These specific geometry-parameter classes are supported:

GeomParamLine3 : a polyline containing a set of (x,y,z) values, where (x,y) are 
the spatial coordinates and z is a non-geometric value (i.e. head). When 
vertices are inserted, the z-value of the new vertex is interpolated.

GeomParamPoly3  : a polygon containing a set of (x,y,z) values, where (x,y) are 
the spatial coordinates and z is a non-geometric value (i.e. head). When 
vertices are inserted, the z-value of the new vertex is interpolated.

GeomParamPoly2  : a polygon containing a set of (x,y) values, where (x,y) are 
the spatial coordinates.

Version History
11-29-04    lsm   Created
******************************************************************************/
#ifndef GEOM_PARAM_ABC_H
#define GEOM_PARAM_ABC_H

#include <stdio.h>
#include "MyTypes.h"
#include "VertexList.h"
#include "GeometryUtility.h"

#define MY_POLYGON_TYPE (0)
#define MY_LINE_TYPE (1)
#define MY_CIRCLE_TYPE (2)

/******************************************************************************
class GeomParamABC

Abstract base class of a geometry-parameter.
******************************************************************************/
class GeomParamABC
{
   public:
      virtual void   Convert(void) = 0;
      virtual bool   Reorder(void) = 0;
      virtual bool   FixVertices(GeomParamABC * pOther) = 0;
      virtual int    GetValStrSize(void) = 0;
      virtual void   GetValAsStr(UnmoveableString valStr) = 0;
      virtual void   Write(FILE * pFile, int type) = 0;
      virtual UnchangeableString GetName(void) = 0;
      virtual void Destroy(void) = 0;
      virtual void InsertVertex(AugVertexList * pNew) = 0;      
      virtual VertexList * FixVertex(Segment2D * pSeg) = 0;
      virtual VertexList * GetVertexList(int * type) = 0;
}; /* end class GeomParamABC */

/******************************************************************************
class GeomParamLine3

Represents a line geometry having 2 spatial coordinates (x,y) and one 
non-spatial value (z) at each vertex.
******************************************************************************/
class GeomParamLine3 : public GeomParamABC
{
   public:      
      GeomParamLine3(void);
      GeomParamLine3(IroncladString name);
      ~GeomParamLine3(void){ Destroy();}
      void Destroy(void);

      void   Convert(void);
      bool   Reorder(void);
      bool   FixVertices(GeomParamABC * pOther);
      int    GetValStrSize(void);
      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);      
      UnchangeableString GetName(void){ return m_pName;}
      void InsertVertex(AugVertexList * pNew);
      VertexList * GetVertexList(int * type){ *type = MY_LINE_TYPE; return m_pFixed;}

   private:
      VertexList * FixVertex(Segment2D * pSeg);

      StringType m_pName;
      AugVertexList * m_pInit;
      VertexList * m_pFixed;
}; /* end class GeomParamLine3 */

/******************************************************************************
class GeomParamPoly3

Represents a polygon geometry having 2 spatial coordinates (x,y) and one 
non-spatial value (z) at each vertex.
******************************************************************************/
class GeomParamPoly3 : public GeomParamABC
{
   public:      
      GeomParamPoly3(void);
      GeomParamPoly3(IroncladString name);
      ~GeomParamPoly3(void){ Destroy();}
      void Destroy(void);

      void   Convert(void);
      bool   Reorder(void);
      bool   FixVertices(GeomParamABC * pOther);
      int    GetValStrSize(void);
      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      UnchangeableString GetName(void){ return m_pName;}
      void InsertVertex(AugVertexList * pNew);
      VertexList * GetVertexList(int * type){ *type = MY_POLYGON_TYPE; return m_pFixed;}

   private:
      VertexList * FixVertex(Segment2D * pSeg);

      StringType m_pName;
      AugVertexList * m_pInit;
      VertexList * m_pFixed;
}; /* end class GeomParamPoly3 */

/******************************************************************************
class GeomParamPoly2

Represents a polygon geometry having 2 spatial coordinates (x,y) at each vertex.
******************************************************************************/
class GeomParamPoly2 : public GeomParamABC
{
   public:      
      GeomParamPoly2(void);
      GeomParamPoly2(IroncladString name);
      ~GeomParamPoly2(void){ Destroy();}
      void Destroy(void);

      void   Convert(void);
      bool   Reorder(void);
      bool   FixVertices(GeomParamABC * pOther);      
      int    GetValStrSize(void);
      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      UnchangeableString GetName(void){ return m_pName;}
      void InsertVertex(AugVertexList * pNew);
      VertexList * GetVertexList(int * type){ *type = MY_POLYGON_TYPE; return m_pFixed;}

   private:
      VertexList * FixVertex(Segment2D * pSeg);
      StringType m_pName;
      AugVertexList * m_pInit;
      VertexList * m_pFixed;
}; /* end class GeomParamPoly2 */

/******************************************************************************
class GeomParamCirc4

Represents a circle geometry having center at (x,y) radius of r and one 
non-spatial value (z).
******************************************************************************/
class GeomParamCirc4 : public GeomParamABC
{
   public:      
      GeomParamCirc4(void);
      GeomParamCirc4(IroncladString name, AugCircle * pData);
      ~GeomParamCirc4(void){ Destroy();}
      void Destroy(void);

      void   Convert(void);
      bool   Reorder(void){ return true;}
      bool   FixVertices(GeomParamABC * pOther);      
      int    GetValStrSize(void){ return 100;}
      void   GetValAsStr(UnmoveableString valStr);
      void   Write(FILE * pFile, int type);
      UnchangeableString GetName(void){ return m_pName;}
      void InsertVertex(AugVertexList * pNew){ return;}
      VertexList * GetVertexList(int * type);

   private:
      VertexList * FixVertex(Segment2D * pSeg);
      StringType m_pName;
      AugCircle  * m_pInit;
      Circle2D   m_Fixed;
      double     m_Zcur;
}; /* end class GeomParamCirc4 */

#endif /* GEOM_PARAM_ABC_H */
