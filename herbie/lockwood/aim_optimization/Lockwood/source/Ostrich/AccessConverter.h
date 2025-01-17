/******************************************************************************
File     : AccessConverter.h
Author   : Tsu-Wei Webber Chen
Copyright: 2010, Tsu-Wei Webber Chen

AccessConverter deals with Microsoft Access files.

Version History
01-25-10    twc   Created file.
******************************************************************************/
#ifndef ACCESS_CONVERTER_H
#define ACCESS_CONVERTER_H

#include "DatabaseABC.h"
#include "Exception.h"
#include "MyTypes.h"

#ifdef WIN32
#include <string>
using namespace std;
#endif

/******************************************************************************
class AccessConverter

An instance of the DatabaseABC abstract base class.
******************************************************************************/
class AccessConverter : public DatabaseABC
{
   public:
	  AccessConverter(void);
     void Initialize(char * line);
	  bool ReadFromFile(void);
	  ~AccessConverter(void);
	  void Destroy(void);
	  void Convert(void);

     void InsertDbase(DatabaseABC * pNxt);
     DatabaseABC * GetNext(void) { return m_pNxt;}
     bool WriteParameter(char * pName, char * pValue);
     void ReadResponse(void);
     void DeleteASCIIFile(void);

   private:
      DatabaseABC * m_pNxt;
      bool m_bIsEmpty; //true if the converter has not been initialized
      #ifdef WIN32
	      string m_connectionString;
      #else
         char m_connectionString[DEF_STR_SZ]; 
      #endif
	   char m_accessType[DEF_STR_SZ];
      char m_fileName[DEF_STR_SZ];
      char m_table[DEF_STR_SZ];
      char m_keyColumn[DEF_STR_SZ];
      char m_key[DEF_STR_SZ];
      char m_column[DEF_STR_SZ];
      char m_param[DEF_STR_SZ];
      char m_name[DEF_STR_SZ];
	   void CreateConnectionString(void);
}; /* end class AccessConverter */

#endif /* ACCESS_CONVERTER_H */

