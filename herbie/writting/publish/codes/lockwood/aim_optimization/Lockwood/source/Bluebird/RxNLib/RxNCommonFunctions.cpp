
#include "RxNLibraryInclude.h"

//Global Variables 
//************************************************************************
char*   RxNError=NULL;
int     RxNErrorCode=0;
void (*ExitFunction)(char *statement, badcode code);

/******************************************************************************
  Library exit strategy functions-------------------------------------------------
*******************************************************************************/
void   InitializeRxNLibrary(void (*ExternalExitFunct)(char *statement, badcode code)){
  ExitFunction=ExternalExitFunct;
}
void   RxNExitGracefully(char *statement, badcode code){
	ExitFunction(statement, code);
  /*if (statement!=NULL){
    RxNError= new char[strlen(statement)+1];
    RxNError=strcpy(RxNError,statement);
	}
  RxNErrorCode=(int)(code);*/
}
//******************************************************************************
void   RxNExitGracefullyIf(bool condition, char *statement, badcode code){
  if (condition){RxNExitGracefully(statement,code);}
}
//******************************************************************************
/*bool   CheckRxNError(char *statement, int &code){
  if (RxNError!=NULL){
    statement= new char[strlen(RxNError)+1];
    statement=strcpy(statement,RxNError);
	}
	code = RxNErrorCode;
	return (RxNErrorCode!=0);
}*/