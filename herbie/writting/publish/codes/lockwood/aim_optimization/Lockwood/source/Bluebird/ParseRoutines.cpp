#include "GeneralInclude.h"

/************************************************************************
								TokenizeLine
*************************************************************************
  tokenizes a sentence delimited by  spaces, tabs & return characters
parameters:
	FILE is the file from which the line is requested
  out is the array of strings in the line
	numwords is the number of strings in the line
	returns true if file has ended
-------------------------------------------------------------------------*/
bool TokenizeLine(ifstream &FILE, char **out, int &numwords){

	static char wholeline     [MAXCHARINLINE];
	static char *tempwordarray[MAXINPUTITEMS];
  char *p;
	int ct(0),w;
	char delimiters[6];
	delimiters[0]=' '; //space
	delimiters[1]='\t';//tab
	delimiters[2]=','; //comma
	delimiters[3]='\r';//carriage return (ignored in visual c++)
	delimiters[4]='\n';//newline -necc?
	delimiters[5]='\0';//last delimiter

	FILE.getline(wholeline,MAXCHARINLINE);             //get entire line as 1 string
	
	if ((*wholeline)==0){return true;}                 //if line is null, break, return true
	p=strtok(wholeline, delimiters);                   //tokenize first word in line
	while (p){                                         //sift through words, place in temparray, count line length  
    tempwordarray[ct]=p;
	  p=strtok(NULL, delimiters);
    ct++;
	}
  for (w=0; w<ct; w++){                                              //copy temp array of words into out[] 
		if (w>MAXINPUTITEMS){return true;}
		//	ExitGracefully("Tokenizeline::bad allocation",OUT_OF_MEMORY);}
		out[w]=tempwordarray[w];
		//cout<<out[w]<<"|";
	}
	//cout <<endl;
	numwords=ct;
	return false;
}
//******************************************************************************
void   ImproperFormat(char **s, int l){
	cout <<"line"<< l << "is wrong length"<<endl; //FOR NOW
}
//******************************************************************************
void   WriteEllipsisToScreen(const int i,const int N, const int modN){
	if ((N>modN) && (i%(N/modN)==0)) {cout<< ".";}
}