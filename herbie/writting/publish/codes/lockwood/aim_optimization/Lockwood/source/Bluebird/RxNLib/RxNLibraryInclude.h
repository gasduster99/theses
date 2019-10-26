//RxNLibraryInclude.h
#ifndef RXNINCLUDE_H
#define RXNINCLUDE_H

#include "GeneralInclude.h"
#include <math.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <typeinfo>
#include <strstream>
#include <time.h>

//library exit strategy functions-------------------------------------------------
void   InitializeRxNLibrary(void (*ExternalExitFunct)(char *statement, badcode code));
void   RxNExitGracefully(char *statement, badcode code);
void   RxNExitGracefullyIf(bool condition, char *statement, badcode code);

#endif