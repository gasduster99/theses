/******************************************************************************
File      : SuperMuseUtility.h
Author    : L. Shawn Matott
Copyright : 2007, L. Shawn Matott

This file contains a c-style interface to the global SuperMUSE instance.

Version History
07-16-07    lsm   Created.
******************************************************************************/
#ifndef SUPER_MUSE_UTILITY_H
#define SUPER_MUSE_UTILITY_H

#include <stdio.h>
#include "SuperMUSE.h"
#include "ModelABC.h"

extern "C" {
bool IsSuperMUSE(void);
void EnableSuperMUSE(void);
void DisableSuperMUSE(void);
SuperMUSE * GetSuperMusePtr();
void InitSuperMUSE(FILE * pFile, ModelABC * pModel);
void DestroySuperMUSE(void);
void WriteSuperMuseSetupToFile(FILE * pFile);
void CleanSuperMUSE(void);
}

#endif /* SUPER_MUSE_UTILITY_H */

