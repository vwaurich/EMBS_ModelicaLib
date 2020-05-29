#ifndef SID_DATA_H
#define SID_DATA_H


#include "../src/ReadSIDLib/readSID.h"
#include "ModelicaUtilities.h"


typedef struct{
  SID_Data* sid;
} extObjSID;

void* SID_Constructor(char* fileName){
	extObjSID* sid = (extObjSID*)malloc(sizeof(extObjSID));
	ModelicaFormatMessage("SID_Data.h: Create SID object of file %s \n",fileName);
	sid->sid = SIDFileConstructor(fileName);
	return (void *)sid;
}

void SID_Destructor(void* p_eo){
	extObjSID* sid = (extObjSID*)p_eo;
	free(sid);
	}
	
	
int getNumberOfNodes(void* p_eo){
	extObjSID* sid = (extObjSID*)p_eo;
	int i = SIDFile_getNumberOfNodes(sid->sid);
	return i;
}

int getNumberOfModes(void* p_eo){
	extObjSID* sid = (extObjSID*)p_eo;
	int i = SIDFile_getNumberOfModes(sid->sid);
	return i;
}

void getM0ArrforNode(void* p_eo, int nodeIdx, double* m0){
	extObjSID* sid = (extObjSID*)p_eo;
	SIDFile_getM0ForNode(sid->sid, nodeIdx, m0);
}

#endif