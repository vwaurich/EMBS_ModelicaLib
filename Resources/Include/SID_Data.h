#ifndef SID_DATA_H
#define SID_DATA_H


#include "../src/ReadSIDLib/readSID.h"

typedef struct{
  int i1;
  int i2;
  SID_Data sid;
} extObjSID;

void* SID_Constructor(const char* fileName){
	extObjSID* sid = (extObjSID*)malloc(sizeof(extObjSID));
	sid->i1 = 10;
	sid->i2 = 2;
	return (void *)sid;
}

void SID_Destructor(void* p_eo){
	extObjSID* sid = (extObjSID*)p_eo;
	free(sid);
	}
	
	
int getNumberOfNodes(void* p_eo){
	return 8;
}

#endif