#ifndef SID_DATA_H
#define SID_DATA_H


#include "../src/ReadSIDLib/readSID.h"
#include "ModelicaUtilities.h"


typedef struct{
  SID_Data* sid;
} extObjSID;

void* SID_Constructor(char* fileName){
	extObjSID* extObj = (extObjSID*)malloc(sizeof(extObjSID));
	ModelicaFormatMessage("SID_Data.h: Create SID object of file %s \n",fileName);
	extObj->sid = SIDFileConstructor(fileName);
	return (void *)extObj;
}

void SID_Destructor(void* p_eo){
	extObjSID* extObj = (extObjSID*)p_eo;
	//int i=0;
	//for(i=0;i<(extObj->sid->numNodes);i++)
	//{
	//	free(extObj->sid->nodes[i].orig.M0);
	//	free(extObj->sid->nodes[i].orig.M1);
	//	free(extObj->sid->nodes[i].orig.Mn);
	//	free(extObj->sid->nodes[i].phi.M0);
	//	free(extObj->sid->nodes[i].phi.M1);
	//	free(extObj->sid->nodes[i].phi.Mn);
	//	free(extObj->sid->nodes[i].psi.M0);
	//	free(extObj->sid->nodes[i].psi.M1);
	//	free(extObj->sid->nodes[i].psi.Mn);
	//	free(extObj->sid->nodes[i].AP.M0);
	//	free(extObj->sid->nodes[i].AP.M1);
	//	free(extObj->sid->nodes[i].AP.Mn);
	//}
	//free(extObj->sid->mdCM.M0);
	//free(extObj->sid->mdCM.M1);
	//free(extObj->sid->mdCM.Mn);
	//free(extObj->sid->J.M0);
	//free(extObj->sid->J.M1);
	//free(extObj->sid->J.Mn);
	//free(extObj->sid->Ct.M0);
	//free(extObj->sid->Ct.M1);
	//free(extObj->sid->Ct.Mn);
	//free(extObj->sid->Cr.M0);
	//free(extObj->sid->Cr.M1);
	//free(extObj->sid->Cr.Mn);
	//free(extObj->sid->Me.M0);
	//free(extObj->sid->Me.M1);
	//free(extObj->sid->Me.Mn);
	//free(extObj->sid->Gr.M0);
	//free(extObj->sid->Gr.M1);
	//free(extObj->sid->Gr.Mn);
	//free(extObj->sid->Ge.M0);
	//free(extObj->sid->Ge.M1);
	//free(extObj->sid->Ge.Mn);
	//free(extObj->sid->Oe.M0);
	//free(extObj->sid->Oe.M1);
	//free(extObj->sid->Oe.Mn);
	//free(extObj->sid->ksigma.M0);
	//free(extObj->sid->ksigma.M1);
	//free(extObj->sid->ksigma.Mn);
	//free(extObj->sid->Ke.M0);
	//free(extObj->sid->Ke.M1);
	//free(extObj->sid->Ke.Mn);
	//free(extObj->sid->De.M0);
	//free(extObj->sid->De.M1);
	//free(extObj->sid->De.Mn);
	//free(extObj->sid->nodes);
	//free(extObj->sid);
	//free(extObj);
}
	
int getNumberOfNodes(void* p_eo){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->numNodes;
}

int getNumberOfModes(void* p_eo){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->numModes;
}


//Taylor Retrieval

int getOriginOrder(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].orig.order;
}
	
int getOriginNRows(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].orig.nrow;
}
	
int getOriginNCols(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].orig.ncol;
}

int getOriginNq(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].orig.nq;
}

int getOriginNqn(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].orig.nqn;
}

int getOriginStructure(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].orig.structure;
}
	
double getOriginM0(void* p_eo, int nodeIdx, int r, int c, int dimR, int dimC){
	extObjSID* extObj = (extObjSID*)p_eo;
	double* dp =  extObj->sid->nodes[nodeIdx-1].orig.M0;
	int i = (c-1) + (r-1) * dimC;
	//ModelicaFormatMessage("r %d c %d dimR %d dimC %d Idx %d\n",r,c,dimR,dimC, i);
	//ModelicaFormatMessage("node %d %f %f %f\n",nodeIdx, dp[0],dp[1],dp[2]);
	return extObj->sid->nodes[nodeIdx-1].orig.M0[i];
}

double getOriginM1(void* p_eo, int nodeIdx, int r, int q, int c, int dimR, int dimQ, int dimC){
	extObjSID* extObj = (extObjSID*)p_eo;
	double* dp =  extObj->sid->nodes[nodeIdx-1].orig.M1;
	int i = (c-1) + (r-1) * dimC + (q-1) * dimC*dimR;
	return extObj->sid->nodes[nodeIdx].orig.M1[i];
}


#endif