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
//=================================

//origin-taylor
//=================================
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
	return extObj->sid->nodes[nodeIdx-1].orig.M1[i];
}

//phi-taylor
//=================================
int getPhiOrder(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].phi.order;
}
	
int getPhiNRows(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].phi.nrow;
}
	
int getPhiNCols(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].phi.ncol;
}

int getPhiNq(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].phi.nq;
}

int getPhiNqn(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].phi.nqn;
}

int getPhiStructure(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].phi.structure;
}
	
double getPhiM0(void* p_eo, int nodeIdx, int r, int c, int dimR, int dimC){
	extObjSID* extObj = (extObjSID*)p_eo;
	double* dp =  extObj->sid->nodes[nodeIdx-1].phi.M0;
	int i = (c-1) + (r-1) * dimC;
	return extObj->sid->nodes[nodeIdx-1].phi.M0[i];
}

double getPhiM1(void* p_eo, int nodeIdx, int r, int q, int c, int dimR, int dimQ, int dimC){
	extObjSID* extObj = (extObjSID*)p_eo;
	double* dp =  extObj->sid->nodes[nodeIdx-1].phi.M1;
	int i = (c-1) + (r-1) * dimC + (q-1) * dimC*dimR;
	return extObj->sid->nodes[nodeIdx-1].phi.M1[i];
}

//psi-taylor
//=================================
int getPsiOrder(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].psi.order;
}
	
int getPsiNRows(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].psi.nrow;
}
	
int getPsiNCols(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].psi.ncol;
}

int getPsiNq(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].psi.nq;
}

int getPsiNqn(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].psi.nqn;
}

int getPsiStructure(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].psi.structure;
}
	
double getPsiM0(void* p_eo, int nodeIdx, int r, int c, int dimR, int dimC){
	extObjSID* extObj = (extObjSID*)p_eo;
	double* dp =  extObj->sid->nodes[nodeIdx-1].psi.M0;
	int i = (c-1) + (r-1) * dimC;
	return extObj->sid->nodes[nodeIdx-1].psi.M0[i];
}

double getPsiM1(void* p_eo, int nodeIdx, int r, int q, int c, int dimR, int dimQ, int dimC){
	extObjSID* extObj = (extObjSID*)p_eo;
	double* dp =  extObj->sid->nodes[nodeIdx-1].psi.M1;
	int i = (c-1) + (r-1) * dimC + (q-1) * dimC*dimR;
	return extObj->sid->nodes[nodeIdx-1].psi.M1[i];
}


//AP-taylor
//=================================
int getAPOrder(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].AP.order;
}
	
int getAPNRows(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].AP.nrow;
}
	
int getAPNCols(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].AP.ncol;
}

int getAPNq(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].AP.nq;
}

int getAPNqn(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].AP.nqn;
}

int getAPStructure(void* p_eo, int nodeIdx){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->nodes[nodeIdx-1].AP.structure;
}
	
double getAPM0(void* p_eo, int nodeIdx, int r, int c, int dimR, int dimC){
	extObjSID* extObj = (extObjSID*)p_eo;
	double* dp =  extObj->sid->nodes[nodeIdx-1].AP.M0;
	int i = (c-1) + (r-1) * dimC;
	return extObj->sid->nodes[nodeIdx-1].AP.M0[i];
}

double getAPM1(void* p_eo, int nodeIdx, int r, int q, int c, int dimR, int dimQ, int dimC){
	extObjSID* extObj = (extObjSID*)p_eo;
	double* dp =  extObj->sid->nodes[nodeIdx-1].AP.M1;
	int i = (c-1) + (r-1) * dimC + (q-1) * dimC*dimR;
	ModelicaFormatMessage("node %d %f %f %f %f %f %f %f %f %f\n",nodeIdx, dp[0],dp[1],dp[2],dp[4],dp[5],dp[6],dp[7],dp[8]);
	ModelicaFormatMessage("r %d q %d c %d dimR %d dimQ %d dimC %d Idx %d\n",r,q,c,dimR,dimQ,dimC, i);

	
	return extObj->sid->nodes[nodeIdx-1].AP.M1[i];
}



//All remaining taylors
//=================================

taylor* getTaylorPtr(void* p_eo, const char* taylorName){
	extObjSID* extObj = (extObjSID*)p_eo;
	taylor* t;
	if(!strcmp(taylorName,"mdCM")){
		t = &extObj->sid->mdCM;
	}
	else if(!strcmp(taylorName,"J")){
		t = &extObj->sid->J;
	}
	return t;
}

int getTaylorOrder(void* p_eo, const char* taylorName){
	taylor* t = getTaylorPtr(p_eo,taylorName);
	return t->order;
}

int getTaylorNrow(void* p_eo, const char* taylorName){
	taylor* t = getTaylorPtr(p_eo,taylorName);
	return t->nrow;
}

int getTaylorNcol(void* p_eo, const char* taylorName){
	taylor* t = getTaylorPtr(p_eo,taylorName);
	return t->ncol;
}

int getTaylorNq(void* p_eo, const char* taylorName){
	taylor* t = getTaylorPtr(p_eo,taylorName);
	return t->nq;
}

int getTaylorNqn(void* p_eo, const char* taylorName){
	taylor* t = getTaylorPtr(p_eo,taylorName);
	return t->nqn;
}

int getTaylorStructure(void* p_eo, const char* taylorName){
	taylor* t = getTaylorPtr(p_eo,taylorName);
	return t->structure;
}

double getTaylorM0(void* p_eo, const char* taylorName, int r, int c, int dimR, int dimC){
	taylor* t = getTaylorPtr(p_eo, taylorName);
	int i = (c-1) + (r-1) * dimC;
	return t->M0[i];
}

double getTaylorM1(void* p_eo, const char* taylorName, int r, int q, int c, int dimR, int dimQ, int dimC){
	taylor* t = getTaylorPtr(p_eo, taylorName);
	int i = (c-1) + (r-1) * dimC + (q-1) * dimC*dimR;
	return t->M1[i];
}

#endif