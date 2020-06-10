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
	extObj->sid = (extObjSID*)SIDFileConstructor(fileName);
	
	ModelicaFormatMessage("mdCM  [%d,%d, %d]  at %p\n",extObj->sid->mdCM.nrow,extObj->sid->mdCM.ncol,extObj->sid->mdCM.nq,(void*)&extObj->sid->mdCM);
	ModelicaFormatMessage("J  [%d,%d, %d]  at %p\n",extObj->sid->J.nrow,extObj->sid->J.ncol,extObj->sid->J.nq,(void*)&extObj->sid->J);
	ModelicaFormatMessage("Ct  [%d,%d, %d]  at %p\n",extObj->sid->Ct.nrow,extObj->sid->Ct.ncol,extObj->sid->Ct.nq,(void*)&extObj->sid->Ct);
	ModelicaFormatMessage("Cr  [%d,%d, %d]  at %p\n",extObj->sid->Cr.nrow,extObj->sid->Cr.ncol,extObj->sid->Cr.nq,(void*)&extObj->sid->Cr);
	ModelicaFormatMessage("Me  [%d,%d, %d]  at %p\n",extObj->sid->Me.nrow,extObj->sid->Me.ncol,extObj->sid->Me.nq,(void*)&extObj->sid->Me);
	ModelicaFormatMessage("Gr  [%d,%d, %d]  at %p\n",extObj->sid->Gr.nrow,extObj->sid->Gr.ncol,extObj->sid->Gr.nq,(void*)&extObj->sid->Gr);
	ModelicaFormatMessage("Ge  [%d,%d, %d]  at %p\n",extObj->sid->Ge.nrow,extObj->sid->Ge.ncol,extObj->sid->Ge.nq,(void*)&extObj->sid->Ge);
	ModelicaFormatMessage("Oe  [%d,%d, %d]  at %p\n",extObj->sid->Oe.nrow,extObj->sid->Oe.ncol,extObj->sid->Oe.nq,(void*)&extObj->sid->Oe);
	ModelicaFormatMessage("ksigma  [%d,%d, %d]  at %p\n",extObj->sid->ksigma.nrow,extObj->sid->ksigma.ncol,extObj->sid->ksigma.nq,(void*)&extObj->sid->ksigma);
	ModelicaFormatMessage("Ke  [%d,%d, %d]  at %p\n",extObj->sid->Ke.nrow,extObj->sid->Ke.ncol,extObj->sid->Ke.nq,(void*)&extObj->sid->Ke);
	ModelicaFormatMessage("De  [%d,%d, %d]  at %p\n",extObj->sid->De.nrow,extObj->sid->De.ncol,extObj->sid->De.nq,(void*)&extObj->sid->De);
	
	return (void *)extObj;
}

void SID_Destructor(void* p_eo){
	extObjSID* extObj = (extObjSID*)p_eo;
	int i=0;
	for(i=0;i<(extObj->sid->numNodes);i++)
	{
		free(extObj->sid->nodes[i].orig.M0);
		free(extObj->sid->nodes[i].orig.M1);
		free(extObj->sid->nodes[i].orig.Mn);
		free(extObj->sid->nodes[i].phi.M0);
		free(extObj->sid->nodes[i].phi.M1);
		free(extObj->sid->nodes[i].phi.Mn);
		free(extObj->sid->nodes[i].psi.M0);
		free(extObj->sid->nodes[i].psi.M1);
		free(extObj->sid->nodes[i].psi.Mn);
		free(extObj->sid->nodes[i].AP.M0);
		free(extObj->sid->nodes[i].AP.M1);
		free(extObj->sid->nodes[i].AP.Mn);
	}
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

double getMass(void* p_eo){
	extObjSID* extObj = (extObjSID*)p_eo;
	return extObj->sid->modeStruct.mass;
}



//Taylor Retrieval
//=================================

taylor* getTaylorByName(SID_Data* sid, const char* name)
{
	//ModelicaFormatMessage("getTaylorByName %s\n",name);
	if(!strcmp(name,"mdCM")){
		return &sid->mdCM;
	}
	else if (!strcmp(name,"J")){
	    return &sid->J;
	}
	else if (!strcmp(name,"Ct")){
		return &sid->Ct;
	}
	else if (!strcmp(name,"Cr")){
		return &sid->Cr;
	}
	else if (!strcmp(name,"Me")){
		return &sid->Me;
	}
	else if (!strcmp(name,"Gr")){
		return &sid->Gr;	
	}
	else if (!strcmp(name,"Ge")){
		return &sid->Ge;
	}
	else if (!strcmp(name,"Oe")){
		return &sid->Oe;
	}
	else if (!strcmp(name,"ksigma")){
		return &sid->ksigma;
	}
	else if (!strcmp(name,"Ke")){
		return &sid->Ke;
	}
	else if (!strcmp(name,"De")){
		return &sid->De;
	}
	else{
	ModelicaFormatMessage("getTaylorByName: the given name is not defined: %s\n",name);
	return &sid->mdCM;
	}
}

taylor* getNodeTaylorByName(node* node, const char* name)
{
	//ModelicaFormatMessage("getNodeTaylorByName %s\n",name);
	if(!strcmp(name,"origin")){
		return &node->orig;
	}
	else if (!strcmp(name,"phi")){
	    return &node->phi;
	}
	else if (!strcmp(name,"psi")){
		return &node->psi;
	}
	else if (!strcmp(name,"AP")){
		return &node->AP;
	}
	else{
	ModelicaFormatMessage("getNodeTaylorByName: the given name is not defined: %s\n",name);
	return &node->orig;
	}
}

//M0 for modal objects
void getM0(void* p_eo, const char* taylorName, double* m0, size_t nr, size_t nc){

	extObjSID* extObj = (extObjSID*)p_eo;
	taylor* t = getTaylorByName((extObj->sid), taylorName);
	ModelicaFormatMessage("getM0 %s [%d, %d]  at %p\n",taylorName,nr,nc,t);

	if(nr!=t->nrow || nc!=t->ncol){
		ModelicaFormatMessage(" getM0: %s the given dimensions [%d, %d] are not equal to the stored matrix dimension [%d %d] %p\n",taylorName, nr, nc, t->nrow,t->ncol,(void*)t);
	}
	else{
	  //memcpy(m0,t->M0,sizeof(double)*nr*nc);
	}
}
//M1 for modal objects
void getM1(void* p_eo, const char* taylorName, double* m1, size_t nr, size_t nq, size_t nc){
	extObjSID* extObj = (extObjSID*)p_eo;
	ModelicaFormatMessage("getM1 %s [%d, %d,  %d]\n",taylorName,nr,nq,nc);
	taylor* t = getTaylorByName((extObj->sid), taylorName);

	if(nr!=t->nrow || nc!=t->ncol || nq!=t->nq){
		ModelicaFormatMessage(" getM1: the given dimensions [%d, %d, %d] are not equal to the stored matrix dimension [%d, %d, %d]\n",nr,nq,nc, t->nrow,t->nq,t->ncol);
	}
	else{
	  memcpy(m1,t->M1,sizeof(double)*nr*nc*nq);
	}
}

//M0 for node object
void getM0Node(void* p_eo, const char* taylorName, int nodeIdx, double* m0, size_t nr, size_t nc){

	ModelicaFormatMessage("getM0Node: Get  for node %s [%d]\n",taylorName,nodeIdx);

	extObjSID* extObj = (extObjSID*)p_eo;
	if (nodeIdx >extObj->sid->numNodes){
		ModelicaFormatMessage(" getM0Node: the given nodeidx  [%d] is bigger than the number of nodes [%d]\n", nodeIdx, extObj->sid->numNodes);
	}

	node n = extObj->sid->nodes[nodeIdx-1];
	
	taylor* t = getNodeTaylorByName(&n, taylorName);

	if(nr!=t->nrow || nc!=t->ncol){
		ModelicaFormatMessage(" getM0Node: the given dimensions [%d, %d] are not equal to the stored matrix dimension [%d %d]\n", nr, nc, t->nrow,t->ncol);
	}
	else{
	  memcpy(m0,t->M0,sizeof(double)*nr*nc);
	}
	
}
//mdCM
void getM1Node(void* p_eo, const char* taylorName, int nodeIdx, double* m1, size_t nr, size_t nq, size_t nc){

	extObjSID* extObj = (extObjSID*)p_eo;
	if (nodeIdx > extObj->sid->numNodes){
		ModelicaFormatMessage(" getM0Node: the given nodeidx  [%d] is bigger than the number of nodes [%d]\n", nodeIdx, extObj->sid->numNodes);
	}

	node n = extObj->sid->nodes[nodeIdx-1];
	taylor* t = getNodeTaylorByName(&n, taylorName);
			
	if(nr!=t->nrow || nc!=t->ncol || nq!=t->nq){
		ModelicaFormatMessage(" getM1: the given dimensions [%d, %d, %d] are not equal to the stored matrix dimension [%d, %d, %d]\n",nr,nq,nc, t->nrow,t->nq,t->ncol);
	}
	else{
	  memcpy(m1,t->M1,sizeof(double)*nr*nc*nq);
	}

}


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


#endif