#pragma once

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "ModelicaUtilities.h"



typedef struct
{
	char ielastq[20]; //name of the mode shape
	double freq;
} modeShape;


typedef struct
{
	double mass; //mass of body
	//int nelastq; //number of mode shapes
	//modeShape* ielastq;
	//int currModeIdx;
} refmod;

typedef struct
{
	int order; // order of expansion
	int nrow; //number of rows
	int ncol; // number of columns
	int nq; //dimension of q
	int nqn;//dimension of higher order expansion
	int structure;// 0 =zero matrix, 1=diagonal, 2=symmetric, 3=no particular, 4 = unit matrix
	double* M0;// zero order matrix, nrow X ncol
	double* M1;// first order matrix, nrow X nq X ncol
	double* Mn;// higher order matrix
} taylor; // = class taylor

typedef struct
{
	int id;
	char rFrame[8];
	taylor orig;
	taylor phi;
	taylor psi;
	taylor AP;
} node;

typedef struct {
	int numNodes;
	int numModes;
	int useGeoStiffness;
	refmod modeStruct;
	node* nodes;
	taylor mdCM;
	taylor J;
	taylor Ct;
	taylor Cr;
	taylor Me;
	taylor Gr;
	taylor Ge;
	taylor Oe;
	taylor ksigma;
	taylor Ke;
	taylor De;
	//int currNodeIdx;
} SID_Data;



int getNextInteger(char* line, int startToken, int* value) {
	int len = strlen(line);
	int t = 0;
	int i = startToken;
	int found = 0;
	char intTokens[6] = "";
	char c;
	while (i < len) {
		c = line[i];
		while (isdigit(c)) {
			strncat(intTokens, &c, 1);
			t++;
			i++;
			found = 1;
			c = line[i];
		}
		if (found) {
			break;
		}
		i++;
	}
	int v = atoi(&intTokens);
	*value = v;
	if (found)
		return i;
	else
		return -1;
}

int getNextReal(char* line, int startToken, double* value) {
	int len = strlen(line);
	int i = startToken;
	int found = 0;
	char baseToken[20] = "";
	char exponentToken[20] = "";
	char c;
	int seq = 0;
	while (i < len) {
		c = line[i];
		//start real retrieval
		if (isdigit(c) || c == 45) {
			seq = 1;
		}
		while ((0 < seq) && (seq< 3)) {
			//some digit
			if (isdigit(c) || c == 45 || c == 46) {
				//post comma
				if (seq == 1) {
					strncat(baseToken, &c, 1);
				}
				else if (seq == 2) {
					strncat(exponentToken, &c, 1);
				}
			}
			//exponent
			if (c == 'D' && seq == 1) {
				seq = 2;

			}
			if (c == '\n') {
				break;
			}
			i++;
			found = 1;
			c = line[i];
		}
		i++;
	}
	double base = atof(&baseToken);
	double exponent = atof(&exponentToken);
	double d = base * pow(10, exponent);
	*value = d;
	if (found)
		return i;
	else
		return -1;
}


void parseFirstLine(char* line, SID_Data* sid) {
	int startToken = 0;
	int* numNodes = calloc(sizeof(int), 1);
	int* numModes = calloc(sizeof(int), 1);
	startToken = getNextInteger(line, startToken, numNodes);
	startToken = getNextInteger(line, startToken, numModes);
	sid->numNodes = numNodes[0];
	sid->numModes = numModes[0];
	free(numNodes);
	free(numModes);
}

void parseRefMod(char* line, int bufferSize, FILE* file, SID_Data* sid) {
	while (fgets(line, bufferSize, file) != NULL) {
		if (strstr(line, " mass ")) {
			int startToken = 0;
			double mass = 0;
			int found = getNextReal(line, startToken, &mass);
			sid->modeStruct.mass = mass;
		}
		else if (strstr(line, "end refmod"))
		{
			return;
		}
	}
	printf("Did not found end of refmod\n");
}

taylor initTaylor(int order, int nr, int ncol, int nq, int nqn, int structure) {
	taylor t = { .order = 0,.nrow = nr,.ncol = ncol,.nq = nq,.nqn = nqn,.structure = structure,.M0 = NULL,.M1 = NULL,.Mn = NULL };
	t.M0 = calloc(t.nrow*t.ncol, sizeof(double));
	t.M1 = calloc(t.nrow*t.nq*t.ncol, sizeof(double));
	if (t.nqn != 0)
		t.Mn = calloc(t.nrow*t.nqn*t.ncol, sizeof(double));
	return t;
}

int matrix2Index(int r, int c, int numR, int numC)
{
	return (r - 1) * numC + c - 1;
}

int matrix3Index(int r, int q, int c, int numR, int numQ, int numC)
{
	return (q - 1) * numC * numR + (r - 1)*numC + c - 1;
}



taylor parseTaylor(char* line, int bufferSize, FILE* file) {
	int order = 0;
	int nrow = 0;
	int ncol = 0;
	int nq = 0;
	int nqn = 0;
	int structure = 0;
	taylor t;
	int gotObjects = 0;
	while (fgets(line, bufferSize, file) != NULL) {
		if (gotObjects == 6) {
			t = initTaylor(order, nrow, ncol, nq, nqn, structure);
			gotObjects = 7;
		}

		if (strstr(line, "order ")) {
			getNextInteger(line, 0, &order);
			gotObjects++;
		}
		else if (strstr(line, "nrow ")) {
			getNextInteger(line, 0, &nrow);
			gotObjects++;
		}
		else if (strstr(line, "ncol ")) {
			getNextInteger(line, 0, &ncol);
			gotObjects++;
		}
		else if (strstr(line, "nqn ")) {
			getNextInteger(line, 0, &nqn);
			gotObjects++;
		}
		else if (strstr(line, "nq ")) {
			getNextInteger(line, 0, &nq);
			gotObjects++;
		}
		else if (strstr(line, "structure ")) {
			getNextInteger(line, 0, &structure);
			gotObjects++;
		}
		else if (strstr(line, "m0")) {
			int idx1 = 0;
			int idx2 = 0;
			int tokenIdx = 0;
			double val = 0.0;
			tokenIdx = getNextInteger(line, tokenIdx, &idx1);// this is the 0 in m0
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx1);
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx2);
			getNextReal(line, tokenIdx, &val);
			int i = matrix2Index(idx1, idx2, nrow, ncol);
			//printf("set M0 i %d  = %f\n",i,val);
			t.M0[i] = val;
			//printM_2(t.M0, nrow, ncol);
		}
		else if (strstr(line, "m1")) {
			int idx1 = 0;
			int idx2 = 0;
			int idx3 = 0;
			int tokenIdx = 0;
			double val = 0.0;
			tokenIdx = getNextInteger(line, tokenIdx, &idx1);// this is the 0 in m0
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx1);
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx2);
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx3);
			getNextReal(line, tokenIdx, &val);
			int i = matrix3Index(idx1, idx2, idx3, nrow, nq, ncol);
			//printf("set M1 i %d  = %f\n", i, val);
			t.M1[i] = val;
		}
		else if (strstr(line, "end "))
		{
			//printTaylor(t);
			return t;
		}
	}
	printf("Did not found end of taylor\n");
	return t;
}

void parseNode(char* line, int bufferSize, FILE* file, SID_Data* sid, int nodeIdx) {
	taylor t;
	while (fgets(line, bufferSize, file) != NULL) {
		if (strstr(line, "origin")) {
			t = parseTaylor(line, bufferSize, file, sid);
			sid->nodes[nodeIdx].orig = t;
		}
		if (strstr(line, "phi")) {
			t = parseTaylor(line, bufferSize, file, sid);
			sid->nodes[nodeIdx].phi = t;
		}
		if (strstr(line, "psi")) {
			t = parseTaylor(line, bufferSize, file, sid);
			sid->nodes[nodeIdx].psi = t;
		}
		if (strstr(line, "AP")) {
			t = parseTaylor(line, bufferSize, file, sid);
			sid->nodes[nodeIdx].AP = t;
		}
		else if (strstr(line, "end node"))
		{
			return;
		}
	}
	printf("Did not found end of node\n");
}

void parseNodes(char* line, int bufferSize, FILE* file, SID_Data* sid) {
	int nodeIdx = 0;
	while (fgets(line, bufferSize, file) != NULL) {
		if (strstr(line, "new node")) {
			nodeIdx++;
			int* nodeID = calloc(sizeof(int), 1);
			getNextInteger(line, 0, nodeID);
			sid->nodes[nodeIdx].id = nodeID[0];
			parseNode(line, bufferSize, file, sid, nodeIdx);
		}
		else if (strstr(line, "end node"))
		{
			if (nodeIdx == sid->numNodes) {
			}
		}
		else if (strstr(line, "end frame")) {
			return;
		}
	}
	printf("Did not found end of frame\n");
}


void SIDFileDestructor_C(void* p_sid) {
	SID_Data* sid = (SID_Data*)p_sid;
}

void* SIDFileConstructor_C(const char* fileName)
{
	FILE *sidFile;
	if ((sidFile = fopen(fileName, "r")) == NULL) {
		printf("Could not open %s\n", fileName);
	}

	SID_Data* sid = calloc(1, sizeof(SID_Data));
	printf("Allocated SID struct of size %d\n", sizeof(SID_Data));

	//read line by line ( http://openbook.rheinwerk-verlag.de/c_von_a_bis_z/016_c_ein_ausgabe_funktionen_016.htm#mjcea47bd6d32a4a8f51be329a672845d7 )
	int bufferSize = 512;
	char *line = malloc(sizeof(char*) * 512);

	int init = 0;
	taylor t;
	while (fgets(line, bufferSize, sidFile) != NULL) {
		if (!init) {
			parseFirstLine(line, sid);
			init = 1;
			sid->nodes = calloc(sizeof(node), sid->numNodes);
		}
		if (strstr(line, "refmod")) {
			parseRefMod(line, bufferSize, sidFile, sid);
		}
		else if (strstr(line, "frame")) {
			parseNodes(line, bufferSize, sidFile, sid);
		}
		else if (strstr(line, "mdCM ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->mdCM = t;
		}
		else if (strstr(line, "J ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->J = t;
		}
		else if (strstr(line, "Ct ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->Ct = t;
		}
		else if (strstr(line, "Cr ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->Cr = t;
		}
		else if (strstr(line, "Me ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->Me = t;
		}
		else if (strstr(line, "Gr ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->Gr = t;
		}
		else if (strstr(line, "Ge ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->Ge = t;
		}
		else if (strstr(line, "Oe ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->Oe = t;
		}
		else if (strstr(line, "ksigma ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->ksigma = t;
		}
		else if (strstr(line, "Ke ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->Ke = t;
		}
		else if (strstr(line, "De ")) {
			t = parseTaylor(line, bufferSize, sidFile, sid);
			sid->De = t;
		}

		//puts(line);
	}

	return (void*)sid;
}


//retrieve basic data
//=================================

int getNumberOfNodes(void* p_sid) {
	SID_Data* sid = (SID_Data*)p_sid;
	return sid->numNodes;
}

int getNumberOfModes(void* p_sid) {
	SID_Data* sid = (SID_Data*)p_sid;
	return sid->numModes;
}

//All remaining taylors
//=================================

taylor* getTaylorPtr(void* p_sid, const char* taylorName) {
	SID_Data* sid = (SID_Data*)p_sid;
	taylor* t;
	if (!strcmp(taylorName, "mdCM")) {
		t = &sid->mdCM;
	}
	else if (!strcmp(taylorName, "J")) {
		t = &sid->J;
	}
	return t;
}

int getTaylorOrder(void* p_sid, const char* taylorName) {
	taylor* t = getTaylorPtr(p_sid, taylorName);
	return t->order;
}

int getTaylorNrow(void* p_sid, const char* taylorName) {
	taylor* t = getTaylorPtr(p_sid, taylorName);
	return t->nrow;
}

int getTaylorNcol(void* p_sid, const char* taylorName) {
	taylor* t = getTaylorPtr(p_sid, taylorName);
	return t->ncol;
}

int getTaylorNq(void* p_sid, const char* taylorName) {
	taylor* t = getTaylorPtr(p_sid, taylorName);
	return t->nq;
}

int getTaylorNqn(void* p_sid, const char* taylorName) {
	taylor* t = getTaylorPtr(p_sid, taylorName);
	return t->nqn;
}

int getTaylorStructure(void* p_sid, const char* taylorName) {
	taylor* t = getTaylorPtr(p_sid, taylorName);
	return t->structure;
}

double getTaylorM0(void* p_sid, const char* taylorName, int r, int c, int dimR, int dimC) {
	taylor* t = getTaylorPtr(p_sid, taylorName);
	int i = (c - 1) + (r - 1) * dimC;
	return t->M0[i];
}

double getTaylorM1(void* p_sid, const char* taylorName, int r, int q, int c, int dimR, int dimQ, int dimC) {
	taylor* t = getTaylorPtr(p_sid, taylorName);
	int i = (c - 1) + (r - 1) * dimC + (q - 1) * dimC*dimR;
	return t->M1[i];
}




//printing
//=================================

void printM_2(double* M, int dimRow, int dimCol) {
	for (int r = 1; r <= dimRow; r++) {
		for (int c = 1; c <= dimCol; c++) {
			int i = matrix2Index(r, c, dimRow, dimCol);
			printf("M[%d,%d]:%d = %f\n", r, c, i, M[i]);
		}
		printf("\n");
	}
}

void printM_3(double* M, int dimRow, int dimQ, int dimCol) {
	for (int r = 1; r <= dimRow; r++) {
		for (int c = 1; c <= dimCol; c++) {
			for (int q = 1; q <= dimQ; q++) {
				int i = matrix3Index(r, q, c, dimRow, dimQ, dimCol);
				printf("M[%d,%d,%d]:%d = %f\n", r, q, c, i, M[i]);
			}
		}
	}
}

void printTaylor(taylor t) {
	printf("\t order=%d\n", t.order);
	printf("\tnrow=%d  ncol=%d  nq=%d  nqn=%d structure=%d \n", t.nrow, t.ncol, t.nq, t.nqn, t.structure);
	printf("\t______________M0__________________\n");
	printM_2(t.M0, t.nrow, t.ncol);
	printf("\t______________M1__________________\n");
	printM_3(t.M1, t.nrow, t.nq, t.ncol);
}