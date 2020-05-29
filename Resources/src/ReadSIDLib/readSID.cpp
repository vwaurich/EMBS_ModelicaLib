#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include "readSID.h"

int getM2Idx(int r, int c, int dimR, int dimC) {
	return c + r * dimC;
}

int getM3Idx(int r, int q, int c, int dimR, int dimQ, int dimC) {
	return c + r * dimC + q * dimC*dimR;
}


void printM_2(double* M, int dimRow, int dimCol) {
	for (int r = 0; r < dimRow; r++) {
		for (int c = 0; c < dimCol; c++) {
			std::cout << M[getM2Idx(r,c,dimRow,dimCol)] << " \t ";
		}
		std::cout << std::endl;
	}
}


void printM_3(double* M, int dimRow, int dimNQ, int dimCol) {
	for (int r = 0; r < dimRow; r++) {
		for (int c = 0; c < dimCol; c++) {
			for (int nq = 0; nq < dimNQ; nq++) {
				std::cout << "M[" << r+1 << ", " << nq+1 << ", " << c+1 << "] = " << M[getM3Idx(r,nq,c,dimRow,dimNQ,dimCol)] << std::endl;
			}
		}
	}
}

void printTaylor(taylor t) {
	std::cout << "\t" << "order:" << t.order << ")" << std::endl;
	std::cout << "\t" << "ncol:" << t.ncol << " nrow:" << t.nrow << " nqn:" << t.nq << " n_nqn" << t.nqn << std::endl;
	std::cout << "\t" << "______________M0__________________" << std::endl;
	printM_2(t.M0, t.nrow, t.ncol);
	std::cout << "\t" << "______________M1__________________" << std::endl;
	printM_3(t.M1, t.nrow, t.nq, t.ncol);
}


int allocateSIDstructure(std::string line, SID_Data* sid)
{
	std::stringstream ss(line);
	std::string numNodeStr, numModeStr, geoStiffStr, str;
	bool useGeoStifffness = false;
	ss >> numNodeStr;
	ss >> numModeStr;
	while (ss >> str)
	{
		if (str == "Geo_Stiff=") {
			ss >> str;
			if (str == "yes") {
				useGeoStifffness = true;
			}
		}
	}
	sid->numNodes = std::stoi(numNodeStr);
	sid->numModes = std::stoi(numModeStr);
	sid->useGeoStiffness = useGeoStifffness;
	sid->modeStruct.mass = 0;
	sid->modeStruct.nelastq = sid->numNodes;
	sid->modeStruct.ielastq = (modeShape*)calloc(sid->numModes, sizeof(modeShape));

	sid->nodes = (node*)calloc(sid->numNodes, sizeof(node));
	sid->currNodeIdx = -1;//when encountering a new node, increase to 0
	return 0;
}


double stringDouble_SID(std::string s) {
	double base = 0.0;
	double decExponent = 0.0;
	double ret = 0.0;
	int pos = 0;
	pos = s.find("D+", 0);
	if(pos != std::string::npos)
	{
		base = std::stod(s.substr(0, pos));
		decExponent = std::stod(s.substr(pos + 2, s.length()));
		ret = base * pow(10, decExponent);
	}
	else {
		pos = s.find("D-", 0);
		if (pos != std::string::npos)
		{
			base = std::stod(s.substr(0, pos));
			decExponent = std::stod(s.substr(pos + 2, s.length()));
			ret = base * pow(10, -decExponent);
		}
	}

	return ret;
}

std::string stringBetween(std::string s, std::string start, std::string stop) {
	int pos0 = s.find(start.c_str(), 0);
	int pos1 = s.find(stop.c_str(), pos0);
	return s.substr(pos0+1, pos1 - pos0-1);
}

void getIndeces2(std::string s, int* i1, int* i2) {
	int pos = s.find("(", 0);
	std::string idxStr;
	std::string subs = s.substr(pos + 1, 8);
	std::stringstream ss(subs);
	ss >> idxStr;
	*i1 = std::stoi(idxStr);
	ss >> idxStr;
	*i2 = std::stoi(idxStr);
}

void getIndeces3(std::string s, int* i1, int* i2, int* i3){
	int pos = s.find("(", 0);
	std::string idxStr;
	std::string subs = s.substr(pos+1, 8);
	std::stringstream ss(subs);
	ss >> idxStr;
	*i1 = std::stoi(idxStr);
	ss >> idxStr;
	*i2 = std::stoi(idxStr);
	ss >> idxStr;
	*i3 = std::stoi(idxStr);
}

bool findString(std::string s, std::string findStr)
{
	bool b=false;
	try {
		int p = (int)s.find(findStr.c_str(), 0);
		b = std::string::npos != p;
	}
	catch (std::exception e)
	{
		return b;
	}
	return b;
}

int parseRefMod(std::ifstream* sidFile, SID_Data* sid) {
	std::string line;
	std::stringstream ss(line);
	std::string str;
	int pos = 0;
	while (getline(*sidFile, line)) {
		if (findString(line,"mass")){
			str = stringBetween(line, "=", "/n");
			sid->modeStruct.mass = stringDouble_SID(str);
		}
		else if (findString(line, "nelastq")) {
			str = stringBetween(line, "=", "/n");
			sid->modeStruct.nelastq = std::stoi(str);
			sid->modeStruct.currModeIdx = 0;
		}
		else if (findString(line, "ielastq")) {
			int arrIdx = std::stoi(stringBetween(line, "(", ")"));
			std::string modeIdx = stringBetween(line, "=", ":");
			double freq = std::stod(stringBetween(line, ":", "Hz"));
			modeShape modeShape;
			modeShape.freq = freq;
			memcpy(modeShape.ielastq, modeIdx.c_str(), modeIdx.length());
			sid->modeStruct.ielastq[arrIdx - 1] = modeShape;
			sid->modeStruct.currModeIdx++;
		}
		else if (line.find("end refmod")) {
			return 0;
		}
	}
	return 0;
}


void setMatrix2D(double value, double* mat, int r, int c, int numR, int numC)
{
	mat[getM2Idx(r,c,numR,numC)] = value;
}

void setMatrix3D(double value, double* mat, int r, int nq, int c, int numR, int numNQ, int numC)
{
	//traverse the cols,rows and q
	mat[getM3Idx(r,nq,c,numR,numNQ,numC)] = value;
}

int parseTaylor(std::ifstream* sidFile, taylor* tayl) {
	std::string line;
	std::stringstream ss(line);
	std::string str;
	int pos = 0;
	while (getline(*sidFile, line))
	{
		if (findString(line, "order"))
		{
			tayl->order = std::stoi(stringBetween(line, "=", "/n"));
		}
		else if (findString(line, "nrow"))
		{
			tayl->nrow = std::stoi(stringBetween(line, "=", "/n"));
		}
		else if (findString(line, "ncol"))
		{
			tayl->ncol = std::stoi(stringBetween(line, "=", "/n"));
		}
		else if (findString(line, "nq "))//Attention, the space seperates the token from nqn
		{
			tayl->nq = std::stoi(stringBetween(line, "=", "/n"));
		}
		else if (findString(line, "nqn "))
		{
			tayl->nqn = std::stoi(stringBetween(line, "=", "/n"));
		}
		else if (findString(line, "structure"))
		{
			tayl->structure = std::stoi(stringBetween(line, "=", "/n"));
			//after structure, lets allocate the matrices
			tayl->M0 = (double*)calloc(tayl->ncol *
				tayl->nrow, sizeof(double));
			tayl->M1 = (double*)calloc(tayl->ncol *
				tayl->nrow *
				tayl->nq, sizeof(double));
			tayl->Mn = (double*)calloc(tayl->ncol *
				tayl->nrow *
				tayl->nqn, sizeof(double));
		}
		else if (findString(line, "m0"))
		{
			int idxr=0;
			int idxc=0;
			getIndeces2(line, &idxr, &idxc);
			double val = stringDouble_SID(stringBetween(line, "=", "/n"));
			setMatrix2D(val, tayl->M0, idxr-1, idxc-1, tayl->nrow, tayl->ncol);
		}
		else if (findString(line, "m1"))
		{
			int idxr = 0;
			int idxnq = 0;
			int idxc = 0;
			getIndeces3(line, &idxr, &idxnq, &idxc);
			double val = stringDouble_SID(stringBetween(line, "=", "/n"));
			setMatrix3D(val, tayl->M1, idxr-1, idxnq-1, idxc-1, tayl->nrow, tayl->nq, tayl->ncol);
		}
		else if (findString(line, "mn"))
		{
			int idxr = std::stoi(stringBetween(line, "(", ","));
			int idxnqn = std::stoi(stringBetween(line, ",", ","));
			int idxc = std::stoi(stringBetween(line, ",", ")"));
			double val = stringDouble_SID(stringBetween(line, "=", "/n"));
			setMatrix3D(val, tayl->Mn, idxr-1, idxnqn-1, idxc-1, tayl->nrow, tayl->nqn, tayl->ncol);
		}
		else if (findString(line, "end")) {
			return 0;
		}
	}
	return 0;
}

int parseNode(std::ifstream* sidFile, SID_Data* sid) {
	std::string line;
	std::stringstream ss(line);
	std::string str;
	int pos = 0;
	while (getline(*sidFile, line))
	{
		if (findString(line, "rframe"))
		{
			str = stringBetween(line, "=", "/n");
			memcpy(sid->nodes[sid->currNodeIdx].rFrame, str.c_str(), 8);
		}
		else if (findString(line, "origin"))
		{
			parseTaylor(sidFile, &sid->nodes[sid->currNodeIdx].orig);
		}
		else if (findString(line, "phi"))
		{
			parseTaylor(sidFile, &sid->nodes[sid->currNodeIdx].phi);
		}
		else if (findString(line, "psi"))
		{
			parseTaylor(sidFile, &sid->nodes[sid->currNodeIdx].psi);
		}
		else if (findString(line, "AP"))
		{
			parseTaylor(sidFile, &sid->nodes[sid->currNodeIdx].AP);
		}	
		else if (findString(line, "end node")) {
			return 0;
		}
	}
	return 0;
}


int parseFrame(std::ifstream* sidFile, SID_Data* sid) {
	std::string line;
	std::stringstream ss(line);
	std::string str;
	int pos = 0;
	while (getline(*sidFile, line)) {
		if (findString(line, "new node")) {
			sid->currNodeIdx++;
			sid->nodes[sid->currNodeIdx].id = std::stoi(stringBetween(line,"=", "/n"));
			parseNode(sidFile, sid);
		}
		else if (line.find("end frame")) {
			return 0;
		}
	}
	return 0;
}


int parseModal(std::ifstream* sidFile, SID_Data* sid) {
	std::string line;
	std::stringstream ss(line);
	std::string str;
	int pos = 0;
	while (getline(*sidFile, line)) 
	{
		if (findString(line, "refmod"))
		{
			parseRefMod(sidFile, sid);
		}
		else if (findString(line, "frame"))
		{
			parseFrame(sidFile, sid);
		}
		else if (findString(line, "mdCM"))
		{
			parseTaylor(sidFile, &sid->mdCM);
		}
		else if (findString(line, "J"))
		{
			parseTaylor(sidFile, &sid->J);
		}
		else if (findString(line, "Ct"))
		{
			parseTaylor(sidFile, &sid->Ct);
		}
		else if (findString(line, "Cr"))
		{
			parseTaylor(sidFile, &sid->Cr);
		}
		else if (findString(line, "Me"))
		{
			parseTaylor(sidFile, &sid->Me);
		}
		else if (findString(line, "Gr"))
		{
			parseTaylor(sidFile, &sid->Gr);
		}
		else if (findString(line, "Ge"))
		{
			parseTaylor(sidFile, &sid->Ge);
		}
		else if (findString(line, "Oe"))
		{
			parseTaylor(sidFile, &sid->Oe);
		}
		else if (findString(line, "ksigma"))
		{
			parseTaylor(sidFile, &sid->ksigma);
		}
		else if (findString(line, "Ke"))
		{
			parseTaylor(sidFile, &sid->Ke);
		}
		else if (findString(line, "De"))
		{
			parseTaylor(sidFile, &sid->De);
		}
		else if (findString(line,"end modal")) {
			return 0;
		}
	}
	return 0;
}

int parsePart(std::ifstream* sidFile, SID_Data* sid) {
	std::string line;
	std::stringstream ss(line);
	std::string str;
	int pos = 0;
	while (getline(*sidFile, line)) {
		if (findString(line, "new modal"))
		{
			parseModal(sidFile, sid);
		}
		else if(findString(line, "end part"))	{
			return 0;
		}
	}
	return 0;
}

int processLine(std::ifstream* sidFile, SID_Data* sid) {
	std::string line;
	std::stringstream ss(line);
	std::string str;
	int pos = 0;
	while (getline(*sidFile, line)) {
		if (findString(line, "part"))
		{
			parsePart(sidFile, sid);
		}

	}
	return 0;
}

void* SIDFileConstructor(char* fileName)
{
	std::string line;
	std::ifstream sidFile(fileName);

	SID_Data sid;
	if (sidFile.is_open())
	{
		getline(sidFile, line);
		//first line holds information about number of nodes and number of modes
		allocateSIDstructure(line, &sid);
		//process the rest of the file;
		processLine(&sidFile, &sid);
		sidFile.close();
	}
	else { std::cout << "Unable to open file"<<std::endl; }
	return (void*)&sid;
}
