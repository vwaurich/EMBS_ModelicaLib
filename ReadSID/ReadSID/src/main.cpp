#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include "../include/SIDtypes.h"

static int lineNumber;

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
	std::cout << "numNodes " << numNodeStr << std::endl;

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
	double decExponent = 0.0f;
	int pos = 0;
	pos = s.find("D+", pos);
	base = std::stod(s.substr(0, pos));
	decExponent = std::stod(s.substr(pos + 2, s.length()));
	return base * pow(10, decExponent);
}

std::string stringBetween(std::string s, std::string start, std::string stop) {
	int pos0 = s.find(start.c_str(), 0);
	int pos1 = s.find(stop.c_str(), pos0);
	return s.substr(pos0+1, pos1 - pos0-1);
}

bool findString(std::string s, std::string findStr, int* pos)
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
		if (findString(line,"mass",&pos)){
			str = stringBetween(line, "=", "/n");
			sid->modeStruct.mass = stringDouble_SID(str);
		}
		else if (findString(line, "nelastq", &pos)) {
			str = stringBetween(line, "=", "/n");
			sid->modeStruct.nelastq = std::stoi(str);
			sid->modeStruct.currModeIdx = 0;
		}
		else if (findString(line, "ielastq", &pos)) {
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
	mat[(r - 1) * numC + c] = value;
}

void setMatrix3D(double value, double* mat, int r, int nq, int c, int numR, int numNQ, int numC)
{
	//traverse the rows, than the qn dim than the cols
	mat[(r - 1) * numC * numNQ + (nq-1)*numC + c] = value;
}

int parseOrigin(std::ifstream* sidFile, SID_Data* sid) {
	std::string line;
	std::stringstream ss(line);
	std::string str;
	int pos = 0;
	while (getline(*sidFile, line))
	{
		try {
			if (findString(line, "order", &pos))
			{
				sid->nodes[sid->currNodeIdx].orig.order = std::stoi(stringBetween(line, "=", "/n"));
			}
			else if (findString(line, "nrow", &pos))
			{
				sid->nodes[sid->currNodeIdx].orig.nrow = std::stoi(stringBetween(line, "=", "/n"));
			}
			else if (findString(line, "ncol", &pos))
			{
				sid->nodes[sid->currNodeIdx].orig.ncol = std::stoi(stringBetween(line, "=", "/n"));
			}
			else if (findString(line, "nq", &pos))
			{
				sid->nodes[sid->currNodeIdx].orig.nq = std::stoi(stringBetween(line, "=", "/n"));
			}
			else if (findString(line, "nqn", &pos))
			{
				sid->nodes[sid->currNodeIdx].orig.nqn = std::stoi(stringBetween(line, "=", "/n"));
			}
			else if (findString(line, "structure", &pos))
			{
				sid->nodes[sid->currNodeIdx].orig.structure = std::stoi(stringBetween(line, "=", "/n"));
				//after structure, lets allocate the matrices
				sid->nodes[sid->currNodeIdx].orig.M0 = (double*)calloc(sid->nodes[sid->currNodeIdx].orig.ncol *
					sid->nodes[sid->currNodeIdx].orig.nrow, sizeof(double));
				sid->nodes[sid->currNodeIdx].orig.M1 = (double*)calloc(sid->nodes[sid->currNodeIdx].orig.ncol *
					sid->nodes[sid->currNodeIdx].orig.nrow *
					sid->nodes[sid->currNodeIdx].orig.nq, sizeof(double));
				sid->nodes[sid->currNodeIdx].orig.Mn = (double*)calloc(sid->nodes[sid->currNodeIdx].orig.ncol *
					sid->nodes[sid->currNodeIdx].orig.nrow *
					sid->nodes[sid->currNodeIdx].orig.nqn, sizeof(double));
			}
			else if (findString(line, "m0", &pos))
			{
				int idxr = std::stoi(stringBetween(line, "(", ","));
				int idxc = std::stoi(stringBetween(line, ",", ")"));
				double val = stringDouble_SID(stringBetween(line, "=", "/n"));
				setMatrix2D(val, sid->nodes[sid->currNodeIdx].orig.M0, idxr, idxc, sid->nodes[sid->currNodeIdx].orig.nrow, sid->nodes[sid->currNodeIdx].orig.ncol);
			}
			else if (findString(line, "m1", &pos))
			{
				std::string idxRs = stringBetween(line, "(", ",");
				std::string idxNQs = stringBetween(line, ",", ",");
				std::string idxCs = stringBetween(line, ",", ")");

				int idxr = std::stoi(idxRs);
				int idxnq = std::stoi(idxNQs);
				int idxc = std::stoi(idxCs);
				double val = stringDouble_SID(stringBetween(line, "=", "/n"));
				setMatrix3D(val, sid->nodes[sid->currNodeIdx].orig.M1, idxr, idxnq, idxc, sid->nodes[sid->currNodeIdx].orig.nrow, sid->nodes[sid->currNodeIdx].orig.nq, sid->nodes[sid->currNodeIdx].orig.ncol);
			}
			else if (findString(line, "mn", &pos))
			{
				int idxr = std::stoi(stringBetween(line, "(", ","));
				int idxnqn = std::stoi(stringBetween(line, ",", ","));
				int idxc = std::stoi(stringBetween(line, ",", ")"));
				double val = stringDouble_SID(stringBetween(line, "=", "/n"));
				setMatrix3D(val, sid->nodes[sid->currNodeIdx].orig.Mn, idxr, idxnqn, idxc, sid->nodes[sid->currNodeIdx].orig.nrow, sid->nodes[sid->currNodeIdx].orig.nqn, sid->nodes[sid->currNodeIdx].orig.ncol);
			}
			else if (findString(line, "end origin", &pos)) {
				return 0;
			}
		}
		catch (std::exception e)
		{
			std::cout << "exception" << std::endl;
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
		if (findString(line, "rframe", &pos))
		{
			str = stringBetween(line, "=", "/n");
			memcpy(sid->nodes[sid->currNodeIdx].rFrame, str.c_str(), 8);
		}
		else if (findString(line, "origin", &pos))
		{
			parseOrigin(sidFile, sid);
		}
		else if (findString(line, "phi", &pos))
		{
		}
		else if (findString(line, "psi", &pos))
		{
		}
		else if (findString(line, "AP", &pos))
		{
		}	
		else if (findString(line, "end node", &pos)) {
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
		if (findString(line, "new node", &pos)) {
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
		if (findString(line, "refmod", &pos))
		{
			parseRefMod(sidFile, sid);
		}
		else if (findString(line, "frame", &pos))
		{
			parseFrame(sidFile, sid);
		}
		else if (findString(line,"end modal",&pos)) {
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
		if (findString(line, "new modal", &pos))
		{
			parseModal(sidFile, sid);
		}
		else if(findString(line, "end part", &pos))	{
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
		if (findString(line, "part", &pos))
		{
			parsePart(sidFile, sid);
		}
	}
	return 0;
}

int main(int arg_count, char* arg_vec[]) {
	std::string line;
	std::ifstream sidFile("C:/Program Files/Dymola 2020x/Modelica/Library/FlexibleBodies 2.3.0/Resources/Data/cartopPragV32.SID_FEM");

	lineNumber = 1;
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
	return 0;
}
