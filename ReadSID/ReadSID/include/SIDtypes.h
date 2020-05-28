
typedef struct
{
	char ielastq[20]; //name of the mode shape
	double freq;
} modeShape;


typedef struct
{
	double mass; //mass of body
	int nelastq; //number of mode shapes
	modeShape* ielastq;
	int currModeIdx;
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
} origin; // = class taylor

typedef struct
{
	int id;
	char rFrame[8];
	origin orig;
} node;

typedef struct {
	int numNodes;
	int numModes;
	bool useGeoStiffness;
	refmod modeStruct;
	node* nodes;
	int currNodeIdx;
} SID_Data;

