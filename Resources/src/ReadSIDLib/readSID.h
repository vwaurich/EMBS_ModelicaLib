
#ifdef EXPORTLIB
   #define LIB_API __declspec(dllexport)
#else
   #define LIB_API  __declspec(dllimport)
#endif


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
	int currNodeIdx;
} SID_Data;



#ifdef __cplusplus
extern "C"
{
#endif

	LIB_API void* SIDFileConstructor(char* fileName);
	LIB_API int getM3Idx(int r, int q, int c, int dimR, int dimQ, int dimC);//zero based
	LIB_API int getM2Idx(int r, int c, int dimR, int dimC);//zero based

#ifdef __cplusplus
}
#endif


