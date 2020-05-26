#include <stdio.h>
#include <stdlib.h>

//==========Typedefs=================

typedef enum { FALSE = 0, TRUE } boolean;

typedef struct
{
	int eigenmode;
	double freq;
} mode;

typedef struct
{
	int id;
} node;

typedef struct  {
	int numNodes;
	int numModes;
	boolean useGeoStiffness;
	double mass;
	mode* modes;
	node* nodes;
} SID_Data;


//==========Functions=================

SID_Data initSIDstruct()
{
	SID_Data sid;
	sid.numNodes = 0;
	sid.numModes = 0;
	sid.useGeoStiffness = FALSE;
	return sid;
}


boolean isDigit(char* str)
{
	for (int i = 0; i<strlen(str); i++)
		if (!isdigit(str[i]))
		{
			return FALSE;
		}
	return TRUE;
}


int stringInt(char* str)
{
	return (int)strtol(str, (char **)NULL, 10);
}


int allocateStructure(char* line, SID_Data* sid)
{
	char * pch;
	pch = strtok(line, " ,");


	int i = 1;
	while (pch != NULL)
	{
		//numNodes and numModes
		if (isDigit(pch)){
			if (i == 1)	{
				sid->numNodes = stringInt(pch);

				i++;
			}
			else if (i == 2){
				sid->numModes = stringInt(pch);
				i++;
			}
		}
		//geometric stiffness
		if (!strcmp(pch, "Geo_Stiff=")){
			pch = strtok(NULL, " ,.-");
			if (!strcmp(pch, "yes")){
				sid->useGeoStiffness = TRUE;
			}
		}
		pch = strtok(NULL, " ,.-");//new starting location for previous token search
	}
	// allocate arrays in SID
	sid->modes = calloc(sid->numModes, sizeof(mode));
	sid->nodes = calloc(sid->numNodes, sizeof(node));

	return sid;
}

int processLine(char* line, SID_Data* sid)
{

}


int initStructure(char* fileName)
{
	SID_Data sid = initSIDstruct();
	boolean allocatedStructure = FALSE;

	//open file
	FILE *fp = fopen(fileName, "r");
	if (fp == NULL) {
		perror("Unable to open file!");
		exit(1);
	}

	// get single line as a string
	char chunk[256];
	size_t lineLength = sizeof(chunk);
	char *line = (char*)malloc(lineLength);
	if (line == NULL) {
		perror("Unable to allocate memory for the line buffer.");
		exit(1);
	}
	line[0] = '\0';
	while (fgets(chunk, sizeof(chunk), fp) != NULL) {
		// Resize the line buffer if necessary
		size_t len_used = strlen(line);
		size_t chunk_used = strlen(chunk);

		if (lineLength - len_used < chunk_used) {
			lineLength *= 2;
			if ((line = (char*)realloc(line, lineLength)) == NULL) {
				perror("Unable to reallocate memory for the line buffer.");
				free(line);
				exit(1);
			}
		}
		// Copy the chunk to the end of the line buffer
		strncpy(line + len_used, chunk, lineLength - len_used);
		len_used += chunk_used;

		// Check if line contains '/n', if yes process the line of text
		if (line[len_used - 1] == '\n') {
			//process line
			if (!allocatedStructure)
			{
				allocateStructure(line, &sid);
				allocatedStructure = TRUE;
			}
			else
			{
				processLine(line, &sid);
			}
			line[0] = '/0';
		}
	} //while 
	fclose(fp);
	free(line);
	return 0;
}


int main() {
	initStructure("C:/Program Files/Dymola 2020x/Modelica/Library/FlexibleBodies 2.3.0/Resources/Data/cartopPragV32.SID_FEM");
 }