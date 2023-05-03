// File: graphMain.cpp
// Author: Andrew Westlund
// Date: 4/21/23
// Pupose: To receive a file with lines of edges, store it into a graph, square the graph,
// all for incremental updates of the squared graph, and perform primitives on the squared graph.

#include <cstdio>
#include <ctype.h>			// Provides isdigit
#include <iostream>			// Provides cout
#include <fstream>
#include <cstring>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include "graph.h"
#include "lluAVL.h"

using namespace std;

// Find the true difference between two times
double msDiffTime(struct timespec start, struct timespec stop) {
    long seconds = stop.tv_sec - start.tv_sec;
    long ns = stop.tv_nsec - start.tv_nsec;

    if(ns < 0) {
    	seconds--;
        ns += 1000000000;
	}
    return (1000.0 * ((double)seconds + (double)ns/(double)(1000000000)));
}

UINT16 getTypeID(UINT8 *typeArr, vID v) {
	return ((UINT16)(typeArr[v.id]));
}

// Preconditions: v.id must be smaller than numOfVertices (size of typeArr)
bool setTypeID(UINT8 *typeArr, vID v, UINT8 typeID) {
	bool rc = false;
	// If the type has not been set yet, set it
	if (getTypeID(typeArr, v) == 0) {
		typeArr[v.id] = typeID;
		rc = true;
	}
	return rc;
}

// Extend the graph class
class bvGraph : public graph {
	public:
		// Keep the same constructor and destructor
		bvGraph(UINT64 numOfVertices) : graph(numOfVertices) {}
		~bvGraph() {graph::~graph();}
		
		// Returns the biggest connection between any two entities of the graph
		UINT32 getMaxCooc(UINT64 numOfEntities, UINT64 &maxEnt1, UINT64 &maxEnt2) {
			vID z;
			Edge *eNList;
			UINT32 maxCooc = 0;
			UINT32 outDegree;
			UINT64 maxA, maxB;
			// For every entity on the graph
			for (z.id = 0; z.id < numOfEntities; z.id++) {
				eNList = getVneighbors(z);
				outDegree = getVoutDegree(z);
				// For every entity that entity shares a connection with
				for (UINT32 c = 0; c < outDegree; c++) {
					if (eNList[c].weight > maxCooc) {
						maxCooc = eNList[c].weight;
						maxA = z.id;
						maxB = eNList[c].vid.id;
					}
				}
			}
			maxEnt1 = maxA;
			maxEnt2 = maxB;
			return maxCooc;
		}
		
		// Expand will find the entities with the greatest coocs with the given entity that fit the type
		void expand(vID u, UINT16 typeID, Edge * expandArray, UINT16 expandCapacity, UINT16 &expandSize, UINT8 * eTypeArray) {
			// Get the entities that U shares a paper with
			UINT32 posU = 0;
			UINT32 outDegreeU = getVoutDegree(u);
			Edge *uNeighbors = getVneighbors(u);
			// cout << "outDegreeU: " << outDegreeU << endl;
			
			UINT16 arraySize = 0;
			// For every cooc entity of U
			while (posU < outDegreeU) {
				// If the type of the current entity is what we are looking for
				if (getTypeID(eTypeArray, uNeighbors[posU].vid) == typeID) {
					UINT64 currID = uNeighbors[posU].vid.id;
					UINT32 currExpandWeight = uNeighbors[posU].weight;
					// cout << "Current cooc ID: " << currEdge.vid.id << " - Current cooc weight: " << currEdge.weight << endl;
					int i;
					UINT64 tempID;
					UINT32 tempWeight;
					// Find if the entity's cooc is big enough to take a spot on the array and shift the others down;
					// If the array is full, drop the entity with the weakest cooc
					for (i = 0; i < arraySize; i++) {
						if (expandArray[i].weight < currExpandWeight) {
							tempID = expandArray[i].vid.id;
							tempWeight = expandArray[i].weight;
							expandArray[i].vid.id = currID;
							expandArray[i].weight = currExpandWeight;
							currID = tempID;
							currExpandWeight = tempWeight;
						}
					}
					// If the array is not full, increase the size and don't drop an entity
					if (i < expandCapacity) {
						expandArray[i].vid.id = currID;
						expandArray[i].weight = currExpandWeight;
						arraySize++;
					}
						
					// for (int j = 0; j < arraySize; j++) {
						// cout << "Position " << j << ": Entity ID: " << coocArray[j].vid.id << " Cooc Weight: " << coocArray[j].weight << endl;
					// }
				}
				posU++;
			}
			expandSize = arraySize;
		}
		
		// Given a master entity and an array of entities, connect will put all the entities from the connectEIDArray that 
		// fit the type and have a cooc with the master entity on the connectArray.
		void connect(vID u, UINT16 typeID, UINT32 * connectEIDArray, Edge * connectArray, UINT16 connectCapacity, UINT16 &connectSize, UINT8 * eTypeArray) {
			UINT32 posList = 0;
			UINT32 arraySize = 0;
			vID eID;
			
			// Get the array of entities that U shares a paper with and its size
			UINT32 outDegreeU = getVoutDegree(u);
			Edge *uNeighbors = getVneighbors(u);
			
			// For every entity on the given array
			for (; posList < connectCapacity; posList++) {
				eID.id = connectEIDArray[posList];
				if (getTypeID(eTypeArray, eID) == typeID) {	// Check if the entity matches the type
					
					// Binary Search to find if the current entity coocs the master entity and its cooc
					// Limits are there to prevent overflow
					UINT32 lowerBound, upperBound, lowerBoundLimit, upperBoundLimit, mid;
					lowerBound = 0;
					lowerBoundLimit = pow(2, (8 * sizeof(outDegreeU))) - 1;
					// Set the upper bound to start as the last possible position of the neighbor list
					upperBoundLimit = 0;
					upperBound = outDegreeU - 1;

					// fprintf(stderr, "Lower Bound: %u, Upper Bound: %u\n", lowerBound, upperBound);
					if (outDegreeU > 0) {
						// While the indicated id has not been found and the entire list has not been searched
						while (lowerBound <= upperBound) {
							// Set mid to be inbetween the lower and upper bounds
							mid = lowerBound + ((upperBound - lowerBound) / 2);
						
							if ((uNeighbors[mid]).vid.id > eID.id) {
								if (mid == upperBoundLimit) {
									break;
								}
								upperBound = mid - 1;
							} else if ((uNeighbors[mid]).vid.id < eID.id) {
								if (mid == lowerBoundLimit) {
									break;
								}
								lowerBound = mid + 1;
							} else {
								// If the two entities are connected, put it on the array
								connectArray[arraySize].vid.id = (uNeighbors[mid]).vid.id;
								connectArray[arraySize].weight = (uNeighbors[mid]).weight;
								arraySize++;
								break;
							}
						}
					}
				}
			}
			connectSize = arraySize;
		}
		
		// Get the array of the entities with the greatest shared cooc to both u and v that fits a given type
		void bridge(vID u, vID v, UINT16 typeID, Edge * coocArray, UINT16 coocCapacity, UINT16 &coocSize, UINT8 * eTypeArray) {	
			// Get the list of entities that u and v point to
			UINT32 posU = 0;
			UINT32 outDegreeU = getVoutDegree(u);
			Edge *uNeighbors = getVneighbors(u);
	
			UINT32 posV = 0;
			UINT32 outDegreeV = getVoutDegree(v);
			Edge *vNeighbors = getVneighbors(v);

			UINT16 arraySize = 0;
			// Initialize the values of the array
			// for (int a = 0; a < coocCapacity; a++) {
				// coocArray[a].vid.id = 0;
			// }
	
			// For each entity that u or v point to (entities that share a common parent paper with u or v)
			while ((posU < outDegreeU) && (posV < outDegreeV)) {
				// If u and v's current neighbors do not match, keep going on the entity that is behind
				if (uNeighbors[posU].vid.id < vNeighbors[posV].vid.id) {
					posU++;
				} else if (uNeighbors[posU].vid.id > vNeighbors[posV].vid.id) {
					posV++;
				} else {	// If u and v share a common neighbor
					// If the type of the current entity is what we are looking for
					if (getTypeID(eTypeArray, uNeighbors[posU].vid) == typeID) {
						UINT64 currID = uNeighbors[posU].vid.id;	// eID that is pointed by both u and v
						UINT32 currCoocWeight = uNeighbors[posU].weight + vNeighbors[posV].weight;
						// UINT32 currCoocWeight = uNeighbors[posU].weight;
						// UINT32 currCoocWeight = vNeighbors[posV].weight;
				
						int i;
						UINT64 tempID;
						UINT32 tempWeight;
						// Find if the current cooced entity deserves a spot on the array and shift everything down;
						// If the array is full, the weakest cooc entity is dropped
						for (i = 0; i < arraySize; i++) {
							if (coocArray[i].weight < currCoocWeight) {
								tempID = coocArray[i].vid.id;
								tempWeight = coocArray[i].weight;
								coocArray[i].vid.id = currID;
								coocArray[i].weight = currCoocWeight;
								currID = tempID;
								currCoocWeight = tempWeight;
							}
						}
						// If the array is not yet full, increase the size
						if (i < coocCapacity) {
							coocArray[i].vid.id = currID;
							coocArray[i].weight = currCoocWeight;
							arraySize++;
						}
					}
					posU++;
					posV++;
				}
			}
			coocSize = arraySize;
		}
		
		// Get an array of papers that point both to U and V
		void bibliography(vID u, vID v, UINT64 * biblioArray, UINT32 biblioCapacity, UINT32 &biblioSize) {
			// Get the arrays of the papers that point to U and V and their sizes
			UINT32 posU = 0;
			UINT32 inDegreeU = getVinDegree(u);
			Edge *uBackNeighbors = getVbackNeighbors(u);
	
			UINT32 posV = 0;
			UINT32 inDegreeV = getVinDegree(v);
			Edge *vBackNeighbors = getVbackNeighbors(v);
			
			UINT32 arraySize = 0;
			
			// For every paper that points to U or V
			while ((posU < inDegreeU) && (posV < inDegreeV)) {
				// cout << "uBackNeighbors[posU].vid.id: " << uBackNeighbors[posU].vid.id << " - vBackNeighbors[posV].vid.id: " << vBackNeighbors[posV].vid.id << endl;
				
				// If u and v's current papers do not match, make the entity with the smaller paper move forward
				if (uBackNeighbors[posU].vid.id < vBackNeighbors[posV].vid.id) {
					posU++;
				} else if (uBackNeighbors[posU].vid.id > vBackNeighbors[posV].vid.id) {
					posV++;
				} else {	// If u and v share a common paper
					// Append the paper to the end of the array
					biblioArray[arraySize] = uBackNeighbors[posU].vid.id;
					// cout << "biblioArray[" << arraySize << "]: " << biblioArray[arraySize];
					
					// Once the paper is full, break out of the loop
					if (arraySize >= biblioCapacity) {
						break;
					}
					arraySize++;
					posU++;
					posV++;
				}
			}
			// return coocArray;
			biblioSize = arraySize;
		}
		
		friend class bvAVL;
};

// Extends AVL Tree structure that stores unsigned long longs
class bvAVL : public lluAVL {
	private:
		// Preconditions: Each node stores an id to the paper that the eID is an edge to
		// Postconditions: Each "edge" on the tree will be transferred to the graph in numerical order of the paper ID
		void addGEdges(bvGraph *bvg, vID eID, node *p) {
			if (p) {
				addGEdges(bvg, eID, p->left);
				vID pID;
				pID.id = p->val;
				bvg->addEdge2(pID, eID);
				addGEdges(bvg, eID, p->right);
				// Clear the tree as you go so the tree is available for the next eTreeArray
				delete p;
			}
		}
		
		// Preconditions: Each node stores an id to an entity that the eID1 entity has a shared paper with
		// Postconditions: Each G2 edge on the tree will be transferred to the graph in numerical order of the entity IDs.
		// For every entity node that eID1 shares a paper with, that entity has a tree where eID1 is a node.
		void addG2Edges(bvGraph *bvg, vID eID1, node *p) {
			if (p) {
				addG2Edges(bvg, eID1, p->left);
				vID eID2;
				eID2.id = p->val;
				// eID2 will not have a backedge pointing to eID1 since eID2 will also have an edge to eID1
				bvg->addEdge2Forward(eID1, eID2, p->weight);
				addG2Edges(bvg, eID1, p->right);
				// Clear the tree as you go so the tree is available for the next eTreeArray
				delete p;
			}
		}
		
	public:
		// Same constructor and destructor as parent class
		bvAVL() : lluAVL(){}
		~bvAVL() {lluAVL::~lluAVL();}
		
		// Preconditions: Each node stores an id to the paper that the eID is an edge to
		// Postconditions: Each "edge" on the tree will be transferred to the graph in numerical order of the paper ID
		void addGEdges(bvGraph *bvg, vID eID) {
			// Make the space for back-edges to be appended without reallocation.
			// bvg->createV(eID, 0, tCount);
			bvg->createVBackward(eID, tCount);
			// Call recursive function to add G edges
			addGEdges(bvg, eID, root);
			// Clear the tree afterwards for the next time eTreeArray is used
			root = NULL;
			tCount = 0;
		}
		
		// Preconditions: Each node stores an id to an entity that the eID1 entity has a shared paper with
		// Postconditions: Each G2 edge on the tree will be transferred to the graph in numerical order of the entity IDs.
		// For every entity node that eID1 shares a paper with, that entity has a tree where eID1 is a node.
		void addG2Edges(bvGraph *bvg, vID eID1) {
			// Make the space for edges to be appended without reallocation.
			bvg->createVForward(eID1, tCount);
			// Call recursive function to add G2 edges
			addG2Edges(bvg, eID1, root);
			// Clear the tree afterwards to free unnecessary space
			root = NULL;
			tCount = 0;
		}
};


int main(int argc, char **argv) {
	
	cout << "Begin!" << endl;
   
// --------- Graph and generater variable defaults ---------------------

	INT16 opt;						// For reading in user input   
	
	// Stores whether a file and/or incremental file is provided, as well as if those files are binary
	bool dumpToFile = false, hasFile = false, isBin = false, hasInc = false, incBin = false;
	// Stores which primitives the user wants to run
	bool doExpand = false, doConnect = false, doBridge = false, doBiblio = false;
	// Stores arrays of data for each run of a primitive so that multiple can be run at one execution
	UINT16 *expandTypeID;
	UINT32 *expandEID;
	UINT16 expandIterations;
	UINT16 *connectTypeID;
	UINT16 *connectEIDArraySize;
	UINT32 *connectEID;
	UINT32 **connectEIDArray;
	UINT16 connectIterations;
	UINT16 *bridgeTypeID;
	UINT32 *bridgeEID1, *bridgeEID2;
	UINT16 bridgeIterations;
	UINT32 *biblioEID1, *biblioEID2, *biblioCapacity;
	UINT16 biblioIterations;
	
	vID x;							// Temp vIDs
	vID z;							
	UINT16 typeID;					// Temp entity type ID

	INT16 allocation = 16;			// How much to grow neighbor lists by when full
	UINT64 vNum = 32000000;			// Num of total vertices, |Papers| + |Entities|
	UINT64 numOfEntities = 2000000; // Num of entities there can be; paper ids will start where entity ids end
	// Allows us to save space (Uses the max entity ID found from a separate search + 1)
	UINT64 actualEntities = 1110736 + 1;		// 1110736 is the biggest entity used in data
	// UINT64 actualEntities = 6;
	UINT64 numOfPapers = vNum - numOfEntities;
	UINT64 numOfEdges = 1500000000;	// Will be accurately reassigned if binary file
	UINT64 count = 0;				// Used to find how many iterations have passed if COUNT is defined
	// Assume that incremental file will not go over 1000 edges
	UINT64 numOfIncEdges = 1000;	// Will be accurately reassigned if binary file
	FILE *myBinFileIn, *myIncBinFileIn;
	// For times
	struct timespec t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;

	char *file = NULL;
	char *incFile = NULL;
	
	opterr = 0;

// -------------------- Getop Flags----------------------------------

	// while ((opt = getopt(argc, argv, "hf:be:cr:")) != -1) {	// For in-command args
	while ((opt = getopt(argc, argv, "hf:bi:ne:c:r:p:")) != -1) {		// For cin args
		switch (opt) {
			case 'h':
				cout << "Command Flag Options:" << endl;
				cout << "-f [filename]: Read in filename and produce a graph based on the paperID-entityID pairs established in the file." << endl;
				cout << "-b: Read in the file given as a binary file." << endl;
				
				cout << "-i [filename]: Read in incremental data file and update the graph based on the new paperID-entityID pairs established in the file." << endl;
				cout << "-n: Read in the incremental file as a binary file." << endl;
				
				// cout << "-e [EntityID,EntityTypeID]: Perform the expand primitive on graph to find the passed entity's ten greatest connections of the given type." << endl;
				// cout << "\tEntity1ID will fit entity IDs up to 2^32 - 1, and the EntityTypeID will fit entity type IDs up to 2 ^ 16 - 1." << endl;
				cout << "-e [numOfExpands]: The user will be prompted for an Entity ID and an Entity ID Type. The expand primitive will be performed on the graph to" << endl;
				cout << "\tfind the given entity's ten greatest connections to other entities of the given type." << endl;
				cout << "\tEntity1ID will fit entity IDs up to 2^32 - 1, and the EntityTypeID will fit entity type IDs up to 2^16 - 1." << endl;
				
				cout << "-c [numOfConnects]: The user will be promted for an Entity ID, an EntityID Type and a list of Entity IDs of a size specified by the user." << endl;
				cout << "\tThe connect primitive will be performed on the graph to see which entities of the passed list are of the given type and have" << endl;
				cout << "\ta cooc with the initially passed entity." << endl;
				cout << "\tAll EntityIDs can fit entity IDs up to 2^32 - 1, and the EntityTypeID will fit entity type IDs up to 2 ^ 16 - 1." << endl;
				
				// cout << "-r [Entity1ID,Entity2ID,EntityTypeID]: Perform the bridge primitive on graph to find the ten largest coocs shared between the two passed entitis of the given type." << endl;
				// cout << "\tEntity1ID and Entity2ID will both fit entity IDs up to 2^32 - 1, and the EntityTypeID will fit entity type IDs up to 2 ^ 16 - 1." << endl;
				cout << "-r [numOfBridges]: The user will be prompted for two Entity IDs and an Entity ID Type. The bridge primitive will be performed on the graph to" << endl;
				cout << "\tfind the ten entities of the given type with the largest coocs shared between the two passed entities." << endl;
				cout << "\tEntityID1 and EntityID2 will both fit entity IDs up to 2^32 - 1, and the EntityTypeID will fit entity type IDs up to 2 ^ 16 - 1." << endl;
				
				cout << "-p [numOfBibliographies]: The user will be prompted for two Entity IDs and a Size for the bibliography. The bibliography primitive will be performed on the graph to" << endl;
				cout << "\tfind the first Size many papers that the two passed Entity IDs have in common." << endl;
				cout << "\tEntityID1, EntityID2, and Size will all fit numbers from 0 up to 2^32 - 1." << endl;
				abort();
			case 'f':
				hasFile = true;
				file = optarg;
				break;
			case 'b':
			{
				if (!hasFile) {
					fprintf(stderr, "File flag must precede binary flag.");
					abort();
				}
				isBin = true;
				ifstream binFile(file, ios_base::in | ios_base::binary);
				binFile.seekg(0, ios::end);		// Go to the final byte
				// 4 bytes per unsigned int, and 4 ints per edge
				numOfEdges = binFile.tellg() / 16;
				binFile.seekg(0);				// Go back to the first byte of the file
				myBinFileIn = fopen(file, "rb");
				break;
			}

			case 'i':
				if (!hasFile) {
					fprintf(stderr, "File flag must precede increment flag.");
					abort();
				}
				hasInc = true;
				incFile = optarg;
				break;
			case 'n':
			{
				if (!hasInc) {
					fprintf(stderr, "Increment flag must precede increment binary flag.");
					abort();
				}
				incBin = true;
				ifstream incBinFile(incFile, ios_base::in | ios_base::binary);
				incBinFile.seekg(0, ios::end);		// Go to the final byte
				// 4 bytes per unsigned int, and 4 ints per edge
				numOfIncEdges = incBinFile.tellg() / 16;
				incBinFile.seekg(0);				// Go back to the first byte of the file
				myIncBinFileIn = fopen(incFile, "rb");
				break;
			}			
			
			case 'e':
			{
				if (!hasFile) {
					fprintf(stderr, "File flag must precede expand flag.");
					abort();
				}
				// Get how many times the user wants to perform the expand primitive
				expandIterations = atoi(optarg);
				doExpand = true;
				
				// Read in arguments as a comma-separated list at execution command
				/*
				// cout << "sizeof(optarg): " << sizeof(optarg) << " - " << strlen(optarg) << endl;
				char *argument = optarg;
				// cout << "sizeof(argument): " << sizeof(argument) << " - " << strlen(argument) << endl;
				// for (int i = 0; i < sizeof(optarg); i++) {
					// cout << argument[i] << endl;
					// cout << (UINT32)argument[i] << endl;
				// }
				char **numStrings = new char*[2];
				UINT32 *nums = new UINT32[2];
				int commaIndex = 1;		// The 0th char cannot be a comma
				int numStart = 0, numPos;
				int numOfCommas = 0;
				// int argLength = (sizeof(argument) / sizeof(argument[0]));
				int argLength = strlen(argument);
				
				for (; commaIndex < argLength; commaIndex++) {
					if (argument[commaIndex] == ',') {
						numPos = 0;
						numStrings[numOfCommas] = new char[commaIndex - numStart];
						while (numStart < commaIndex) {
							numStrings[numOfCommas][numPos] = argument[numStart];
							numStart++;
							numPos++;
						}
						nums[numOfCommas] = atoi(numStrings[numOfCommas]);
						numStart = commaIndex + 1;
						numOfCommas++;
						if (numOfCommas >= 1) {
							break;
						}
					}
				}
				if (commaIndex == argLength) {
					fprintf(stderr, "Value-list cannot end with a comma.");
					abort();
				}
				if (numOfCommas < 1) {
					fprintf(stderr, "Not enough comma-separated values. Three are necessary.");
					abort();
				} else {
					// while (((UINT32)argument[argLength - 1]) == 0) {
						// argLength--;
					// }
					// cerr << "Num Start: " << numStart << " - argLength: " << argLength << endl;
					numPos = 0;
					numStrings[numOfCommas] = new char[argLength - numStart];
					while (numStart < argLength) {
						numStrings[numOfCommas][numPos] = argument[numStart];
						numStart++;
						numPos++;
					}
					nums[numOfCommas] = atoi(numStrings[numOfCommas]);
				}
				expandEID = nums[0];
				expandTypeID = nums[1];
				*/
				break;
			}
			case 'c':
			{
				if (!hasFile) {
					fprintf(stderr, "File flag must precede connect flag.");
					abort();
				}
				// Get how many times the user wants to perform the connect primitive
				connectIterations = atoi(optarg);
				doConnect = true;
				break;
			}
			case 'r':
			{
				if (!hasFile) {
					fprintf(stderr, "File flag must precede bridge flag.");
					abort();
				}
				// Get how many times the user wants to perform the bridge primitive
				bridgeIterations = atoi(optarg);
				doBridge = true;
				
				// Read in arguments as a comma-separated list at execution command
				/*
				// cout << "sizeof(optarg): " << sizeof(optarg) << " - " << strlen(optarg) << endl;
				char *argument = optarg;
				// cout << "sizeof(argument): " << sizeof(argument) << " - " << strlen(argument) << endl;
				// for (int i = 0; i < sizeof(optarg); i++) {
					// cout << argument[i] << endl;
					// cout << (UINT32)argument[i] << endl;
				// }
				char **numStrings = new char*[3];
				UINT32 *nums = new UINT32[3];
				int commaIndex = 1;		// The 0th char cannot be a comma
				int numStart = 0, numPos;
				int numOfCommas = 0;
				// int argLength = (sizeof(argument) / sizeof(argument[0]));
				int argLength = strlen(argument);
				
				for (; commaIndex < argLength; commaIndex++) {
					if (argument[commaIndex] == ',') {
						numPos = 0;
						numStrings[numOfCommas] = new char[commaIndex - numStart];
						while (numStart < commaIndex) {
							numStrings[numOfCommas][numPos] = argument[numStart];
							numStart++;
							numPos++;
						}
						nums[numOfCommas] = atoi(numStrings[numOfCommas]);
						numStart = commaIndex + 1;
						numOfCommas++;
						if (numOfCommas >= 2) {
							break;
						}
					}
				}
				if (commaIndex == argLength) {
					fprintf(stderr, "Value-list cannot end with a comma.");
					abort();
				}
				if (numOfCommas < 2) {
					fprintf(stderr, "Not enough comma-separated values. Three are necessary.");
					abort();
				} else {
					// while (((UINT32)argument[argLength - 1]) == 0) {
						// argLength--;
					// }
					// cerr << "Num Start: " << numStart << " - argLength: " << argLength << endl;
					numPos = 0;
					numStrings[numOfCommas] = new char[argLength - numStart];
					while (numStart < argLength) {
						numStrings[numOfCommas][numPos] = argument[numStart];
						numStart++;
						numPos++;
					}
					nums[numOfCommas] = atoi(numStrings[numOfCommas]);
				}
				bridgeEID1 = nums[0];
				bridgeEID2 = nums[1];
				bridgeTypeID = (UINT16)nums[2];
				*/
				// cout << bridgeEID1 << " - " << bridgeEID2 << " - " << bridgeTypeID << endl;
				break;
			}
			case 'p':
			{
				if (!hasFile) {
					fprintf(stderr, "File flag must precede bibliography flag.");
					abort();
				}
				// Get how many times the user wants to perform the bibliography primitive
				biblioIterations = atoi(optarg);
				doBiblio = true;
				break;
			}
			default:
				fprintf(stderr, "Unrecognized option.\n");
				abort();
		}
	}
	for (int i = optind; i < argc; i++) {
		fprintf(stderr, "Non-option argument %s\n", argv[i]); 
	}
   
   if(!hasFile) {		// If no input flag given, tell them they need to do something
		cerr << "-f [filename] option required." << endl;	// For now, if there is no input, break
		abort();
	}
	
	// Read In Values
	if (doExpand) {
		// Each index i of the arrays will store the information for the ith performance 
		expandEID = new UINT32[expandIterations];		// Stores the Entity IDs to expand
		expandTypeID = new UINT16[expandIterations];	// Stores the Entity Type ID to compare against

		// Ask for the input as many times as the user has requested in the execute command args
		for (int i = 0; i < expandIterations; i++) {
			cout << "Expand Entity ID: ";
			cin >> expandEID[i];
			cout << "Expand Entity Type ID: ";
			cin >> expandTypeID[i];
			cout << endl;
		}
	}
	
	if (doConnect) {
		// Each index i of the arrays will store the information for the ith performance 
		connectEID = new UINT32[connectIterations];				// Stores the Entity IDs to connect
		connectTypeID = new UINT16[connectIterations];			// Stores the Entity Type ID to compare against
		connectEIDArraySize = new UINT16[connectIterations];	// Stores how many entities from the array connect with the given entity
		connectEIDArray = new UINT32*[connectIterations];		// Stores arrays of entity IDs to check for a connection
		
		// Ask for the input as many times as the user has requested in the execute command args
		for (int i = 0; i < connectIterations; i++) {
			cout << "Connect Entity ID: ";
			cin >> connectEID[i];
			cout << "Connect Entity Type ID: ";
			cin >> connectTypeID[i];
			cout << "Number of Entity IDs are in the list: ";
			cin >> connectEIDArraySize[i];
			
			// Set the ith entity array to hold the specified number of entities
			connectEIDArray[i] = new UINT32[connectEIDArraySize[i]];
			// Enter entity ids one at a time to be put onto the array
			for (int j = 0; j < connectEIDArraySize[i]; j++) {
				cout << "Enter Entity ID " << j << ": ";
				cin >> connectEIDArray[i][j];
			}
			cout << endl;
		}
	}
	
	if (doBridge) {
		bridgeEID1 = new UINT32[bridgeIterations];		// Stores the Entity IDs to bridge
		bridgeEID2 = new UINT32[bridgeIterations];		// Stores the other Entity IDs to bridge
		bridgeTypeID = new UINT16[bridgeIterations];	// Stores the Entity Type ID to compare against
		
		// Ask for the input as many times as the user has requested in the execute command args
		for (int i = 0; i < bridgeIterations; i++) {
			cout << "Bridge Entity ID 1: ";
			cin >> bridgeEID1[i];
			cout << "Bridge Entity ID 2: ";
			cin >> bridgeEID2[i];
			cout << "Bridge Entity Type ID: ";
			cin >> bridgeTypeID[i];
			cout << endl;
		}
	}
	if (doBiblio) {
		biblioEID1 = new UINT32[biblioIterations];		// Stores the Entity IDs to bibliography
		biblioEID2 = new UINT32[biblioIterations];		// Stores the other Entity IDs to bibliography
		biblioCapacity = new UINT32[biblioIterations];	// Stores the number of papers to find
		
		// Ask for the input as many times as the user has requested in the execute command args
		for (int i = 0; i < biblioIterations; i++) {
			cout << "Bibliography Entity ID 1: ";
			cin >> biblioEID1[i];
			cout << "Bibliography Entity ID 2: ";
			cin >> biblioEID2[i];
			cout << "Maximum Number of Papers to Put in the Bibliography: ";
			cin >> biblioCapacity[i];
			cout << endl;
		}
	}
	
// --------------------- Init Graph ----------------------------------------------------
	clock_gettime(CLOCK_MONOTONIC, &t0);

	bvGraph *myGraph = new bvGraph(vNum);
	myGraph->setAllocation(allocation);
	
	clock_gettime(CLOCK_MONOTONIC, &t1);
	fprintf(stdout, "The time taken to allocate the graph: %6.3fms.\n", msDiffTime(t0, t1));
	cout << "Number of edges: " << numOfEdges << endl;
	
// ----------------------- Create Helper Arrays ------------------------------------------------	

	// typeIDArray[entity e's ID] will hold e's first given entity type ID
	UINT8 *typeIDArray = new UINT8[actualEntities];
	
	// eAVLArray[entity e's ID] will hold all the papers that point to entity e 
	bvAVL *eAVLArray = new bvAVL[actualEntities];
	
	// Stores the edges that are read in so they may be processed in parallel;
	// The array[index that is a multiple of 4] is a paper, then the publication year, then the entity type ID,
	// then the entity the paper points to
	UINT32 *edgeArray = new UINT32[numOfEdges * 4];
	
	// The index of each paper's count will be pID - numOfEntities for 0-based indexing;
	// each value will have the number of entities the paper references
	UINT16 *paperCountArray = new UINT16[numOfPapers];
	
// ----------------------- Prepare for Read-in ------------------------------------------------	
	
	UINT64 tempIndex = 0;
	// Set both arrays at the same time
	// for (; tempIndex < numOfEntities; tempIndex++) {
	for (; tempIndex < actualEntities; tempIndex++) {
		typeIDArray[tempIndex] = 0;
		paperCountArray[tempIndex] = 0;
	}
	// Here, tempIndex == numOfEntities, so just fill out the rest
	for (; tempIndex < numOfPapers; tempIndex++) {
		paperCountArray[tempIndex] = 0;
	}
	
	// Set up parallel arrays
	UINT16 numOfThreads;		// Stores how many threads there are
	// entityRange[myThread] will give the smallest entity that myThread can modify and entityRange[myThread+1] - 1 
	// will give the largest entity myThread can modify
	UINT32 *entityRange;
	UINT64 *edgeArrayIndex;		// Will store each threads current index to the edge array
	#ifdef PAR
	{
		#pragma omp parallel
		{
			if (!omp_get_thread_num()) {	// Only set it up once (if it's the 0th thread)
				numOfThreads = omp_get_num_threads();
				
				// We need spare spot for final thread's range
				entityRange = new UINT32[numOfThreads + 1];
				edgeArrayIndex = new UINT64[numOfThreads];
				
				// Stores how many entities each thread can modify
				UINT32 div = actualEntities / (numOfThreads);
				// Stores how many threads get one extra entity
				UINT32 rem = actualEntities % (numOfThreads);
				
				UINT32 currTop = 0;
				UINT16 thread = 0;
				// For each thread
				for (; thread < numOfThreads; thread++) {
					// Set the thread's lowest entity to be currTop
					entityRange[thread] = currTop;
					currTop += div;
					// Give the first rem many threads the extra labor
					if (thread < rem) {
						currTop++;
					}
					// Set every thread to start looking at the beginning of the edgeArray
					edgeArrayIndex[thread] = 0;
				}
				// thread == num_threads
				// Sets the last thread's upper entity to be the biggest entity in the file
				entityRange[thread] = actualEntities;
			}
		}
	}
	#endif
	
	// Test to verify the entity ranges and indices of each thread
	// for (int i = 0; i < numOfThreads; i++) {
		// cerr << "entityRange[" << i << "]: " << entityRange[i] << endl;
		// cerr << "edgeArrayIndex[" << i << "]: " << edgeArrayIndex[i] << endl;
	// }
	// cerr << "entityRange[" << numOfThreads << "]: " << entityRange[numOfThreads] << endl;
	
	clock_gettime(CLOCK_MONOTONIC, &t2);
	fprintf(stdout, "The time taken to prepare for read: %6.3fms.\n", msDiffTime(t1, t2));
	
// ----------------------- Begin Read ------------------------------------------------	
	
	vID pID;						// vID for paper
	vID eID;						// vID for entity
	// In ASCII file, each edge's paperID is read into a, year into b, typeID into c, and entityID into e
	UINT64 a, b, c, d;
	// After the file is done being read, readIndex * 4 == number of edges in the file
	UINT64 readIndex = 0;			// readIndex / 4 determines which edge the reader is currently at
	UINT32 currPaper, currEntity;	// Stores the ID of the current edge's paper and entity
	UINT16 myThread;				// Stores each thread's rank
	
	cout << "\nStart of Read" << endl;
	
// ----------------------- ASCII Read ------------------------------------------------		

	if (!isBin) {	// Assume file to be ASCII unless binary is specified
		fstream myfile(file, ios_base::in);
		while (myfile >> a >> b >> c >> d) {
			// Put the paper's ID in the edge array
			edgeArray[readIndex] = a + numOfEntities;
			// Mark the paper to have one more entity it references
			paperCountArray[a]++;	
			// Put year and type ID into the array
			edgeArray[readIndex + 1] = b;
			edgeArray[readIndex + 2] = c;
			// Put the entity in the array
			edgeArray[readIndex + 3] = d;
			readIndex += 4;
			
			// Set the entity's type if it is not yet set
			eID.id = d;
			if (eID.id > actualEntities) {
				fprintf(stderr, "Number of actual entities declared: %llu - Entity of added edge: %llu", actualEntities, eID.id);
				abort();
			}
			setTypeID(typeIDArray, eID, c);
		}
		myfile.close();
	} else {	// If file specified as binary
	
// ---------------------- Binary Read ------------------------------------------------		

		// Directly transfer the file into the edge array
		fread(edgeArray, sizeof(unsigned int), (numOfEdges * 4), myBinFileIn);
		fclose(myBinFileIn);
		// For every edge
		UINT64 readEdgeTimes4;
		for (UINT64 readEdge = 0; readEdge < numOfEdges; readEdge++) {
			readEdgeTimes4 = readEdge * 4;
			// Mark the paper to have one more entity it references
			paperCountArray[edgeArray[readEdgeTimes4]]++;
			// Adjust the paperID to match the namespace
			edgeArray[readEdgeTimes4] += numOfEntities;
			
			// Set the entity's type if it is not yet set
			eID.id = edgeArray[(readEdgeTimes4) + 3];
			if (eID.id > actualEntities) {
				fprintf(stderr, "Number of actual entities declared: %llu - Entity of added edge: %llu", actualEntities, eID.id);
				abort();
			}
			setTypeID(typeIDArray, eID, edgeArray[(readEdgeTimes4) + 2]);
		}
		readIndex = numOfEdges * 4;
	}
	cout << "End of Read" << endl;
	clock_gettime(CLOCK_MONOTONIC, &t3);
	fprintf(stdout, "The time taken to read: %6.3fms.\n", msDiffTime(t2, t3));
	cout << "Begin Inserting Edges into the Tree to the Graph" << endl;
	
	struct timespec tTree1, tTree2;		// For the times of a certain number of tree insertions
	#ifdef PAR
	{
		#pragma omp parallel private(myThread, currEntity, currPaper, tTree1, tTree2)
		{
			clock_gettime(CLOCK_MONOTONIC, &tTree1);
			myThread = omp_get_thread_num();
			count = 0;				// To see how many edges have been added if COUNT is defined
			
			// For chomping-at-the-bit strategy, the + constant was to ensure the reader was a full edge ahead
			while (edgeArrayIndex[myThread] + 3 < readIndex) {
				#ifdef COUNT 
				{
				if ((count % 100000000 == 0) && (myThread == 0)) {
					clock_gettime(CLOCK_MONOTONIC, &tTree2);
					fprintf(stdout, "The time taken to insert %llu edges into the AVLs: %6.3fms.\n", count, msDiffTime(tTree1, tTree2));
					tTree1 = tTree2;
				}
				count++;
				}
				#endif
				
				// edgeArrayIndex[myThread] gets the index of the paper of the edge that the thread is currently looking at;
				// add 3 to that, and edgeArray[that sum] will the entity the thread is currently looking at
				currEntity = edgeArray[edgeArrayIndex[myThread] + 3];
				
				// Check if the entity fits within the thread's entity range
				if ((currEntity >= entityRange[myThread]) && (currEntity < entityRange[myThread + 1])) {
				// if ((edgeArray[edgeArrayIndex[myThread] + 3] >= entityRange[myThread]) && (edgeArray[edgeArrayIndex[myThread] + 3] < entityRange[myThread + 1])) {
					// Gets the paper of the edge that myThread is currently at in edgeArray
					currPaper = edgeArray[edgeArrayIndex[myThread]];
					
					// The current entity's AVL tree will store the IDs of all the papers that point to the entity
					eAVLArray[currEntity].insertV(currPaper);
					
					// Debugging section
					// #pragma omp critical
					// {
						// cerr << "edgeArrayIndex[" << myThread << "]: " << edgeArrayIndex[myThread] << endl;
						// cerr << "edgeArrayIndex[" << myThread << "] + 3: " << edgeArrayIndex[myThread] + 3 << endl;
						// cerr << "edgeArray[" << edgeArrayIndex[myThread] << "]: " << edgeArray[edgeArrayIndex[myThread]] << endl;
						// cerr << "edgeArray[" << edgeArrayIndex[myThread] + 3 << "]: " << edgeArray[edgeArrayIndex[myThread] + 3] << endl;
						// cerr << "currEntity: " << currEntity << endl;
						// cerr << "entityRange[" << myThread << "]: " << entityRange[myThread] << endl;
						// cerr << "entityRange[" << myThread + 1 << "]: " << entityRange[myThread + 1] << endl;
						// cerr << "(edgeArray[edgeArrayIndex[myThread] + 3] >= entityRange[myThread]): " << (edgeArray[edgeArrayIndex[myThread] + 3] >= entityRange[myThread]) << endl;
						// cerr << "(edgeArray[edgeArrayIndex[myThread] + 3] < entityRange[myThread + 1]): " << (edgeArray[edgeArrayIndex[myThread] + 3] < entityRange[myThread + 1]) << endl;
						// cerr << ((edgeArray[edgeArrayIndex[myThread] + 3] >= entityRange[myThread]) && (edgeArray[edgeArrayIndex[myThread] + 3] < entityRange[myThread + 1])) << endl;
						// cerr << "eAVLArray[" << edgeArray[edgeArrayIndex[myThread] + 3] << "].insertV(" << edgeArray[edgeArrayIndex[myThread]] << ")" << endl;
						// cerr << endl;
					// }
				}
				// Move myThread to look at the next edge
				edgeArrayIndex[myThread] += 4;
			}
		}
	}
	#else
	{
		clock_gettime(CLOCK_MONOTONIC, &tTree1);
		UINT64 currEdgeIndex = 0;
		while (currEdgeIndex + 3 < readIndex) {
			#ifdef COUNT
			{
				if (count % 1000000 == 0) {
					// treeClock2 = clock();
					clock_gettime(CLOCK_MONOTONIC, &tTree2);
					// cerr << "Iteration " << iteration << " for thread 0." << endl;
					fprintf(stdout, "The time taken to insert %llu edges into the AVLs: %6.3fms.\n", count, msDiffTime(tTree1, tTree2));
					tTree1 = tTree2;
				}
				count++;
			}
			#endif
			
			// edgeArray[currEdgeIndex + 3] gets the current edge's entity ID
			// Set the current edge's AVL tree to have the paper that points to the entity
			eAVLArray[edgeArray[currEdgeIndex + 3]].insertV(edgeArray[currEdgeIndex]);
			currEdgeIndex += 4;
		}
	}
	#endif
	
	clock_gettime(CLOCK_MONOTONIC, &t4);
	fprintf(stdout, "End of Inserting Edges\nThe time taken to insert edges: %6.3fms.\n", msDiffTime(t3, t4));
	cout << "Begin Transferring Edges from the Tree to the Graph" << endl;
	
	int useless;	// Allows for a buffer, during which one can access how much memory is being used by the process
	#ifdef MEM
	{
		cout << "Check for memory of AVL trees full. (Enter 0) "; 
		cin >> useless;
	}
	#endif
	
	delete [] edgeArray;
	
	UINT32 currPapersNumOfEntities;
	// Set the size of each paper's neighborList to the number of entities the paper references
	for (pID.id = numOfEntities; pID.id < vNum; pID.id++) {
		currPapersNumOfEntities = paperCountArray[(pID.id - numOfEntities)];
		if (currPapersNumOfEntities > 0) {
			// myGraph->createV(pID, currPapersNumOfEntities, 0);
			myGraph->createVForward(pID, currPapersNumOfEntities);
		}
	}
	delete [] paperCountArray;
	
	// Go through each entity's AVL tree and add the edges in an in-order traversal
	for (eID.id = 0; eID.id < actualEntities; eID.id++) {
		if (eAVLArray[eID.id].count() > 0) {
			// cerr << eID.id << endl;
			// eAVLArray[eID.id].printIt();
			// cerr << eAVLArray[eID.id].count();
			// cerr << endl;
			
			eAVLArray[eID.id].addGEdges(myGraph, eID);
		}
	}
	
	// cerr << "G's unused space: " << myGraph->getUnusedSpace() << endl;
	
	cout << "End of Edge Transfer" << endl;
	clock_gettime(CLOCK_MONOTONIC, &t5);
	fprintf(stdout, "The time taken to transfer the edges: %6.3fms.\n", msDiffTime(t4, t5));
	fprintf(stdout, "The time taken to populate G from the start of read: %6.3fms.\n", msDiffTime(t2, t5));
	cout << "Total degree of G: " << myGraph->getTotalDegree() << endl << endl;
	
	#ifdef MEM
	{
		cout << "Check for memory of G. (Enter 0) "; 
		cin >> useless;
	}
	#endif
	
	// myGraph->printDebugGraph();
	// myGraph->printGraph();
	
// ------------------------------ Square The Graph --------------------------------------------	
	
	cout << "Start of Calculating G^2" << endl;
	Edge *neighborList;			// Stores a paper's edges; each edge will have an entity ID unique to the list
	UINT32 currOutDegree;		// Stores the number of entities the paper points to
	
	struct timespec tG2Tree1, tG2Tree2;
	UINT64 entityID1, entityID2;	// Pair of entities that a paper points to; add an edge between them for G2
	
	// Parallel G2 Insertion to AVLs
	#ifdef PAR
	{
	#pragma omp parallel private(myThread, pID, tG2Tree1, tG2Tree2, neighborList, currOutDegree, count, entityID1, entityID2)
	{
		clock_gettime(CLOCK_MONOTONIC, &tG2Tree1);
		myThread = omp_get_thread_num();
		count = 0;
		
		for (pID.id = numOfEntities; pID.id < vNum; pID.id++) {	// for each paper
			// Get the edges that paper points has and how many there are
			neighborList = myGraph->getVneighbors(pID);
			currOutDegree = myGraph->getVoutDegree(pID);
			if (currOutDegree > 0) {
				// #ifdef DEBUG
				// {
					// #pragma omp critical
					// {
						// cerr << "Thread: " << myThread << " - Paper: " << z.id << " - Current Out Degree: " << currOutDegree << endl;
						// myGraph->printVneighbors(z);
					// }
				// }
				// #endif
				
				// Have each pair of entities pointed to by the same paper also point to each other.
				// For each edge of the paper (except the last)
				for (UINT32 i = 0; i < currOutDegree - 1; i++) {
					entityID1 = neighborList[i].vid.id;
					// For each subsequent edge
					for (UINT32 j = i + 1; j < currOutDegree; j++) {
						
						#ifdef COUNT
						{
							if ((count % 100000000 == 0) && (myThread == 0)) {
								clock_gettime(CLOCK_MONOTONIC, &tG2Tree2);
								fprintf(stderr, "The time taken to insert %llu edges to AVLS on thread 0: %6.3fms.\n", count, msDiffTime(tG2Tree1, tG2Tree2));
								tG2Tree1 = tG2Tree2;
							}
							// count++;
						}
						#endif
						
						// eAVLArray[neighborList[i].vid.id].insertV(neighborList[j].vid.id);
						// eAVLArray[neighborList[j].vid.id].insertV(neighborList[i].vid.id);
						entityID2 = neighborList[j].vid.id;
						
						// Check if Entity 1 is modifiable by myThread
						if ((entityID1 >= entityRange[myThread]) && (entityID1 < entityRange[myThread + 1])) {
							eAVLArray[entityID1].insertV(entityID2);
							#ifdef COUNT
							{
								count++;
							}
							#endif
						}
						// Check if Entity 2 is modifiable by myThread
						if ((entityID2 >= entityRange[myThread]) && (entityID2 < entityRange[myThread + 1])) {
							eAVLArray[entityID2].insertV(entityID1);
							#ifdef COUNT
							{
								count++;
							}
							#endif
						}
					}
				}

				// #ifdef DEBUG
				// {
					// fprintf(stderr, "New Degree: %u\n\t", currOutDegree);
					// myGraph->printVneighbors(z);
				// }
				// #endif
			}
		}
	}
	}
	#else
	{
		clock_gettime(CLOCK_MONOTONIC, &tG2Tree1);
		count = 0;
		
		for (pID.id = numOfEntities; pID.id < vNum; pID.id++) {	// for each paper
			// Get the edges that paper points has and how many there are
			neighborList = myGraph->getVneighbors(pID);
			currOutDegree = myGraph->getVoutDegree(pID);
			if (currOutDegree > 0) {
				// #ifdef DEBUG
				// {
					// #pragma omp critical
					// {
						// cerr << "Paper: " << z.id << " - Current Out Degree: " << currOutDegree << endl;
						// myGraph->printVneighbors(z);
					// }
				// }
				// #endif
				
				// Have each pair of entities pointed to by the same paper also point to each other.
				// For each edge of the paper (except the last)
				for (UINT32 i = 0; i < currOutDegree - 1; i++) {
					entityID1 = neighborList[i].vid.id;
					// For each subsequent edge
					for (UINT32 j = i + 1; j < currOutDegree; j++) {
						
						#ifdef COUNT
						{
							if ((count % 100000000 == 0) && (myThread == 0)) {
								clock_gettime(CLOCK_MONOTONIC, &tG2Tree2);
								fprintf(stderr, "The time taken to insert %llu edges to AVLS on thread 0: %6.3fms.\n", count, msDiffTime(tG2Tree1, tG2Tree2));
								tG2Tree1 = tG2Tree2;
							}
						}
						#endif
						
						entityID2 = neighborList[j].vid.id;
						eAVLArray[entityID1].insertV(entityID2);
						eAVLArray[entityID2].insertV(entityID1);
						#ifdef COUNT
						{
							count += 2;
						}
						#endif
					}
				}

				// #ifdef DEBUG
				// {
					// fprintf(stderr, "New Degree: %u\n\t", currOutDegree);
					// myGraph->printVneighbors(z);
				// }
				// #endif
			}
		}
	}
	#endif
	
	cout << "Done Inserting G2 into the Trees." << endl;
	clock_gettime(CLOCK_MONOTONIC, &t6);
	fprintf(stdout, "The time taken to insert G2's edges into the trees: %6.3fms.\n", msDiffTime(t5, t6));
	cout << "Begin Transferring to G2." << endl;
	
	#ifdef MEM
	{
		cout << "Check for memory of G + AVL trees of G2. (Enter 0) "; 
		cin >> useless;
	}
	#endif
	
	// Serial G2 Transfer from trees to graph
	for (eID.id = 0; eID.id < actualEntities; eID.id++) {
		if (eAVLArray[eID.id].count() > 0) {
			// cerr << "\nEntity: " << eID.id << endl;
			// eAVLArray[eID.id].printIt();
			eAVLArray[eID.id].addG2Edges(myGraph, eID);
		}
	}
	
	// Parallel Transfer of G2's Trees to Graph
	// Currently not working, perhaps due to totalDegree increment being a race condition
	/*
	#pragma omp parallel private(eID, myThread)
	{
		myThread = omp_get_thread_num();
		for (eID.id = entityRange[myThread]; eID.id < entityRange[myThread + 1]; eID.id++) {
			// totalG2Edges += eAVLArray[eID.id].count();
			
			if (eAVLArray[eID.id].count() > 0) {
				// cerr << "\nEntity: " << eID.id << endl;
				// eAVLArray[eID.id].printIt();
				eAVLArray[eID.id].addG2Edges(myGraph, eID);
			}
		}
	}
	*/

	//cerr << "G2: " << endl;
	// myGraph->printDebugGraph();
	// myGraph->printGraphWeights();
	cout << "End of Populating G2" << endl;
	// cerr << "The total number of edges to expect on the graph is " << (totalG2Edges + myGraph->getTotalDegree()) << "." << endl;
    clock_gettime(CLOCK_MONOTONIC, &t7);
	fprintf(stdout, "The time taken to transfer G2's edges: %6.3fms.\n", msDiffTime(t6, t7));
	fprintf(stdout, "The time taken to populate G2: %6.3fms.\n", msDiffTime(t5, t7));
	cout << "Total Degree of G2: " << myGraph->getTotalDegree() << endl << endl;
	
	delete [] eAVLArray;
	
	#ifdef MEM
	{
		cout << "Check for memory of G2. (Enter 0) "; 
		cin >> useless;
	}
	#endif

	// Get all entities and their types
	/*
	UINT64 currTypeID;
	for (z.id = 0; z.id < maxEntity; z.id++) {
		currTypeID = myGraph->getTypeID(z);
		if (currTypeID != 0) {
			cout << "Type of entity " << z.id << ": " << currTypeID << endl;
		}
	}
	*/
	
// ---------------------- Incrementally Update the Graph -------------------------------------	

	if (hasInc) {
		readIndex = 0;
		UINT32 *incEdgeArray = new UINT32[numOfIncEdges * 4];
		
		if (!incBin) {
			// ASCII Read
			fstream myIncFile(incFile, ios_base::in);
			while (myIncFile >> a >> b >> c >> d) {
				// cout << "a: " << a << " - b: " << b << " - c: " << c << " - d: " << d << endl;
				
				// Put the paper in the edge array
				incEdgeArray[readIndex] = a + numOfEntities;
				// paperCountArray[a]++;	// Unnecessary for inc
				// Put the year and type into the array
				incEdgeArray[readIndex + 1] = b;
				incEdgeArray[readIndex + 2] = c;
				// Put the entity in the array
				incEdgeArray[readIndex + 3] = d;
				
				// If the edge's type ID hasn't been set yet, set it
				eID.id = d;
				if (eID.id > actualEntities) {
					fprintf(stderr, "Number of actual entities declared: %llu - Entity of added edge: %llu", actualEntities, eID.id);
					abort();
				}
				setTypeID(typeIDArray, eID, c);
				
				readIndex += 4;
			}
			myIncFile.close();
		} else {
			// Binary Read
			// Scan in the incremental binary file directly into incEdgeArray
			fread(incEdgeArray, sizeof(unsigned int), (numOfIncEdges * 4), myIncBinFileIn);
			fclose(myIncBinFileIn);
			
			UINT64 incEdgeTimes4 = 0;
			// For every edge
			for (UINT64 incEdge = 0; incEdge < numOfIncEdges; incEdge++) {
				incEdgeTimes4 = incEdge * 4;
				// fprintf(stderr, "Edge: %d", readEdge);
				// paperCountArray[incEdgeArray[incEdgeTimes4]]++;	// Not necessary

				// Add the namespace to each paper ID in the edge array
				incEdgeArray[incEdgeTimes4] += numOfEntities;
				
				// If the edge's type ID hasn't been set yet, set it
				eID.id = incEdgeArray[incEdgeTimes4 + 3];
				if (eID.id > actualEntities) {
					fprintf(stderr, "Number of actual entities declared: %llu - Entity of added edge: %llu", actualEntities, eID.id);
					abort();
				}
				setTypeID(typeIDArray, eID, incEdgeArray[incEdgeTimes4 + 2]);
			}
			readIndex = numOfIncEdges * 4;
		}
		
		// End of Incremental Read
		clock_gettime(CLOCK_MONOTONIC, &t8);
		fprintf(stdout, "The time taken to read the incremental file: %6.3fms.\n", msDiffTime(t7, t8));
		
		// Sanity check for the incremental edge array
		// for (int i = 0; i < readIndex / 4; i++) {
			// cout << incEdgeArray[i * 4] << "\t" << incEdgeArray[i * 4 + 1] << "\t" << incEdgeArray[i * 4 + 2] << "\t" << incEdgeArray[i * 4+3] << endl;
		// }
		
		// Incrementally add every new edge to the graph sequentially;
		// For every edge on the incremental edge array
		for (UINT64 currIncEdge = 0; currIncEdge < readIndex / 4; currIncEdge++) {
			pID.id = incEdgeArray[currIncEdge * 4];
			eID.id = incEdgeArray[(currIncEdge * 4) + 3];
			// cerr << "pID.id: " << pID.id << " - eID.id: " << eID.id << endl;
			
			// If an edge from the paper to the entity did not already exist;
			// Add an edge from the paper to the entity
			if ((pID.id != 0) && (myGraph->addEdge(pID, eID))) {
				// Find each edge that the paper points to 
				Edge *oldEntities = myGraph->getVneighbors(pID);
				UINT32 pOutDegree = myGraph->getVoutDegree(pID);
				
				// Increment the cooc between the new entity and all the paper's existing entities
				for (UINT32 currEnt = 0; currEnt < pOutDegree; currEnt++) {
					// Since the new entity has been added to the paper's neighbor list, avoid adding an edge to itself
					if (oldEntities[currEnt].vid.id != eID.id) {
						myGraph->addG2Edge(eID, oldEntities[currEnt].vid, 1, true);
					}
				}
			}
		}
		
		delete [] incEdgeArray;
		
		// End of Incremental Read
		clock_gettime(CLOCK_MONOTONIC, &t9);
		fprintf(stdout, "The time taken to add the incremental edges: %6.3fms.\n", msDiffTime(t8, t9));
		fprintf(stdout, "The total time taken for incremental data: %6.3fms.\n\n", msDiffTime(t7, t9));
		
		// myGraph->printGraphWeights();
	}
	
	
// ------------------------------ Find Coocs --------------------------------------------	
	
	// Find biggest cooc
	/*
	UINT64 maxCoocEntityID1, maxCoocEntityID2;
	UINT32 maxCooc = myGraph->getMaxCooc(numOfEntities, maxCoocEntity1, maxCoocEntity1);
	cout << "The max cooc is between entity " << maxCoocEntityID1 << " and entity " << maxCoocEntityID2 << " with a cooc weight of: " << maxCooc << endl;
	*/
	
	// Old way of getting max cooc
	/*
	UINT16 maxCooc = 0;
	Edge *eNList;
	UINT16 outDegree;
	UINT64 maxA;
	UINT64 maxB;
	for (z.id = 0; z.id < maxEntity; z.id++) {
		eNList = myGraph->getNList(z);
		outDegree = myGraph->getVoutDegree(z);
		for (UINT16 c = 0; c < outDegree; c++) {
			if (eNList[c].weight > maxCooc) {
				maxCooc = eNList[c].weight;
				maxA = z.id;
				maxB = eNList[c].vid.id;
			}
		}
	}
	cout << "Biggest cooc: " << maxCooc << " edges between entity " << maxA << " and entity " << maxB << endl;
	*/
	
	clock_gettime(CLOCK_MONOTONIC, &t10);
	// fprintf(stdout, "The total time taken to find the max cooc: %6.3fms.\n\n", msDiffTime(t9, t10));
	if (doExpand || doConnect || doBridge || doBiblio) {
		cout << "Begin Primitives." << endl;
	}
	
// ------------------------------ Perform Expand Primitive --------------------------------------------	
	
	if (doExpand) {
		cout << "Begin Expands." << endl;
		struct timespec tExpandStart, tExpandStop;
		clock_gettime(CLOCK_MONOTONIC, &tExpandStart);
		
		// Stores an array for the top coocs for each entity's expansion
		Edge **expandArray = new Edge*[expandIterations];
		// Stores the size of each entity's expansion (max can be expandCapacity)
		UINT16 *expandSize = new UINT16[expandIterations];
		
		#ifdef PARPRIM
		{
		// For every requested expansion
		#pragma omp parallel for private(x, typeID)
		for (int a = 0; a < expandIterations; a++) {
			// Get the ath arguments
			x.id = expandEID[a];
			typeID = expandTypeID[a];	
			UINT16 expandCapacity = 10;
			expandArray[a] = new Edge[expandCapacity];
			// The expansion array for the ath entity given is stored in the ath expandArray, and the 
			// number of elements in it is stored in the ath expandSize element
			myGraph->expand(x, typeID, expandArray[a], expandCapacity, expandSize[a], typeIDArray);
		}
		// After the expansions, go through and print out the results of each expansion serially
		for (int i = 0; i < expandIterations; i++) {
			cout << "Entity: " << expandEID[i] << " - Type ID: " << expandTypeID[i] << endl;
			for (int j = 0; j < expandSize[i]; j++) {
				cout << "Position " << j << ": Entity ID: " << expandArray[i][j].vid.id << " Expansion Weight: " << expandArray[i][j].weight << endl;
			}
			cout << endl;
		}
		// Print the total time taken since the expansions were done in parallel
		clock_gettime(CLOCK_MONOTONIC, &tExpandStop);
		fprintf(stdout, "The time taken to find the expands: %6.3fms.\n\n", msDiffTime(tExpandStart, tExpandStop));
		}
		#else
		{
		// For every requested expansion
		for (int i = 0; i < expandIterations; i++) {
			// Get the ith arguments 
			x.id = expandEID[i];
			typeID = expandTypeID[i];
			UINT16 expandCapacity = 10;
			expandArray[i] = new Edge[expandCapacity];
			
			// The expansion array for the ath entity given is stored in the ath expandArray, and the 
			// number of elements in it is stored in the ath expandSize element
			myGraph->expand(x, typeID, expandArray[i], expandCapacity, expandSize[i], typeIDArray);
			
			// Print out the results and time for each individual expansion
			cout << "Entity: " << expandEID[i] << " - Type ID: " << expandTypeID[i] << endl;
			for (int j = 0; j < expandSize[i]; j++) {
				cout << "Position " << j << ": Entity ID: " << expandArray[i][j].vid.id << " Expansion Weight: " << expandArray[i][j].weight << endl;
			}
			clock_gettime(CLOCK_MONOTONIC, &tExpandStop);
			fprintf(stdout, "Expand Time: %6.3fms.\n\n", msDiffTime(tExpandStart, tExpandStop));
			tExpandStart = tExpandStop;
		}	
		}
		#endif
	}
	
// ------------------------------ Perform Connect Primitive --------------------------------------------	

	if (doConnect) {
		cout << "Begin Connects." << endl;
		struct timespec tConnectStart, tConnectStop;
		clock_gettime(CLOCK_MONOTONIC, &tConnectStart);
		
		// Stores an array for the entities given that successfully connect with the given comparison entity
		Edge **connectArray = new Edge*[connectIterations];
		// Stores the size of each entity's number of successful connections (max can be connectCapacity)
		UINT16 *connectSize = new UINT16[connectIterations];
		
		#ifdef PARPRIM
		{
		// For each requested connection
		#pragma omp parallel for private(x, typeID)
		for (int a = 0; a < connectIterations; a++) {
			// Get the ath arguments
			x.id = connectEID[a];
			typeID = connectTypeID[a];
			// Number of entities the user put in the array of entities to check
			UINT16 connectCapacity = connectEIDArraySize[a];
			connectArray[a] = new Edge[connectCapacity];
			// The connection array for the ath entity given is stored in the ath connectArray, the array of entities the user 
			// requested to check are in the ath connectEIDArray element and the number of elements in it is stored in the ath 
			// connectSize element
			myGraph->connect(x, typeID, connectEIDArray[a], connectArray[a], connectCapacity, connectSize[a], typeIDArray);
		}
		// For each requested connection, print the results serially
		for (int i = 0; i < connectIterations; i++) {
			cout << "Entity: " << connectEID[i] << " - Type ID: " << connectTypeID[i] << endl;
			for (int j = 0; j < connectSize[i]; j++) {
				cout << "Position " << j << ": Entity ID: " << connectArray[i][j].vid.id << " Connection Weight: " << connectArray[i][j].weight << endl;
			}
			cout << endl;
		}
		// Print the total time taken for the connects since they were done in parallel
		clock_gettime(CLOCK_MONOTONIC, &tConnectStop);
		fprintf(stdout, "The time taken to find the connects: %6.3fms.\n\n", msDiffTime(tConnectStart, tConnectStop));
		}
		#else
		{
		// For every requested connection
		for (int i = 0; i < connectIterations; i++) {
			// Get the ith connect arguments
			x.id = connectEID[i];
			typeID = connectTypeID[i];
			// Number of entities the user wants to check for a connection with the ith comparison entity
			UINT16 connectCapacity = connectEIDArraySize[i];
			connectArray[i] = new Edge[connectCapacity];
			// The connection array for the ath entity given is stored in the ith connectArray, the array of entities the user 
			// requested to check are in the ith connectEIDArray element and the number of elements in it is stored in the ith 
			// connectSize element
			myGraph->connect(x, typeID, connectEIDArray[i], connectArray[i], connectCapacity, connectSize[i], typeIDArray);
			
			// Print out each connection's output and time after it happens
			cout << "Entity: " << connectEID[i] << " - Type ID: " << connectTypeID[i] << endl;
			for (int j = 0; j < connectSize[i]; j++) {
				cout << "Position " << j << ": Entity ID: " << connectArray[i][j].vid.id << " Connection Weight: " << connectArray[i][j].weight << endl;
			}
			clock_gettime(CLOCK_MONOTONIC, &tConnectStop);
			fprintf(stdout, "Connect time: %6.3fms.\n\n", msDiffTime(tConnectStart, tConnectStop));
			tConnectStart = tConnectStop;
		}
		}
		#endif
	}
	
// ------------------------------ Perform Bridge Primitive --------------------------------------------	
	
	// General Bridge Primitive
	if (doBridge) {
		cout << "Begin Bridges." << endl;
		struct timespec tBridgeStart, tBridgeStop;
		clock_gettime(CLOCK_MONOTONIC, &tBridgeStart);
		
		// Stores an array for the top shared coocs for each pair of given entities
		Edge **coocArray = new Edge*[bridgeIterations];
		// Stores the size of each entity's bridge (max can be bridgeCapacity)
		UINT16 *coocSize = new UINT16[bridgeIterations];
		
		#ifdef PARPRIM
		{
		// For every request bridge
		#pragma omp parallel for private(x, z, typeID)
		for (int a = 0; a < bridgeIterations; a++) {
			// Get the ath bridge arguments
			x.id = bridgeEID1[a];
			z.id = bridgeEID2[a];
			typeID = bridgeTypeID[a];
			UINT16 coocCapacity = 10;
			coocArray[a] = new Edge[coocCapacity];
			// The bridge array for the ath pair of entities given is stored in the ath coocArray array, and the 
			// number of elements in it is stored in the ath coocSize element
			myGraph->bridge(x, z, typeID, coocArray[a], coocCapacity, coocSize[a], typeIDArray);
		}
		// For finding U's cooc weight and V's cooc weight with the coocArray entities
		/*
		UINT32 uIndex, vIndex;
		Edge *uNeighbors, *vNeighbors;
		*/
		// For every requested bridge, print out the results serially
		for (int i = 0; i < bridgeIterations; i++) {
			cout << "Entity One: " << bridgeEID1[i] << " - Entity Two: " << bridgeEID2[i] << " - Type ID: " << bridgeTypeID[i] << endl;
			for (int j = 0; j < coocSize[i]; j++) {
				cout << "Position " << j << ": Entity ID: " << coocArray[i][j].vid.id << " Cooc Weight: " << coocArray[i][j].weight << endl;
				// For finding U's cooc weight and V's cooc weight with the coocArray entities
				/*
				x.id = bridgeEID1[i];
				z.id = bridgeEID2[i];
				myGraph->isEdge(x, coocArray[i][j].vid, uIndex);
				myGraph->isEdge(z, coocArray[i][j].vid, vIndex);
				// cout << "myGraph->isEdge(x, coocArray[i][j].vid, uIndex): " << myGraph->isEdge(x, coocArray[i][j].vid, uIndex) << endl;
				// cout << "myGraph->isEdge(z, coocArray[i][j].vid, vIndex): " << myGraph->isEdge(z, coocArray[i][j].vid, vIndex) << endl;
				// cout << "uIndex: " << uIndex << " - vIndex: " << vIndex << endl;
				uNeighbors = myGraph->getVneighbors(x);
				vNeighbors = myGraph->getVneighbors(z);
				// cout << "uNeighbors[uIndex].vid.id: " << uNeighbors[uIndex].vid.id << "vNeighbors[vIndex].vid.id: " << vNeighbors[vIndex].vid.id << endl;
				cout << "Cooc Weight with U: " << uNeighbors[uIndex].weight << " - Cooc Weight with V: " << vNeighbors[vIndex].weight << endl;
				*/
			}
			cout << endl;
		}
		// Print out the total time taken, since they were all done in parallel
		clock_gettime(CLOCK_MONOTONIC, &tBridgeStop);
		fprintf(stdout, "The time taken to find the bridges: %6.3fms.\n\n", msDiffTime(tBridgeStart, tBridgeStop));
		}
		#else
		{
		// For every requested bridge
		for (int i = 0; i < bridgeIterations; i++) {
			// Get the ith arguments
			x.id = bridgeEID1[i];
			z.id = bridgeEID2[i];
			typeID = bridgeTypeID[i];
			UINT16 coocCapacity = 10;
			coocArray[i] = new Edge[coocCapacity];
			
			// The bridge array for the ath pair of entities given is stored in the ath coocArray array, and the 
			// number of elements in it is stored in the ath coocSize element
			myGraph->bridge(x, z, typeID, coocArray[i], coocCapacity, coocSize[i], typeIDArray);
			
			// For finding U's cooc weight and V's cooc weight with the coocArray entities
			/*
			UINT32 uIndex, vIndex;
			Edge *uNeighbors, *vNeighbors;
			*/
			// Print the individual bridge's results and time after it happens
			cout << "Entity One: " << bridgeEID1[i] << " - Entity Two: " << bridgeEID2[i] << " - Type ID: " << bridgeTypeID[i] << endl;
			for (int j = 0; j < coocSize[i]; j++) {
				cout << "Position " << j << ": Entity ID: " << coocArray[i][j].vid.id << " Cooc Weight: " << coocArray[i][j].weight << endl;
				// For finding U's cooc weight and V's cooc weight with the coocArray entities
				/*
				x.id = bridgeEID1[i];
				z.id = bridgeEID2[i];
				myGraph->isEdge(x, coocArray[i][j].vid, uIndex);
				myGraph->isEdge(z, coocArray[i][j].vid, vIndex);
				cout << "myGraph->isEdge(x, coocArray[i][j].vid, uIndex): " << myGraph->isEdge(x, coocArray[i][j].vid, uIndex) << endl;
				cout << "myGraph->isEdge(z, coocArray[i][j].vid, vIndex): " << myGraph->isEdge(z, coocArray[i][j].vid, vIndex) << endl;
				cout << "uIndex: " << uIndex << " - vIndex: " << vIndex << endl;
				uNeighbors = myGraph->getVneighbors(x);
				vNeighbors = myGraph->getVneighbors(z);
				cout << "uNeighbors[uIndex].vid.id: " << uNeighbors[uIndex].vid.id << "vNeighbors[vIndex].vid.id: " << vNeighbors[vIndex].vid.id << endl;
				cout << "Cooc Weight with U: " << uNeighbors[uIndex].weight << " - Cooc Weight with V: " << vNeighbors[vIndex].weight << endl;
				*/
			}
			clock_gettime(CLOCK_MONOTONIC, &tBridgeStop);
			fprintf(stdout, "Bridge time: %6.3fms.\n\n", msDiffTime(tBridgeStart, tBridgeStop));
			tBridgeStart = tBridgeStop;
		}
		}
		#endif
	}
	
	// Bridge for baby test
	/*
	x.id = 2;
	z.id = 4;
	UINT16 typeID = 11;
	UINT16 coocCapacity = 10;
	UINT16 coocSize;
	Edge coocArray[coocCapacity];
	// myGraph->bridge(coocArray, coocCapacity, x, z, typeID);
	myGraph->bridge(x, z, typeID, coocArray, coocCapacity, coocSize, typeIDArray);
	for (int i = 0; i < coocSize; i++) {
		cout << "Position " << i << ": Entity ID: " << coocArray[i].vid.id << " Cooc Weight: " << coocArray[i].weight << endl;
	}
	*/
		
	// Bridge for biggest cooc in 100k
	/*
	x.id = 522584;
	z.id = 522585;
	UINT16 typeID1 = getTypeID(typeIDArray, x);
	UINT16 typeID2 = getTypeID(typeIDArray, z);
	UINT16 coocCapacity1 = 10, coocCapacity2 = 10;
	UINT16 coocSize1, coocSize2;
	Edge coocArray1[coocCapacity1];
	// myGraph->bridge(coocArray1, arrayCapacity1, x, z, typeID1);
	myGraph->bridge(x, z, typeID1, coocArray1, coocCapacity1, coocSize1, typeIDArray);
	cout << "Entity One: " << x.id << " - Entity Two: " << z.id << " - Type ID: " << typeID1 << endl;
	for (int i = 0; i < coocSize1; i++) {
		cout << "Position " << i << ": Entity ID: " << coocArray1[i].vid.id << " Cooc Weight: " << coocArray1[i].weight << endl;
	}
	cout << "\nEntity One: " << x.id << " - Entity Two: " << z.id << " - Type ID: " << typeID2 << endl;
	Edge coocArray2[coocCapacity2];
	// myGraph->bridge(coocArray2, arrayCapacity2, x, z, typeID2);
	myGraph->bridge(x, z, typeID2, coocArray2, coocCapacity2, coocSize2, typeIDArray);
	for (int i = 0; i < coocSize2; i++) {
		cout << "Position " << i << ": Entity ID: " << coocArray2[i].vid.id << " Cooc Weight: " << coocArray2[i].weight << endl;
	}
	clock_gettime(CLOCK_MONOTONIC, &t8);
	fprintf(stderr, "The time taken to find the bridge of %llu and %llu of type %d: %6.3fms.\n", x.id, z.id, typeID1, msDiffTime(t7, t8));
	*/
	
	// Find biggest bridge in G2
	/*
	UINT16 arrayCapacity = 1;
	Edge coocArray[arrayCapacity];
	UINT16 maxBridge = 0;
	UINT64 maxE1, maxE2, maxSharedE;
	for (x.id = 0; x.id < maxEntity - 1; x.id++) {
		for (z.id = x.id + 1; z.id < maxEntity; z.id++) {
			myGraph->bridge(coocArray, arrayCapacity, x, z);
			if (coocArray[0].vid.id != -1) {
				if (coocArray[0].weight > maxBridge) {
					maxSharedE = coocArray[0].vid.id;
					maxmaxBridge = coocArray[0].weight;
					maxE1 = x.id;
					maxE2 = z.id;
				}
			}
		}
	}
	cout << "Max Shared Entity: " << maxSharedE << " - Max Cooc Weight: " << maxBridge << " - Max E1: " << maxE1 << " - Max E2: " << maxE2 << endl;
	*/
	
// ------------------------------ Perform Bibliography Primitive --------------------------------------------	
	
	if (doBiblio) {
		cout << "Begin Bibliographies." << endl;
		struct timespec tBiblioStart, tBiblioStop;
		clock_gettime(CLOCK_MONOTONIC, &tBiblioStart);
		
		// Stores an array for the top shared papers for each pair of entities's bibliography
		UINT64 **biblioArray = new UINT64*[biblioIterations];
		// Stores the size of each pair of entities' biblioigraphies (max can be its biblioCapacity argument)
		UINT32 *biblioSize = new UINT32[biblioIterations];
		
		#ifdef PARPRIM
		{
		// For every request bibliography
		#pragma omp parallel for private(x, z)
		for (int a = 0; a < biblioIterations; a++) {
			// Get the ath arguments
			x.id = biblioEID1[a];
			z.id = biblioEID2[a];
			// biblioCapacity[a] stores how many papers the user wants on the list
			// To make it fit every paper link, instead of asking for how many papers the user wants,
			// you can set it to be the max number of papers that reference entity 1 or 2 (max inDegree
			// of the two entities)
			biblioArray[a] = new UINT64[biblioCapacity[a]];
			
			// The bibliography array for the ath pair of entities given is stored in the ath biblioArray array, 
			// and the number of elements in it is stored in the ath biblioSize element
			myGraph->bibliography(x, z, biblioArray[a], biblioCapacity[a], biblioSize[a]);
		}
		
		// Print out the results of each requested bibliography serially
		for (int i = 0; i < biblioIterations; i++) {
			cout << "Entity 1: " << biblioEID1[i] << " Entity 2: " << biblioEID2[i] << endl;
			for (int j = 0; j < biblioSize[i]; j++) {
				cout << "Position " << j << ": Paper: " << biblioArray[i][j] << endl;
			}
			cout << endl;
		}
		// Print out the time for all of the bibliographies since they were all done in parallel
		clock_gettime(CLOCK_MONOTONIC, &tBiblioStop);
		fprintf(stdout, "The time taken to find the bibliographies: %6.3fms.\n\n", msDiffTime(tBiblioStart, tBiblioStop));
		}
		#else
		{
		// For each requested bibliography
		for (int i = 0; i < biblioIterations; i++) {
			// Get the ith arguments
			x.id = biblioEID1[i];
			z.id = biblioEID2[i];
			// biblioCapacity[a] stores how many papers the user wants on the list
			// To make it fit every paper link, instead of asking for how many papers the user wants,
			// you can set it to be the max number of papers that reference entity 1 or 2 (max inDegree
			// of the two entities)
			biblioArray[i] = new UINT64[biblioCapacity[i]];
			
			// The bibliography array for the ath pair of entities given is stored in the ath biblioArray array, 
			// and the number of elements in it is stored in the ath biblioSize element
			myGraph->bibliography(x, z, biblioArray[i], biblioCapacity[i], biblioSize[i]);
		
			// Print out the bibliography's results and time after it is found
			cout << "Entity 1: " << biblioEID1[i] << " Entity 2: " << biblioEID2[i] << endl;
			for (int j = 0; j < biblioSize[i]; j++) {
				cout << "Position " << j << ": Paper: " << biblioArray[i][j] << endl;
			}
			clock_gettime(CLOCK_MONOTONIC, &tBiblioStop);
			fprintf(stdout, "The time taken to find the bibliographies: %6.3fms.\n\n", msDiffTime(tBiblioStart, tBiblioStop));
			tBiblioStart = tBiblioStop;
		}
		}
		#endif
	}

	
	// End of Primitives
	clock_gettime(CLOCK_MONOTONIC, &t11);
	if (doExpand || doConnect || doBridge || doBiblio) {
		fprintf(stdout, "The time taken to perform the primitives: %6.3fms.\n", msDiffTime(t10, t11));
	}

// -------------------------------- Dump to File  --------------------------------------------
	
	// Incomplete: After G2 is found, dump the graph into a file that can be quickly scanned in so that G2 is reset
	if (dumpToFile) {
		myGraph->dumpGraph();
	}	
	
// -------------------------------- Make Hisotgram --------------------------------------------
	
	// By Michael VanDusen
	// Shows how many vertices have how many edges, as well as the maximum number of edges one vertex has
    /*
    int max = 0;
    int maxID = 0;
    int maxV = myGraph->getVertexCount();
    UINT32 *outDegreeHist = UINT32 int[11];
    for (int i = 0; i < 11 ;i++){
        outDegreeHist[i] = 0;
    }
    // 0   1-20
    // 1   21-50
    // 2   51-100
    // 3   101-200
    // 4   201-400
    // 5   401-800
    // 6   800-1000
    // 7   1001-1500
    // 8   1501-2000
    // 9   2001-3185
    for(z.id = 0; z.id < maxV; z.id++) {
        if (myGraph->getVoutDegree(z) == 0) {
			outDegreeHist[10]++;
		}
		else if(myGraph->getVoutDegree(z) < 21 && myGraph->getVoutDegree(z) > 0){
            outDegreeHist[0]++;
        }
        else if(myGraph->getVoutDegree(z) < 51){
            outDegreeHist[1]++;
        }
        else if(myGraph->getVoutDegree(z) < 101){
            outDegreeHist[2]++;
        }
        else if(myGraph->getVoutDegree(z) < 201){
            outDegreeHist[3]++;
        }
        else if(myGraph->getVoutDegree(z) < 401){
            outDegreeHist[4]++;
        }
        else if(myGraph->getVoutDegree(z) < 801){
            outDegreeHist[5]++;
        }
        else if(myGraph->getVoutDegree(z) < 1001){
            outDegreeHist[6]++;
        }
        else if(myGraph->getVoutDegree(z) < 1501){
            outDegreeHist[7]++;
        }
        else if(myGraph->getVoutDegree(z) < 2001){
            outDegreeHist[8]++;
        }
        else{
            outDegreeHist[9]++;
        }
    }
    clock_t clock4 = clock(); // after Histogram
    for (z.id = 0; z.id < maxV; z.id++){		// find the vertex with the most eges outdegree
        if(max < myGraph->getVoutDegree(z)){
            max = myGraph->getVoutDegree(z);
            maxID = z.id;
        }    
    }
	*/
    
// --------------------- Print Histogram -----------------------------------------------

	/*
    cout << "Out Degree Histogram" << endl << endl;
	cout << "0             " << outDegreeHist[10] << endl;
    cout << "1 - 20        " << outDegreeHist[0] << endl;
    cout << "21 - 50       " << outDegreeHist[1] << endl;
    cout << "51 - 100      " << outDegreeHist[2] << endl;
    cout << "101 - 200     " << outDegreeHist[3] << endl;
    cout << "201 - 400     " << outDegreeHist[4] << endl;
    cout << "401 - 800     " << outDegreeHist[5] << endl;
    cout << "801 - 1000    " << outDegreeHist[6] << endl;
    cout << "1001 - 1500   " << outDegreeHist[7] << endl;
    cout << "1501 - 2000   " << outDegreeHist[8] << endl;
    cout << "2001+         " << outDegreeHist[9] << endl;
	
    // cout << "Vertex " << maxID << " has the largest out degrees which is: " << max << endl;
	// print out the vertex with most edges (optional fill free to comment out)
	// z.id = maxID;
    // myGraph->printVneighbors(z);	
	*/

// ----------------------------- Print Clocks ----------------------------------------------
	
	
	// ttime = (float)(clock1 - clock0) / (float)CLOCKS_PER_SEC;
    // cout << "The time taken to allocate the graph: " << ttime <<endl;
	// ttime = (float)(clock2 - clock1) / (float)CLOCKS_PER_SEC;
    // cerr << "The time taken to prepare for the read: " << ttime << endl;
	// ttime = (float)(clockEndRead - clock2) / (float)CLOCKS_PER_SEC;
	// cerr << "The time taken to read: " << ttime << endl;
    // ttime = (float)(clock3 - clockEndRead) / (float)CLOCKS_PER_SEC;
    // cerr << "The time taken to insert p to e edges in the AVLs: " << ttime << endl;
    // ttime = (float)(clock4 - clock3) / (float)CLOCKS_PER_SEC;
    // cerr << "The time taken to transfer edges from the tree to the graph: " << ttime << endl;
	// ttime = (float)(clock4 - clock2) / (float)CLOCKS_PER_SEC;
    // cerr << "The time taken to populate the graph from the start of reading in: " << ttime << endl;
	// ttime = (float)(clock5 - clock4) / (float)CLOCKS_PER_SEC;
    // cerr << "The time taken to populate G2: " << ttime <<endl;
    // ttime = (float)(clock4 - clock3) / (float)CLOCKS_PER_SEC;
    // cout << "The time taken to make Histogram: " << ttime <<endl;
    // ttime = (float)(clock5 - clock4) / (float)CLOCKS_PER_SEC;
    // cout << "The time taken to find max: " << ttime <<endl;
    // ttime = (float)(clock()-clock0) / (float)CLOCKS_PER_SEC;
	// cout << "The time taken for the total program: " << ttime << endl;
	
	struct timespec tEnd;
	clock_gettime(CLOCK_MONOTONIC, &tEnd);
	fprintf(stdout, "\nThe total time taken: %6.3fms.\n", msDiffTime(t0, tEnd));
	// cout << "The unused space is " << myGraph->getUnusedSpace() << endl;
	// myGraph->printGraphWeights();
	
	#ifdef MEM
	{
		cout << "Final check for memory. (Enter 0) "; 
		cin >> useless;
	}
	#endif
	cout << "End of Main" << endl;
	
    return 0;
}
