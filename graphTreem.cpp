// File: mvgraphm.cpp
// Author: Michael VanDusen
// Date: 2/22/22
// Pupose: To test the graph program made by the ORU Research team

#include <cstdio>
#include <ctype.h>			// Provides isdigit
#include <iostream>			// Provides cout
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include "graph32.h"
#include "longAVL.h"

using namespace std;

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

class bvGraph : public graph {
	public:
		bvGraph(UINT64 numOfVertices) : graph(numOfVertices) {}
		~bvGraph() {graph::~graph();}
		
		UINT32 getMaxCooc(UINT64 numOfEntities) {
			vID z;
			Edge *eNList;
			UINT32 maxCooc = 0;
			UINT32 outDegree;
			UINT64 maxA, maxB;
			for (z.id = 0; z.id < numOfEntities; z.id++) {
				eNList = getVneighbors(z);
				outDegree = getVoutDegree(z);
				for (UINT32 c = 0; c < outDegree; c++) {
					if (eNList[c].weight > maxCooc) {
						maxCooc = eNList[c].weight;
						maxA = z.id;
						maxB = eNList[c].vid.id;
					}
				}
			}
			return maxCooc;
		}
		friend class bvAVL;
};

class bvAVL : public longAVL {
	private:
		void addAllEdges(bvGraph *bvg, vID eID, node *p) {
			if (p) {
				addAllEdges(bvg, eID, p->left);
				vID pID;
				pID.id = p->val;
				bvg->addEdge2(pID, eID);
				addAllEdges(bvg, eID, p->right);
			}
		}
		
	public:
		bvAVL() : longAVL(){}
		~bvAVL() {longAVL::~longAVL();}
		
		void addAllEdges(bvGraph *bvg, vID eID) {
			bvg->createV(eID, 0, tCount);
			addAllEdges(bvg, eID, root);
		}
};

int main(int argc, char **argv) {
	
	cout << "Begin!" << endl;
    clock_t clock0 = clock();	//Start of the program
   
// --------- Graph and generater variable defaults ---------------------
   
	bool dumpToFile = false;
	bool hasFile = false;
	
	vID x;							// Temp vIDs
	vID z;							
	INT16 opt;						// For reading in user input
	INT16 allocation = 16;			// How much to grow neighbor lists by
	UINT64 vNum = 32000000;			// Num of total vertices, |Papers| + |Entities|
	UINT64 numOfEntities = 2000000; // Num of entities there can be
	UINT64 numOfPapers = vNum - numOfEntities;

	char *file = NULL;
	
	opterr = 0;

// -------------------- Getop Flags----------------------------------

	while ((opt = getopt(argc, argv, "f:")) != -1) {
		switch (opt) {
			case 'h':
				cout << "Command Flag Options:" << endl;
				cout << "-f [filename]: Read in filename and produce a graph based on the paperID-entityID pairs established in the file." << endl;
				abort();
			case 'f':
				hasFile = true;
				file = optarg;
				break;
			default:
				fprintf(stderr, "Unrecognized option.\n");
				abort();
		}
	}
	for (int i = optind; i < argc; i++) {
		printf ("Non-option argument %s\n", argv[i]); 
	}
   
   if(!hasFile) {		// If no input flag given, tell them they need to do something
		cout << "-f [filename] option required." << endl;	// For now, if there is no input, break
		abort();
	}
	
// --------------------- Init Graph ----------------------------------------------------
   
    // graph *myGraph = new graph(vNum);
	bvGraph *myGraph = new bvGraph(vNum);
	myGraph->setAllocation(allocation);
	
// ----------------------- Create Type Array ------------------------------------------------	

	// UINT64 typeArrSize = ((UINT64)(ceil(numOfEntities / 8)) + 1) * 8;
	// UINT8 *typeIDArray = new UINT8[typeArrSize];
	UINT8 *typeIDArray = new UINT8[numOfEntities];
	bvAVL *eAVLArray = new bvAVL[numOfEntities];
	
	// The index of each paper's count will be pID - numOfEntities
	UINT16 *paperCountArray = new UINT16[numOfPapers];
	
	UINT64 t = 0;
	// Set both arrays at the same time
	// #pragma omp parallel for 
	// for (UINT64 t = 0; t < typeArrSize; t++) {
	for (; t < numOfEntities; t++) {
		typeIDArray[t] = 0;
		paperCountArray[t] = 0;
	}
	// t will == numOfEntities, so just fill out the rest
	for (; t < numOfPapers; t++) {
		paperCountArray[t] = 0;
	}
	
	clock_t clock1 = clock();	//Start of read
    float ttime = (float)(clock1 - clock0) / (float)CLOCKS_PER_SEC;
    cout << "The time taken to allocate the graph: " << ttime << endl;
	
// ----------------------- Populate the Graph ------------------------------------------------	
	
	cout << "\nStart of Read" << endl;
	
	// G Maker
	vID pID;
	vID eID;
	UINT64 a, b, c, d;
	UINT8 type;
	UINT64 twoToTheEighth = (UINT64)(pow(2, 8));
	UINT64 count = 0, limit = 80;
	fstream myfile(file, ios_base::in);
	
	while (myfile >> a >> b >> c >> d) {
		 // if (count < limit) {
			pID.id = a + numOfEntities;
			
			eID.id = d;
			
			#ifdef COUNT
				if (count % 10000 == 0) {
					fprintf(stderr, "\nCount: %llu, Paper: %llu, Entity: %llu.\n", count, pID.id, eID.id);
				}
			#endif
			
			eAVLArray[eID.id].insertV(pID.id);
			// myGraph->addEdge(pID, eID);
			// a == pID.id - numOfEntities; allows 0-based indexing of the count array
			paperCountArray[a]++;
			setTypeID(typeIDArray, eID, c);
			
			// if (count == limit - 1) {
				// myGraph->printDebugGraph();
				// break;
				// fprintf(stderr, "Paper capacity: %d, Paper back capacity: %d.\n Entity capacity: %d, Entity back capacity: %d\nPaper Out Degree: %d, Paper In Degree: %d.\n Entity Out Degree: %d, Entity In Degree: %d\n", 
					// myGraph->getCapacity(pID), myGraph->getBackCapacity(pID), myGraph->getCapacity(eID), myGraph->getBackCapacity(eID),
					// myGraph->getVoutDegree(pID), myGraph->getVinDegree(pID), myGraph->getVoutDegree(eID), myGraph->getVinDegree(eID));
				// myGraph->printVneighbors(pID);
				// myGraph->printVbackNeighbors(eID);
			// }
			count++;
		// } else {
			// break;
		// }
	}
	
	// Set the size of each paper's neighborList
	for (pID.id = numOfEntities; pID.id < vNum; pID.id++) {
		myGraph->createV(pID, paperCountArray[(pID.id - numOfEntities)], 0);
	}
	
	for (eID.id = 0; eID.id < numOfEntities; eID.id++) {
		eAVLArray[eID.id].addAllEdges(myGraph, eID);
	}
	
	// x.id = 361056;
	// cerr << "\nData of " << x.id << ":\nOut Degree: " << myGraph->getVoutDegree(x) << "\tOut Capacity: " << myGraph->getVcapacity(x) << endl;
	// cerr << "In Degree: " << myGraph->getVinDegree(x) << "\tIn Capacity: " << myGraph->getVbackCapacity(x) << endl;
	cerr << "Total degree of G: " << myGraph->getTotalDegree() << endl;
	cerr << "End of Read" << endl;
	clock_t clock2 = clock();		//Graph has been populated
	
	
    ttime = (float)(clock2 - clock1) / (float)CLOCKS_PER_SEC;
    cout << "The time taken to populate the graph: " << ttime << endl;
	
// ------------------------------ Square The Graph --------------------------------------------	
	// Parallel G2 Maker
	/*
	int numOfThreads;
	#pragma omp parallel
	{
		numOfThreads = omp_get_num_threads();
	}
	*/
	
	// Calculate G^2
	cout << "\nStart of Calculating G^2" << endl;
	Edge *neighborList;
	UINT32 currOutDegree;
	
	for (z.id = numOfEntities; z.id < vNum; z.id++) {	// for each paper
		// Get the entities that paper points to
		neighborList = myGraph->getVneighbors(z);
		currOutDegree = myGraph->getVoutDegree(z);
		
		if (currOutDegree > 0) {
			
			// #ifdef DEBUG
			// {
				// fprintf(stderr, "Old Paper: %llu; Degree: %u\n\t", z.id, currOutDegree);
				// myGraph->printVneighbors(z);
			// }
			// #endif
			// Have each pair of entities pointed to by the same paper also point to each other
			for (UINT32 i = 0; i < currOutDegree - 1; i++) {
				for (UINT32 j = i + 1; j < currOutDegree; j++) {
					myGraph->addEdge(neighborList[i].vid, neighborList[j].vid, 1, true);
					myGraph->addEdge(neighborList[j].vid, neighborList[i].vid, 1, true);
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
	clock_t clock3 = clock();// before the histogram calculations
	cout << "End of Populating G2" << endl;
	
    ttime = (float)(clock3 - clock2) / (float)CLOCKS_PER_SEC;
    cout << "The time taken to populate G2: " << ttime <<endl;

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
	
	// cout << "Total degree of G2: " << myGraph->getTotalDegree() << endl;
	
// ------------------------------ Find Coocs --------------------------------------------	
	
	// Find biggest cooc
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
	
// ------------------------------ Use Bridge Primitive --------------------------------------------	
	
	// Bridge for baby test
	/*
	x.id = 2;
	z.id = 4;
	UINT64 typeID = 11;
	UINT16 arrayCapacity = 10;
	Edge coocArray[arrayCapacity];
	myGraph->bridge(coocArray, arrayCapacity, x, z, typeID);
	for (int i = 0; i < arrayCapacity; i++) {
		if (coocArray[i].vid.id == -1) {
			break;
		}
		cout << "Position " << i << ": Entity ID: " << coocArray[i].vid.id << " Cooc Weight: " << coocArray[i].weight << endl;
	}
	*/
	
	// Bridge for biggest cooc in 100k
	/*
	x.id = 522584;
	z.id = 522585;
	UINT64 typeID1 = myGraph->getTypeID(x);
	UINT64 typeID2 = myGraph->getTypeID(z);
	UINT16 arrayCapacity = 10;
	Edge coocArray[arrayCapacity];
	myGraph->bridge(coocArray, arrayCapacity, x, z, typeID1);
	cout << "Entity One: " << x.id << " - Entity Two: " << z.id << " - Type ID: " << typeID1 << endl;
	for (int i = 0; i < arrayCapacity; i++) {
		if (coocArray[i].vid.id == -1) {
			break;
		}
		cout << "Position " << i << ": Entity ID: " << coocArray[i].vid.id << " Cooc Weight: " << coocArray[i].weight << endl;
	}
	cout << "\nEntity One: " << x.id << " - Entity Two: " << z.id << " - Type ID: " << typeID2 << endl;
	Edge coocArray2[arrayCapacity];
	myGraph->bridge(coocArray2, arrayCapacity, x, z, typeID2);
	for (int i = 0; i < arrayCapacity; i++) {
		if (coocArray2[i].vid.id == -1) {
			break;
		}
		cout << "Position " << i << ": Entity ID: " << coocArray2[i].vid.id << " Cooc Weight: " << coocArray2[i].weight << endl;
	}
	*/
	
	// Find biggest pair of coocs
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
	
	
	if (dumpToFile) {
		myGraph->dumpGraph();
	}
	
// ------------------------------ Stats Finder --------------------------------------------	

	// z.id = 1 + numOfEntities;
	// myGraph->printVneighbors(z);
	// myGraph->printVBackNeighbors(z);
	// z.id = 959996;
	// myGraph->printVneighbors(z);
	// myGraph->printVBackNeighbors(z);

	// for your use
	// for (z.id = 0; z.id < vNum; z.id++) {
		// myGraph->printVneighborsWeights(z);
		// myGraph->printVBackNeighbors(z);
	// }
	/*
	for (z.id = 0; z.id <= maxEntity; z.id++) {
		myGraph->printVneighborsWeights(z);
		myGraph->printVBackNeighbors(z);
	}
	for (z.id = numOfEntities; z.id <= maxPaper; z.id++) {
		myGraph->printVneighbors(z);
		myGraph->printVBackNeighbors(z);
	}
	*/
	
	// vID x;
	// x.id = 438432;
	// z.id = 623492;
	// UINT8 weightSum = 0;
	// Edge *aList = myGraph->getBackNList(x);
	// UINT16 inA = myGraph->getVinDegree(x);
	// Edge *bList = myGraph->getBackNList(z);
	// UINT16 inB = myGraph->getVinDegree(z);
	// cout << "Shared papers between " << x.id << " and " << z.id << ": ";
	// for (int i = 0; i < inA; i++) {
		// for (int j = 0; j < inB; j++) {
			// if (aList[i].vid.id == bList[j].vid.id) {
				// cout << aList[i].vid.id << " ";
				// weightSum++;
			// }
		// }
	// }
	// cout << "\nTotal weight between " << x.id << " and " << z.id << ": " << weightSum;
	
	
	// for (z.id = 1 + numOfEntities; z.id < 4 + numOfEntities; z.id++) {
		// cout << "Paper " << (z.id - numOfEntities) << ": \n";
		// myGraph->printVneighbors(z);
		// myGraph->printVneighborsWeights(z);
		// myGraph->printVBackNeighbors(z);
		// Edge *edges = myGraph->getNList(z);
		// UINT16 outDeg = myGraph->getVoutDegree(z);
		// for (int a = 0; a < outDeg; a++) {
			// cout << "Entity " << edges[a].vid.id << ": \n";
			// myGraph->printVneighbors(edges[a].vid);
			// myGraph->printVneighborsWeights(edges[a].vid);
			// myGraph->printVBackNeighbors(edges[a].vid);
		// }
	// }
	// for (z.id = 1; z.id < 6; z.id++) {
		// myGraph->printVneighbors(z);
		// myGraph->printVBackNeighbors(z);
	// }
	
	
    // for (z.id = 0; z.id < print10;z.id++){
        // myGraph->printVneighbors(z);
    // }
	
// -------------------------------- Make Hisotgram --------------------------------------------
	
    /*
    int max = 0;
    int maxID = 0;
    int maxV = myGraph->getVertexCount();
    unsigned int *outDegreeHist = new unsigned int[11];
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
   
    clock_t clock5 = clock(); //after max
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
	*/

// ----------------------------- print clocks ----------------------------------------------
	cout << "\nTotal number of edges is " << myGraph->getTotalDegree() << endl;
    // cout << "Vertex " << maxID << " has the largest out degrees which is: " << max << endl;
	
	// print out the vertex with most edges (optional fill free to comment out)
	// z.id = maxID;
    // myGraph->printVneighbors(z);	
	
    // float ttime = (float)(clock1 - clock0) / (float)CLOCKS_PER_SEC;
    // cout << "\nThe time taken to allocate the graph: " << ttime <<endl;
    // ttime = (float)(clock2 - clock1) / (float)CLOCKS_PER_SEC;
    // cout << "The time taken to populate the graph: " << ttime <<endl;
    // ttime = (float)(clock3 - clock2) / (float)CLOCKS_PER_SEC;
    // cout << "The time taken to populate G2: " << ttime <<endl;
    // ttime = (float)(clock4 - clock3) / (float)CLOCKS_PER_SEC;
    // cout << "The time taken to make Histogram: " << ttime <<endl;
    // ttime = (float)(clock5 - clock4) / (float)CLOCKS_PER_SEC;
    // cout << "The time taken to find max: " << ttime <<endl;

	cout << "\nThe unused space is " << myGraph->getUnusedSpace() << endl;
	cout << "End of Main" << endl;
	
    return 0;
}
