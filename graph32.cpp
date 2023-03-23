// File: graph.cpp
// Author: Dr. Stephen Wheat, Austin Lehman, Samuel Udall, Michael Vandusen, Andrew Westlund
// Date: 1/31/22
// Pupose: This program will recieve a graph and will create a list of all the vertices
// and will list all the vertices adjacent to each vertex.
// The both the list and the lists of adjacent vertices will be dynamically allocated
// and will be sorted least to greatest.

#include <cstdio>
#include <iostream>			// Provides cout
#include <cassert>
#include <cmath>			// Provides roof
#include "graph32.h"			// Provides header file containing definitions

using namespace std;

// -------------------------------------- Neighbors -----------------------------------------

// Constructor for neighbors object
Neighbors::Neighbors() {
	nList = NULL;
	degree = 0;
	capacity = 0;
	myPos = 0;
}

// Destructor for neighbors object
Neighbors::~Neighbors() {
	delete [] nList;
}

// WILL COMPLETELY DELETE OLD NEIGHBOR
void Neighbors::setNeighbors(UINT32 newCapacity) {
	Edge *newList = new Edge[newCapacity];
	delete [] nList;
	nList = newList;
	degree = 0;
	capacity = newCapacity;
	myPos = 0;
}

void Neighbors::addNeighbor0(Edge e) {
	UINT32 j;
	for (j = degree; j > 0 ; j--) {
		// moves array entries greater than vID n to the right to insert n in order
		if (nList[j-1].vid.id > e.vid.id) {
			nList[j] = nList[j-1];
		} else {
			break;
		}
	}
	
	nList[j] = e;		// Add it to its correct spot
	// degree++;
}

bool Neighbors::isNeighbor(vID v, UINT32 &index) {
	#ifdef DEBUG
		fprintf(stderr, "Begin Neighbors::isNeighbor\n");
	#endif
	bool isANeighbor= false;	// Equals false if v is not a neighbor
	// Limits are there since if upperBound and mid were 0 and the upperBound were
	// to decrement, it would go to 2^32 - 1, keeping the loop going. Likewise (but
	// far less likely) for lowerBound's limit in the opposite direction.
	UINT32 lowerBound, upperBound, lowerBoundLimit, upperBoundLimit, mid;
	lowerBound = 0;
	lowerBoundLimit = pow(2, (8 * sizeof(degree))) - 1;
	// Set the upper bound to start as the last possible position of the neighbor list
	upperBoundLimit = 0;
	upperBound = degree - 1;
	#ifdef DEBUG
		fprintf(stderr, "isNeighbor beginning variables:\n");
		fprintf(stderr, "\tlowerBound: %u; upperBound: %u; degree: %u\n\t", lowerBound, upperBound, degree);
	
		for (UINT32 a = 0; a < degree; a++) {
			fprintf(stderr, "Array[%u]: %llu", a, nList[a].vid.id);
		}
	#endif

	// fprintf(stderr, "Lower Bound: %u, Upper Bound: %u\n", lowerBound, upperBound);
	if (degree > 0) {
		// While the indicated id has not been found and the entire list has not been searched
		while (lowerBound <= upperBound) {
			#ifdef DEBUG
				fprintf(stderr, "\tlowerBound: %u; upperBound: %u\n", lowerBound, upperBound);
			#endif
			// Set mid to be inbetween the lower and upper bounds
			
			mid = lowerBound + ((upperBound - lowerBound) / 2);
			#ifdef DEBUG
				fprintf(stderr, "\tmid: %u\n", mid);
			#endif
			// fprintf(stderr, "Mid: %lld", mid);
		
			if ((nList[mid]).vid.id > v.id) {
				if (mid == upperBoundLimit) {
					break;
				}
				upperBound = mid - 1;
			} else if ((nList[mid]).vid.id < v.id) {
				if (mid == lowerBoundLimit) {
					break;
				}
				lowerBound = mid + 1;
			} else {
				isANeighbor= true;
				break;
			}
		}
	}
	
	// If the vertex with the indicated ID is not pointed to by the current vertex, return -1
	// if (!isANeighbor) {
		// mid = 0;
	// }
	
	index = mid;
	#ifdef DEBUG
		fprintf(stderr, "End Neighbors::isNeighbor\n");
	#endif
    return isANeighbor;
}

bool Neighbors::addNeighbor(Edge e, INT16 allocation) {
	#ifdef DEBUG
		fprintf(stderr, "Begin Neighbors::addNeighbor\n");
	#endif
	bool rc = false;
	UINT32 factor = 2;
	UINT32 index;
	bool isANeighbor = isNeighbor(e.vid, index);
	// fprintf(stderr, "Before isNeighbor.\n");
	if (!isANeighbor) {		// if the neighbor does not already exist.
		// fprintf(stderr, "Not neighbor.\n");
		if (degree == capacity) {	// if the list needs to be expanded
			try {
				UINT32 degOverFac = degree / factor;
				UINT32 incSize = (allocation > degOverFac) ? allocation : degOverFac;
				UINT32 i = 0;
				Edge *newList = new Edge[capacity + incSize];
				
				for (; i < capacity; i++) {
					if (nList[i].vid.id > e.vid.id) {
						break;
					}
					newList[i] = nList[i];				// copy array over to the new bigger array
				}
				
				newList[i] = e;
				i++;
				
				for (; i < capacity + 1; i++) {
					newList[i] = nList[i - 1];
				}
				
				delete [] nList;						// delete old array.
				nList = newList;
				capacity += incSize;
			} catch (bad_alloc& badAlloc) {
				cerr << "Not enough memory: " << badAlloc.what() << endl;
			}
		} else {
			addNeighbor0(e);
		}
		degree++;
	}
	#ifdef DEBUG
		fprintf(stderr, "End Neighbors::addNeighbor\n");
	#endif
	return !isANeighbor;
}

// Meant to be used to simply append a value to the end of the list; will reallocate if needed,
// but createV should allocate the amount of space necessary for it to not be needed;
// No safeguard from adding an edge that is already there.
void Neighbors::addNeighbor2(Edge e, INT16 allocation) {
	UINT32 factor = 2;
	if (degree == capacity) {	// if the list needs to be expanded
		try {
			UINT32 degOverFac = degree / factor;
			UINT32 incSize = (allocation > degOverFac) ? allocation : degOverFac;
			UINT32 i = 0;
			Edge *newList = new Edge[capacity + incSize];
				
			for (; i < degree; i++) {
				newList[i] = nList[i];				// copy array over to the new bigger array
			}
				
			newList[i] = e;						// i == degree, appends e to end	
			delete [] nList;						// delete old array.
			nList = newList;
			capacity += incSize;
		} catch (bad_alloc& badAlloc) {
			cerr << "Not enough memory: " << badAlloc.what() << endl;
		}
	} else {
		nList[degree] = e;
	}
	degree++;
}

bool Neighbors::delNeighbor(Edge e) {
    bool rc;
	UINT32 index;
    
	rc = isNeighbor(e.vid, index);
    if (rc) {	// If the edge exists
		
		// Create new array with decremented capacity.
		Edge *temp = new Edge[capacity - 1];
		
		for (UINT32 i = 0; i < index; i++) {
			temp[i] = nList[i];
		}

		for (UINT32 i = index; i < capacity - 1; i++) {
			temp[i] = nList[i + 1];
		}

		delete [] nList;
		nList = temp;
	
		capacity--;
		degree--;
	}

    return rc;
}

UINT32 Neighbors::getWeight(vID v) {
	UINT32 index;
	bool isANeighbor = isNeighbor(v, index);
	UINT32 weight = 0;
	if (isANeighbor) {
		weight = nList[index].weight;
	}
	return weight;
}

void Neighbors::setWeight(vID v, UINT32 w) {		// set the weight of an edge
	UINT32 index;
	bool isANeighbor= isNeighbor(v, index);
	UINT32 weight = 0;
	if (isANeighbor) {
		nList[index].weight = w;
	}
}

void Neighbors::incPos() { 	// my position within myself tracking
	myPos++;
}

void Neighbors::decPos() {
	myPos--;
}

UINT32 Neighbors::getDegree() const {
	return degree;
}

UINT32 Neighbors::getCapacity() const {
	return capacity;
}

UINT32 Neighbors::getMyPos() const {
	return myPos;
}

void Neighbors::setMyPos(UINT32 newPos) {
	myPos = newPos;
}

Edge * Neighbors::getNList() const {
	return nList;
}

void Neighbors::printNeighbors() const {
    for (UINT32 j = 0; j < degree; j++) {	// print each value in nList
		cout << (nList[j].vid.id) << " ";
    }
    cout << endl;
}

void Neighbors::printNeighborsWeights() const {
    for (UINT32 j = 0; j < degree; j++) {	// print the weights
		cout << (nList[j].weight) << " ";
    }
    cout << endl;
}

void Neighbors::clearNeighbors() {
	if (nList) {
		delete [] nList;
	}
	nList = NULL;
	degree = 0;
	capacity = 0;
	myPos = 0;
}

// -------------------------------------- Vertex -----------------------------------------

// Constructor for vertex object
vertex::vertex() {
	// forward = new Neighbors();
	// forward.Neighbors();
	// backward = new Neighbors();
	// backward.Neighbors();
}

// Destructor for vertex object
vertex::~vertex() {
	// clearNeighbors();
    forward.~Neighbors();
	// clearBackNeighbors();
	backward.~Neighbors();
}

void vertex::setVertex(UINT32 capacity, UINT32 backCapacity) {
	forward.setNeighbors(capacity);
	backward.setNeighbors(backCapacity);
}

bool vertex::isNeighbor(vID v, UINT32 &index) {
	return forward.isNeighbor(v, index);
}

bool vertex::isBackNeighbor(vID v, UINT32 &index) {
	return backward.isNeighbor(v, index);
}

bool vertex::addNeighbor(Edge e, INT16 allocation) {
	#ifdef DEBUG
		fprintf(stderr, "Do vertex::addNeighbor\n");
	#endif
	return forward.addNeighbor(e, allocation);
}

void vertex::addNeighbor2(Edge e, INT16 allocation) {
	forward.addNeighbor2(e, allocation);
}

bool vertex::addBackNeighbor(Edge e, INT16 allocation) {
	#ifdef DEBUG
		fprintf(stderr, "Do vertex::addBackNeighbor\n");
	#endif
	return backward.addNeighbor(e, allocation);
}

void vertex::addBackNeighbor2(Edge e, INT16 allocation) {
	backward.addNeighbor2(e, allocation);
}

bool vertex::delNeighbor(Edge e) {
	return forward.delNeighbor(e);
}

bool vertex::delBackNeighbor(Edge e) {
	return backward.delNeighbor(e);
}

UINT32 vertex::getWeight(vID v) {
	return forward.getWeight(v);
}

void vertex::setWeight(vID v, UINT32 w) {		// set the weight of an edge
	forward.setWeight(v, w);
}

UINT32 vertex::getBackWeight(vID v) {
	return backward.getWeight(v);
}

void vertex::setBackWeight(vID v, UINT32 w) {		// set the weight of an edge
	backward.setWeight(v, w);
}

UINT32 vertex::getMyPos() const {
	return forward.getMyPos();
}

void vertex::setMyPos(UINT32 newPos) {
	forward.setMyPos(newPos);
}

void vertex::incPos() { 	// my position within myself tracking
	forward.incPos();
}

void vertex::decPos() {
	forward.decPos();
}

UINT32 vertex::getMyBackPos() const {
	return backward.getMyPos();
}

void vertex::setMyBackPos(UINT32 newBackPos) {
	backward.setMyPos(newBackPos);
}

void vertex::incBackPos() {
	backward.incPos();
}

void vertex::decBackPos() {
	backward.decPos();
}

UINT32 vertex::getOutDegree() const {
    return forward.getDegree();
}

UINT32 vertex::getInDegree() const {
    return backward.getDegree();
}

UINT32 vertex::getCapacity() const {
	return forward.getCapacity();
}

UINT32 vertex::getBackCapacity() const {
	return backward.getCapacity();
}

Edge * vertex::getNeighbors() const {
	return forward.getNList();
}

Edge * vertex::getBackNeighbors() const {
	return backward.getNList();
}

void vertex::printNeighbors() const {
    forward.printNeighbors();
}

void vertex::printNeighborsWeights() const {
    forward.printNeighborsWeights();
}

void vertex::printBackNeighbors() const {
    backward.printNeighbors();
}

void vertex::printBackNeighborsWeights() const {
    backward.printNeighborsWeights();
}

// Reset the list of neighbors
void vertex::clearNeighbors() {
    forward.clearNeighbors();
}

// Reset the list of back-neighbors
void vertex::clearBackNeighbors() {
	backward.clearNeighbors();
}

// -------------------------------------- Graph -----------------------------------------

// Constructor for graph
graph::graph(UINT64 numberOfVertices) {
	vertexCount = numberOfVertices;
    vList = new vertex[vertexCount];
    totalDegree = 0;
	allocation = 1;
}

// Destructor for graph
graph::~graph() {
    clear();
    delete [] vList;
}

// Private
Edge graph::getEdge(vID v, UINT32 weight) {
	Edge e;
	e.vid = v;
	e.weight = weight;
	// cout << "e.vid: " << e.vid.id << ", e.weight: " << e.weight << endl;
	return e;
}

void graph::setAllocation(INT16 a) {
	if (a < 1) {
		allocation = 1;
	} else if (a > 16) {
		allocation = 16;
	} else {
		allocation = a;
	}
}

void graph::createV(vID v, UINT32 capacity, UINT32 backCapacity) {
	assert(v.id < vertexCount);
	vList[v.id].setVertex(capacity, backCapacity);
}

bool graph::isEdge(vID u, vID v, UINT32 weight) {
    // If vertex u points to vertex v, return true
	UINT32 index;
	return vList[u.id].isNeighbor(v, index);
}

bool graph::isBackEdge(vID u, vID v, UINT32 weight) {
    // If vertex u points to vertex v, return true
	UINT32 index;
    return vList[u.id].isBackNeighbor(v, index);
}

bool graph::addEdge(vID u, vID v, UINT32 weight, bool inc) {
	assert(u.id < vertexCount);
	assert(v.id < vertexCount);
	#ifdef DEBUG
		fprintf(stderr, "Begin graph::addEdge\n");
	#endif
    bool rc;
	Edge e = getEdge(v, weight);
	rc = vList[u.id].addNeighbor(e, allocation);
	
	if (rc) {
		// If the edge is successfully added, increase vertex v's inDegree by 1 and add u to v's back-neighbor list
		Edge f = getEdge(u, weight);
		vList[v.id].addBackNeighbor(f, allocation);
		totalDegree++;

		if (u.id > v.id) {
			vList[u.id].incPos();
		} else if (v.id > u.id){
			vList[v.id].incBackPos();
		}
	} else {		// if the edge exists
		if (inc) {	// if the user wants to increase the weight
			UINT32 currWeight = vList[u.id].getWeight(v);
			UINT32 currBackWeight = vList[v.id].getBackWeight(v);
			vList[u.id].setWeight(v, (currWeight + weight));
			vList[v.id].setBackWeight(u, (currBackWeight + weight));
		// CONSIDER NOT SETTING IT IF INC IS FALSE AND MAKING RC = ISNEIGHBOR || INC
		} else {	// otherwise, set the weight
			vList[u.id].setWeight(v, weight);
			vList[v.id].setBackWeight(u, weight);
		}
	}
	// Will only return if a new edge has been added, not if it has been set or incremented
	#ifdef DEBUG
		fprintf(stderr, "End graph::addEdge\n");
	#endif
    return rc;
}

void graph::addEdge2(vID u, vID v, UINT32 weight) {
	assert(u.id < vertexCount);
	assert(v.id < vertexCount);
	Edge e = getEdge(v, weight), f = getEdge(u, weight);
	vList[u.id].addNeighbor2(e, allocation);
	vList[v.id].addNeighbor2(f, allocation);
	totalDegree++;
}

bool graph::delEdge(vID u, vID v, UINT32 weight) {
    bool rc;
	Edge e = getEdge(v, weight);
    rc = vList[u.id].delNeighbor(e);

    if (rc) {
		// If the edge is successfully deleted, decrease vertex v's inDegree by 1
		Edge f = getEdge(u, weight);
		vList[v.id].delBackNeighbor(f);
		totalDegree--;
		
		if (u.id > v.id) {
			vList[u.id].decPos();
		} else if (v.id > u.id) {
			vList[v.id].decBackPos(); 
		}
    }

    return rc;
}



void graph::dumpGraph() {
	Edge *E;
	UINT64 *A;
	int n = 0;
	FILE * pfile;
    pfile = fopen("dumpOutput.txt","w");
	for (UINT64 i = 0; i < vertexCount; i++) {
		E = vList[i].getNeighbors();
		n = vList[i].getOutDegree();
		for (UINT32 j = 0; j < n; j++) {
			// fprintf(pfile, "%u %u \n",i,E[j].vid.id);
		}
	}
	
    fclose(pfile);
	
}

UINT32 graph::getVoutDegree(vID v) const {
	assert(v.id < vertexCount);
    return (vList[v.id].getOutDegree());
}

UINT32 graph::getVinDegree(vID v) const {
	assert(v.id < vertexCount);
    return (vList[v.id].getInDegree());
}

UINT32 graph::getVmyPos(vID v) const {
	assert(v.id < vertexCount);
	return (vList[v.id].getMyPos());
}

void graph::setVmyPos(vID v, UINT32 newPos) {
	assert(v.id < vertexCount);
	vList[v.id].setMyPos(newPos);
}

UINT32 graph::getVmyBackPos(vID v) const {
	assert(v.id < vertexCount);
	return (vList[v.id].getMyBackPos());
}

void graph::setVmyBackPos(vID v, UINT32 newBackPos) {
	assert(v.id < vertexCount);
	vList[v.id].setMyBackPos(newBackPos);
}

UINT64 graph::getVertexCount() const {
    // Return the number of vertices there are in the entire adjacency list
    return vertexCount;
}

Edge * graph::getVneighbors(vID v) {
	assert(v.id < vertexCount);
	return vList[v.id].getNeighbors();
}

Edge * graph::getVbackNeighbors(vID v) {
	assert(v.id < vertexCount);
	return vList[v.id].getBackNeighbors();
}

UINT64 graph::getTotalDegree() const {
    // Return |E| - this only works for a directed graph
    return totalDegree;
}

// DELETE LATER - DEBUGGING / INFO PURPOSES
UINT64 graph::getUnusedSpace() const {
	UINT64 space = 0;
	for (int i = 0; i < vertexCount; i++) {
		space += (vList[i].getCapacity() - vList[i].getOutDegree());
	}
	return space;
}

// DELETE LATER - DEBUGGING / INFO PURPOSES
UINT32 graph::getVcapacity(vID v) const {
	return vList[v.id].getCapacity();
}

// DELETE LATER - DEBUGGING / INFO PURPOSES
UINT32 graph::getVbackCapacity(vID v) const {
	return vList[v.id].getBackCapacity();
}

void graph::printVneighbors(vID v) {
    // Print out the neighbors for vertex v
    cout << "Vertex " << v.id << "'s neighbors are: ";
    vList[v.id].printNeighbors();
}

void graph::printVbackNeighbors(vID v) {
    // Print out the back neighbors for vertex v
    cout << "Vertex " << v.id << "'s backedge neighbors are: ";
    // vList[v.id].printBackNeighbors(numOfEntities);
    vList[v.id].printBackNeighbors();
}

void graph::printVneighborsWeights(vID v) {
    // Print out the neighbors for vertex v
    cout << "Vertex " << v.id << "'s weights are: ";
    vList[v.id].printNeighborsWeights();
}

void graph::printGraph() {
    // Print out the ID of each vertex that the indicated vertex points to
    for (UINT64 i = 0; i < vertexCount; i++) {
		cout << "Vertex " << i << "'s neighbors are: ";
		vList[i].printNeighbors();
		cout << "Vertex " << i << "'s weights are:   ";
		vList[i].printNeighborsWeights();
		cout << "Vertex " << i << "'s back-neighbors are: ";
		vList[i].printBackNeighbors();
		cout << "Vertex " << i << "'s back-weights are:   ";
		vList[i].printBackNeighborsWeights();	
		cout << endl;
    }
}

void graph::printDebugGraph() {
	vID v;
	UINT32 outDegree, inDegree;
	UINT64 sumOfOutDegrees = 0, sumOfInDegrees = 0;
	for (v.id = 0; v.id < vertexCount - 1; v.id++) {
		outDegree = getVoutDegree(v);
		inDegree = getVinDegree(v);
		if (outDegree > 0 || inDegree > 0) {
			fprintf(stderr, "Data for vertex %llu:\nOutDegree: %u\tCapacity: %u\n", v.id, outDegree, getVcapacity(v));
			// printNeighbors();
			fprintf(stderr, "InDegree: %u\tBackCapacity: %u\n", inDegree, getVbackCapacity(v));
			sumOfOutDegrees += outDegree;
			sumOfInDegrees += inDegree;
		}
	}
	fprintf(stderr, "Total number of out-degrees: %llu\nTotal number of in-degrees: %llu\n", sumOfOutDegrees, sumOfInDegrees);
}

// Clear the neighbor list of a particular vertex
void graph::clearNlist(vID v) {
	Edge *clearList = getVneighbors(v);
	UINT64 currID;
	Edge currEdge;
	
	// Decrement all of the vertices affected
	for (UINT32 i = 0; i < vList[v.id].getOutDegree(); i++) {
		currID = clearList[i].vid.id;
		currEdge = getEdge(v, 1);
		vList[currID].delBackNeighbor(currEdge);
	}
	// Decrease total degree by the appropriate amount
	totalDegree = totalDegree - vList[v.id].getOutDegree();
	vList[v.id].clearNeighbors();
}

// Public
void graph::clear() {
    for (UINT64 i = 0; i < vertexCount; i++) {
		vList[i].clearNeighbors();
		vList[i].clearBackNeighbors();
    }
    totalDegree = 0;
}
