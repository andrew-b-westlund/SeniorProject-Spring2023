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
#include "graph.h"			// Provides header file containing definitions

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
// Useful when you know how many edges the neighbor is going to have; paired with addNeighbor2
void Neighbors::setNeighbors(UINT32 newCapacity) {
	Edge *newList = new Edge[newCapacity];
	delete [] nList;
	nList = newList;
	degree = 0;
	capacity = newCapacity;
	myPos = 0;
}

// O(n/2) == O(n) where n is the number of edges on the neighbors list
void Neighbors::addNeighbor0(Edge e) {
	UINT32 j;
	for (j = degree; j > 0 ; j--) {
		// moves array entries greater than vID n to the right to insert the edge where it belongs
		if (nList[j-1].vid.id > e.vid.id) {
			nList[j] = nList[j-1];
		} else {
			break;
		}
	}
	
	nList[j] = e;		// Add it to its correct spot
}

// Preconditions: The edge list is in-order
// Postcondition: A binary search is performed to find if there is an edge pointing to v;
// if there is, return true and set index to v's position on the edge array; otherwise, return false
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
	// Set the upper bound to start as the last position of the edge list
	upperBoundLimit = 0;
	upperBound = degree - 1;
	#ifdef DEBUG
		fprintf(stderr, "isNeighbor beginning variables:\n");
		fprintf(stderr, "\tlowerBound: %u; upperBound: %u; degree: %u\n\t", lowerBound, upperBound, degree);
		
		// Print every edge
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
			
			// Set mid to be in between the lower and upper bounds
			mid = lowerBound + ((upperBound - lowerBound) / 2);
			#ifdef DEBUG
				fprintf(stderr, "\tmid: %u\n", mid);
			#endif
			// fprintf(stderr, "Mid: %lld", mid);
		
			if ((nList[mid]).vid.id > v.id) {
				// Lower end of list has been reached, prevent overflow
				if (mid == upperBoundLimit) {
					break;
				}
				upperBound = mid - 1;
			} else if ((nList[mid]).vid.id < v.id) {
				// Upper end of list has been reached, prevent overflow
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
	
	// Index is set whether or not a neighbor is found
	index = mid;
	#ifdef DEBUG
		fprintf(stderr, "End Neighbors::isNeighbor\n");
	#endif
    return isANeighbor;
}

// Preconditions: The edge list is in order
// Postconditions: If the edge is not on the array yet, insert the edge into where it belongs and return true;
// otherwise, return false
bool Neighbors::addNeighbor(Edge e, INT16 allocation) {
	#ifdef DEBUG
		fprintf(stderr, "Begin Neighbors::addNeighbor\n");
	#endif
	bool rc = false;
	// If the edge list ends up being massive, we'll want to allocate more memory than the 
	// set allocation size to save time on reallocation
	UINT32 factor = 2;
	UINT32 index;	// Could potentially be used, but serves to fulfill method params
	bool isANeighbor = isNeighbor(e.vid, index);
	
	if (!isANeighbor) {				// If the neighbor does not already exist.
		if (degree == capacity) {	// If the list needs to grow
			try {
				UINT32 degOverFac = degree / factor;
				UINT32 incSize = (allocation > degOverFac) ? allocation : degOverFac;
				UINT32 i = 0;
				Edge *newList = new Edge[capacity + incSize];
				
				// Copy every edge over to the bigger array until i is the index of where the edge to be added 
				for (; i < capacity; i++) {
					if (nList[i].vid.id > e.vid.id) {
						break;
					}
					newList[i] = nList[i];
				}
				
				// Put the added edge where it belongs
				newList[i] = e;
				i++;
				
				// Copy over the rest of the edges to the bigger array
				for (; i < capacity + 1; i++) {
					newList[i] = nList[i - 1];
				}
				
				delete [] nList;
				nList = newList;
				capacity += incSize;
			} catch (bad_alloc& badAlloc) {
				cerr << "Not enough memory: " << badAlloc.what() << endl;
			}
		} else {	// If the array has space
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
// No safeguard from adding an edge that is already there or adding edges out of order.
void Neighbors::addNeighbor2(Edge e, INT16 allocation) {
	UINT32 factor = 2;
	if (degree == capacity) {	// If the list needs to expand
		try {
			UINT32 degOverFac = degree / factor;
			UINT32 incSize = (allocation > degOverFac) ? allocation : degOverFac;
			UINT32 i = 0;
			Edge *newList = new Edge[capacity + incSize];
				
			// Copy every edge over to the bigger array
			for (; i < degree; i++) {
				newList[i] = nList[i];
			}
				
			newList[i] = e;						// i == degree, appends e to end	
			delete [] nList;					// delete old array.
			nList = newList;
			capacity += incSize;
		} catch (bad_alloc& badAlloc) {
			cerr << "Not enough memory: " << badAlloc.what() << endl;
		}
	} else {	// If there is enough space, append the edge
		nList[degree] = e;
	}
	degree++;
}

// Preconditions: The edge array is in order
// Postconditions: If the passed edge is on the array, is removed from the list,
// the array's size and degree are decremented, and true is returned;
// otherwise, return false
// Could be improved to only create a new array of a new size when degree and capacity 
// have a certain difference in size and not change capacity unless that occurs
bool Neighbors::delNeighbor(Edge e) {
	UINT32 index;
	bool isANeighbor = isNeighbor(e.vid, index);
    if (isANeighbor) {	// If the edge exists
		
		// Create new array with decremented capacity.
		Edge *temp = new Edge[capacity - 1];
		
		// Copy every edge less than the given edge over
		for (UINT32 i = 0; i < index; i++) {
			temp[i] = nList[i];
		}

		// Copy every edge greater than the given edge over
		for (UINT32 i = index; i < capacity - 1; i++) {
			temp[i] = nList[i + 1];
		}

		delete [] nList;
		nList = temp;
	
		capacity--;
		degree--;
	}

    return isANeighbor;
}

// Preconditions: None
// Postconditions: If an edge with the passed ID exists, return the edge's weight;
// otherwise, return 0
UINT32 Neighbors::getWeight(vID v) {
	UINT32 index;
	bool isANeighbor = isNeighbor(v, index);
	UINT32 weight = 0;
	if (isANeighbor) {
		weight = nList[index].weight;
	}
	return weight;
}

// Preconditions: None
// Postconditions: If an edge with the passed ID exists, set the edge's weight to
// be the passed weight
void Neighbors::setWeight(vID v, UINT32 w) {		// set the weight of an edge
	UINT32 index;
	bool isANeighbor= isNeighbor(v, index);
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
    for (UINT32 j = 0; j < degree; j++) {	// Print each edge
		cout << (nList[j].vid.id) << " ";
    }
    cout << endl;
}

void Neighbors::printNeighborsWeights() const {
    for (UINT32 j = 0; j < degree; j++) {	// Print each edge and its weight
		cout << "vID: " << (nList[j].vid.id) << " - weight: " << (nList[j].weight) << "; ";
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
	// Nothing needs to happen, the forward and backward objects already exist
}

// Destructor for vertex object
vertex::~vertex() {
	// clearNeighbors();
    forward.~Neighbors();
	// clearBackNeighbors();
	backward.~Neighbors();
}

// Preconditions: None
// Postconditions: Forward will be set to have an empty edge list with
// the passed number of places available
void vertex::setVertexForward(UINT32 capacity) {
	forward.setNeighbors(capacity);
}

// Preconditions: None
// Postconditions: Backward will be set to have an empty edge list with
// the passed number of places available
void vertex::setVertexBackward(UINT32 backCapacity) {
	backward.setNeighbors(backCapacity);
}

// Preconditions: The vertex's edge array is in order
// Postconditions: If the vertex has an edge pointing to the specified
// ID, index is set to that edge's position and true is returned;
// otherwise, index is still set but false is returned
bool vertex::isNeighbor(vID v, UINT32 &index) {
	return forward.isNeighbor(v, index);
}

// Preconditions: The vertex's back-edge array is in order
// Postconditions: If the vertex is pointed to by the specified
// ID, index is set to that back-edge's position and true is returned;
// otherwise, index is still set but false is returned
bool vertex::isBackNeighbor(vID v, UINT32 &index) {
	return backward.isNeighbor(v, index);
}

// Preconditions: The vertex's edge array is in order
// Postconditions: If an edge from the vertex to e's vID does not exist,
// add the edge to the vertex's edge array (if it's full, increase its size); 
// otherwise, return false
bool vertex::addNeighbor(Edge e, INT16 allocation) {
	#ifdef DEBUG
		fprintf(stderr, "Do vertex::addNeighbor\n");
	#endif
	return forward.addNeighbor(e, allocation);
}

// Meant to be used to simply append an edge to the end of the vertex's edge array; will reallocate if needed,
// but createV should allocate the amount of space necessary for it to not be needed;
// No safeguard from adding an edge that is already there or adding edges out of order.
void vertex::addNeighbor2(Edge e, INT16 allocation) {
	forward.addNeighbor2(e, allocation);
}

// Preconditions: The vertex's back-edge array is in order
// Postconditions: If a back-edge from the vertex to e's vID does not exist,
// add the back-edge to the vertex's back-edge array (if it's full, increase its size); 
// otherwise, return false
bool vertex::addBackNeighbor(Edge e, INT16 allocation) {
	#ifdef DEBUG
		fprintf(stderr, "Do vertex::addBackNeighbor\n");
	#endif
	return backward.addNeighbor(e, allocation);
}

// Meant to be used to simply append a back-edge to the end of the vertex's back-edge array; will reallocate if needed,
// but createV should allocate the amount of space necessary for it to not be needed;
// No safeguard from adding a back-edge that is already there or adding back-edges out of order.
void vertex::addBackNeighbor2(Edge e, INT16 allocation) {
	backward.addNeighbor2(e, allocation);
}

// Preconditions: The edge array is in order
// Postconditions: If the passed edge exists, it is removed from the vertex's edge array,
// the array's capacity and degree are decremented, and true is returned;
// otherwise, return false
bool vertex::delNeighbor(Edge e) {
	return forward.delNeighbor(e);
}

// Preconditions: The back-edge array is in order
// Postconditions: If the passed back-edge exists, it is removed from the vertex's back-edge array,
// the back-edge array's capacity and degree are decremented, and true is returned;
// otherwise, return false
bool vertex::delBackNeighbor(Edge e) {
	return backward.delNeighbor(e);
}

// Returns the weight of the vertex's edge to v; if the edge doesn't exist, 0 is returned
UINT32 vertex::getWeight(vID v) {
	return forward.getWeight(v);
}

// Set the weight of the vertex's edge to v to be the passed weight if the edge exists
void vertex::setWeight(vID v, UINT32 w) {
	forward.setWeight(v, w);
}

// Returns the weight of the vertex's back-edge to v; if the back-edge doesn't exist, 0 is returned
UINT32 vertex::getBackWeight(vID v) {
	return backward.getWeight(v);
}

// Set the weight of the vertex's back-edge to v to be the passed weight if the back-edge exists
void vertex::setBackWeight(vID v, UINT32 w) {
	backward.setWeight(v, w);
}

// The vertex's hypothetical position on its own edge array is returned
UINT32 vertex::getMyPos() const {
	return forward.getMyPos();
}

// The vertex's hypothetical position on its own edge array is set to the passed value
void vertex::setMyPos(UINT32 newPos) {
	forward.setMyPos(newPos);
}

// The vertex's hypothetical position on its own edge array is incremented
void vertex::incPos() {
	forward.incPos();
}

// The vertex's hypothetical position on its own edge array is decremented
void vertex::decPos() {
	forward.decPos();
}

// The vertex's hypothetical position on its own back-edge array is returned
UINT32 vertex::getMyBackPos() const {
	return backward.getMyPos();
}

// The vertex's hypothetical position on its own back-edge array is set to the passed value
void vertex::setMyBackPos(UINT32 newBackPos) {
	backward.setMyPos(newBackPos);
}

// The vertex's hypothetical position on its own back-edge array is incremented
void vertex::incBackPos() {
	backward.incPos();
}

// The vertex's hypothetical position on its own back-edge array is decremented
void vertex::decBackPos() {
	backward.decPos();
}

// The number of edges the vertex has is returned
UINT32 vertex::getOutDegree() const {
    return forward.getDegree();
}

// The number of back-edges the vertex has is returned
UINT32 vertex::getInDegree() const {
    return backward.getDegree();
}

// Not so sure if we want this to be public
UINT32 vertex::getCapacity() const {
	return forward.getCapacity();
}

// Not so sure if we want this to be public
UINT32 vertex::getBackCapacity() const {
	return backward.getCapacity();
}

// Returns the vertex's edge array
Edge * vertex::getNeighbors() const {
	return forward.getNList();
}

// Returns the vertex's back-edge array
Edge * vertex::getBackNeighbors() const {
	return backward.getNList();
}

// Prints the IDs of the vertex's edges
void vertex::printNeighbors() const {
    forward.printNeighbors();
}

// Prints the IDs of the vertex's edges and the edges' weights
void vertex::printNeighborsWeights() const {
    forward.printNeighborsWeights();
}

// Prints the IDs of the vertex's back-edges
void vertex::printBackNeighbors() const {
    backward.printNeighbors();
}

// Prints the IDs of the vertex's back-edges and the back-edges' weights
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

// Private, used to have an edge to pass to other functions
Edge graph::getEdge(vID v, UINT32 weight) {
	Edge e;
	e.vid = v;
	e.weight = weight;
	// cout << "e.vid: " << e.vid.id << ", e.weight: " << e.weight << endl;
	return e;
}

// Sets how much the edge/back-edge arrays will grow by when full
void graph::setAllocation(INT16 a) {
	if (a < 1) {
		allocation = 1;
	} else if (a > 16) {
		allocation = 16;
	} else {
		allocation = a;
	}
}

// Sets v's vertex's edge array to be empty and have capacity empty places
// Sets v's vertex's back-edge array to be empty and have backCapacity empty places
void graph::createV(vID v, UINT32 capacity, UINT32 backCapacity) {
	createVForward(v, capacity);
	createVBackward(v, backCapacity);
}

// Sets v's vertex's edge array to be empty and have capacity empty places
void graph::createVForward(vID v, UINT32 capacity) {
	assert(v.id < vertexCount);
	vList[v.id].setVertexForward(capacity);
}

// Sets v's vertex's back-edge array to be empty and have backCapacity empty places
void graph::createVBackward(vID v, UINT32 backCapacity) {
	assert(v.id < vertexCount);
	vList[v.id].setVertexBackward(backCapacity);
}

// If vertex u points to vertex v, return true
bool graph::isEdge(vID u, vID v) {
	UINT32 index;
	return vList[u.id].isNeighbor(v, index);
}

// If vertex u points to vertex v, return true and set index to the edge's position in u's edge array
bool graph::isEdge(vID u, vID v, UINT32 &index) {
	return vList[u.id].isNeighbor(v, index);
}

// If vertex u is pointed to by vertex v, return true
bool graph::isBackEdge(vID u, vID v) {
	UINT32 index;
    return vList[u.id].isBackNeighbor(v, index);
}

// If vertex u is pointed to by vertex v, return true and set index to the back-edge's position in u's back-edge array
bool graph::isBackEdge(vID u, vID v, UINT32 &index) {
    return vList[u.id].isBackNeighbor(v, index);
}

// Preconditions: vertex(u) and vertex(v) must be existing vertices
// Postconditions: If inc is false, then if u has an edge to v, the weight from u to v will be overwritten by the  
// new weight; if u does not have an edge to v, then an edge with the passed weight from u to v will be added,
// and a similar back-edge from v to u will be added.
// If inc is true, then if an edge already exists, the edge will increase by the passed weight.
// If the edge already existed, false will be returned; otherwise, return true.
bool graph::addEdge(vID u, vID v, UINT32 weight, bool inc) {
	assert(u.id < vertexCount);
	assert(v.id < vertexCount);
	#ifdef DEBUG
		fprintf(stderr, "Begin graph::addEdge\n");
	#endif
	
	Edge e = getEdge(v, weight);
	bool rc = vList[u.id].addNeighbor(e, allocation);
	
	if (rc) {
		// If the edge is successfully added, add u to v's back-neighbor list
		Edge f = getEdge(u, weight);
		vList[v.id].addBackNeighbor(f, allocation);
		totalDegree++;

		// Change the vertices's (back-)position within its own (back-)edge array as necessary
		if (u.id > v.id) {
			vList[u.id].incPos();
		} else if (v.id > u.id){
			vList[v.id].incBackPos();
		}
		
	} else {		// If the edge already exists
		if (inc) {	// If the user wants to increase the weight
			UINT32 currWeight = vList[u.id].getWeight(v);
			UINT32 currBackWeight = vList[v.id].getBackWeight(u);
			vList[u.id].setWeight(v, (currWeight + weight));
			vList[v.id].setBackWeight(u, (currBackWeight + weight));
		// CONSIDER NOT SETTING IT IF INC IS FALSE AND MAKING RC = ISNEIGHBOR || INC
		} else {	// otherwise, set the weight
			vList[u.id].setWeight(v, weight);
			vList[v.id].setBackWeight(u, weight);
		}
	}
	// Will only return true if a new edge has been added, not if it has been set or increased
	#ifdef DEBUG
		fprintf(stderr, "End graph::addEdge\n");
	#endif
    return rc;
}

// Preconditions: vertex(u) and vertex(v) must be existing vertices
// Postconditions: If a forward edge is added from u to v, a forward edge will also be added from v to u 
// without adding any backedges; if an edge from u to v already exists, then if inc is true, then the weights
// of the forward edges from u to v and vice versa will be increased by the weight amount and false is returned;
// if an edge from u to v already exists and if inc is false, then the weights of the forward edges will be set to weight.
bool graph::addG2Edge(vID u, vID v, UINT32 weight, bool inc) {
	assert(u.id < vertexCount);
	assert(v.id < vertexCount);
	#ifdef DEBUG
		fprintf(stderr, "Begin graph::addG2Edge\n");
	#endif
	
	Edge e = getEdge(v, weight);
	bool rc = vList[u.id].addNeighbor(e, allocation);
	
	if (rc) {
		// If the edge is successfully added, add u to v's neighbor list
		Edge f = getEdge(u, weight);
		vList[v.id].addNeighbor(f, allocation);
		totalDegree += 2;

		// Change the vertices's position within its own edge array as necessary
		if (u.id > v.id) {
			vList[u.id].incPos();
		} else if (v.id > u.id){
			vList[v.id].incPos();
		}
		
	} else {		// If the edge from u to v already exists
		if (inc) {	// If the user wants to increase the weight, increase it for both forward edges
			UINT32 currWeightU = vList[u.id].getWeight(v);
			UINT32 currWeightV = vList[v.id].getWeight(u);
			vList[u.id].setWeight(v, (currWeightU + weight));
			vList[v.id].setWeight(u, (currWeightV + weight));
		// CONSIDER NOT SETTING IT IF INC IS FALSE AND MAKING RC = ISNEIGHBOR || INC
		} else {	// otherwise, set the weights
			vList[u.id].setWeight(v, weight);
			vList[v.id].setWeight(u, weight);
		}
	}
	// Will only return if the new edges has been added, not if they have been set or increased
	#ifdef DEBUG
		fprintf(stderr, "End graph::addG2Edge\n");
	#endif
    return rc;
}

// Append an edge from u to v and back-edge from v to u to the appropriate (back-)neighbor lists
void graph::addEdge2(vID u, vID v, UINT32 weight) {
	// assert(u.id < vertexCount);
	// assert(v.id < vertexCount);
	// Edge e = getEdge(v, weight), f = getEdge(u, weight);
	// vList[u.id].addNeighbor2(e, allocation);
	// vList[v.id].addBackNeighbor2(f, allocation);
	// totalDegree++;
	addEdge2Forward(u, v, weight);
	addEdge2Backward(u, v, weight);
}

// Append an edge from u to v with the passed weight at the end of u's neighbor list
void graph::addEdge2Forward(vID u, vID v, UINT32 weight) {
	assert(u.id < vertexCount);
	assert(v.id < vertexCount);
	// cout << u.id << " forward" << endl;
	Edge e = getEdge(v, weight);
	vList[u.id].addNeighbor2(e, allocation);
	totalDegree++;
}

// Append a back-edge from v to u with the passed weight at the end of v's back-neighbor list
void graph::addEdge2Backward(vID u, vID v, UINT32 weight) {
	assert(u.id < vertexCount);
	assert(v.id < vertexCount);
	// cout << v.id << " backward" << endl;
	Edge f = getEdge(u, weight);
	vList[v.id].addBackNeighbor2(f, allocation);
}

// Deletes the edge from u to v and back-edge from v to u if the edges exist
bool graph::delEdge(vID u, vID v, UINT32 weight) {
	Edge e = getEdge(v, weight);
    bool rc = vList[u.id].delNeighbor(e);	// Returns if the edge existed and deletes it

    if (rc) {
		// If the edge existed and was successfully deleted, delete the back-edge from v's back-neighbor list
		Edge f = getEdge(u, weight);
		vList[v.id].delBackNeighbor(f);
		totalDegree--;
		
		// Change the vertices's position within its own (back-)edge array as necessary
		if (u.id > v.id) {
			vList[u.id].decPos();
		} else if (v.id > u.id) {
			vList[v.id].decBackPos(); 
		}
    }
    return rc;
}

// INCOMPLETE: I propose it should put the graph into a binary file, where the first number
// indicates the number of vertices, and every vertex will start with a pair of numbers indicating 
// the number of forward and backward edges the vertex has, followed by the list of forward edge IDs 
// and weights followed by the list of backward edge IDs and weights.
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

// Returns the number of edges v has
UINT32 graph::getVoutDegree(vID v) const {
	assert(v.id < vertexCount);
    return (vList[v.id].getOutDegree());
}

// Returns the number of back-edges v has
UINT32 graph::getVinDegree(vID v) const {
	assert(v.id < vertexCount);
    return (vList[v.id].getInDegree());
}

// Returns v's hypothetical position within its own neighbor list
UINT32 graph::getVmyPos(vID v) const {
	assert(v.id < vertexCount);
	return (vList[v.id].getMyPos());
}

// Sets v's hypothetical position within its own neighbor list
void graph::setVmyPos(vID v, UINT32 newPos) {
	assert(v.id < vertexCount);
	vList[v.id].setMyPos(newPos);
}

// DELETE LATER - DEBUGGING / INFO PURPOSES
UINT32 graph::getVcapacity(vID v) const {
	return vList[v.id].getCapacity();
}

// DELETE LATER - DEBUGGING / INFO PURPOSES
UINT32 graph::getVbackCapacity(vID v) const {
	return vList[v.id].getBackCapacity();
}

// Returns v's hypothetical position within its own back-neighbor list
UINT32 graph::getVmyBackPos(vID v) const {
	assert(v.id < vertexCount);
	return (vList[v.id].getMyBackPos());
}

// Sets v's hypothetical position within its own back-neighbor list
void graph::setVmyBackPos(vID v, UINT32 newBackPos) {
	assert(v.id < vertexCount);
	vList[v.id].setMyBackPos(newBackPos);
}

// Return v's neighbor list
Edge * graph::getVneighbors(vID v) {
	assert(v.id < vertexCount);
	return vList[v.id].getNeighbors();
}

// Return v's back-neighbor list
Edge * graph::getVbackNeighbors(vID v) {
	assert(v.id < vertexCount);
	return vList[v.id].getBackNeighbors();
}

// Get the number of vertices in the graph
UINT64 graph::getVertexCount() const {
    return vertexCount;
}

// Return |E|
UINT64 graph::getTotalDegree() const {
    return totalDegree;
}

// DELETE LATER - DEBUGGING / INFO PURPOSES
UINT64 graph::getUnusedSpace() const {
	UINT64 space = 0;
	for (int i = 0; i < vertexCount; i++) {
		// If the capacity does not equal the number of edges on the graph
		// if (vList[i].getCapacity() - vList[i].getOutDegree()) {
			// cerr << "Vertex: " << i << endl;
			// cerr << vList[i].getCapacity() << " - " << vList[i].getOutDegree() << endl << endl;
		// }
		space += (vList[i].getCapacity() - vList[i].getOutDegree());
	}
	return space;
}

// Print out the vertices that vertex v points to
void graph::printVneighbors(vID v) {
    cout << "Vertex " << v.id << "'s neighbors are: ";
    vList[v.id].printNeighbors();
}

// Print out the vertices that point to vertex v
void graph::printVbackNeighbors(vID v) {
    cout << "Vertex " << v.id << "'s backedge neighbors are: ";
    // vList[v.id].printBackNeighbors(numOfEntities);
    vList[v.id].printBackNeighbors();
}

// Print out the vertices that vertex v points to and those edges' weights
void graph::printVneighborsWeights(vID v) {
    // Print out the neighbors for vertex v
    cout << "Vertex " << v.id << "'s weights are: ";
    vList[v.id].printNeighborsWeights();
}

// Print out the vertices that point to vertex v and those edges' weights
void graph::printVbackNeighborsWeights(vID v) {
    // Print out the neighbors for vertex v
    cout << "Vertex " << v.id << "'s backedge weights are: ";
    vList[v.id].printBackNeighborsWeights();
}

// Print out the ID of each vertex that the indicated vertex points to and that points to it
void graph::printGraph() {
    for (UINT64 i = 0; i < vertexCount; i++) {
		if (vList[i].getInDegree() > 0 || vList[i].getOutDegree() > 0) {
			cout << "Vertex " << i << "'s neighbors are: ";
			vList[i].printNeighbors();
			cout << "Vertex " << i << "'s back-neighbors are: ";
			vList[i].printBackNeighbors();
			cout << endl;
		}
    }
}

// Print out the ID of each vertex that the indicated vertex points to and that points to it,
// as well as all the edges' weights
void graph::printGraphWeights() {
    for (UINT64 i = 0; i < vertexCount; i++) {
		if (vList[i].getInDegree() > 0 || vList[i].getOutDegree() > 0) {
			// cout << "Vertex " << i << "'s neighbors are: ";
			// vList[i].printNeighbors();
			cout << "Vertex " << i << "'s weights are:   ";
			vList[i].printNeighborsWeights();
			// cout << "Vertex " << i << "'s back-neighbors are: ";
			// vList[i].printBackNeighbors();
			cout << "Vertex " << i << "'s back-weights are:   ";
			vList[i].printBackNeighborsWeights();	
			cout << endl;
		}
    }
}

// DELETE LATER - DEBUGGING / INFO PURPOSES
void graph::printDebugGraph() {
	vID v;
	UINT32 outDegree, inDegree, capacity, backCapacity;
	UINT64 sumOfOutDegrees = 0, sumOfInDegrees = 0;
	for (v.id = 0; v.id < vertexCount - 1; v.id++) {
		outDegree = getVoutDegree(v);
		inDegree = getVinDegree(v);
		capacity = getVcapacity(v);
		backCapacity = getVbackCapacity(v);
		if (outDegree > 0 || inDegree > 0 || capacity > 0 || backCapacity > 0) {
			// fprintf(stderr, "Data for vertex %llu:\nOutDegree: %u\tCapacity: %u\n", v.id, outDegree, capacity);
			// printNeighbors();
			// fprintf(stderr, "InDegree: %u\tBackCapacity: %u\n", inDegree, backCapacity);
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
	
	// Remove the passed vertex from all of the vertex's neighbors' back neighbor lists
	for (UINT32 i = 0; i < vList[v.id].getOutDegree(); i++) {
		currID = clearList[i].vid.id;
		currEdge = getEdge(v, 1);
		vList[currID].delBackNeighbor(currEdge);
	}
	// Decrease total degree by the appropriate amount and clear the vertex's neighbor list (back-neighbor list persists)
	totalDegree = totalDegree - vList[v.id].getOutDegree();
	vList[v.id].clearNeighbors();
}

// Clear every vertex on the graph
void graph::clear() {
    for (UINT64 i = 0; i < vertexCount; i++) {
		vList[i].clearNeighbors();
		vList[i].clearBackNeighbors();
    }
    totalDegree = 0;
}
