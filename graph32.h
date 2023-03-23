// File: graph.h
// Author: Dr. Stephen Wheat, Austin Lehman, Samuel Udall, Michael Vandusen, Andrew Westlund
// Date: 1/31/22
// Pupose: This program will recieve a graph and will create a list of all the vertices
// and will list all the vertices adjacent to each vertex.
// The both the list and the lists of adjacent vertices will be dynamically allocated
// and will be sorted least to greatest. 

#ifndef __GRAPH_H
#define __GRAPH_H

// Create typedefs for variables
typedef unsigned long long UINT64;
typedef unsigned int UINT32;
typedef unsigned short UINT16;
typedef unsigned char UINT8;

typedef long long INT64;
typedef int INT32;
typedef short INT16;
typedef char INT8;

// -------------------------------------- VID -----------------------------------------

class vID {
    public:
	UINT64 id;
};

// -------------------------------------- Edge -----------------------------------------

class Edge {
	public:
	vID vid;
	UINT32 weight;
	// padding to round out memory to multiples of 64
	UINT32 pad;
};

// -------------------------------------- Neighbors -----------------------------------------

class Neighbors {
	private:
	Edge *nList;		// An array of edge objects (neighbor list for a given index).
	UINT32 degree;		// The count of egdes.
	UINT32 capacity;	// The capacity of the vertex array.
	UINT32 myPos;		// Identifies where this vertex should be in the nList if it were there.
	UINT32 pad;			// Pads to even multiple of 64
	
	// Preconditions: vID v is not already on the adjacency list
	// Postconditions: Will insert the vID to the adjacency list in ascending order
	void addNeighbor0(Edge e);
	
	public:
	// Constructor
	Neighbors();
	
	// Destructor
	~Neighbors();
	
	// WILL COMPLETELY DELETE OLD NEIGHBOR
	// Preconditions: None
	// Postconditions: Will set Neighbors' nList to be the size of the given capacity, will set Neighbors' capacity
	// to the passed capacity, and will set degree and myPos to 0; essentially a constructor for a specific size.
	void setNeighbors(UINT32 newCapacity);
	
	// CONSIDER MAKING INDEX OPTIONAL (OVERLOAD THE FUNCTION)
	// Preconditions: None. 
	// Postconditions: Will return whether v is in the nList.
	// Index will be turned into the index of v in the nList, 0 if it is not found.
	bool isNeighbor(vID v, UINT32 &index);
	
	// Preconditions: None
	// Postconditions: Will check if the vertex is already on the vertex's neighbor list;
	// if the vertex is not already on it, it then checks if there is free space on the list, if there is no free space,
	// it will increase the list and add the edge. Returns whether the neighbor was added.
	bool addNeighbor(Edge e, INT16 allocation);
	
	// Preconditions: None
	// Postconditions: Appends a value to the end of the list; will reallocate if needed,
	// but createV should allocate the amount of space necessary for it to not be needed;
	// no safeguard from adding an edge that is already there.
	void addNeighbor2(Edge e, INT16 allocation);
	
	// CONSIDER ADDING AN ALLOCATION PARAM TO DECREASE CAPACITY BY WHEN DEGREE + ALLOCATION <= CAPACITY
	// Preconditions: None
	// Postconditions: Will return true if the vertex is on the vertex's adjacency list and then 
	// removes it; if it is not on the list, will return false.
	bool delNeighbor(Edge e);
	
	// Preconditions: None.
	// Postconditions: If v is not on nList return 0;
	// else, return the weight of the edge with v.
	UINT32 getWeight(vID v);
	
	// Preconditions: None.
	// Postconditions: If an edge with v is on the nList, change the weight of the edge.
	void setWeight(vID v, UINT32 weight);
	
	// Preconditions: None.
	// Postconditions: Increments myPos.
	void incPos();
	
	// Preconditions: None.
	// Postconditions: Decrements myPos.
	void decPos();
	
	// Preconditions: None
	// Postconditions: Returns the nList's size.
	UINT32 getDegree() const;
	
	// Preconditions: None
    // Postconditions: Returns the capacity of the nList.
	UINT32 getCapacity() const;
	
	// Preconditions: None
    // Postconditions: Returns the vertex's position on the nList.
	UINT32 getMyPos() const;
	
	// Preconditions: None
	// Postconditions: Will set the Neighbors' myPos to the passed position
	void setMyPos(UINT32 newPos);
	
	// Preconditions: None
	// Postconditions: Returns the nList.
	Edge * getNList() const;
	
	// Preconditions: None.
	// Postconditions: Will print the list of adjacent vertices.
	void printNeighbors() const;
	
	// Preconditions: None.
	// Postconditions: Will print the list of adjacent vertices weights.
	void printNeighborsWeights() const;
	
	// Preconditions: None.
	// Postconditions: Clears the nList and resets the stored values.
	void clearNeighbors();
};

// -------------------------------------- Vertex -----------------------------------------

class vertex {
    private:
	Neighbors forward;
	Neighbors backward;
	// Neighbors *forward;
	// Neighbors *backward;
	
    public:
	// Constructor
	vertex();
	
	// Destructor
	~vertex();
	
	// 
	// 
	void setVertex(UINT32 capacity, UINT32 backCapacity);
	
	// Preconditions: None.
	// Postconditions: Will return true if this vertex has an edge to v and will set index to v's position;
	// else, false is returned.
	bool isNeighbor(vID v, UINT32 &index);
	
	// Preconditions: None.
	// Postconditions: Will return true if this vertex has a back-edge to v and will set index to v's position;
	// else, false is returned.
	bool isBackNeighbor(vID v, UINT32 &index);

	// Preconditions: None
	// Postconditions: Will check if e's vertex is already on the vertex's edge list;
	// if the vertex is not already on it, it then checks if there is free space on the list; if there is no free space,
	// it will increase the capacity by alloc, then add e as a neighbor. If e is already there, return false.
	bool addNeighbor(Edge e, INT16 alloc);
	
	// Preconditions: None.
	// Postconditions: If enough space is available, the edge will be appended to the vertex's neighbor list in constant time;
	// else, it will be appended to the list in linear time.
	void addNeighbor2(Edge e, INT16 allocation);
	
	// Preconditions: None
	// Postconditions: Will check if e's vertex is already on the vertex's back-edge list;
	// if the vertex is not already on it, it then checks if there is free space on the list; if there is no free space,
	// it will increase the capacity by alloc, then add e as a back-neighbor. If e is already there, return false.
	bool addBackNeighbor(Edge e, INT16 a);
	
	// Preconditions: None.
	// Postconditions: If enough space is available, the edge will be appended to the vertex's back-neighbor list in constant time;
	// else, it will be appended to the list in linear time.
	void addBackNeighbor2(Edge e, INT16 allocation);
	
	// Preconditions: None
	// Postconditions: Will return true if e's vertex is on the vertex's neighbor list and then 
	// removes it; if it is not on the list, will return false.
	bool delNeighbor(Edge e);
	
	// Preconditions: None
	// Postconditions: Will return true if e's vertex is on the vertex's backedge neighbor list and then 
	// removes it; if it is not on the list, will return false.
	bool delBackNeighbor(Edge e);
	
	// Preconditions: None.
	// Postconditions: If v is not on the neighbor list, return 0;
	// else, return the weight of the edge with v.
	UINT32 getWeight(vID v);
	
	// Preconditions: None.
	// Postconditions: If v is a neighbor, return true and change the weight of the edge;
	// else, return false.
	void setWeight(vID v, UINT32 weight);
	
	// Preconditions: None.
	// Postconditions: If v is not on back-neighbor list, return 0;
	// else, return the weight of the edge with v.
	UINT32 getBackWeight(vID v);
	
	// Preconditions: None.
	// Postconditions: If v is a back-neighbor, return true and change the weight of the edge;
	// else, return false.
	void setBackWeight(vID v, UINT32 weight);
	
	// Preconditions: None
	// Postconditions: Returns the vertex's myPos
	UINT32 getMyPos() const;
	
	// Preconditions: None
	// Postconditions: Will set the vertex's myPos to the passed position
	void setMyPos(UINT32 newPos);
	
	// Preconditions: None.
	// Postconditions: Increments myPos.
	void incPos();
	
	// Preconditions: None.
	// Postconditions: Decrements myPos.
	void decPos();
	
	// Preconditions: None
	// Postconditions: Returns the vertex's myBackPos
	UINT32 getMyBackPos() const;
	
	// Preconditions: None
	// Postconditions: Will set the vertex's myBackPos to the passed position
	void setMyBackPos(UINT32 newBackPos);
	
	// Preconditions: None.
	// Postconditions: Increments myBackPos.
	void incBackPos();
	
	// Preconditions: None.
	// Postconditions: Decrements myBackPos.
	void decBackPos();
	
	// Preconditions: None
	// Postconditions: Returns the number of this vertex's neighbors.
	UINT32 getOutDegree() const;
	
	// Preconditions: None
	// Postconditions: Returns the number of vertices that this vertex is a neighbor to.
	UINT32 getInDegree() const;
	
	// DELETE LATER - DEBUGGING / INFO PURPOSES
	// Preconditions: None
    // Postconditions: Returns the capacity of the neighbor list.
	UINT32 getCapacity() const;
	
	// DELETE LATER - DEBUGGING / INFO PURPOSES
	// Preconditions: None
    // Postconditions: Returns the capacity of the back-neighbor list.
	UINT32 getBackCapacity() const;
	
	// Preconditions: None.
	// Postconditions: Returns the list of neighbors.
	Edge * getNeighbors() const;
	
	// Preconditions: None.
	// Postconditions: Returns the list of back-neighbors.
	Edge * getBackNeighbors() const;
	
	// Preconditions: None.
	// Postconditions: Will print the list of adjacent vertices.
	void printNeighbors() const;
	
	// Preconditions: None.
	// Postconditions: Will print the list of the neighbor's weights.
	void printNeighborsWeights() const;
	
	// Preconditions: None.
	// Postconditions: Will print the list of backedge-adjacent vertices 
	void printBackNeighbors() const;
	
	// Preconditions: None.
	// Postconditions: Will print the list of the back-neighbor's weights.
	void printBackNeighborsWeights() const;
	
	// Preconditions: None.
	// Postconditions: Clears the list of neighbors and resets relevant values.
	void clearNeighbors();
	
	// Preconditions: None.
	// Postconditions: Clears the list of back-neighbors and resets relevant values.
	void clearBackNeighbors();
};

// -------------------------------------- Graph -----------------------------------------

class graph {
    private:
	vertex *vList;		// An array of vertex objects.
	UINT64 vertexCount;	// Stores how many vertices there are in the graph
	UINT64 totalDegree;	// Stores how many edges there are in the graph
	INT16 allocation;   // The value of how much the vertex capacity increments by
	UINT16 pad1;		// padding to round out memory to multiples of 8
	UINT32 pad2;
	
	// Preconditions: None.
	// Postconditions: Returns an edge that has vID v and wieght 'weight'
	Edge getEdge(vID v, UINT32 weight);
    
	public:
	// Constructor
	graph(UINT64 numberOfVertices);
	
	// Destructor
	~graph();
	
	// Preconditions: None.
	// Postconditions: Will set allocation to a if a is between 1 and 16; if under 1, will be set to 1
	// and if over 16 will be set to 16
	void setAllocation(INT16 a);
	
	// Preconditions: v.id is less than vNum
	// Postconditions: neighborLists of size capacity will have a one-time instantiation (reduces linear copying over)
	void createV(vID v, UINT32 capacity, UINT32 backCapacity);
	
	// Preconditions: vertex(u) and vertex(v) must be existing vertices on the graph.
	// Postconditions: If there is an edge going from vertex(u) to vertex(v), return true; otherwise, return false.
	bool isEdge(vID u, vID v, UINT32 weight = 1);
	
	// Preconditions: vertex(u) and vertex(v) must be existing vertices on the graph.
	// Postconditions: If there is a backedge going from vertex(u) to vertex(v), return true; otherwise, return false.
	bool isBackEdge(vID u, vID v, UINT32 weight = 1);

	// Preconditions: vertex(u) and vertex(v) must be existing vertices
	// Postconditions: If inc is false, then if u has an edge to v, false will be returned, but if 
	// u does not have an edge to v, then an edge with the passed weight from u to v will be added,
	// a similar back-edge from v to u will be added, and true will be returned.
	// If inc is true, then even if an edge already exists, the edge will increase by the passed weight
	// and true will be returned.
	bool addEdge(vID u, vID v, UINT32 weight = 1, bool inc = false);
	
	// Preconditions: vertex(u) and vertex(v) must be existing vertices
	// Postconditions: Essentially the instant addEdge; an edge from u to v will be added quickly from
	// vertex(u) to vertex(v) without seeing if it already exists and without finding the index where it
	// belongs; throws an edge right at the end of the u's neighborList and v's backNeighborList in constant
	// time, unless more memory must be allocated.
	// Most efficient if used after createV function sets the size of u and v's (back)capacities and (back)neighborLists.
	void addEdge2(vID u, vID v, UINT32 weight = 1);
	
	// CONSIDER ADDING bool dec = false AS A PARAM
	// Preconditions: v(u) and v(v) must be existing vertices
	// Postconditions: If the edge already exists, the edge that points from vertex(u) to 
	// vertex(v) is removed, the , and true is returned.
	// Otherwise, return false.
	bool delEdge(vID u, vID v, UINT32 weight = 1);
	
	// CONSIDER ASKING ABOUT ADDING A SETEDGE(u, v) AND GETEDGE(u, v) FUNCTION
	
	//
	//
	void dumpGraph();

	// Preconditions: None.
	// Postconditions: Returns the vertex outdegree
	UINT32 getVoutDegree(vID v) const;
	
	// Preconditions: None.
	// Postconditions: Returns the vertex indegree
	UINT32 getVinDegree(vID v) const;
	
	// DELETE LATER - DEBUGGING / INFO PURPOSES
	UINT32 getVcapacity(vID v) const;
	
	// DELETE LATER - DEBUGGING / INFO PURPOSES
	UINT32 getVbackCapacity(vID v) const;
	
	// Preconditions: The passed vertex's id is withing the graph's indexing bounds.
	// Postconditions: Returns the passed vertex's myPos.
	UINT32 getVmyPos(vID v) const;
	
	// Preconditions: The passed vertex's id is withing the graph's indexing bounds.
	// Postconditions: The passed vertex's myPos is set to the passed position.
	void setVmyPos(vID v, UINT32 newPos);
	
	// Preconditions: The passed vertex's id is withing the graph's indexing bounds.
	// Postconditions: Returns the passed vertex's myBackPos.
	UINT32 getVmyBackPos(vID v) const;
	
	// Preconditions: The passed vertex's id is withing the graph's indexing bounds.
	// Postconditions: The passed vertex's myBackPos is set to the passed position.
	void setVmyBackPos(vID v, UINT32 newBackPos);

	// Preconditions: v must be the id of a vertex on the graph.
	// Postconditions: The neighbor list of vertex(v) is returned	
	Edge * getVneighbors(vID v);
	
	// Preconditions: v must be the id of a vertex on the graph.
	// Postconditions: The backedge neighbor list of vertex(v) is returned
	Edge * getVbackNeighbors(vID v);

	// Preconditions: None
	// Postconditions: Returns |V|
	UINT64 getVertexCount() const;

	// Preconditions: None
	// Postconditions: Returns |E|
	UINT64 getTotalDegree() const;
	
	// DELETE LATER - DEBUGGING / INFO PURPOSES
	UINT64 getUnusedSpace() const;
	
	// Preconditions: None.
	// Postconditions: Will print the list of all vertices and their adjacent vertices.
	void printGraph();
	
	// Preconditions: vNum must be the id of a vertex on the graph.
	// Postconditions: The vertex's neighbors are printed.
	void printVneighbors(vID v);
	
	// Preconditions: v must have an id of a vertex on the graph.
	// Postconditions: The vertex's backedge neighbors are printed.
	void printVbackNeighbors(vID v);
	
	// Preconditions: vNum must be the id of a vertex on the graph.
	// Postconditions: The vertex's neighbors weights are printed.
	void printVneighborsWeights(vID v);
	
	// Preconditions: vNum must be the id of a vertex on the graph.
	// Postconditions: The vertex's neighbors weights are printed.
	void printVbackNeighborsWeights(vID v);

	// Preconditions: None.
	// Postconditions: clears the graph and resets all values of all nLists.
	void clear();
	
	// Preconditions: None.
	// Postconditions: Clears the nList and resets the stored values of vID v.
    void clearNlist(vID v);
	
	// DELETE LATER - DEBUGGING / INFO PURPOSES
	void printDebugGraph();
};
#endif
