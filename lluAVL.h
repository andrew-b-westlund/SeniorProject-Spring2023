// File: p5.h (lluAVL)
// Author: Andrew Westlund
// Date: 3/17/2021
// Pupose: Build skills on C++ class creation and integration while
// implementing an AVL tree and understanding its structure and functions.
// Implements booleans, recursion, reference nod pointer variables,
// tree rotations, height calculations, and AVL tree functions and members.

#ifndef P4_H
#define P1_H

// Create typedefs for variables
typedef unsigned long long UINT64;
typedef unsigned int UINT32;
typedef unsigned short UINT16;
typedef unsigned char UINT8;

typedef long long INT64;
typedef int INT32;
typedef short INT16;
typedef char INT8;

class lluAVL {
	private:
		// Class provides nodes to the tree
		class node {
			public:
				node(UINT64 v); // Constructor
				UINT64 val;	 // Stores the value of the node
				UINT32 height;	 // Stores the height of the node
				UINT32 weight;	 // Stores how many times the node has been inserted
				node *left;	 // Points to the node to the left
				node *right; // Points to the node to the right
		};
		
		node *root;	// Points to the root node of the tree
		
		UINT32 tCount; // Stores how many values are in the tree

		// Preconditions: Node cannot be NULL
		// Postconditions: Returns the lowest value in the tree
		// to which p points
		UINT64 findMin(node *p) const;

		// Preconditions: None
		// Postconditions: Inserts val into the tree and returns
		// true if val isn't already in the tree; otherwise false is 
		// returned and the tree remains unchanged
		bool insertV(node *&p, UINT64 v);

		// Preconditions: None
		// Postconditions: Removes val from the tree and returns true
		// if successful; otherwise, false is returned and the tree
		// remains unchanged
		bool deleteV(node *&p, UINT64 v);

		// Preconditions: None
		// Postconditions: Returns true if val is found in the tree;
		// otherwise, returns false
		bool isIn(node *p, UINT64 v) const;

		// Preconditions: None
		// Postconditions: The tree's values are printed to the output
		// in ascending order along with their heights
		void printIt(node *p) const;

		// Preconditions: None
		// Postconditions: Removes all nodes from the tree, making it empty
		void clear(node *p);

		// Preconditions: None 
		// Postconditions: Balances the subtree of node p
		void bal(node *&p);

		// Preconditions: None
        // Postconditions: Returns the height of node p
        UINT32 height(node *p);

		// Preconditions: None
        // Postconditions: Calculates the height of node p
        UINT32 calcHeight(node *p);

		// Preconditions: None
        // Postconditions: Returns the higher value of the two
        UINT32 max(UINT32 a, UINT32 b);
		
		// Preconditions: None
        // Postconditions: Rotates the subtree of p to the right
        void rotateRight(node *&p);

		// Preconditions: None
        // Postconditions: Rotates the subtree p to the left
        void rotateLeft(node *&p);

	public:
		// Constructor
		lluAVL();

		// Destructor
		~lluAVL();	

		// Preconditions: None
		// Postconditions: Inserts val into the tree and returns
		// true if val isn't already in the tree; otherwise false is 
		// returned and the tree remains unchanged
		bool insertV(UINT64 v);

		// Preconditions: None
		// Postconditions: Removes val from the tree and returns true
		// if successful; otherwise, false is returned and the tree
		// remains unchanged
		bool deleteV(UINT64 v);

		// Preconditions: None
		// Postconditions: Returns true if val is found in the tree;
		// otherwise, returns false
		bool isIn(UINT64 v) const;

		// Preconditions: None
		// Postconditions: The tree's values are printed to the output
		// in ascending order along with their heights
		void printIt() const;

		// Preconditions: None
		// Postconditions: Returns the number of values in the tree
		UINT32 count() const;

		// Preconditions: None
		// Postconditions: Removes all nodes from the tree, making it empty
		void clear();
		
		friend class bvAVL;
};
#endif
