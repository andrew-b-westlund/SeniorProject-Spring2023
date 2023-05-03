// File: p5.cpp (lluAVL)
// Author: Andrew Westlund
// Date: 3/17/2021
// Pupose: Build skills on C++ class creation and integration while
// implementing an AVL tree and understanding its structure and functions.
// Implements booleans, recursion, reference nod pointer variables,
// tree rotations, height calculations, and AVL tree functions and members.

#include <iostream>    // Provides cout
#include <iomanip>     // Provides setw function for setting output width
#include <cstdlib>     // Provides EXIT_SUCCESS
#include <cassert>     // Provides assert function
#include "lluAVL.h"        // Provides header file containing definitions

using namespace std;   // Allows all standard library items to be used

lluAVL::node::node(UINT64 v) {
	val = v;		// Stores the value of the node
	height = 1;		// Stores the height of the node
	weight = 1;
	left = NULL;	// Points to the node to the left
	right = NULL;	// Points to the node to the right
}

UINT64 lluAVL::findMin(node *p) const {
	UINT64 rc = p->val;  // Return the current node's data
	
	if (p->left) {	// If the node has a smaller child
		// Call the function of the smaller child
		rc = findMin(p->left);
	}
	return rc;
}

void lluAVL::bal(node *&p) {
	if (p) {	// If p is not NULL
		// If p's left child's height is two bigger than p's right child's height
		if (height(p->left) - height(p->right) == 2) {

			// If p's left child's left child is taller than p's left child's
			// right child
			if (height((p->left)->left) > height((p->left)->right)) {
				// If it's left-left heavy
				rotateRight(p);	// Rotate the tree of p to the right

			// If p's left child's left child is shorter than or equal to p's
			// left child's right child
			} else {
				// If it's left-right heavy
				// Rotate p's left child's tree to the left
				rotateLeft(p->left);
				rotateRight(p);	// Rotate the changed tree of p to the right
			}

		// If p's left child's height is two smaller than p's right child's height
		} else if (height(p->left) - height(p->right) == -2) {

			// If p's right child's left child is taller than or equal to p's
			// right child's right child
			if (height((p->right)->left) >= height((p->right)->right)) {	
				// If it's right-left heavy
				// Rotate p's right node to to the right
				rotateRight(p->right);
				rotateLeft(p);	// Rotate the changed tree of p to the left

			// If p's right child's left child is shorter than p's right child's
			// right child
			} else {
				// If it's right-right heavy
				rotateLeft(p);	// Rotate the tree of p to the left
			}
		}
		p->height = calcHeight(p); // Reset the height of p
	}
}

UINT32 lluAVL::height(node *p) {
	UINT32 rc = 0;	// Return 0 if p is NULL

	if (p) {
		rc = p->height;	// If p exists, return the height of p
	}

	return rc;
}

UINT32 lluAVL::calcHeight(node *p) {
	UINT32 rc = 0;	// Return 0 if p is NULL

	if (p) {	// If p exists
		// Return the greater height of p's children plus 1
		rc = max(height(p->left), height(p->right)) + 1;
	}

	return rc;
}

UINT32 lluAVL::max(UINT32 a, UINT32 b) {
	// If a is greater than b, return a; otherwise, return b
	return (a > b ? a : b);
}

void lluAVL::rotateRight(node *&p) {
	node *p2 = p->left;			// Set p2 to be p's left child
	// Set p's left child to be its left child's right child
	p->left = p2->right; 
	// Set p to be its former left child's right child
	p2->right = p;
	p->height = calcHeight(p);	// Recalculate p's height
	p2->height = calcHeight(p2);	// Recalculate p2's height
	p = p2;	// Set p2 to be in the position its parent was formerly in
}

void lluAVL::rotateLeft(node *&p) {
	node *p2 = p->right;			// Set p2 to be p's right child
	// Set p's right child to be its right child's left child
	p->right = p2->left;
	// Set p to be its former right child's left child
	p2->left = p;
	p->height = calcHeight(p);	// Recalculate p's height
	p2->height = calcHeight(p2);	// Recalculate p2's height
	p =  p2; // Set p2 to be in the position its parent was formerly in
}

lluAVL::lluAVL() {
	root = NULL;  // Stores the root node of the lluAVL
	tCount = 0;	  // Stores the number of nodes or values in the lluAVL
}

lluAVL::~lluAVL() {
	clear();	  // Calls clear function to get rid of the tree
}

bool lluAVL::insertV(UINT64 v) {
	// Calls private recursive insertV function
	return insertV(root, v);
}

bool lluAVL::insertV(node *&p, UINT64 v) {
	bool rc;	// Stores return code
	
	if (p) {	// If node isn't NULL
	
		// If the current node's value if too small
		if (p->val < v) {
			// Call the function of the bigger child
			rc = insertV(p->right, v);
		
		// If the current node's value is too big
		} else if (p->val > v) {
			// Call the function of the smaller child
			rc = insertV(p->left, v);
		
		// If a node has the value
		} else {
			// New for bvGraph: if the node already exists, increment it
			p->weight++;
			// Do nothing and return false
			rc = false;
		}
	
	// If the node with value v doesn't exist yet
	} else {
		// Create a new node with value v in the right spot
		p = new node(v);
		tCount++;		 // Increment the count of the tree
		rc = true;		 // Return true
	}
	
	// Makes the tree balanced if not already balanced and fixes heights
	bal(p);
	return rc;
}

bool lluAVL::deleteV(UINT64 v) {
	// Calls private recursive deleteV function
	return deleteV(root, v); 
}

bool lluAVL::deleteV(node *&p, UINT64 v) {
	bool rc = false;
	node *p2;	// Points to a temporary node
	UINT64 v2;		// Stores a temporary int value
	
	if (p) {	// If the current node exists
		
		// If the current node is too small
		if (p->val < v) {
			// Call the function of the bigger child
			rc = deleteV(p->right, v);
		
		// If the current node is too big
		} else if (p->val > v) {
			// Call the function of the smaller child
			rc = deleteV(p->left, v);
		
		// If the current node matches int v
		} else {
			
			// If there are no children
			if (!(p->left || p->right)) {
				delete p;
				p = NULL;	// Set p to NULL
				tCount--;	// Decrement the size of count
			
			// Else if there is only one child
			} else if (!(p->left && p->right) && 
			  ( p->left || p->right)) {
				
				// If the one child is smaller
				if (p->left) {
					// Set p2 to smaller child
					p2 = p->left;
				
				// If the one child is bigger
				} else {
					// Set p2 to bigger child
					p2 = p->right;
				}
				delete p;
				p = p2;	// Replace node with child
				tCount--;	// Decrement the size of count
			
			// Else if two children
			} else {
				// Set v2 to the smallest value of the right subtree
				v2 = findMin(p->right);
				// Change node's value to smallest right subtree number
				p->val = v2;
				// Delete node with the value that replaced old node 
				deleteV(p->right, v2);
			}
			// Return true if value was deleted
			rc = true;
		}
	}
	
	bal(p);	// Makes the tree balanced if not already
	return rc;
}

bool lluAVL::isIn(UINT64 v) const {
	// Calls private recursive isIn function
	return isIn(root, v);
}

bool lluAVL::isIn(node *p, UINT64 v) const {
	bool rc = false;	// If v is not found, return false
	
	if (p) {	// If p is not NULL
		
		// If the node matches v
		if (p->val == v) {
			rc = true;	// Return true
		
		// Otherwise if node's value is too big
		} else if (p->val > v) {
			// Recurse using smaller child
			rc = isIn(p->left, v);
		
		// Otherwise if node's value is too big
		} else {
			// Recurse using bigger child
			rc = isIn(p->right, v);
		}	
	}
	return rc;
}

void lluAVL::printIt() const {
	// Calls private recursive printIt function
		printIt(root);
}

void lluAVL::printIt(node *p) const {
	if (p) {	// If node is not NULL
		// Print smaller child nodes
		printIt(p->left);
		// Print out current node's value
		cout << p->val << " : height = " << p->height << " - weight = " << p->weight << endl;
		// Print bigger child nodes
		printIt(p->right);
	}
}

UINT32 lluAVL::count() const {
	return tCount;	// Returns the number of nodes or values in the lluAVL
}

void lluAVL::clear() {
	// Calls private recursive clear function
	clear(root);
	// Resets class variables
	root = NULL;
	tCount = 0;
}

void lluAVL::clear(node *p) {
	if (p) {	// If node is not NULL
		// Clear smaller child nodes
		clear(p->left);
		// Clear bigger child nodes
		clear(p->right);
		// Delete current node
		delete p;
	}
}
