# SeniorProject-Spring2023
Andrew's Graph Creating and Squaring Project in Collaboration with BioVista

This custom graph framework that Andrew has been involved with since his research project in Spring 2022 is meant to create large, sparse graphs. 

Biovista’s Project Prodigy is an AI engine that processes a large graph of published papers to develop diagnoses and proposed treatment plans for complex illnesses.
However, BioVista's current framework cannot incrementally update the data, but it must instead recreate the graph from scratch with all the new data included, a process that takes approximately three weeks.
Having the ability to rapidly add incremental data would allow a potential diagnosis or proposed treatment plan to incorporate the latest science, which is what Andrew's senior project aims to accomplish.

Input: When executing the program, a -f flag followed by the input file's location should be included. This file should be formatted as a repeated sequence of PaperID PublicationYear EntityTypeID EntityID.

Output: Primitives will be added later for the use of the calculations of BioVista.

Updates/Blockers: Andrew and his advisor, Dr. Stephen Wheat, have decided that for this implementation of the graph framework, the insertion of an edge takes too long since some entities are referenced by more than 2^24 papers. Each insertion has a time complexity of O(n), where n is the number of edges that the vertex being added to has. Andrew is thus working on using an AVL tree to insert edges with a time of O(lg(n)), where each edge is put into the graph with a complexity of O(1).

Compilation: g++ graphTreem.cpp graph32.cpp lluAVL.cpp -O2 -o graph

Execution: ./graph -f inputFileLocation
