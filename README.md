# SeniorProject-Spring2023
Andrew's Graph Building, Squaring, Incrementing, and Analyzing Project in Collaboration with BioVista

This custom graph framework that Andrew has been involved with since his research project in Spring 2022 is meant to create large, sparse graphs. 

Biovista’s Project Prodigy is an AI engine that processes a large graph of published papers to develop diagnoses and proposed treatment plans for complex illnesses.
However, BioVista's current framework cannot incrementally update the data, but it must instead recreate the graph from scratch with all the new data included, a process that takes approximately ten days.

Having the ability to rapidly add incremental data would allow a potential diagnosis or proposed treatment plan to incorporate the latest science, which is what Andrew's senior project aims to accomplish.

Input: When executing the program, a -f flag followed by the input file's location should be included (babyTest.txt). This file should be formatted as a repeated sequence of PaperID PublicationYear EntityTypeID EntityID. 
If the file is in binary format, the -b flag should be included. 

A file with new edges can be added by including the -i flag followed by the file. The edges should be in the same format as the first file's: PaperID PublicationYear EntityTypeID EntityID. 
If this file is binary, you can add the -n flag.

There are four different primitives that are described in the final report. Each primitive's flag is meant to be followed by how many times that primitive is to be performed.

-e is the expand's primitive flag.

-c is the connect's primitive flag.

-r is the bridge's primitive flag.

-p is the bibliography's primitive flag.

Each primitive will prompt the user for further input.

Output: The output of the primitives, the outdegree of the graph, and the time duration of each stage is reported.

Serial Compilation: g++ -O2 graphMain.cpp graph.cpp lluAVL.cpp -o graph

Multithreaded Compilation: g++ -O2 -fopenmp -D PAR graphMain.cpp graph.cpp lluAVL.cpp -o graph

Execution: ./graph -f inputFileLocation (-b) (-i incrementalFileLocation (-n) ) (-e #OfTimes) (-c #OfTimes) (-r #OfTimes) (-p #OfTimes)
