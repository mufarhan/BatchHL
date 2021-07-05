# BatchHL

This is the implementation of the paper "BatchHL-Answering Distance Queries on Batch-Dynamic Networks at scale", currently under review in SIGMOD. 

The format of the dataset text file is: 
vertex_u deg_u v1 ... vn, where v1 to vn are neighbors of u. Note that the vertex id are from 0 to (V-1), where V is the number of vertices. There are no self loops in the graph, i.e., no edge from any vertex to itself. 

To see the accepted format for datasets, batch updates and query pairs, you may refer to the Sample folder. After the test inputs are ready, please us the following commands to test BatchHL.

=====================================

1 - Compile source files using the following command:
g++ -O3 -std=c++11 -pthread BatchHL.cpp -o run

=====================================

2 - Construct Labelling:
./run update_labelling @1 @2 @3
@1: name of the dataset
@2: number of landmarks
@3: file to store labelling

Example:
./index construct_labelling skitter.txt 20 skitter

=====================================

3 - Update Labelling:
./run update_labelling @1 @2 @3 @4 @5 @6
@1: name of the dataset
@2: number of landmarks
@3: file to load the labelling from
@4: batch file containing updates
@5: method parameter (0 to run BHL+ or 1 to run BHL)
@6: parallelism parameter (0/1)

Example:
./index update_labelling skitter.txt 20 skitter batch.txt 0 0

=====================================

4 - Perform distance queries
./run query-dis @1 @2 @3 @4 @5
@1: name of the dataset
@2: number of landmarks
@3: file to load the labelling from
@4: file containing query pairs
@5: file to write query results

Example:
./index query_labelling skitter.txt 20 skitter query_pairs.txt results.txt









