# BatchHL

This is an implementation of BatchHL method which can reflect graph changes in batches to answer shortest-path distance queries efficiently over very large batch-dynamic networks.

To see the accepted format for datasets, batch updates and query pairs, you may refer to the Sample folder. After the test inputs are ready, please us the following commands to test BatchHL.

=====================================

#1 - Compile source files using the following command:<br/>
g++ -O3 -std=c++11 main.cpp -o run

=====================================

#2 - Construct Labelling:<br/>
./run construct_labelling @1 @2 @3<br/>
@1: name of the dataset<br/>
@2: number of landmarks<br/>
@3: file to store labelling

Example:<br/>
./run construct_labelling graph.txt 20 graph_labelling

=====================================

#3 - Update Labelling:<br/>
./run update_labelling @1 @2 @3 @4 @5<br/>
@1: name of the dataset<br/>
@2: number of landmarks<br/>
@3: file to load the labelling from<br/>
@4: batch file containing updates<br/>
@5: method parameter (0 to run BHL+ or 1 to run BHL)<br/>

Example:<br/>
./run update_labelling graph.txt 20 graph_labelling batch.txt 0

=====================================

#4 - Perform distance queries<br/>
./run query-dis @1 @2 @3 @4 @5<br/>
@1: name of the dataset<br/>
@2: number of landmarks<br/>
@3: file to load the labelling from<br/>
@4: file containing query pairs<br/>
@5: file to write query results<br/>

Example:<br/>
./run query_labelling graph.txt 20 graph_labelling query_pairs.txt query_results.txt

# References
* Muhammad Farhan, Qing Wang and Henning Koehler, **[BatchHL: Answering Distance Queries on Batch-Dynamic Networks at Scale](https://dl.acm.org/doi/abs/10.1145/3514221.3517883)**. SIGMOD 2022.


