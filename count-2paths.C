#include "parseCommandLine.h"
#include "graphIO.h"
#include "parallel.h"
#include "gettime.h"

//Takes as input a file in SNAP format
//(http://snap.stanford.edu/data/index.html).
//Currently assumes an undirected graph where each undirected edge
//appears once as a pair (u,v). It should not appear in both
//directions.
int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-batch <batchSize> -totalEdges <totalEdges> <input SNAP file>");
  char* iFile = P.getArgument(0);
  long batchSize = P.getOptionLongValue("-batch",10000);
  long totalEdges = P.getOptionLongValue("-totalEdges",1000000);
  edgeArray<uintT> G = readSNAP<uintT>(iFile);
  //minimum of totalEdges and number of edges in the graph
  totalEdges = min(totalEdges,(long)G.nonZeros);
  cout << "starting timer\n";
  timer t;
  t.start();
  long n = max(G.numRows,G.numCols); //number of vertices
  uintT* vertices = newA(uintT,n);
  for(long i=0;i<n;i++) vertices[i] = 0;
  long numBatches = 1+(totalEdges-1)/batchSize;
  cout << "number of batches = " << numBatches << endl;
  long count;
  for(long i=0;i<numBatches;i++) {
    for(long j=i*batchSize;j<min((long)(i+1)*batchSize,(long)totalEdges);j++) {
      vertices[G.E[j].u]++; vertices[G.E[j].v]++;
    }
    //output number of 2-paths. currently it is not efficient
    //because it loops through all vertices. can be made more
    //efficient by keeping track of vertices with degree > 1 and
    //looping over those only
    count = 0;
    for(long k=0;k<n;k++) if(vertices[k] > 1) count += vertices[k]*(vertices[k]-1);
    //cout << "count so far = " << count << endl;
  }
  cout << "total count = " << count << endl;
  t.reportTotal("total time");
}
