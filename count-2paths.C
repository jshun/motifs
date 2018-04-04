#include "parseCommandLine.h"
#include "graphIO.h"
#include "parallel.h"
#include "gettime.h"

//Takes as input a file in SNAP format
//(http://snap.stanford.edu/data/index.html).
//Currently assumes a directed graph where each directed edge
//appears once as a pair (u,v).
int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-batch <batchSize> -totalEdges <totalEdges> <input SNAP file>");
  char* iFile = P.getArgument(0);
  long batchSize = P.getOptionLongValue("-batch",10000);
  long totalEdges = P.getOptionLongValue("-totalEdges",1000000);
  edgeArray<uintT> G = readSNAP<uintT>(iFile);

  cout << "number of edges read = " << (long) G.nonZeros << endl;
  //minimum of totalEdges and number of edges in the graph
  totalEdges = min(totalEdges,(long)G.nonZeros);
  cout << "starting timer\n";
  timer t;
  t.start();
  long n = max(G.numRows,G.numCols); //number of vertices

  // keep track of both out and in-degrees for each vertex.
  uintT* out_degrees = newA(uintT,n);
  uintT* in_degrees = newA(uintT,n);
  for(long i=0;i<n;i++) {
    out_degrees[i] = 0;
    in_degrees[i] = 0;
  }
  long numBatches = 1+(totalEdges-1)/batchSize;
  long count;

  for(long i=0;i<numBatches;i++) {
    for(long j=i*batchSize;j<min((long)(i+1)*batchSize,(long)totalEdges);j++) {
      out_degrees[G.E[j].u]++;
      in_degrees[G.E[j].v]++;
    }
    //output number of 2-paths. currently it is not efficient
    //because it loops through all vertices. can be made more
    //efficient by keeping track of vertices with degree > 1 and
    //looping over those only
    count = 0;
    for(long k=0;k<n;k++) {
      if(in_degrees[k] * out_degrees[k] >= 1) {
        count += in_degrees[k] * out_degrees[k];
      }
    }
  }
  cout << "total count = " << count << endl;
  t.reportTotal("total time");
}
