#include "parseCommandLine.h"
#include "graphIO.h"
#include "parallel.h"
#include "gettime.h"
#include <map>
#include <list>

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
  uintT* outDegrees = newA(uintT,n);
  uintT* inDegrees = newA(uintT,n);
  for(long i=0;i<n;i++) {
    outDegrees[i] = 0;
    inDegrees[i] = 0;
  }

  // keep track of adjacency list for each vertex.
  // NOTE: currently using vanilla STL.
  map<uintT, list<uintT>> inEdges = {};
  map<uintT, list<uintT>> outEdges = {};

  long numBatches = 1+(totalEdges-1)/batchSize;
  long count;
  long listCount;

  for(long i=0;i<numBatches;i++) {
    for(long j=i*batchSize;j<min((long)(i+1)*batchSize,(long)totalEdges);j++) {
      auto src = G.E[j].u;
      auto dst = G.E[j].v;

      outDegrees[src]++;
      inDegrees[dst]++;

      // outgoing edges
      if (!outEdges.insert(make_pair(src, list<uintT>{dst})).second) {
        outEdges[src].push_back(dst);
        //cout << "outEdges[" << src << "] size so far = " << outEdges[src].size() << endl;
      } else {
        outEdges[src] = {dst};
      }

      // incoming edges
      if (!inEdges.insert(make_pair(dst, list<uintT>{src})).second) {
        inEdges[dst].push_back(src);
        //cout << "inEdges[" << dst << "] size so far = " << inEdges[dst].size() << endl;
      } else {
        inEdges[dst] = {src};
      }

    }
    //output number of 2-paths. currently it is not efficient
    //because it loops through all vertices. can be made more
    //efficient by keeping track of vertices with degree > 1 and
    //looping over those only
    count = 0;
    listCount = 0;
    for(long k=0;k<n;k++) {

      // O(1) count
      if(inDegrees[k] * outDegrees[k] >= 1) {
        count += inDegrees[k] * outDegrees[k];
      }

      // Actual listing of edges
      if (inEdges.find(k) != inEdges.end() && outEdges.find(k) != outEdges.end()) {
        for (auto itSrc = inEdges[k].begin(); itSrc != inEdges[k].end(); ++itSrc) {
          for (auto itDst = outEdges[k].begin(); itDst != outEdges[k].end(); ++itDst) {
            //cout << "2-hop edge: " << *itSrc << "->" << *itDst << endl;
            // TODO: finish this up by storing them on an edgeArray
            listCount++;
          }
        }
      }
    }
  }
  cout << "total count = " << count << endl;
  cout << "total count via listing (should match above) = " << listCount << endl;
  t.reportTotal("total time");
}
