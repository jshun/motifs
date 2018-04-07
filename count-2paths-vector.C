#include "parseCommandLine.h"
#include "graphIO.h"
#include "parallel.h"
#include "gettime.h"
#include <map>
#include <list>

struct memoryPool {
  long total, used;
  uintT* pool; 
  memoryPool(long size) { 
    total = size; used = 0; pool = newA(uintT,size); 
  }
  long allocate(long size) {
    if(used + size > total) {
      //get more memory from OS
      cout << used << " " << size << " " << total << " get more memory\n";
      exit(0);
    }
    long start = used;
    used += size;
    return start;
  }
  void del() { free(pool); }
};

struct myVector {
  long maxSize, end; uintT* A;
  void init() { maxSize = 1; end = 0; A = newA(uintT,maxSize);}
  void resize() { 
    uintT* B = newA(uintT,2*maxSize);
    for(long i=0; i<maxSize; i++) {
      B[i] = A[i];
    }
    A = B;
    maxSize *= 2;    
  }
  void add(uintT v) { 
    if(end == maxSize) resize();
    A[end++] = v;
  }
  uintT size() { return end; }
  uintT get(uintT i) { return A[i]; }
};

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

  //memoryPool M(1 << 30);
  myVector* inEdges = newA(myVector,n);
  myVector* outEdges = newA(myVector,n);
  for(long i=0;i<n;i++) { inEdges[i].init(); }
  for(long i=0;i<n;i++) { outEdges[i].init(); }

  long numBatches = 1+(totalEdges-1)/batchSize;
  long listCount;

  for(long i=0;i<numBatches;i++) {
    for(long j=i*batchSize;j<min((long)(i+1)*batchSize,(long)totalEdges);j++) {
      uintT src = G.E[j].u;
      uintT dst = G.E[j].v;
      outEdges[src].add(dst);
      inEdges[dst].add(src);
    }
    //output number of 2-paths. currently it is not efficient
    //because it loops through all vertices. can be made more
    //efficient by keeping track of vertices with degree > 1 and
    //looping over those only
    listCount = 0;
    for(long k=0;k<n;k++) {
      for(long g=0;g<inEdges[k].size();g++) {
      	for(long h=0;h<outEdges[k].size();h++) {
      	  //check to force a read of edge list. should never return true.
      	  //if(inEdges[k].get(g) > n || outEdges[k].get(h) > n) cout << "oops\n";
      	  listCount++;
      	}
      }
    }
  }

  cout << "total count via listing = " << listCount << endl;
  t.reportTotal("total time");
  
  free(inEdges); free(outEdges);
}
