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
  long maxSize, start, end; uintT* A; memoryPool* M;
  void init(memoryPool* _M) { M = _M; maxSize = 16; start = M->allocate(maxSize); end = start;}
  void resize() { 
    long newStart = M->allocate(2*maxSize);
    for(long i=0; i<maxSize; i++) {
      M->pool[newStart+i] = M->pool[start+i]; 
    }
    start = newStart;
    end = start + maxSize;
    maxSize *= 2;    
  }
  void add(uintT v) { 
    if((end-start) == maxSize) resize();
    M->pool[end++] = v;
  }
  uintT size() { return end-start; }
  uintT get(uintT i) { return M->pool[start+i]; }
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

  // keep track of both out and in-degrees for each vertex.
  uintT* outDegrees = newA(uintT,n);
  uintT* inDegrees = newA(uintT,n);
  for(long i=0;i<n;i++) {
    outDegrees[i] = 0;
    inDegrees[i] = 0;
  }

  // keep track of adjacency list for each vertex.
  // NOTE: currently using vanilla STL.
  memoryPool M(1 << 30);
  myVector* inEdges = newA(myVector,n);
  myVector* outEdges = newA(myVector,n);
  for(long i=0;i<n;i++) { inEdges[i].init(&M); }
  for(long i=0;i<n;i++) { outEdges[i].init(&M); }

  long numBatches = 1+(totalEdges-1)/batchSize;
  long count;
  long listCount;

  for(long i=0;i<numBatches;i++) {
    for(long j=i*batchSize;j<min((long)(i+1)*batchSize,(long)totalEdges);j++) {
      uintT src = G.E[j].u;
      uintT dst = G.E[j].v;
      outEdges[src].add(dst);
      inEdges[dst].add(src);
      outDegrees[src]++;
      inDegrees[dst]++;

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

      
      for(long g=0;g<inEdges[k].size();g++) {
      	for(long h=0;h<outEdges[k].size();h++) {
      	  //check to force a read of edge list. should never return true.
      	  //if(inEdges[k].get(g) > n || outEdges[k].get(h) > n) cout << "oops\n";
      	  listCount++;
      	}
      }

    }
  }
  cout << "total count = " << count << endl;
  cout << "total count via listing (should match above) = " << listCount << endl;
  t.reportTotal("total time");

  
  M.del();
  free(inEdges); free(outEdges);
  free(inDegrees); free(outDegrees);
}
