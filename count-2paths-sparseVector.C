#include "gettime.h"
#include "graphIO.h"
#include "parallel.h"
#include "parseCommandLine.h"
#include "myVector.h"
#include "sparseSet.h"

// Takes as input a file in SNAP format
//(http://snap.stanford.edu/data/index.html).
// Currently assumes a directed graph where each directed edge
// appears once as a pair (u,v).
int parallel_main(int argc, char *argv[]) {
  commandLine P(
      argc, argv,
      "-batch <batchSize> -totalEdges <totalEdges> <input SNAP file>");
  char *iFile = P.getArgument(0);
  long batchSize = P.getOptionLongValue("-batch", 10000);
  long totalEdges = P.getOptionLongValue("-totalEdges", 1000000);
  edgeArray<uintT> G = readSNAP<uintT>(iFile);

  // minimum of totalEdges and number of edges in the graph
  totalEdges = min(totalEdges, (long)G.nonZeros);
  cout << "starting timer\n";
  timer t;
  t.start();
  long n = max(G.numRows, G.numCols); // number of vertices

  // memoryPool M(1 << 30);
  myVector *inEdges = newA(myVector, n);
  myVector *outEdges = newA(myVector, n);
 
  for (long i = 0; i < n; i++) {
    inEdges[i].init();
    outEdges[i].init();
  }

  long numBatches = 1 + (totalEdges - 1) / batchSize;
  long listCount = 0;

  sparseSet batchInEdges = sparseSet(batchSize,1);
  sparseSet batchOutEdges = sparseSet(batchSize,1);

  for (long i = 0; i < numBatches; i++) {
    batchInEdges.clearA();
    batchOutEdges.clearA();

    // edges seen in this batch
    for (long j = i * batchSize;
         j < min((long)(i + 1) * batchSize, (long)totalEdges); j++) {
      uintT src = G.E[j].u;
      uintT dst = G.E[j].v;

      batchOutEdges.insert(src,dst);      
      batchInEdges.insert(dst,src);
    }

    //cout << batchInEdges.count() << " " << batchOutEdges.count() << endl;
    //cout << batchInEdges.m << " " << batchOutEdges.m << endl;
    
    _seq<kvPair> In = batchInEdges.entries();
    _seq<kvPair> Out = batchOutEdges.entries();

    //cout << "entries\n";
    //cout << In.n << " " << Out.n << endl;

    for(long k = 0; k < In.n; k++) {
      uintE v = In.A[k].first;
      myVector* vIn = In.A[k].second;
      //intersect vertex v from batch with in-neighbors with its
      //out-neighbors from batch
      myVector* vOut = batchOutEdges.find(v);
      if(vOut != NULL) {
	for(long g = 0; g < vIn->size(); g++) {
	  for(long h = 0; h < vOut->size(); h++) {
	    listCount++;
	  }
	}
      }

     
      //intersect vertex v from batch with in-neighbors with its
      //out-neighbors from existing graph
      for(long h = 0; h < outEdges[v].size(); h++) {
	for(long g = 0; g < vIn->size(); g++) {
	  listCount++;
	}
	
      }

    }

    //cout << "loop over in\n";

    for(long k = 0; k < Out.n; k++) {
      uintE v = Out.A[k].first;
      myVector* vOut = In.A[k].second;
      //intersect vertex v from batch with out-neighbors with its
      //in-neighbors from existing graph
      for(long h = 0; h < inEdges[v].size(); h++) {
	for(long g = 0; g < vOut->size(); g++) {
	  listCount++;
	}
      }
    }

    //cout << "loop over out\n";

    //add batch edges to graph
    for(long k = 0; k < In.n; k++) {
      uintE v = In.A[k].first;
      myVector* vIn = In.A[k].second;
      for(long g = 0; g < vIn->size(); g++) {
	inEdges[v].add(vIn->get(g));
      }
    }

    for(long k = 0; k < Out.n; k++) {
      uintE v = Out.A[k].first;
      myVector* vOut = Out.A[k].second;
      for(long g = 0; g < vOut->size(); g++) {
	outEdges[v].add(vOut->get(g));
      }
    }

    // long inSizeSoFar = 0, outSizeSoFar = 0;
    // for(long k=0;k<n;k++){ inSizeSoFar += inEdges[k].size(); outSizeSoFar += outEdges[k].size();}
    // cout << inSizeSoFar << " " << outSizeSoFar << endl;

    //cout << "add to graph\n";

    In.del(); Out.del();
    // for (long j = i * batchSize;
    //      j < min((long)(i + 1) * batchSize, (long)totalEdges); j++) {
    //   uintT src = G.E[j].u;
    //   uintT dst = G.E[j].v;
    //   outEdges[src].add(dst);
    //   inEdges[dst].add(src);
    // }

    // free(batchInEdges);
    // free(batchOutEdges);
  }

  cout << "total count via listing = " << listCount << endl;
  t.reportTotal("total time");

  batchInEdges.del(); batchOutEdges.del();
  free(inEdges);
  free(outEdges);
}
