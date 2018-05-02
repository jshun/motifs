#include "gettime.h"
#include "graphIO.h"
#include "myVector.h"
#include "parallel.h"
#include "parseCommandLine.h"
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

  myVector *inEdges = newA(myVector, n);
  myVector *outEdges = newA(myVector, n);

  for (long i = 0; i < n; i++) {
    inEdges[i].init();
    outEdges[i].init();
  }

  long numBatches = 1 + (totalEdges - 1) / batchSize;
  long listCount = 0;
  long listCount3Hop = 0;

  // sparseSets to store the in/out edges in a batch
  sparseSet batchInEdges = sparseSet(batchSize, 1);
  sparseSet batchOutEdges = sparseSet(batchSize, 1);

  for (long i = 0; i < numBatches; i++) {
    // clear sparseSets
    batchInEdges.clearA();
    batchOutEdges.clearA();

    // add edges in this batch to sparseSet
    for (long j = i * batchSize;
         j < min((long)(i + 1) * batchSize, (long)totalEdges); j++) {
      uintT src = G.E[j].u;
      uintT dst = G.E[j].v;
      batchOutEdges.insert(src, dst);
      batchInEdges.insert(dst, src);
    }

    // extract the entries from the sparseSets.  In.A is an array of
    // pairs (a,b) where a is the vertex id and b is a pointer to its
    // myVector of in-edges.  Out.A is an array of pairs (a,b) where a
    // is the vertex id and b is a pointer to its myVector of
    // out-edges.
    _seq<kvPair> In = batchInEdges.entries();
    _seq<kvPair> Out = batchOutEdges.entries();

    // loop through all vertices in batch with in-neighbors
    for (long k = 0; k < In.n; k++) {
      uintE v = In.A[k].first;
      myVector *vIn = In.A[k].second;
      // intersect vertex v's in-neighbors with its
      // out-neighbors from batch
      myVector *vOut = batchOutEdges.find(v);
      if (vOut != NULL) {
        for (long g = 0; g < vIn->size(); g++) {
          for (long h = 0; h < vOut->size(); h++) {
            listCount++;

            // 3rd hop for every out neighbor
            uintE v2 = vOut->get(h);
            myVector *vOut2 = batchOutEdges.find(v2);
            if (vOut2 != NULL) {
              for (long h2 = 0; h2 < vOut2->size(); h2++) {
                listCount3Hop++;
              }
            }
          }
        }
      }

      // intersect vertex v's in-neighbors with its
      // out-neighbors from existing graph
      for (long h = 0; h < outEdges[v].size(); h++) {
        for (long g = 0; g < vIn->size(); g++) {
          listCount++;

          // and for every in-neighbor, we go an extra level to get 3-edge
          // paths.
          uintE v2 = vIn->get(g);
          myVector *vIn2 = batchInEdges.find(v2);
          if (vIn2 != NULL) {
            for (long g2 = 0; g2 < vIn2->size(); g2++) {
              listCount3Hop++;
            }
          }
        }
      }
    }

    for (long k = 0; k < Out.n; k++) {
      uintE v = Out.A[k].first;
      myVector *vOut = Out.A[k].second;
      // intersect vertex v's out-neighbors with its
      // in-neighbors from existing graph
      for (long h = 0; h < inEdges[v].size(); h++) {
        for (long g = 0; g < vOut->size(); g++) {
          listCount++;

          // and for every out-neighbor, we go an extra level to get 3-edge
          // paths.
          uintE v2 = vOut->get(g);
          myVector *vOut2 = batchOutEdges.find(v2);
          if (vOut2 != NULL) {
            for (long g2 = 0; g2 < vOut2->size(); g2++) {
              listCount3Hop++;
            }
          }
        }
      }
    }

    // add batch edges to graph
    for (long k = 0; k < In.n; k++) {
      uintE v = In.A[k].first;
      myVector *vIn = In.A[k].second;
      for (long g = 0; g < vIn->size(); g++) {
        inEdges[v].add(vIn->get(g));
      }
    }

    for (long k = 0; k < Out.n; k++) {
      uintE v = Out.A[k].first;
      myVector *vOut = Out.A[k].second;
      for (long g = 0; g < vOut->size(); g++) {
        outEdges[v].add(vOut->get(g));
      }
    }

    In.del();
    Out.del();

    // another way to add batch edges to graph. cache locality seems to be
    // worse.
    // for (long j = i * batchSize;
    //      j < min((long)(i + 1) * batchSize, (long)totalEdges); j++) {
    //   uintT src = G.E[j].u;
    //   uintT dst = G.E[j].v;
    //   outEdges[src].add(dst);
    //   inEdges[dst].add(src);
    // }
  }

  cout << "total count via listing (2 hops) = " << listCount << endl;
  cout << "total count via listing (3 hops) = " << listCount3Hop << endl;
  t.reportTotal("total time");

  // check answer
  long checkCount = 0;
  for (long i = 0; i < n; i++) {
    checkCount += outEdges[i].size() * inEdges[i].size();
  }
  cout << "expected count = " << checkCount << endl;

  batchInEdges.del();
  batchOutEdges.del();
  free(inEdges);
  free(outEdges);
}
