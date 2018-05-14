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

  // edges that we have already batch processed
  myVector *processedInEdges = newA(myVector, n);
  myVector *processedOutEdges = newA(myVector, n);

  for (long i = 0; i < n; i++) {
    processedInEdges[i].init();
    processedOutEdges[i].init();
  }

  long numBatches = 1 + (totalEdges - 1) / batchSize;
  long listCount2Hop = 0;
  long listCount3Hop = 0;

  // sparseSets to store the in/out edges in a batch
  sparseSet batchInEdges = sparseSet(batchSize, 1);
  sparseSet batchOutEdges = sparseSet(batchSize, 1);

  for (long i = 0; i < numBatches; i++) {
    // clear sparseSets: these will hold the new batch of edges
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

    // extract the entries from the sparseSets. BIn.A is an array of
    // pairs (a,b) where a is the vertex id and b is a pointer to its
    // myVector of in-edges. BOut.A is an array of pairs (a,b) where a
    // is the vertex id and b is a pointer to its myVector of
    // out-edges.
    _seq<kvPair> BIn = batchInEdges.entries();
    _seq<kvPair> BOut = batchOutEdges.entries();

    // In the case of 3 hop directed paths, these are all the configurations
    // that yield a 3 hop directed path:
    //
    // in -> out -> out sequences:
    // BatchIn     -> BatchOut     -> BatchOut
    // BatchIn     -> BatchOut     -> ProcessedOut
    // BatchIn     -> ProcessedOut -> BatchOut
    // BatchIn     -> ProcessedOut -> ProcessedOut
    // ProcessedIn -> BatchOut     -> BatchOut
    // ProcessedIn -> BatchOut     -> ProcessedOut
    // ProcessedIn -> ProcessedOut -> BatchOut
    // ProcessedIn -> ProcessedOut -> ProcessedOut
    //
    // out -> in -> in sequences:
    // BatchOut     -> BatchIn     -> BatchIn
    // BatchOut     -> BatchIn     -> ProcessedIn
    // BatchOut     -> ProcessedIn -> BatchIn
    // BatchOut     -> ProcessedIn -> ProcessedIn
    // ProcessedOut -> BatchIn     -> BatchIn
    // ProcessedOut -> BatchIn     -> ProcessedIn
    // ProcessedOut -> ProcessedIn -> BatchIn
    // ProcessedOut -> ProcessedIn -> ProcessedIn

    // loop through all vertices in batch with in-neighbors
    for (long k = 0; k < BIn.n; k++) {

      uintE v = BIn.A[k].first;
      myVector *bIn = BIn.A[k].second;
      // intersect vertex v's in-neighbors with its
      // out-neighbors from batch
      myVector *bOut = batchOutEdges.find(v);
      if (bOut != NULL) {
        // 1st hop: BatchIn
        for (long g = 0; g < bIn->size(); g++) {
          // 2nd hop: BatchOut
          for (long h = 0; h < bOut->size(); h++) {
            // 2-hop config: BatchIn -> BatchOut
            listCount2Hop++;

            // 3rd hop: BatchOut
            uintE v2 = bOut->get(h);
            myVector *bOut2 = batchOutEdges.find(v2);
            if (bOut2 != NULL) {
              for (long h2 = 0; h2 < bOut2->size(); h2++) {
                // 3-hop config: BatchIn -> BatchOut -> BatchOut
                listCount3Hop++;
              }
            }

            // 3rd hop: ProcessedOut
            for (long h2 = 0; h2 < processedOutEdges[v2].size(); h2++) {
              // 3-hop config: BatchIn -> BatchOut -> ProcessedOut
              listCount3Hop++;
            }
          } // end 2nd hop

          // 2nd hop: ProcessedOut
          for (long h = 0; h < processedOutEdges[v].size(); h++) {
            // 3rd hop: BatchOut
            uintE v2 = bOut->get(h);
            myVector *bOut2 = batchOutEdges.find(v2);
            if (bOut2 != NULL) {
              for (long h2 = 0; h2 < bOut2->size(); h2++) {
                // 3-hop config: BatchIn -> ProcessedOut -> BatchOut
                listCount3Hop++;
              }
            }

            // 3rd hop: ProcessedOut
            for (long h2 = 0; h2 < processedOutEdges[v2].size(); h2++) {
              // 3-hop config: BatchIn -> ProcessedOut -> ProcessedOut
              listCount3Hop++;
            }
          } // end 2nd hop
        }
      }

      // intersect vertex v's in-neighbors with its
      // out-neighbors from existing graph
      //
      // 1st hop: ProcessedOut
      for (long h = 0; h < processedOutEdges[v].size(); h++) {
        // 2nd hop: BatchIn
        for (long g = 0; g < bIn->size(); g++) {
          // 2-hop config: ProcessedOut -> BatchIn
          listCount2Hop++;

          // 3rd hop: BatchIn
          uintE v2 = bIn->get(g);
          myVector *bIn2 = batchInEdges.find(v2);
          if (bIn2 != NULL) {
            for (long g2 = 0; g2 < bIn2->size(); g2++) {
              // 3-hop config: ProcessedOut -> BatchIn -> BatchIn
              listCount3Hop++;
            }
          }

          // 3rd hop: ProcessedIn
          for (long h2 = 0; h2 < processedInEdges[v2].size(); h2++) {
            // 3-hop config: ProcessedOut -> BatchIn -> ProcessedIn
            listCount3Hop++;
          }
        } // end 2nd hop: BatchIn

        // 2nd hop: ProcessedIn
        for (long g = 0; g < processedInEdges[v].size(); g++) {
          // 3rd hop: BatchIn
          uintE v2 = bIn->get(g);
          myVector *bIn2 = batchInEdges.find(v2);
          if (bIn2 != NULL) {
            for (long g2 = 0; g2 < bIn2->size(); g2++) {
              // 3-hop config: ProcessedOut -> ProcessedIn -> BatchIn
              listCount3Hop++;
            }
          }
        } // end 2nd hop: ProcessedIn

      } // end 1st hop: ProcessedOut

    } // end all vertices with in-neighbors

    // loop through all vertices in batch with out-neighbors
    for (long k = 0; k < BOut.n; k++) {
      uintE v = BOut.A[k].first;
      myVector *bOut = BOut.A[k].second;
      myVector *bIn = BIn.A[k].second; // used in later configs

      // intersect vertex v's out-neighbors with its
      // in-neighbors from existing graph
      //
      // 1st hop: ProcessedIn
      for (long h = 0; h < processedInEdges[v].size(); h++) {
        // 2nd hop: BatchOut
        for (long g = 0; g < bOut->size(); g++) {
          // 2-hop config: ProcessedIn -> BatchOut
          listCount2Hop++;

          // 3rd hop: BatchOut
          uintE v2 = bOut->get(g);
          myVector *bOut2 = batchOutEdges.find(v2);
          if (bOut2 != NULL) {
            for (long g2 = 0; g2 < bOut2->size(); g2++) {
              // 3-hop config: ProcessedIn -> BatchOut -> BatchOut
              listCount3Hop++;
            }
          }

          // 3rd hop: ProcessedOut
          for (long h2 = 0; h2 < processedOutEdges[v2].size(); h2++) {
            // 3-hop config: ProcessedIn -> BatchOut -> ProcessedOut
            listCount3Hop++;
          }
        } // end 2nd hop: BatchOut

        // 2nd hop: ProcessedOut
        for (long g = 0; g < processedOutEdges[v].size(); g++) {

          // 3rd hop: BatchOut
          uintE v2 = bOut->get(g);
          myVector *bOut2 = batchOutEdges.find(v2);
          if (bOut2 != NULL) {
            for (long g2 = 0; g2 < bOut2->size(); g2++) {
              // 3-hop config: ProcessedIn -> ProcessedOut -> BatchOut
              listCount3Hop++;
            }
          }

        } // end 2nd hop: ProcessedOut

      } // end 1st hop

      // 1st hop: BatchOut
      for (long g = 0; g < bOut->size(); g++) {
        // 2nd hop: BatchIn
        if (bIn != NULL) {
          for (long h = 0; h < bIn->size(); h++) {
            // 3rd hop: ProcessedIn
            uintE v2 = bIn->get(h);
            for (long h2 = 0; h2 < processedInEdges[v2].size(); h2++) {
              // 3-hop config: BatchOut -> BatchIn -> ProcessedIn
              listCount3Hop++;
            }
          } // end 2nd hop: BatchIn
        }

        // 2nd hop: ProcessedIn
        for (long h = 0; h < processedInEdges[v].size(); h++) {
          uintE v2 = bIn->get(h);
          myVector *bIn2 = batchInEdges.find(v2);
          if (bIn2 != NULL) {
            for (long h2 = 0; h2 < bIn2->size(); h2++) {
              // 3-hop config: BatchOut -> ProcessedIn -> BatchIn
              listCount3Hop++;
            }

            for (long h2 = 0; h2 < processedInEdges[v2].size(); h2++) {
              // 3-hop config: BatchOut -> ProcessedIn -> ProcessedIn
              listCount3Hop++;
            }
          }
        } // end 2nd hop: ProcessedIn

      } // end 1st hop: BatchOut

    } // end all vertices with out-neighbors

    // add batch edges to graph, as these are now processed
    for (long k = 0; k < BIn.n; k++) {
      uintE v = BIn.A[k].first;
      myVector *bIn = BIn.A[k].second;
      for (long g = 0; g < bIn->size(); g++) {
        processedInEdges[v].add(bIn->get(g));
      }
    }

    // same as above, but for outgoing edges
    for (long k = 0; k < BOut.n; k++) {
      uintE v = BOut.A[k].first;
      myVector *bOut = BOut.A[k].second;
      for (long g = 0; g < bOut->size(); g++) {
        processedOutEdges[v].add(bOut->get(g));
      }
    }

    BIn.del();
    BOut.del();

    // another way to add batch edges to graph. cache locality seems to be
    // worse.
    // for (long j = i * batchSize;
    //      j < min((long)(i + 1) * batchSize, (long)totalEdges); j++) {
    //   uintT src = G.E[j].u;
    //   uintT dst = G.E[j].v;
    //   processedOutEdges[src].add(dst);
    //   processedInEdges[dst].add(src);
    // }
  }

  cout << "total count via listing (2 hops) = " << listCount2Hop << endl;
  cout << "total count via listing (3 hops) = " << listCount3Hop << endl;
  t.reportTotal("total time");

  // check answer
  long checkCount = 0;
  for (long i = 0; i < n; i++) {
    checkCount += processedOutEdges[i].size() * processedInEdges[i].size();
  }
  cout << "expected count = " << checkCount << endl;

  batchInEdges.del();
  batchOutEdges.del();
  free(processedInEdges);
  free(processedOutEdges);
}
