#include "gettime.h"
#include "graphIO.h"
#include "myVector.h"
#include "parallel.h"
#include "parseCommandLine.h"
#include "sparseSet.h"
#include <math.h>
#include <vector> // used to calculate median for "view size" estimates

// median (used for "view size" estimates)
long median(vector<long> &v) {
  size_t n = v.size() / 2;
  nth_element(v.begin(), v.begin() + n, v.end());
  return v[n];
}

// n choose k in O(k) time
long binomialCoeff(long n, long k) {
  long coeff = 1;

  if (k > n - k) {
    k = n - k;
  }

  for (long i = 0; i < k; ++i) {
    coeff *= (n - i);
    coeff /= (i + 1);
  }

  return coeff;
}

// the fill of a network is the proportion of edges to the total number of
// possible edges:
// http://konect.uni-koblenz.de/statistics/fill
double fill(long numNodes, long numEdges) {
  return (double)numEdges / ((double)numNodes * ((double)numNodes - 1));
}

// from Julian: a lower bound estimate on number of edges in k-hop directed
// spanner may be (n choose (k+1)) * (m/(n choose 2))^k
long lowerBoundEstimate(long numNodes, long numEdges, long k) {
  // 2 different ways of calculating the number of possible k-hop paths:
  long numPossibleKHopPaths = binomialCoeff(numNodes, k + 1); // n choose (k+1)
  long numPossibleKHopPaths2 = 1;
  for (long i = 0; i < k + 1; i++) {
    numPossibleKHopPaths2 *= numNodes - i;
  }

  // 2 different ways of calculating the number of possible edges:
  long numPossibleEdges = binomialCoeff(numNodes, 2);
  long numPossibleEdges2 = numNodes * (numNodes - 1);

  //  double num1HopPaths = (double)numEdges / (double)numPossibleEdges;
  double num1HopPaths = (double)numEdges / (double)numPossibleEdges2;

  //  long estimate = (long)((double)numPossibleKHopPaths * pow(num1HopPaths,
  //  k));
  long estimate = (long)((double)numPossibleKHopPaths2 * pow(num1HopPaths, k));

  return estimate;
}

// from Julian: an estimate which is unlikely to be either upper or lower bound
// for number of edges in directed k-hop spanner:
long notSureIfUpperOrLowerEstimate(long numNodes, double avgDegree, long k) {
  return (long)(numNodes * pow(avgDegree, k));
}

// same as above, but taking into account the fill.
long notSureIfUpperOrLowerEstimateWithFill(long numNodes, long numEdges,
                                           double avgDegree, long k) {
  return (long)(numNodes * pow(avgDegree, k) * fill(numNodes, numEdges));
}

// from Julian: an estimate which is likely to be a very loose upper bound
// for number of edges in directed k-hop spanner:
long looseUpperBound(long numNodes, long maxDegree, long k) {
  return (long)(pow(numNodes * (double)maxDegree, k));
}

// similar to above, but trying to incorporate fill into it.
long looseUpperBoundWithFillFactor(long numNodes, long numEdges, long maxDegree,
                                   long k) {
  // if we assume G is a DAG, then fill is:
  double upperBound = looseUpperBound(numNodes, maxDegree, k);

  return (long)(upperBound * fill(numNodes, numEdges));
}

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
  long n = max(G.numRows, G.numCols); // number of vertices
  cout << "G.numRows = " << G.numRows << ", G.numCols = " << G.numCols << endl;

  // stats
  double avgDegree = (double)(2 * totalEdges) / (double)n;
  cout << "#nodes = " << n << ", #edges = " << totalEdges
       << ", avg degree = " << avgDegree << endl;

  // start!
  cout << "starting timer\n";
  timer t;
  t.start();

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
        for (long hop1 = 0; hop1 < bIn->size(); hop1++) {
          // 2nd hop: BatchOut
          for (long hop2 = 0; hop2 < bOut->size(); hop2++) {
            // 2-hop config: BatchIn -> BatchOut
            listCount2Hop++;

            // 3rd hop: BatchOut
            uintE v2 = bOut->get(hop2);
            myVector *bOut2 = batchOutEdges.find(v2);
            if (bOut2 != NULL) {
              for (long hop3 = 0; hop3 < bOut2->size(); hop3++) {
// 3-hop config: BatchIn -> BatchOut -> BatchOut
#ifdef NDEBUG
                cout << "BIn -> BOut -> BOut: " << v << " -> " << v2 << " -> "
                     << bOut2->get(hop3) << endl;
#endif
                listCount3Hop++;
              }
            }

            // 3rd hop: ProcessedOut
            for (long hop3 = 0; hop3 < processedOutEdges[v2].size(); hop3++) {
// 3-hop config: BatchIn -> BatchOut -> ProcessedOut
#ifdef NDEBUG
              cout << "BIn -> BOut -> POut: " << v << " -> " << v2 << " -> "
                   << processedOutEdges[v2].get(hop3) << endl;
#endif
              listCount3Hop++;
            }
          } // end 2nd hop

          // 2nd hop: ProcessedOut
          for (long hop2 = 0; hop2 < processedOutEdges[v].size(); hop2++) {
            // 3rd hop: BatchOut
            uintE v2 = bOut->get(hop2);
            myVector *bOut2 = batchOutEdges.find(v2);
            if (bOut2 != NULL) {
              for (long hop3 = 0; hop3 < bOut2->size(); hop3++) {
// 3-hop config: BatchIn -> ProcessedOut -> BatchOut
#ifdef NDEBUG
                cout << "BIn -> POut -> BOut: " << v << " -> " << v2 << " -> "
                     << bOut2->get(hop3) << endl;
#endif
                listCount3Hop++;
              }
            }

            // 3rd hop: ProcessedOut
            for (long hop3 = 0; hop3 < processedOutEdges[v2].size(); hop3++) {
#ifdef NDEBUG
              // 3-hop config: BatchIn -> ProcessedOut -> ProcessedOut
              cout << "BIn -> POut -> POut: " << v << " -> " << v2 << " -> "
                   << processedOutEdges[v2].get(hop3) << endl;
#endif
              listCount3Hop++;
            }
          } // end 2nd hop
        }
      }

      // intersect vertex v's in-neighbors with its
      // out-neighbors from existing graph
      //
      // 1st hop: ProcessedOut
      for (long hop1 = 0; hop1 < processedOutEdges[v].size(); hop1++) {
        // 2nd hop: BatchIn
        for (long hop2 = 0; hop2 < bIn->size(); hop2++) {
          // 2-hop config: ProcessedOut -> BatchIn
          listCount2Hop++;

          // 3rd hop: BatchIn
          uintE v2 = bIn->get(hop2);
          myVector *bIn2 = batchInEdges.find(v2);
          if (bIn2 != NULL) {
            for (long hop3 = 0; hop3 < bIn2->size(); hop3++) {
// 3-hop config: ProcessedOut -> BatchIn -> BatchIn
#ifdef NDEBUG
              cout << "POut -> BIn -> BIn: " << v << " -> " << v2 << " -> "
                   << bIn2->get(hop3) << endl;
#endif
              listCount3Hop++;
            }
          }

          // 3rd hop: ProcessedIn
          for (long hop3 = 0; hop3 < processedInEdges[v2].size(); hop3++) {
// 3-hop config: ProcessedOut -> BatchIn -> ProcessedIn
#ifdef NDEBUG
            cout << "POut -> BIn -> PIn: " << v << " -> " << v2 << " -> "
                 << processedInEdges[v2].get(hop3) << endl;
#endif
            listCount3Hop++;
          }
        } // end 2nd hop: BatchIn

        // 2nd hop: ProcessedIn
        for (long hop2 = 0; hop2 < processedInEdges[v].size(); hop2++) {
          // 3rd hop: BatchIn
          uintE v2 = bIn->get(hop2);
          myVector *bIn2 = batchInEdges.find(v2);
          if (bIn2 != NULL) {
            for (long hop3 = 0; hop3 < bIn2->size(); hop3++) {
// 3-hop config: ProcessedOut -> ProcessedIn -> BatchIn
#ifdef NDEBUG
              cout << "POut -> PIn -> BIn: " << v << " -> " << v2 << " -> "
                   << bIn2->get(hop3) << endl;
#endif
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
      for (long hop1 = 0; hop1 < processedInEdges[v].size(); hop1++) {
        // 2nd hop: BatchOut
        for (long hop2 = 0; hop2 < bOut->size(); hop2++) {
          // 2-hop config: ProcessedIn -> BatchOut
          listCount2Hop++;

          // 3rd hop: BatchOut
          uintE v2 = bOut->get(hop2);
          myVector *bOut2 = batchOutEdges.find(v2);
          if (bOut2 != NULL) {
            for (long hop3 = 0; hop3 < bOut2->size(); hop3++) {
// 3-hop config: ProcessedIn -> BatchOut -> BatchOut
#ifdef NDEBUG
              cout << "PIn -> BOut -> BOut: " << v << " -> " << v2 << " -> "
                   << bOut2->get(hop3) << endl;
#endif
              listCount3Hop++;
            }
          }

          // 3rd hop: ProcessedOut
          for (long hop3 = 0; hop3 < processedOutEdges[v2].size(); hop3++) {
// 3-hop config: ProcessedIn -> BatchOut -> ProcessedOut
#ifdef NDEBUG
            cout << "PIn -> BOut -> POut: " << v << " -> " << v2 << " -> "
                 << processedOutEdges[v2].get(hop3) << endl;
#endif
            listCount3Hop++;
          }
        } // end 2nd hop: BatchOut

        // 2nd hop: ProcessedOut
        for (long hop2 = 0; hop2 < processedOutEdges[v].size(); hop2++) {

          // 3rd hop: BatchOut
          uintE v2 = bOut->get(hop2);
          myVector *bOut2 = batchOutEdges.find(v2);
          if (bOut2 != NULL) {
            for (long hop3 = 0; hop3 < bOut2->size(); hop3++) {
// 3-hop config: ProcessedIn -> ProcessedOut -> BatchOut
#ifdef NDEBUG
              cout << "PIn -> POut -> BOut: " << v << " -> " << v2 << " -> "
                   << bOut2->get(hop3) << endl;
#endif
              listCount3Hop++;
            }
          }

        } // end 2nd hop: ProcessedOut

      } // end 1st hop

      // 1st hop: BatchOut
      for (long hop1 = 0; hop1 < bOut->size(); hop1++) {

        // 2nd hop: ProcessedOut
        for (long hop2 = 0; hop2 < processedOutEdges[v].size(); hop2++) {
          uintE v2 = processedOutEdges[v].get(hop2);

          for (long hop3 = 0; hop3 < processedOutEdges[v2].size(); hop3++) {
// 3-hop config: BatchOut -> ProcessedOut -> ProcessedOut
#ifdef NDEBUG
            cout << "BOut -> POut -> POut: " << v << " -> " << v2 << " -> "
                 << processedOutEdges[v2].get(hop3) << endl;
#endif
            listCount3Hop++;
          }
        } // end 2nd hop: ProcessedOut

        // 2nd hop: BatchIn
        if (bIn != NULL) {
          for (long hop2 = 0; hop2 < bIn->size(); hop2++) {
            // 3rd hop: ProcessedIn
            uintE v2 = bIn->get(hop2);
            for (long hop3 = 0; hop3 < processedInEdges[v2].size(); hop3++) {
// 3-hop config: BatchOut -> BatchIn -> ProcessedIn
#ifdef NDEBUG
              cout << "BOut -> BIn -> PIn: " << v << " -> " << v2 << " -> "
                   << processedInEdges[v2].get(hop3) << endl;
#endif
              listCount3Hop++;
            }
          } // end 2nd hop: BatchIn
        }

        // 2nd hop: ProcessedIn
        for (long hop2 = 0; hop2 < processedInEdges[v].size(); hop2++) {
          uintE v2 = bIn->get(hop2);
          myVector *bIn2 = batchInEdges.find(v2);
          if (bIn2 != NULL) {
            for (long hop3 = 0; hop3 < bIn2->size(); hop3++) {
// 3-hop config: BatchOut -> ProcessedIn -> BatchIn
#ifdef NDEBUG
              cout << "BOut -> PIn -> BIn: " << v << " -> " << v2 << " -> "
                   << bIn2->get(hop3) << endl;
#endif
              listCount3Hop++;
            }

            for (long hop3 = 0; hop3 < processedInEdges[v2].size(); hop3++) {
// 3-hop config: BatchOut -> ProcessedIn -> ProcessedIn
#ifdef NDEBUG
              cout << "BOut -> PIn -> PIn: " << v << " -> " << v2 << " -> "
                   << processedInEdges[v2].get(hop3) << endl;
#endif
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

  cout << "total #paths via listing (2 hops) = " << listCount2Hop << endl;
  cout << "total #paths via listing (3 hops) = " << listCount3Hop << endl;
  t.reportTotal("total time");

  // check answer
  long checkCount = 0;
  for (long i = 0; i < n; i++) {
    checkCount += processedOutEdges[i].size() * processedInEdges[i].size();
  }
  cout << "expected #paths count (2 hops) = " << checkCount << endl;

  // max degree is in used in one of the spanner size estimates:
  long maxDegree = 0;
  long maxOutDegree = 0;
  vector<long> degrees;
  vector<long> outDegrees;
  for (long i = 0; i < processedInEdges->size(); i++) {
    long degree = processedInEdges[i].size() + processedOutEdges[i].size();
    maxDegree = max(maxDegree, degree);
    maxOutDegree = max(maxOutDegree, (long)processedOutEdges[i].size());
    degrees.push_back(degree);
    outDegrees.push_back(processedOutEdges[i].size());
  }
  long medianDegree = median(degrees);
  long medianOutDegree = median(outDegrees);
  cout << "maxDegree = " << maxDegree << endl;
  cout << "medianDegree = " << medianDegree << endl;
  cout << "medianOutDegree = " << medianOutDegree << endl;

  // different estimates for size in number of edges:
  // 2-hop:
  cout << endl
       << "actual #edges in 2-hop spanner = " << listCount2Hop * 2 << endl;
  cout << "1. estimated #edges (lower bound) = "
       << lowerBoundEstimate(n, totalEdges, 2) << endl;
  cout << "2a. estimated #edges (neither lower nor upper) = "
       << notSureIfUpperOrLowerEstimate(n, avgDegree, 2) << endl;
  cout << "2b. estimated #edges (neither lower nor upper w/ medianDegree) = "
       << notSureIfUpperOrLowerEstimate(n, medianDegree, 2) << endl;
  cout << "2c. estimated #edges (neither lower nor upper w/ medianOutDegree) = "
       << notSureIfUpperOrLowerEstimate(n, medianOutDegree, 2) << endl;
  cout << "estimated #edges (neither lower nor upper w/ fill) = "
       << notSureIfUpperOrLowerEstimateWithFill(n, totalEdges, avgDegree, 2)
       << endl;
  cout << "estimated #edges (neither lower nor upper w/ fill and medianDegree) "
          "= "
       << notSureIfUpperOrLowerEstimateWithFill(n, totalEdges, medianDegree, 2)
       << endl;
  cout << "estimated #edges (neither lower nor upper w/ fill and "
          "medianOutDegree) = "
       << notSureIfUpperOrLowerEstimateWithFill(n, totalEdges, medianOutDegree,
                                                2)
       << endl;
  cout << "3a. estimated #edges (loose upper) = "
       << looseUpperBound(n, maxDegree, 2) << endl;
  cout << "estimated #edges (loose upper w/ fill) = "
       << looseUpperBoundWithFillFactor(n, totalEdges, maxDegree, 2) << endl;
  cout << "estimated #edges on 2-hop path (from paper) = "
       << 1.0 / 2.0 * 2.0 * avgDegree * (double)n << endl;
  cout << "estimated #edges (tighter loose upper) = "
       << looseUpperBound(n, medianDegree, 2) << endl;
  cout << "3f. estimated #edges (tighter loose upper w/ fill) = "
       << looseUpperBoundWithFillFactor(n, totalEdges, medianDegree, 2) << endl;

  // 3-hop:
  cout << endl
       << "actual #edges in 3-hop spanner = " << listCount3Hop * 3 << endl;
  cout << "1. estimated #edges (lower bound) = "
       << lowerBoundEstimate(n, totalEdges, 3) << endl;
  cout << "2a. estimated #edges (neither lower nor upper) = "
       << notSureIfUpperOrLowerEstimate(n, avgDegree, 3) << endl;
  cout << "2b. estimated #edges (neither lower nor upper w/ medianDegree) = "
       << notSureIfUpperOrLowerEstimate(n, medianDegree, 3) << endl;
  cout << "2c. estimated #edges (neither lower nor upper w/ medianOutDegree) = "
       << notSureIfUpperOrLowerEstimate(n, medianOutDegree, 3) << endl;
  cout << "estimated #edges (neither lower nor upper w/ fill) = "
       << notSureIfUpperOrLowerEstimateWithFill(n, totalEdges, avgDegree, 3)
       << endl;
  cout << "estimated #edges (neither lower nor upper w/ fill and medianDegree) "
          "= "
       << notSureIfUpperOrLowerEstimateWithFill(n, totalEdges, medianDegree, 3)
       << endl;
  cout << "estimated #edges (neither lower nor upper w/ fill and "
          "medianOutDegree) = "
       << notSureIfUpperOrLowerEstimateWithFill(n, totalEdges, medianOutDegree,
                                                3)
       << endl;
  cout << "3a. estimated #edges (loose upper) = "
       << looseUpperBound(n, maxDegree, 3) << endl;
  cout << "estimated #edges (loose upper out deg only) = "
       << looseUpperBound(n, maxOutDegree, 3) << endl;
  cout << "estimated #edges (loose upper out deg only using fill) = "
       << looseUpperBoundWithFillFactor(n, totalEdges, maxOutDegree, 3) << endl;
  cout << "estimated #edges on 3-hop path (from paper) = "
       << 1.0 / 3.0 * 3.0 * avgDegree * (double)n << endl;
  cout << "estimated #edges (tighter loose upper) = "
       << looseUpperBound(n, medianDegree, 3) << endl;
  cout << "3f. estimated #edges (tighter loose upper w/ fill) = "
       << looseUpperBoundWithFillFactor(n, totalEdges, medianDegree, 3) << endl;

  batchInEdges.del();
  batchOutEdges.del();
  free(processedInEdges);
  free(processedOutEdges);
}
