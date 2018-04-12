#include "gettime.h"
#include "graphIO.h"
#include "parallel.h"
#include "parseCommandLine.h"
#include <list>
#include <map>

struct memoryPool {
  long total, used;
  uintT *pool;
  memoryPool(long size) {
    total = size;
    used = 0;
    pool = newA(uintT, size);
  }
  long allocate(long size) {
    if (used + size > total) {
      // get more memory from OS
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
  long maxSize, end;
  uintT *A;
  void init() {
    maxSize = 1;
    end = 0;
    A = newA(uintT, maxSize);
  }
  void resize() {
    uintT *B = newA(uintT, 2 * maxSize);
    for (long i = 0; i < maxSize; i++) {
      B[i] = A[i];
    }
    // free(A); //should free
    A = B;
    maxSize *= 2;
  }
  void add(uintT v) {
    if (end == maxSize)
      resize();
    A[end++] = v;
  }
  uintT size() { return end; }
  uintT get(uintT i) { return A[i]; }
};

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
  uintT *degrees = newA(uintT, n);  // so we can skip nodes with degree < 2
  for (long i = 0; i < n; i++) {
    inEdges[i].init();
    outEdges[i].init();
    degrees[i] = 0;
  }

  long numBatches = 1 + (totalEdges - 1) / batchSize;
  long listCount = 0;

  for (long i = 0; i < numBatches; i++) {
    myVector *batchInEdges = newA(myVector, n);
    myVector *batchOutEdges = newA(myVector, n);
    // needed because malloc call from newA doesn't call default ctor.
    // we may consider adding "placement new" call to newA:
    // https://stackoverflow.com/questions/2995099/malloc-and-constructors
    for (long i = 0; i < n; i++) { batchInEdges[i].A = NULL; }
    for (long i = 0; i < n; i++) { batchOutEdges[i].A = NULL; }

    // edges seen in this batch
    for (long j = i * batchSize;
         j < min((long)(i + 1) * batchSize, (long)totalEdges); j++) {
      uintT src = G.E[j].u;
      uintT dst = G.E[j].v;

      // init if this is the first time we get an out-edge for this vertex
      if (batchOutEdges[src].A == NULL) {
        batchOutEdges[src].init();
      }
      batchOutEdges[src].add(dst);

      // init if this is the first time we get an in-edge for this vertex
      if (batchInEdges[dst].A == NULL) {
        batchInEdges[dst].init();
      }
      batchInEdges[dst].add(src);

      degrees[src]++;
      degrees[dst]++;
    }

    for (long k = 0; k < n; k++) {
      // skip vertices with degree smaller than 2, as it isn't possible to have
      // it in the center of a 2-hop path
      if (degrees[k] < 2) {
        continue;
      }

      // new incoming edges generate 2 hop paths with new outgoing
      // edges and with existing outgoing edges
      if (batchInEdges[k].A != NULL) {
        for (long g = 0; g < batchInEdges[k].size(); g++) {

          if (batchOutEdges[k].A != NULL) {
            for (long h = 0; h < batchOutEdges[k].size(); h++) {
              listCount++;
            }
          }

          for (long h = 0; h < outEdges[k].size(); h++) {
            listCount++;
          }
        }
      }

      // new outgoing edges generate 2 hop paths with existing
      // incoming edges
      if (batchOutEdges[k].A != NULL) {
        for (long g = 0; g < batchOutEdges[k].size(); g++) {
          for (long h = 0; h < inEdges[k].size(); h++) {
            listCount++;
          }
        }
      }
    }

    // store edges processed so far for future passes
    for (long j = i * batchSize;
         j < min((long)(i + 1) * batchSize, (long)totalEdges); j++) {
      uintT src = G.E[j].u;
      uintT dst = G.E[j].v;
      outEdges[src].add(dst);
      inEdges[dst].add(src);
    }

    free(batchInEdges);
    free(batchOutEdges);
  }

  cout << "total count via listing = " << listCount << endl;
  t.reportTotal("total time");

  free(inEdges);
  free(outEdges);
}
