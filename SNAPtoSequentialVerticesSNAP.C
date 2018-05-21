#include "gettime.h"
#include "graphIO.h"
#include "myVector.h"
#include "parallel.h"
#include "parseCommandLine.h"
#include "sparseSet.h"
#include <math.h>

void writeToSNAPFile(edgeArray<uintT> *G, uintT numberEdges,
                     sparseSet *newVertexIds, char *fileName) {
  ofstream file(fileName);

  if (!file.is_open()) {
    std::cout << "Unable to open file: " << fileName << std::endl;
    return;
  }

  for (long i = 0; i < numberEdges; i++) {
    uintT src = G->E[i].u;
    uintT dst = G->E[i].v;

    // Assumption: these are necessarily size 1, and not NULL.
    myVector *newSrcId = newVertexIds->find(src);
    myVector *newDstId = newVertexIds->find(dst);

    if (newSrcId == NULL || newDstId == NULL || newSrcId->size() != 1 ||
        newDstId->size() != 1) {
      cout << "Aborting: number of new ids for each vertex should be exactly 1."
           << endl;
      return;
    }

    file << newSrcId->get(0) << " " << newDstId->get(0) << endl;
  }

  file.close();
}

// Takes as input a file in SNAP format
//(http://snap.stanford.edu/data/index.html).
// Currently assumes a directed graph where each directed edge
// appears once as a pair (u,v).
int parallel_main(int argc, char *argv[]) {
  commandLine P(argc, argv, "<input SNAP file> <output SNAP file>");
  char *iFile = P.getArgument(1);
  char *oFile = P.getArgument(0);

  cout << "input: " << iFile << ", output: " << oFile << endl;
  edgeArray<uintT> G = readSNAP<uintT>(iFile);

  long totalEdges = (long)G.nonZeros;
  long n = max(G.numRows, G.numCols); // number of vertices (non-sequential)
  cout << "G.numRows = " << G.numRows << ", G.numCols = " << G.numCols << endl;

  // map of vertexId to sequentialVertexId.  totalEdges is a valid upper bound
  // on the number of unique vertices.
  sparseSet sequentialVertexIds = sparseSet(totalEdges, 1);
  sequentialVertexIds.clearA();

  // scan through all edges, regenerating ids (0-indexed).
  long maxVertexId = 0;
  for (long i = 0; i < totalEdges; i++) {
    uintT src = G.E[i].u;
    uintT dst = G.E[i].v;

    myVector *srcMatches = sequentialVertexIds.find(src);
    if (srcMatches == NULL) {
      sequentialVertexIds.insert(src, maxVertexId);
      maxVertexId++;
    }

    myVector *dstMatches = sequentialVertexIds.find(dst);
    if (dstMatches == NULL) {
      sequentialVertexIds.insert(dst, maxVertexId);
      maxVertexId++;
    }
  }

  // print out all edges using the new ids
  cout << "Writing re-indexed edges to " << oFile << "..." << endl;
  writeToSNAPFile(&G, totalEdges, &sequentialVertexIds, oFile);
  cout << "Done!" << endl;
}
