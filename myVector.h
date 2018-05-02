#pragma once

struct myVector {
  long maxSize, end, numRemoved;
  uintT *A;
  void init() {
    maxSize = 1;
    end = 0;
    numRemoved = 0;
    A = newA(uintT, maxSize);
  }
  void reset() {
    end = 0;
  }
  void resize(long targetSize) {
    uintT *B = newA(uintT, targetSize);
    for (long i = 0; i < maxSize; i++) {
      // skip removed values
      if (A[i] != -1) {
        B[i] = A[i];
      }
    }
    free(A);
    A = B;
    maxSize = targetSize;
  }
  void add(uintT v) {
    if (end == maxSize) {
      resize(maxSize * 2);
    }
    A[end++] = v;
  }
  // Removes element at index i.
  void remove(uintT i) {
    // -1 is dummy value used as deletion marker.
    if (A[i] != -1) {
      A[i] = -1;
      numRemoved++;
    }
    // trigger garbage collection if too many have been removed.
    if (numRemoved > (size() / 2)) {
      resize(maxSize / 2);
    }
  }
  // Removes first occurrence of value v.
  void removeFirstOcurrence(uintT v) {
    long index = -1;
    for (long i = 0; i < end; i++) {
      if (A[i] == v) {
        index = i;
        break;
      }
    }

    if (index != -1) {
      remove(index);
    }
  }
  uintT size() { return end; }
  uintT get(uintT i) { return A[i]; }
};
