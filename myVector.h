#pragma once

struct myVector {
  long maxSize, end;
  uintT *A;
  void init() {
    maxSize = 1;
    end = 0;
    A = newA(uintT, maxSize);
  }
  void reset() {
    end = 0;
  }
  void resize() {
    uintT *B = newA(uintT, 2 * maxSize);
    for (long i = 0; i < maxSize; i++) {
      B[i] = A[i];
    }
    free(A); 
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
