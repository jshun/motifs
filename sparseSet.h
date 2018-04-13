// This code is based on the paper "Phase-Concurrent Hash Tables for 
// Determinism" by Julian Shun and Guy Blelloch from SPAA 2014.
// Copyright (c) 2014 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights (to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#pragma once
#include "parallel.h"
#include "utils.h"
#include "math.h"
#include "myVector.h"
using namespace std;

// returns the log base 2 rounded up (works on ints or longs or unsigned versions)
template <class T>
static int log2RoundUp(T i) {
  int a=0;
  T b=i-1;
  while (b > 0) {b = b >> 1; a++;}
  return a;
}

typedef pair<uintE,myVector*> kvPair;
class sparseSet {
 public:
  uintT m;
  intT mask;
  //kvPair empty;
  myVector* vectors;
  kvPair* TA;
  float loadFactor;

  // needs to be in separate routine due to Cilk bugs
  void clearA() {
    parallel_for (long i=0; i < m; i++) {
      TA[i].first = UINT_E_MAX;
    }
  }

  void initA() {
    parallel_for (long i=0; i < m; i++) {
      TA[i].second = &vectors[i];
      vectors[i].init();
    }
  }

  struct notEmptyF { 
    int operator() (kvPair a) {return a.first != UINT_E_MAX;}};

  inline uintT hashToRange(uintT h) {return h & mask;}
  inline uintT firstIndex(uintT v) {return hashToRange(hashInt(v));}
  inline uintT incrementIndex(uintT h) {return hashToRange(h+1);}

  // Size is the maximum number of values the hash table will hold.
  // Overfilling the table could put it into an infinite loop.
 sparseSet(long size, float _loadFactor) :
  loadFactor(_loadFactor),
    m((uintT) 1 << log2RoundUp((uintT)(_loadFactor*size)+100)),
    mask(m-1),
    TA(newA(kvPair,m)) 
      { 
	vectors = newA(myVector,m);
	initA();
	//empty=make_pair(UINT_E_MAX,NULL);
	//clearA(TA,m,empty); 
      }

  // Deletes the allocated arrays
  void del() {
    free(TA); 
    free(vectors);
  }

  // nondeterministic insert
  bool insert(uintE src, uintE dst) {
    uintT h = firstIndex(src);
    while (1) {
      //kvPair c;
      //c = TA[h];
      if(TA[h].first == UINT_E_MAX && CAS(&TA[h].first,UINT_E_MAX,src)) {
	//cout << "a\n";
  	TA[h].second->reset(); //reset myVector to size 0
  	TA[h].second->add(dst);
	//cout << "b\n";
  	return 1; //return true if value originally didn't exist
      }
      else if (TA[h].first == src) {
  	TA[h].second->add(dst);
  	return 0;
      }
    
      // move to next bucket
      h = incrementIndex(h);
    }
    return 0; // should never get here
  }

  myVector* find(uintE v) {
    uintT h = firstIndex(v);
    kvPair c = TA[h];
    while (1) {
      if (c.first == UINT_E_MAX) return NULL;
      else if (v == c.first)
  	return c.second;
      h = incrementIndex(h);
      c = TA[h];
    }
  }

  template <class F>
  void map(F f){ 
    parallel_for(long i=0;i<m;i++)
      if(TA[i].first != UINT_E_MAX) f(TA[i]);
  }

  template <class F>
  void mapIndex(F f){ 
    parallel_for(long i=0;i<m;i++)
      if(TA[i].first != UINT_E_MAX) f(TA[i],i);
  }


  // returns all the current entries compacted into a sequence
  _seq<kvPair> entries() {
    bool *FL = newA(bool,m);
    parallel_for (long i=0; i < m; i++) 
      FL[i] = (TA[i].first != UINT_E_MAX);
    _seq<kvPair> R = pack((kvPair*)NULL, FL, (uintT) 0, m, sequence::getA<kvPair,uintE>(TA));
    //sequence::pack(TA,(entry*)NULL,FL,m);
    free(FL);
    return R;
  }

  // returns all the current entries satisfying predicate f compacted into a sequence
  template <class F>
  _seq<kvPair> entries(F f) {
    bool *FL = newA(bool,m);
    parallel_for (long i=0; i < m; i++) 
      FL[i] = (TA[i].first != UINT_E_MAX && f(TA[i]));
    _seq<kvPair> R = pack((kvPair*)NULL, FL, (uintT) 0, m, sequence::getA<kvPair,uintE>(TA));
    //sequence::pack(TA,(entry*)NULL,FL,m);
    free(FL);
    return R;
  }

  // returns the number of entries
  intT count() {
    return sequence::mapReduce<intT>(TA,m,addF<intT>(),notEmptyF());
  }

  // prints the current entries along with the index they are stored at
  void print() {
    cout << "vals = ";
    for (long i=0; i < m; i++) {
      if (TA[i].first != UINT_E_MAX)
  	{ cout << "(" << TA[i].first << "," << TA[i].second << ") ";}
    }
    cout << endl;
  }
};

