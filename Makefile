ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

#compilers
# ifdef CILK
# PCC = g++
# PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -O2 -DCILK $(INTT) $(INTE)
# PLFLAGS = -fcilkplus -lcilkrts

ifdef MKLROOT
PCC = icpc
PCFLAGS = -std=c++14 -O3 $(INTT) $(INTE)

# else ifdef OPENMP
# PCC = g++
# PCFLAGS = -std=c++14 -fopenmp -O3 -DOPENMP $(INTT) $(INTE)

else
PCC = g++
PCFLAGS = -std=c++14 -O2 $(INTT) $(INTE)
endif

COMMON = utils.h parseCommandLine.h parallel.h quickSort.h blockRadixSort.h transpose.h gettime.h
LOCAL_COMMON = graphIO.h sparseSet.h myVector.h
GENERATORS = rMatGraph gridGraph randLocalGraph SNAPtoAdj wghSNAPtoAdj adjGraphAddWeights adjToBinary count-2paths count-2paths-vector count-2paths-sparseVector count-3paths-sparseVector SNAPtoSequentialVerticesSNAP streamAppDAGGenerator

.PHONY: all clean
all: $(GENERATORS)

% : %.C $(COMMON) $(LOCAL_COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

clean :
	rm -f *.o $(GENERATORS)

cleansrc :
	make -s clean
	rm -f $(COMMON)
