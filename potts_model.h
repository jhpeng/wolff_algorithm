#ifndef potts_model_h
#define potts_model_h

typedef struct lattice {
    int nbond;
    int nsite;
    int mnspin;
    int* bond2index;
    double* bondst;
    int* nspin;
} lattice;

typedef struct state {
    int nsite;
    int q;
    int* s;
    int* cluster;
    int* weight;
} state;

#endif
