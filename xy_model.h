#ifndef xy_model_h
#define xy_model_h

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
    double* theta;
    int* cluster;
    int* weight;
} state;

#endif
