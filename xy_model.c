#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include <xy_model.h>
#include <union_find.h>

static double PI=3.14159265359;

lattice* create_square_lattice(int lx, int ly) {
    int nsite = lx*ly;
    int nbond = 2*nsite;
    int mnspin = 2;

    lattice* l = (lattice*)malloc(sizeof(lattice));
    l->bond2index = (int*)malloc(sizeof(int)*nbond*mnspin);
    l->bondst = (double*)malloc(sizeof(double)*nbond);
    l->nspin = (int*)malloc(sizeof(int)*nbond);

    int i_bond=0;
    for(int j=0;j<ly;j++) {
        for(int i=0;i<lx;i++) {
            l->bondst[i_bond] = 1.0;
            l->nspin[i_bond]  = 2;
            l->bond2index[i_bond*mnspin+0] = j*lx+i;
            l->bond2index[i_bond*mnspin+1] = j*lx+(i+1)%lx;

            i_bond++;
        }
    }
    for(int j=0;j<ly;j++) {
        for(int i=0;i<lx;i++) {
            l->bondst[i_bond] = 1.0;
            l->nspin[i_bond]  = 2;
            l->bond2index[i_bond*mnspin+0] = j*lx+i;
            l->bond2index[i_bond*mnspin+1] = ((j+1)%ly)*lx+i;

            i_bond++;
        }
    }

    return l;
}

void free_lattice(lattice* l) {
    free(l->bond2index);
    free(l->bondst);
    free(l->nspin);
    free(l);
}

state* create_state(int nsite, gsl_rng* rng) {
    state* s = (state*)malloc(sizeof(state));
    
    s->theta = (double*)malloc(sizeof(double)*nsite);
    s->cluster = (int*)malloc(sizeof(int)*nsite);
    s->weight  = (int*)malloc(sizeof(int)*nsite);

    for(int i=0;i<nsite;i++) {
        s->theta[i] = gsl_rng_uniform_pos(rng)*PI*2;
    }

    return s;
}

void free_state(state* s) {
    free(s->theta);
    free(s->cluster);
    free(s->weight);
    free(s);
}

double per(double theta) {
    if(theta>=0) return theta-(int)(theta*0.5/PI)*2*PI;
    else return theta-(int)(theta*0.5/PI-1)*2*PI;
}

double dot(double theta1, double theta2) {
    return cos(theta2-theta1);
}

double flip(double r, double theta) {
    double rp = per(r+PI);
    double delta = theta-r;

    return per(rp-delta);
}

int main() {
    int Lx = 4;
    int Ly = 4;
    unsigned long int seed = 9733984098;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    lattice* l = create_square_lattice(Lx,Ly);
    state* s = create_state(l->nsite,rng);

    free_lattice(l);
    free_state(s);

    return 0;
}
