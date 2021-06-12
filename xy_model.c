#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include <xy_model.h>
#include <union_find.h>

static double PI=3.14159265359;

lattice* create_square_lattice(int lx, int ly, double beta) {
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
            l->bondst[i_bond] = beta;
            l->nspin[i_bond]  = 2;
            l->bond2index[i_bond*mnspin+0] = j*lx+i;
            l->bond2index[i_bond*mnspin+1] = j*lx+(i+1)%lx;

            i_bond++;
        }
    }
    for(int j=0;j<ly;j++) {
        for(int i=0;i<lx;i++) {
            l->bondst[i_bond] = beta;
            l->nspin[i_bond]  = 2;
            l->bond2index[i_bond*mnspin+0] = j*lx+i;
            l->bond2index[i_bond*mnspin+1] = ((j+1)%ly)*lx+i;

            i_bond++;
        }
    }

    l->nbond = nbond;
    l->nsite = nsite;
    l->mnspin = mnspin;

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

    s->nsite = nsite;

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

void update(state* s, lattice* l, gsl_rng* rng) {
    int nsite  = s->nsite;
    int nbond  = l->nbond;
    int mnspin = l->mnspin;
    
    //initialize the cluster
    for(int i=0;i<nsite;i++) {
        s->cluster[i] = i;
        s->weight[i]  = 1;
    }

    //randomly pick up an spin for cluster
    int index = (int)(gsl_rng_uniform_pos(rng)*nsite);
    s->weight[index] = -1;

    //randomly select a direction for refference
    double r = gsl_rng_uniform_pos(rng)*2*PI;

    //construct the cluster base on Wolff algorithm
    int s1,s2;
    double theta1,theta2,dis,w;
    for(int i=0;i<nbond;i++) {
        s1 = l->bond2index[i*mnspin+0];
        s2 = l->bond2index[i*mnspin+1];

        theta1 = s->theta[s1];
        theta2 = s->theta[s2];

        w = 2*(l->bondst[i])*dot(r,theta1)*dot(r,theta2);

        dis = gsl_rng_uniform_pos(rng);
        if(dis<(1-exp(w))){
            merge(s->cluster,s->weight,s1,s2);
        }
    }

    //flip the cluster
    int* cluster = s->cluster;
    int* weight  = s->weight;
    double theta;
    for(int i=0;i<nsite;i++) {
        if(weight[root(cluster,i)]<0) {
            theta = s->theta[i];
            s->theta[i] = flip(r,theta);
        }
    }
}


int MEASUREMENT_COUNTER=0;
int BLOCK_SIZE=1000;
int NOBS=3;
double* BLOCK_DATA;
void measurement(state* s) {
    //initialize the BLOCK_DATA in the first time
    if(MEASUREMENT_COUNTER==0) {
        BLOCK_DATA = (double*)malloc(sizeof(double)*BLOCK_SIZE*NOBS);
    }

    //measure the observable
    int nsite = s->nsite;
    double tempmx = 0;
    double tempmy = 0;

    for(int i=0;i<nsite;i++) {
        tempmx += cos(s->theta[i]);
        tempmy += sin(s->theta[i]);
    }

    double mag = sqrt(tempmx*tempmx + tempmy*tempmy);
    double mag2 = mag*mag;
    double mag4 = mag*mag*mag*mag;

    //save each measurement in BLOCK_DATA
    BLOCK_DATA[BLOCK_SIZE*0+(MEASUREMENT_COUNTER%BLOCK_SIZE)] = mag;
    BLOCK_DATA[BLOCK_SIZE*1+(MEASUREMENT_COUNTER%BLOCK_SIZE)] = mag2;
    BLOCK_DATA[BLOCK_SIZE*2+(MEASUREMENT_COUNTER%BLOCK_SIZE)] = mag4;

    MEASUREMENT_COUNTER++;

    if((MEASUREMENT_COUNTER%BLOCK_SIZE)==0) {
        mag=0;
        mag2=0;
        mag4=0;
        for(int i=0;i<BLOCK_SIZE;i++) {
            mag  += BLOCK_DATA[BLOCK_SIZE*0+i];
            mag2 += BLOCK_DATA[BLOCK_SIZE*1+i];
            mag4 += BLOCK_DATA[BLOCK_SIZE*2+i];
        }
        mag  = mag/BLOCK_SIZE;
        mag2 = mag2/BLOCK_SIZE;
        mag4 = mag4/BLOCK_SIZE;

        printf("%.12e %.12e %.12e \n",mag,mag2,mag4);
    }
}

int main() {
    int Lx = 64;
    int Ly = 64;
    double Beta = 1.0;
    unsigned long int seed = 9733984098;
    BLOCK_SIZE=2000;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    lattice* l = create_square_lattice(Lx,Ly,Beta);
    state* s = create_state(l->nsite,rng);

    for(int i=0;;i++) {
        update(s,l,rng);
        measurement(s);
    }

    free_lattice(l);
    free_state(s);

    return 0;
}
