#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
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

lattice* create_honeycomb_lattice(int lx, int ly, double beta) {
    int nsite = lx*ly;
    int nbond = (nsite/2)*3;
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
            if((i+j)%2==0) {
                l->bondst[i_bond] = beta;
                l->nspin[i_bond]  = 2;
                l->bond2index[i_bond*mnspin+0] = j*lx+i;
                l->bond2index[i_bond*mnspin+1] = ((j+1)%ly)*lx+i;

                i_bond++;
            }
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

int MAX_NC=10;
int CLUSTER_SIZE=0;
int* LATTICE_MAP=NULL;
int* MAP_COUNTER=NULL;
void update_single_cluster(state* s, lattice* l, gsl_rng* rng) {
    int s1,s2;
    int nsite  = s->nsite;
    int nbond  = l->nbond;
    int mnspin = l->mnspin;

    //create the LATTICE_MAP in the first time
    if(LATTICE_MAP==NULL) {
        LATTICE_MAP = (int*)malloc(sizeof(int)*MAX_NC*nsite);
        MAP_COUNTER = (int*)malloc(sizeof(int)*nsite);

        for(int i=0;i<nsite;i++) MAP_COUNTER[i]=0;

        for(int i_bond=0;i_bond<nbond;i_bond++) {
            for(int i=0;i<(l->nspin[i_bond]);i++) {
                s1 = l->bond2index[i_bond*mnspin+i];

                for(int j=0;j<(l->nspin[i_bond]);j++) {
                    s2 = l->bond2index[i_bond*mnspin+j];

                    if(s1!=s2) {
                        LATTICE_MAP[s1*MAX_NC+MAP_COUNTER[s1]] = s2;
                        MAP_COUNTER[s1]++;
                    }
                }
            }
        }

        //initialize the cluster
        for(int i=0;i<nsite;i++) {
            s->weight[i] = 1;
        }

/*
        //check the LATTICE_MAP
        for(int i=0;i<nsite;i++) {
            printf("%d : ",i);
            for(int j=0;j<MAP_COUNTER[i];j++) {
                printf("%d ",LATTICE_MAP[i*MAX_NC+j]);
            }
            printf("\n");
        }
*/

    }
    
    //randomly pick up an spin for cluster
    int cluster_size=1;
    int index = (int)(gsl_rng_uniform_pos(rng)*nsite);
    s->cluster[0] = index;
    s->weight[index] = 0;
    
    //randomly select a direction for refference
    double r = gsl_rng_uniform_pos(rng)*2*PI;

    //construct the cluster base on Wolff algorithm
    double theta1,theta2,dis,w;
    for(int i=0;i<cluster_size;i++) {
        for(int j=0;j<MAP_COUNTER[i];j++) {
            s1 = s->cluster[i];
            s2 = LATTICE_MAP[s1*MAX_NC+j];

            if(s->weight[s2]) {
                theta1 = s->theta[s1];
                theta2 = s->theta[s2];

                w = -2*(l->bondst[i])*dot(r,theta1)*dot(r,theta2);

                dis = gsl_rng_uniform_pos(rng);
                if(dis<(1-exp(w))){
                    s->cluster[cluster_size] = s2;
                    s->weight[s2] = 0;
                    cluster_size++;
                }
            }
        }
    }
    CLUSTER_SIZE=cluster_size;

    //flip the cluster
    double theta;
    for(int i=0;i<cluster_size;i++) {
        s1 = s->cluster[i];
        theta = s->theta[s1];
        s->theta[s1] = flip(r,theta);
        s->weight[s1] = 1;
    }

    //printf("%d %d \n",index,cluster_size);
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

int main(int argc, char** argv) {
    int Lx = atoi(argv[1]);
    int Ly = atoi(argv[2]);
    double Beta = atof(argv[3]);
    int NBLOCK = atoi(argv[4]);
    BLOCK_SIZE = atoi(argv[5]);
    int THERMAL = atoi(argv[6]);
    unsigned long int seed = atoi(argv[7]);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    //lattice* l = create_square_lattice(Lx,Ly,Beta);
    lattice* l = create_honeycomb_lattice(Lx,Ly,Beta);
    state* s = create_state(l->nsite,rng);

    int counter=0;
    int cumulative_size = 0;
    time_t start = clock();
    time_t end;
    for(int i=0;i<THERMAL;) {
        update_single_cluster(s,l,rng);
        cumulative_size += CLUSTER_SIZE;
        counter++;
        if(cumulative_size>(Lx*Ly)){
            end = clock();
            printf("%d %f\n",counter,(double)(end-start)/CLOCKS_PER_SEC);
            start = clock();
            counter=0;
            cumulative_size=0;
            i++;
        }
    }

    for(int i=0;i<(NBLOCK*BLOCK_SIZE);i++) {
        update_single_cluster(s,l,rng);
        measurement(s);
    }

    free_lattice(l);
    free_state(s);

    return 0;
}
