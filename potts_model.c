#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include <potts_model.h>
#include <union_find.h>

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
            l->bondst[i_bond] = 1.0-exp(-beta);
            l->nspin[i_bond]  = 2;
            l->bond2index[i_bond*mnspin+0] = j*lx+i;
            l->bond2index[i_bond*mnspin+1] = j*lx+(i+1)%lx;

            i_bond++;
        }
    }
    for(int j=0;j<ly;j++) {
        for(int i=0;i<lx;i++) {
            l->bondst[i_bond] = 1.0-exp(-beta);
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
            l->bondst[i_bond] = 1.0-exp(-beta);
            l->nspin[i_bond]  = 2;
            l->bond2index[i_bond*mnspin+0] = j*lx+i;
            l->bond2index[i_bond*mnspin+1] = j*lx+(i+1)%lx;

            i_bond++;
        }
    }
    for(int j=0;j<ly;j++) {
        for(int i=0;i<lx;i++) {
            if((i+j)%2==0) {
                l->bondst[i_bond] = 1.0-exp(-beta);
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

state* create_state(int nsite, int q, gsl_rng* rng) {
    state* conf = (state*)malloc(sizeof(state));
    
    conf->s = (int*)malloc(sizeof(int)*nsite);
    conf->cluster = (int*)malloc(sizeof(int)*nsite);
    conf->weight  = (int*)malloc(sizeof(int)*nsite);

    for(int i=0;i<nsite;i++) {
        conf->s[i] = (int)(gsl_rng_uniform_pos(rng)*q);
    }

    conf->nsite = nsite;
    conf->q = q;

    return conf;
}

void free_state(state* conf) {
    free(conf->s);
    free(conf->cluster);
    free(conf->weight);
    free(conf);
}

int MAX_NC=10;
int CLUSTER_SIZE=0;
int* LATTICE_MAP=NULL;
int* MAP_COUNTER=NULL;
void update_single_cluster(state* conf, lattice* l, double beta, gsl_rng* rng) {
    int s1,s2;
    int nsite  = conf->nsite;
    int qstat  = conf->q;
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
            conf->weight[i] = 1;
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
    conf->cluster[0] = index;
    conf->weight[index] = 0;
    
    //refference s for including to cluster
    double s_ref = conf->s[index];
    double p_add = 1.0-exp(-beta);

    //construct the cluster base on Wolff algorithm
    for(int i=0;i<cluster_size;i++) {
        for(int j=0;j<MAP_COUNTER[i];j++) {
            s1 = conf->cluster[i];
            s2 = LATTICE_MAP[s1*MAX_NC+j];

            if(conf->weight[s2] && (conf->s[s2])==s_ref) {
                if(gsl_rng_uniform_pos(rng)<p_add){
                    conf->cluster[cluster_size] = s2;
                    conf->weight[s2] = 0;
                    cluster_size++;
                }
            }
        }
    }
    CLUSTER_SIZE=cluster_size;

    //flip the cluster
    int check=1;
    int s_new=0;
    while(check) {
        s_new = (int)(gsl_rng_uniform_pos(rng)*qstat);
        if(s_new!=s_ref) check=0;
    }
    for(int i=0;i<cluster_size;i++) {
        s1 = conf->cluster[i];
        conf->s[s1] = s_new;
        conf->weight[s1] = 1;
    }

    //printf("%d %d \n",index,cluster_size);
}

int MEASUREMENT_COUNTER=0;
int BLOCK_SIZE=1000;
int NOBS=4;
double* BLOCK_DATA;
void measurement(state* conf, lattice* l) {
    //initialize the BLOCK_DATA in the first time
    if(MEASUREMENT_COUNTER==0) {
        BLOCK_DATA = (double*)malloc(sizeof(double)*BLOCK_SIZE*NOBS);
    }

    //measure the observable
    int nsite = conf->nsite;
    int mnspin = l->mnspin;
    int nbond = l->nbond;
    double q1 = 0.0;
    double energy = 0.0;

    int i1,i2;
    for(int i_bond=0;i_bond<nbond;i_bond++) {
        i1 = l->bond2index[i_bond*mnspin+0];
        i2 = l->bond2index[i_bond*mnspin+1];
        if((conf->s[i1])==(conf->s[i2])) energy+=1.0;
    }
    energy = -energy/nsite;

    for(int i=0;i<nsite;i++) {
        if(conf->s[i]==0) q1+=1.0;
    }
    q1 = q1/nsite;

    double mag = q1;
    double mag2 = q1*q1;
    double mag4 = q1*q1*q1*q1;

    //save each measurement in BLOCK_DATA
    BLOCK_DATA[BLOCK_SIZE*0+(MEASUREMENT_COUNTER%BLOCK_SIZE)] = mag;
    BLOCK_DATA[BLOCK_SIZE*1+(MEASUREMENT_COUNTER%BLOCK_SIZE)] = mag2;
    BLOCK_DATA[BLOCK_SIZE*2+(MEASUREMENT_COUNTER%BLOCK_SIZE)] = mag4;
    BLOCK_DATA[BLOCK_SIZE*3+(MEASUREMENT_COUNTER%BLOCK_SIZE)] = energy;

    MEASUREMENT_COUNTER++;

    if((MEASUREMENT_COUNTER%BLOCK_SIZE)==0) {
        mag=0;
        mag2=0;
        mag4=0;
        //energy=0;
        for(int i=0;i<BLOCK_SIZE;i++) {
            mag  += BLOCK_DATA[BLOCK_SIZE*0+i];
            mag2 += BLOCK_DATA[BLOCK_SIZE*1+i];
            mag4 += BLOCK_DATA[BLOCK_SIZE*2+i];
            //energy += BLOCK_DATA[BLOCK_SIZE*3+i];
        }
        mag  = mag/BLOCK_SIZE;
        mag2 = mag2/BLOCK_SIZE;
        mag4 = mag4/BLOCK_SIZE;
        //energy = energy/BLOCK_SIZE;

        //printf("%.12e %.12e %.12e %.12e \n",mag,mag2,mag4,energy);
        printf("%.16e \n",energy);
    }
}

int main(int argc, char** argv) {
    int Lx = atoi(argv[1]);
    int Ly = atoi(argv[2]);
    int Qstat = atoi(argv[3]);
    int Ltype = atoi(argv[4]);
    double T = atof(argv[5]);
    double Beta = 1.0/T;
    int NBLOCK = atoi(argv[6]);
    BLOCK_SIZE = atoi(argv[7]);
    int SWEEP = atoi(argv[8]);
    int THERMAL = atoi(argv[9]);
    unsigned long int seed = atoi(argv[10]);

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);

    lattice* l = create_square_lattice(Lx,Ly,Beta);
    if(Ltype==1) { 
        free_lattice(l);
        l = create_honeycomb_lattice(Lx,Ly,Beta);
    }
    state* conf = create_state(l->nsite,Qstat,rng);

    //int counter=0;
    //int cumulative_size = 0;
    //time_t start = clock();
    //time_t end;
    for(int i=0;i<THERMAL;i++) {
        update_single_cluster(conf,l,Beta,rng);
        
        /*cumulative_size += CLUSTER_SIZE;
        counter++;
        if(cumulative_size>(Lx*Ly)){
        //if(0){
            end = clock();
            printf("%d %d %f\n",i,counter,(double)(end-start)/CLOCKS_PER_SEC);
            start = clock();
            counter=0;
            cumulative_size=0;
            i++;
        }*/
    }

    for(int i=0;i<(NBLOCK*BLOCK_SIZE);i++) {
        for(int j=0;j<SWEEP;j++) 
            update_single_cluster(conf,l,Beta,rng);
        measurement(conf,l);
    }

    free_lattice(l);
    free_state(conf);
    gsl_rng_free(rng);
    free(BLOCK_DATA);
    free(LATTICE_MAP);
    free(MAP_COUNTER);

    return 0;
}
