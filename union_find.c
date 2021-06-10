#include <stdlib.h>

int root(int* p, int v) {
    while(p[v]!=v) {
        v = p[v];
    }

    return v;
}

static void compress(int* p, int v, int r) {
    int s;
    while(p[v]!=v)  {
        s = p[v];
        p[v] = r;
        v = s;
    }
}

void merge(int* p, int* w, int va, int vb) {
    int ra = root(p,va);
    int rb = root(p,vb);

    if(ra==rb) return;

    int r0 = ra;
    int r1 = rb;

    if(w[ra]<0 || w[rb]<0) {
        w[ra] = -abs(w[ra]);
        w[rb] = -abs(w[rb]);
        if(w[rb]<w[ra]) {
            r0 = rb;
            r1 = ra;
        }
    } else {
        if(w[ra]<w[rb]) {
            r0 = rb;
            r1 = ra;
        }
    }

    w[r0] = w[r0]+w[r1];
    p[r1] = r0;

    compress(p,va,r0);
    compress(p,vb,r0);
}
