#ifndef ODEINT_H
#define ODEINT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int (*odefunc)(double t, const double *y, double *dydt, void *params);
    int nx;
    void *param;
} odeint_t;

extern int odeintv2(int method, odeint_t *odeint, double dt, double *xt, int nx);

#ifdef __cplusplus
}
#endif
#endif