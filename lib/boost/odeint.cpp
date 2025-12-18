#include "odeint.h"
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <vector>

using namespace std;
using namespace boost::numeric;
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type;

struct ODESystem
{
    ODESystem(odeint_t *ode) {odeint=ode;}
    odeint_t *odeint;

    void operator()(const state_type &x, state_type &dxdt, double t) const
    {
        double *y=(double*)calloc(x.size(),sizeof(double));
        double *dydt=(double*)calloc(x.size(),sizeof(double));

        for (int i=0;i<x.size();i++) y[i]=x[i];
        odeint->odefunc(t,y,dydt,odeint->param);

        for (int i=0;i<dxdt.size();i++) dxdt[i]=dydt[i];
        free(y); free(dydt);
    }
};

static int odeintv2_method_0(odeint_t *odeint, double dt, double *xt, int nx)
{
    state_type stateInOut;
    stateInOut.resize(nx);
    double t=0.0,h=0.0001;
    int type=0;

    if (dt<0.0) h=-h;
    for (int i=0;i<nx;i++) stateInOut[i]=xt[i];

    integrate_adaptive(make_controlled<adaptive_adams_bashforth_moulton<8,state_type>>(1E-10,1E-10),ODESystem(odeint),stateInOut,0.0,dt,h);
    for (int i=0;i<nx;i++) xt[i]=stateInOut[i];
    return 1;
}

static int odeintv2_method_1(odeint_t *odeint, double dt, double *xt, int nx)
{
    state_type stateInOut;
    stateInOut.resize(nx);
    double t=0.0,h=0.0001;
    int type=0;

    if (dt<0.0) h=-h;
    for (int i=0;i<nx;i++) stateInOut[i]=xt[i];

    controlled_adams_bashforth_moulton<adaptive_adams_bashforth_moulton<8,state_type>> stepper;
    integrate_adaptive(stepper,ODESystem(odeint),stateInOut,0.0,dt,h);
    for (int i=0;i<nx;i++) xt[i]=stateInOut[i];
    return 1;
}

static int odeintv2_method_2(odeint_t *odeint, double dt, double *xt, int nx)
{
    state_type stateInOut;
    stateInOut.resize(nx);
    double t=0.0,h=0.0001;
    int type=0;

    if (dt<0.0) h=-h;
    for (int i=0;i<nx;i++) stateInOut[i]=xt[i];

    integrate_adaptive(make_controlled<runge_kutta_fehlberg78<state_type>>(1E-10,1E-10),ODESystem(odeint),stateInOut,0.0,dt,h);
    for (int i=0;i<nx;i++) xt[i]=stateInOut[i];
    return 1;
}

static int odeintv2_method_3(odeint_t *odeint, double dt, double *xt, int nx)
{
    state_type stateInOut;
    stateInOut.resize(nx);
    double t=0.0,h=0.0001;
    int type=0;

    if (dt<0.0) h=-h;
    for (int i=0;i<nx;i++) stateInOut[i]=xt[i];

    integrate_adaptive(make_controlled<runge_kutta_cash_karp54<state_type>>(1E-10,1E-10),ODESystem(odeint),stateInOut,0.0,dt,h);
    for (int i=0;i<nx;i++) xt[i]=stateInOut[i];
    return 1;
}

extern int odeintv2(int method, odeint_t *odeint, double dt, double *xt, int nx)
{
    switch (method) {
        case 0: return odeintv2_method_0(odeint,dt,xt,nx);
        case 1: return odeintv2_method_1(odeint,dt,xt,nx);
        case 2: return odeintv2_method_2(odeint,dt,xt,nx);
        case 3: return odeintv2_method_3(odeint,dt,xt,nx);
        default:
            return odeintv2_method_0(odeint,dt,xt,nx);
    }
}