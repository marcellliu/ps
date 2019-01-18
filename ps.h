#ifndef PS_H
#define PS_H
//#include

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

struct Data{
    icube I;
    imat mask;
};

class ps
{
public:
    ps();
    ~ps();
    bool ini();

private:
    field<cube> I;
    field<mat> Iv;
    vec Phi;
    imat mask;
    mat z;
    mat K;
    mat S;
    mat Dir;
    vec mu;
    int nrows;
    int ncols;
    int ratio;
    int nchannels;
    int nimgs;
    float z0 = 100;
    float fx;
    float fy;
    float x0;
    float y0;
    uint npix;

    void main_func();
    void make_gradient(sp_mat *Dx, sp_mat *Dy);
    void shading_func(vec z, cube tz, mat px_rep, mat py_rep, sp_mat Dx_rep, sp_mat Dy_rep);
    void t_func(vec z, vec u, vec v, cube *T_filed, cube *grad_t);

    // parameters during operation

    field<imat> Omega;
    umat index_matrix;

};

#endif // PS_H
