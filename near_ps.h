#ifndef NEAR_PS_H
#define NEAR_PS_H

#include "lib.h"

class near_ps
{
public:
    near_ps(QJsonObject set, int zi);
    ~near_ps();

    //    void init(cv::Mat image = 0 );
    void init();

private:
    // load
    field<mat> I;
    imat mask;
    mat Phi;
    mat S;
    mat z;
    mat K;
    mat Dir;
    vec mu;
    // return
    mat X;
    mat Y;
    mat Z;
    // config
    int estimator = 0;
    int precond = 0;
    bool shadows = false;
    bool semi_calibrated = true;
    double thr_norm = 1e-2;
    float  lambda = 0.1;
    int maxit = 1;
    // temple
    int nrows;
    int ncols;
    int ratio;
    int z_0;
    int nchannels;
    int nimgs;
    int npixels;
    float fx;
    float fy;
    float x0;
    float y0;
    mat Iv;
    umat index_matrix;
    mat px_rep;
    mat py_rep;
    sp_mat Dx_rep;
    sp_mat Dy_rep;

    field<sp_mat> make_gradient();
    field<cube> t_fun(vec z, vec u_tilde, vec v_tilde);
    mat shading_fun(vec z, cube tz);
    mat r_fun(mat rho, mat shadz, mat II, mat W_idx);
    rowvec J_fun(mat rho, mat shadz, mat II, mat W_idx);
    mat psi_fun(mat x);
    mat chi_fun(mat x);
    mat phi_fun(mat x);
    mat w_fun(mat x);
};

#endif // NEAR_PS_H
