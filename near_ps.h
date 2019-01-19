#ifndef NEAR_PS_H
#define NEAR_PS_H

#include <iostream>

#include <armadillo>
//#include <opencv2/opencv.hpp>
//#include <opencv2/core.hpp>

#include <QString>
#include <QDir>
#include <QFileDialog>
#include <QFile>
#include <QDebug>
#include <QObject>

#include <QJsonObject>
#include <QJsonDocument>
#include <QJsonParseError>
#include <QJsonArray>
#include <QDateTime>

using namespace std;
using namespace arma;

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
    vec Phi;
    mat S;
    mat z0;
    mat K;
    mat Dir;
    vec mu;
    // return
    mat X;
    mat Y;
    mat Z;
    // config
    int estimator = 0;
    bool shadows = false;
    double thr_norm = 1e-2;
    float  lambda = 0.1;
    // temple
    int nrows;
    int ncols;
    int ratio;
    int z_0;
    int nchannels;
    int nimgs;
    int npixel;
    float fx;
    float fy;
    float x0;
    float y0;
    umat index_matrix;

    field<sp_mat> make_gradient();
    field<cube> t_fun(vec z, vec u_tilde, vec v_tilde);
    mat shading_fun(vec z, cube tz, mat px_rep, mat py_rep, sp_mat Dx_rep, sp_mat Dy_rep);
    mat r_fun(mat rho, mat shadz, mat II, mat W_idx);
    mat J_fun();
    mat psi_fun(mat x);
    mat chi_fun(mat x);
    mat phi_fun(mat x);
    mat w_fun(mat x);
};

#endif // NEAR_PS_H
