#include "ps.h"

//****** Data  ******//
//

ps::ps()
{
    nimgs = 8;
    nchannels = 1;
    ratio = 1;

    I = field<cube>(nchannels);
    mask = ones<imat>(3,3);
    mask(1,0) = 0;
    mask(2,0) = 0;
    mask(2,2) = 0;

    I(0) = randu<cube>(3,3,nimgs);

    nrows = I(0).n_rows;
    ncols = I(0).n_cols;

    Phi = vec(1);
    Phi = 1.0;

    S = mat(nimgs,3);
    S.randu();
    Dir = mat(nimgs,3);
    Dir.randu();

    mu = ones<vec>(nimgs);

    ini();
    main_func();
}

ps::~ps()
{

}

bool ps::ini()
{
    K = ones<mat>(3,3);

    if (ratio != 0) {
        uvec ir = regspace<uvec>(0,ratio,nrows-1);
        uvec ic = regspace<uvec>(0,ratio,ncols-1);
        for (int ch=0;ch<nchannels;ch++) {
            for (int i=0;i<nimgs;i++) {
                I(ch).slice(i) = I(ch).slice(i).elem(ir,ic);
            }
        }
        mask = mask.elem(ir, ic);
        K.cols(0,2) = K.cols(0,2)/ratio;
        ncols = I(0,0).n_cols;
        nrows = I(0,0).n_rows;
    }
    z = ones(ncols,nrows)*z0;
    return true;
}

void ps::main_func()
{
    // intrinsics
    fx = K(0,0);
    fy = K(1,1);
    x0 = K(0,2);
    y0 = K(1,2);

    sp_mat Dx, Dy;
    make_gradient(&Dx, &Dy);

    if (Phi.n_rows == 1) {
        Phi = ones<vec>(3)*Phi(0);
    }

    double max_I = 1;

    for (int ch=0;ch<nchannels;ch++) {
        I(ch) = I(ch)/Phi(ch);
        if (max_I<I(ch).max()) max_I = I(ch).max();
    }

    // Scaled pixel units
    uvec u = regspace<uvec>(1,ncols);
    uvec v = regspace<uvec>(1,nrows);
    umat uu = repmat(u,1,nrows).t();
    umat vv = repmat(v,1,ncols);
    mat u_tilde = conv_to<mat>::from(uu)-x0;
    mat v_tilde = conv_to<mat>::from(vv)-y0;

    // Use a bounding box
    uvec imask = find(mask>0);
    npix = imask.n_elem;

    // Some useful variables
    //    mat px_rep = repmat(u_tilde(imask),nimgs,1);
    //    mat py_rep = repmat(v_tilde(imask),nimgs,1);
    //    sp_mat Dx_rep = repmat(Dx,nimgs,1);
    //    sp_mat Dy_rep = repmat(Dy,nimgs,1);

    // Vectorize data
    Iv = field<mat>(nchannels);
    mat It(imask.n_elem,nimgs);
    for (int ch=0;ch<nchannels;ch++) {
        I(ch).reshape(nrows*ncols,nimgs,1);
        It = I(ch);
        It = It(imask,regspace<uvec>(0,nimgs-1));
        Iv(ch) = It;
    }

    // Sort images to remove shadows and highlights
    field<umat> W_idx = field<umat>(nchannels);
    for (int ch=0;ch<nchannels;ch++) {
        W_idx(ch) = zeros<umat>(imask.n_elem,nimgs);
        for (uint i=0;i<imask.n_elem;i++) {
            uvec id = sort_index(Iv(ch).row(i));
            id = id + id*(imask.n_elem-1) + i;
            id = id(regspace<uvec>(2,6));
            W_idx(ch)(id) = ones<uvec>(id.n_elem);
        }
    }

    // Initialize variables
    vec z_tilde;
    vec rho_tilde;
    //    mat XYZ;
    field<mat> rho(nchannels);
    uvec bg = find(mask==0);
    z(bg) = ones<vec>(bg.n_elem)*datum::nan;
    z_tilde = log(z(imask));
    //    sp_mat zz = conv_to<sp_mat>::from(z_tilde);
    for (int ch=0;ch<nchannels;ch++) {
        rho(ch) = ones(nrows,ncols)/max_I;
    }
    ////    XYZ(0) = z*u_tilde;
    ////    XYZ(1) = z*v_tilde;
    ////    XYZ(2) = z;

    vec zx = Dx*z_tilde;
    vec zy = Dy*z_tilde;
    mat Nx = zeros(nrows,ncols);
    mat Ny = zeros(nrows,ncols);
    mat Nz = zeros(nrows,ncols);
    Nx(imask) = fx*zx;
    Ny(imask) = fy*zy;
    Nz(imask) = -u_tilde(imask)%zx-v_tilde(imask)%zy -1;
    mat dz = sqrt(Nx*Nx+Ny*Ny+Nz*Nz);
    mat N(nrows*3,ncols);
    N.submat(span(0,nrows-1),span(0,ncols-1)) = Nx/dz;
    N.submat(span(nrows,2*nrows-1),span(0,ncols-1)) = Ny/dz;
    N.submat(span(2*nrows,3*nrows-1),span(0,ncols-1)) = Nz/dz;
    //    rho_tilde = (rho/dz);

    // Initial energy
    float energy = 0;
    cube T_filed;
    cube grad_t;
    t_fun(z_tilde,u_tilde(imask)/fx,v_tilde(imask)/fy,&T_filed,&grad_t);


    cout <<1<< endl;
}

void ps::make_gradient(sp_mat *Dx, sp_mat *Dy)
{
    imat Omega_padded = zeros<imat>(nrows+2,ncols+2);
    Omega_padded.submat(span(1,nrows),span(1,ncols)) = mask;

    Omega = field<imat>(4);
    // pixel who have bottom neighbour
    Omega(0) = mask%Omega_padded(2,1,size(nrows,ncols));
    // pixel who have top neighbour
    Omega(1) = mask%Omega_padded(0,1,size(nrows,ncols));
    // pixel who have right neighbour
    Omega(2) = mask%Omega_padded(1,2,size(nrows,ncols));
    // pixel who have left neighbour
    Omega(3) = mask%Omega_padded(1,0,size(nrows,ncols));

    uvec imask = find(mask>0);
    index_matrix = zeros<umat>(nrows,ncols);
    for (uint i=0;i<imask.size();i++) {
        index_matrix(imask(i)) = i+1;
    }

    uvec idx_c;
    umat sub;
    uvec indices_center;
    uvec indices_right;
    uvec II;
    uvec JJ;
    vec KK;

    // When there is a neighbor on the bottom : forward differences
    idx_c = find(Omega(0)==1);
    sub = ind2sub(size(mask),idx_c);
    sub.row(0) = sub.row(0)+1;
    indices_center = index_matrix(idx_c)-1;
    indices_right = index_matrix(sub2ind(size(mask),sub))-1;
    II = indices_center;
    JJ = indices_right;
    KK = ones<vec>(size(indices_center));
    II = join_cols(II,indices_center);
    JJ = join_cols(JJ,indices_center);
    KK = join_cols(KK,-KK);
    sp_mat Dup(imask.n_elem,imask.n_elem);
    for (uint i=0;i<KK.n_elem;i++) {
        Dup(II(i),JJ(i)) = KK(i);
    }

    // When there is a neighbor on the top : backword differences
    idx_c = find(Omega(1)==1 && Omega(0)==0);
    sub = ind2sub(size(mask),idx_c);
    sub.row(0) = sub.row(0)-1;
    indices_center = index_matrix(idx_c)-1;
    indices_right = index_matrix(sub2ind(size(mask),sub))-1;
    II = indices_center;
    JJ = indices_right;
    KK = -ones<vec>(size(indices_center));
    II = join_cols(II,indices_center);
    JJ = join_cols(JJ,indices_center);
    KK = join_cols(KK,-KK);
    sp_mat Dum(imask.n_elem,imask.n_elem);
    for (uint i=0;i<KK.n_elem;i++) {
        Dum(II(i),JJ(i)) = KK(i);
    }

    // When there is a neighbor on the right : forward differences
    idx_c = find(Omega(2)==1);
    sub = ind2sub(size(mask),idx_c);
    sub.row(1) = sub.row(1)+1;
    indices_center = index_matrix(idx_c)-1;
    indices_right = index_matrix(sub2ind(size(mask),sub))-1;
    KK = ones<vec>(size(indices_center));
    II = join_cols(indices_center,indices_center);
    JJ = join_cols(indices_right,indices_center);
    KK = join_cols(KK,-KK);
    sp_mat Dvp(imask.n_elem,imask.n_elem);
    for (uint i=0;i<KK.n_elem;i++) {
        Dvp(II(i),JJ(i)) = KK(i);
    }

    // When there is a neighbor on the left : backword differences
    idx_c = find(Omega(3)==1 && Omega(2)==0);
    sub = ind2sub(size(mask),idx_c);
    sub.row(1) = sub.row(1)-1;
    indices_center = index_matrix(idx_c)-1;
    indices_right = index_matrix(sub2ind(size(mask),sub))-1;
    II = indices_center;
    JJ = indices_right;
    KK = -ones<vec>(size(indices_center));
    II = join_cols(II,indices_center);
    JJ = join_cols(JJ,indices_center);
    KK = join_cols(KK,-KK);
    sp_mat Dvm(imask.n_elem,imask.n_elem);
    for (uint i=0;i<KK.n_elem;i++) {
        Dvm(II(i),JJ(i)) = KK(i);
    }

    *Dx = Dvp + Dvm;
    *Dy = Dup + Dum;
}

void ps::shading_func(vec z, cube tz, mat px_rep, mat py_rep, sp_mat Dx_rep, sp_mat Dy_rep)
{

}

void ps::t_func(vec z, vec u, vec v, cube *T, cube *gT)
{
    npix = z.n_elem;
    vec exp_z = exp(z);
    mat XYZ = zeros(npix,3);
    XYZ.col(0) = z;
    XYZ.col(1) = u;
    XYZ.col(2) = v;
    // T_field
    mat a_filed = zeros<mat>(npix,nimgs);
    mat da_filed = zeros<mat>(npix,nimgs);
    cube T_filed= zeros<cube>(npix,nimgs,3);
    cube grad_t = zeros<cube>(npix,nimgs,3);
    for (int i=0;i<nimgs;i++) {
        // Unit lighting filed
        T_filed.slice(0).col(i) = S(i,0)-XYZ.col(0);
        T_filed.slice(1).col(i) = S(i,1)-XYZ.col(1);
        T_filed.slice(2).col(i) = S(i,2)-XYZ.col(2);
        vec normS_i = pow(T_filed.slice(0).col(i),2)+pow(T_filed.slice(1).col(i),2)+pow(T_filed.slice(2).col(i),2);
        normS_i = sqrt(normS_i);
        vec scal_prod = -T_filed.slice(0).col(i)*Dir(i,0)-T_filed.slice(1).col(i)*Dir(i,1)-T_filed.slice(2).col(i)*Dir(i,2);
        a_filed.col(i) = pow(scal_prod,mu(i))/pow(normS_i,3+mu(i));
        da_filed.col(i) = mu(i)*pow(scal_prod,mu(i)-1)%(XYZ.col(0)*Dir(i,0)+XYZ.col(1)*Dir(i,1)+XYZ.col(2)*Dir(i,2))/pow(normS_i,mu(i)+3)-(mu(i)+3)*pow(scal_prod,mu(i))%(-T_filed.slice(0).col(i)%XYZ.col(0)-T_filed.slice(1).col(i)%XYZ.col(1)-T_filed.slice(2).col(i)%XYZ.col(2))/pow(normS_i,mu(i)+5);
        // Final lighting filed
        T_filed.slice(0).col(i) = T_filed.slice(0).col(i)%a_filed.col(i);
        T_filed.slice(1).col(i) = T_filed.slice(1).col(i)%a_filed.col(i);
        T_filed.slice(2).col(i) = T_filed.slice(2).col(i)%a_filed.col(i);
        grad_t.slice(0).col(i) = (-exp_z%u)%(a_filed.col(i)+da_filed.col(i))+S(i,0)*da_filed.col(i);
        grad_t.slice(1).col(i) = (-exp_z%v)%(a_filed.col(i)+da_filed.col(i))+S(i,1)*da_filed.col(i);
        grad_t.slice(2).col(i) = (-exp_z  )%(a_filed.col(i)+da_filed.col(i))+S(i,2)*da_filed.col(i);
    }
    *T = T_filed;
    *gT = grad_t;
}
