#include "near_ps.h"

near_ps::near_ps(QJsonObject set, int zi)
{
    QJsonValue JsonValur;
    QJsonArray JsonArray;
    nimgs = set.value("NL").toInt();
    /* load S data */
    S = zeros(nimgs,3);
    JsonValur = set.value("S");
    JsonArray = JsonValur.toArray();
    for (int i=0;i<nimgs;i++) {
        for (int j=0;j<3;j++)
            S(i,j) = JsonArray.at(i).toArray().at(j).toDouble();
    }
    //    cout << S;
    /* load Dir data */
    Dir = zeros(nimgs,3);
    JsonValur = set.value("Dir");
    JsonArray = JsonValur.toArray();
    for (int i=0;i<nimgs;i++) {
        for (int j=0;j<3;j++)
            Dir(i,j) = JsonArray.at(i).toArray().at(j).toDouble();
    }
    //    cout << Dir;
    /* load Phi data */
    Phi = zeros<vec>(nimgs);
    JsonValur = set.value("Phi");
    JsonArray = JsonValur.toArray();
    for (int i=0;i<nimgs;i++)
        Phi(i) = JsonArray.at(i).toDouble();
    //    cout << Phi;
    /* load K data */
    K = zeros(3,3);
    JsonValur = set.value("K");
    JsonArray = JsonValur.toArray();
    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++)
            K(i,j) = JsonArray.at(i).toArray().at(j).toDouble();
    }
    //    cout << K;
    /* load mu data */
    mu = zeros<vec>(nimgs);
    JsonValur = set.value("mu");
    JsonArray = JsonValur.toArray();
    for (int i=0;i<nimgs;i++)
        mu(i) = JsonArray.at(i).toDouble();
    //    cout << mu;

    nrows = 3;
    ncols = 3;
    nchannels =1;

    I = field<mat>(nimgs,nchannels);
    I(0,0) = randu<mat>(nrows,ncols);
    I(1,0) = randu<mat>(nrows,ncols);
    I(2,0) = randu<mat>(nrows,ncols);
    I(3,0) = randu<mat>(nrows,ncols);
    I(4,0) = randu<mat>(nrows,ncols);
    I(5,0) = randu<mat>(nrows,ncols);
    I(6,0) = randu<mat>(nrows,ncols);
    I(7,0) = randu<mat>(nrows,ncols);

    z_0 = zi;
}

near_ps::~near_ps()
{

}

void near_ps::init()
{
    QDateTime time = QDateTime::currentDateTime();
    int timeT = time.toTime_t();
    qDebug() << timeT;
    //        nrows = image.rows;
    //        ncols = image.cols;
    //        mask = zeros<imat>(nrows,ncols);
    //        for(int r=0;r<nrows;r++){
    //            for(int c=0;c<ncols;c++){
    //                mask(r,c) = int(image.at<uchar>(r,c))/255;
    //            }
    //        }

    //        mask = zeros<imat>(nrows,ncols);
    mask = ones<imat>(nrows,ncols);
    mask(1,0) = 0;
    mask(2,0) = 0;
    mask(2,2) = 0;

    // intrinsics
    fx = K(0,0);
    fy = K(1,1);
    x0 = K(0,2);
    y0 = K(1,2);

    sp_mat Dx, Dy;
    field<sp_mat> DxDy = make_gradient();
    Dx = DxDy(0);
    Dy = DxDy(1);
    DxDy.clear();

    double max_I = 1;

    // Scaled pixel units
    vec u = regspace<vec>(1,ncols);
    vec v = regspace<vec>(1,nrows);
    mat uu = repmat(u,1,nrows).t();
    mat vv = repmat(v,1,ncols);
    mat u_tilde = uu-x0;
    mat v_tilde = vv-y0;

    // Use a bounding box
    uvec imask = find(mask>0);
    if (imask.n_elem<0)
        return;
    npixel = imask.n_elem;
    umat iimask = ind2sub(size(mask),imask);
    int imin = iimask.row(0).min();
    int imax = iimask.row(0).max();
    int jmin = iimask.row(1).min();
    int jmax = iimask.row(1).max();
    u_tilde = u_tilde.submat(imin,jmin,size(imax-imin+1,jmax-jmin+1));
    v_tilde = v_tilde.submat(imin,jmin,size(imax-imin+1,jmax-jmin+1));
    for (uint i=0;i<I.n_rows;i++) {
        I(i,0) =  I(i,0).submat(imin,jmin,size(imax-imin+1,jmax-jmin+1));
    }
    mask = mask.submat(imin,jmin,size(imax-imin+1,jmax-jmin+1));
    imask = find(mask>0);
    nrows = mask.n_rows;
    ncols = mask.n_cols;

    mat px_rep = repmat(u_tilde(imask),1,nimgs);
    mat py_rep = repmat(v_tilde(imask),1,nimgs);
    sp_mat Dx_rep = repmat(Dx,nimgs,1);
    sp_mat Dy_rep = repmat(Dy,nimgs,1);

    // Vectorize data
    mat t=zeros(0,0);
    for (int i=0;i<nimgs;i++) {
        t = join_rows(I(i)(imask),t);
    }
    I.clear();
    I = field<mat>(nchannels);
    I(0) = t;
    t.clear();

    // Sort images to remove shadows and higlights
    field<umat> W_idx = field<umat>(nchannels);
    for (int ch=0;ch<nchannels;ch++) {
        W_idx(ch) = zeros<umat>(imask.n_elem,nimgs);
        for (uint i=0;i<imask.n_elem;i++) {
            uvec id = sort_index(I(ch).row(i));
            id = id + id*(imask.n_elem-1) + i;
            id = id(regspace<uvec>(2,6));
            W_idx(ch)(id) = ones<uvec>(id.n_elem);
        }
    }

    // Initialize variables
    vec z_tilde(npixel);
    cube rho(nrows,ncols,nchannels);
    mat rho_tilde(npixel,nchannels);

    z0 = zeros(nrows,ncols)*datum::nan;
    z0(imask) = ones<vec>(npixel)*z_0;
    z_tilde = log(z0(imask));
    /// channels need to be set
    rho.slice(0) = ones<mat>(nrows,ncols)/max_I;
    X = z0%u_tilde/fx;
    Y = z0%v_tilde/fy;
    Z = z0;

    vec zx = Dx*z_tilde;
    vec zy = Dy*z_tilde;
    mat Nx = zeros(nrows,ncols);
    mat Ny = zeros(nrows,ncols);
    mat Nz = zeros(nrows,ncols);
    Nx(imask) = fx*zx;
    Ny(imask) = fy*zy;
    Nz(imask) = -u_tilde(imask)%zx-v_tilde(imask)%zy - 1;
    mat dz = sqrt(Nx*Nx+Ny*Ny+Nz*Nz);

    for (int ch=0;ch<nchannels;ch++) {
        mat rho_t = rho.slice(ch)/dz;
        rho_tilde.col(ch) = rho_t(imask);
    }

    // Initilize energy
    double energy = 0;
    cube Tz,grad_Tz;
    auto res = t_fun(z_tilde,u_tilde(imask)/fx,v_tilde(imask)/fy);
    Tz = res(0);
    grad_Tz = res(1);
    res.clear();

    mat psi = shading_fun(z_tilde,Tz,px_rep,py_rep,Dx_rep,Dy_rep);
    psi.reshape(npixel,nimgs);

    for (int i=0;i<nchannels;i++) {
//        energy =

//        Ich = I(:,:,ch);
//                Wch_idx = W_idx(:,:,ch);
//                energy = energy+J_fcn(rho_tilde(:,ch),Wch_idx(:).*psi(:),Wch_idx(:).*Ich(:));
    }

    time = QDateTime::currentDateTime();
    timeT = time.toTime_t();
    qDebug() << timeT;
}

field<sp_mat> near_ps::make_gradient()
{
    imat Omega_padded = zeros<imat>(nrows+2,ncols+2);
    Omega_padded.submat(span(1,nrows),span(1,ncols)) = mask;

    imat Omega0,Omega1,Omega2,Omega3;
    // pixel who have bottom neighbour
    Omega0 = mask%Omega_padded(2,1,size(nrows,ncols));
    // pixel who have top neighbour
    Omega1 = mask%Omega_padded(0,1,size(nrows,ncols));
    // pixel who have right neighbour
    Omega2 = mask%Omega_padded(1,2,size(nrows,ncols));
    // pixel who have left neighbour
    Omega3 = mask%Omega_padded(1,0,size(nrows,ncols));

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
    idx_c = find(Omega0==1);
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
    idx_c = find(Omega1==1 && Omega0==0);
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
    idx_c = find(Omega2==1);
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
    idx_c = find(Omega3==1 && Omega2==0);
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

    field<sp_mat> res(2);
    res(0) = Dvp + Dvm;
    res(1) = Dup + Dum;

    return res;
}

field<cube> near_ps::t_fun(vec z,vec u_tilde,vec v_tilde) {
    int npixel = z.n_elem;

    vec exp_z = exp(z);
    vec X = exp_z%u_tilde;
    vec Y = exp_z%v_tilde;
    vec Z = exp_z;

    cube T_field(npixel,nimgs,3);
    mat a_field(npixel,nimgs);
    mat de_field(npixel,nimgs);

    for (int i=0;i<nimgs;i++) {
        T_field.slice(0).col(i) = S(i,0)-X;
        T_field.slice(1).col(i) = S(i,1)-X;
        T_field.slice(2).col(i) = S(i,2)-X;
        vec normS_i = sqrt(pow(T_field.slice(0).col(i),2)+pow(T_field.slice(1).col(i),2)+pow(T_field.slice(2).col(i),2));
        vec scal_prod = -T_field.slice(0).col(i)*Dir(i,0)-T_field.slice(1).col(i)*Dir(i,1)-T_field.slice(2).col(i)*Dir(i,2);
        a_field.col(i) = pow(scal_prod,mu(i))/pow(normS_i,(3+mu(i)));
        de_field.col(i) = (mu(i)*pow(scal_prod,(mu(i)-1))%(X*Dir(i,0)+Y*Dir(i,1)+Z*Dir(i,2)))/pow(normS_i,(mu(i)+3))-(mu(i)+3)*pow(scal_prod,mu(i))%(-T_field.slice(0).col(i)%X-T_field.slice(1).col(i)%Y-T_field.slice(2).col(i)%Z)/pow(normS_i,(mu(i)+5));
        T_field.slice(0).col(i) = T_field.slice(0).col(i)%a_field.col(i);
        T_field.slice(1).col(i) = T_field.slice(1).col(i)%a_field.col(i);
        T_field.slice(2).col(i) = T_field.slice(2).col(i)%a_field.col(i);
    }

    cube grad_t(npixel,nimgs,3);
    grad_t.slice(0) = repmat(-exp_z%u_tilde,1,nimgs)%(a_field+de_field)+repmat(S.col(0).t(),npixel,1)%de_field;
    grad_t.slice(1) = repmat(-exp_z%v_tilde,1,nimgs)%(a_field+de_field)+repmat(S.col(1).t(),npixel,1)%de_field;
    grad_t.slice(2) = repmat(-exp_z,1,nimgs)%(a_field+de_field)+repmat(S.col(2).t(),npixel,1)%de_field;

    field<cube> res(2);
    res(0) = T_field;
    res(1) = grad_t;

    return res;
}

mat near_ps::J_fun() {

}

mat near_ps::shading_fun(vec z, cube tz, mat px_rep, mat py_rep, sp_mat Dx_rep, sp_mat Dy_rep) {
    vec tx = reshape(fx*tz.slice(0)-px_rep%tz.slice(2),npixel*nimgs,1);
    vec ty = reshape(fy*tz.slice(1)-py_rep%tz.slice(2),npixel*nimgs,1);

    umat location(2,npixel*nimgs);
    location.row(0) = regspace<urowvec>(0,npixel*nimgs-1);
    location.row(1) = regspace<urowvec>(0,npixel*nimgs-1);
    sp_mat stx(location,tx);
    sp_mat sty(location,ty);
    stx = stx*Dx_rep;
    sty = sty*Dy_rep;

    mat r = (stx+sty)*z-reshape(tz.slice(2),npixel*nimgs,1);
    return r;
}

mat near_ps::r_fun(mat rho,mat shadz,mat II,mat W_idx) {
    mat res = rho*Phi;
    res.reshape(nimgs*npixel,1);
    res = W_idx%(res%psi_fun(shadz)-II);
    return res;
    //    res = res%psi_fun(shadz)-II;
}

mat near_ps::psi_fun(mat x) {
    if (shadows) {
        uvec i = find(x>0);
        x(i) = zeros<vec>(i.size());
        return x;
    }
    else
        return x;
}

mat near_ps::chi_fun(mat x) {
    if (shadows) {
        mat t = conv_to<mat>::from(x>=0);
        return t;
    }
    else
        return ones<mat>(x.size());
}

mat near_ps::phi_fun(mat x) {
    switch (estimator) {
    case 0:
        return 0.5*x%x;
    case 1:
        return 0.5*log(1+x%x/(lambda*lambda));
    case 2:
    {
        double thr_norm = 1e-2;
        return (pow(abs(x),lambda)/(lambda*pow(thr_norm,lambda-2)));
    }
    case 3:
        return 0.5*lambda*lambda*x%x/(x%x+lambda*lambda);
    case 4:
        return 0.5*lambda*lambda*(1-exp(-x%x/(lambda*lambda)));
    case 5:
        return (1/6)*(abs(x)>=lambda)%(1-pow((1-x%x/(lambda*lambda)),3))*lambda*lambda+(abs(x)<lambda)*lambda*lambda;
    }

}

mat near_ps::w_fun(mat x) {
    switch (estimator) {
    case 0:
        return ones<mat>(x.size());
    case 1:
        return 1/(lambda*lambda+pow(x,2));
    case 2:
    {
        double thr_norm = 1e-2;
        mat t = thr_norm*ones<mat>(x.size());
        return lambda*(t);
    }
    case 3:
        return pow(lambda,4)/pow((x%x+lambda*lambda),2);
    case 4:
        return exp(-x%x/(lambda*lambda));
    case 5:
        return (abs(x)<=lambda)*pow((1-x%x/(lambda*lambda)),2);
    default:
        break;
    }

}
