clear
clc
addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

parpool(15)

load('Y_Gau_Matern_bootstrap.mat')
name = '171_200';
range = 171:200;
R = 200;
beta_fit_all = zeros(R, length(beta_hat));

lb = [-5 -2 -2 -2 0 0 1e-3];
ub = [5 2 2 2 10 Inf Inf];

parfor i = range
    negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r_dist, b_mat, Y(i, :));
    [beta_fit, f_min] = nonsta_Matern_fit(negloglik1, beta_hat, lb, ub, false);
    beta_fit_all(i, :) = beta_fit;
end

save(['Gau_Matern_bootstrap_', name, '.mat'], 'beta_fit_all', 'beta_hat')

delete(gcp)
