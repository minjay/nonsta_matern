clear
clc

addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

load('data_EOF_regr_new.mat')
resid = resid_all(1, :);
load('beta_hat_good_init.mat')

rng(1)

% sampling
n = 4000;
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
[pot_samples, index] = datasample(resid', n, 'Replace', false,...
    'Weights', w);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

% stretching
[x, y, z] = trans_coord(theta_samples*4, phi_samples);

% get distance matrix
r_dist = get_chordal_dist(x, y, z);

% non-stationary variance function
% the first column is all ones
load('ns.mat')
b_mat_full = kron(b_mat, ones(size(theta, 1), 1));
b_mat = b_mat_full(index, :);
b_mat = [ones(n, 1) b_mat];

r = size(b_mat, 2)-1;
lb = [-10*ones(1, r+1) 0 0 1e-3];
ub = [10*ones(1, r+1) 10 Inf Inf];

R = 200;
beta_fit_all = zeros(R, length(beta_hat));

beta = beta_hat(1:end-1);
tau = beta_hat(end);

cov_mat = get_cov_nonsta_Matern(beta, r_dist, b_mat)+tau^2*eye(n);
Y = mvnrnd(zeros(1, n), cov_mat, R);

save('Y_Gau_Matern_bootstrap.mat', 'Y', 'r_dist', 'b_mat', 'beta_hat', 'lb', 'ub')
