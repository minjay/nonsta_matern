clear

addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

load('data_EOF_regr_new.mat')
load('post_samples_exp2.mat')
resid = resid_all(1, :);

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

% rescale the observations
Y = pot_samples/1e3;

beta_init = [beta_hat(1:r+1) 2 10 1e-2];
negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r_dist, b_mat, Y);

lb = [-10*ones(1, r+1) 0 0 1e-3];
ub = [10*ones(1, r+1) 10 Inf Inf];

[beta_hat, f_min] = nonsta_Matern_fit(negloglik1, beta_init, lb, ub, true);

save('beta_hat_good_init.mat', 'beta_hat')

delete(gcp)
