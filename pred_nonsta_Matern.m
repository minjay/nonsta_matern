clear

addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

load('data_EOF_regr_new.mat')
% use new beta_hat from fit_nonsta_Matern_good_init.m
load('beta_hat_good_init.mat')
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

% non-stationary variance function
% the first column is all ones
load('ns.mat')
b_mat_full = kron(b_mat, ones(size(theta, 1), 1));
b_mat = b_mat_full(index, :);
b_mat = [ones(n, 1) b_mat];

m = size(b_mat, 2);

std_vec_est = exp(b_mat*reshape(beta_hat(1:m), m, 1))*1e3;
plot(theta_samples, std_vec_est, '.')

n_r = 1e3;
r_vec = linspace(0, 2, n_r);
corr_vec = zeros(n_r, 1);
nu = beta_hat(m+1);
a = beta_hat(m+2);
for i = 1:n_r
    corr_vec(i) = Matern(r_vec(i), nu, a);
end
% convert chordal distance to great-circle distance
r_vec = asin(r_vec/2)*2;
% plot correlation function
plot(r_vec, corr_vec, 'LineWidth', 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta_vec = theta(:);
phi_vec = phi(:);
% stretching
[x, y, z] = trans_coord(theta_vec*4, phi_vec);

% total number of locations
N = length(x);

% get distance matrix
r = get_chordal_dist(x, y, z);

b_mat = [ones(N, 1) b_mat_full];

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat = get_cov_nonsta_Matern(beta, r, b_mat)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(pot_samples/1e3, n, 1);

SigmaP0 = cov_mat(:, index);
Y_pred_Matern = SigmaP0*tmp*1e3;

Y_err_Matern = resid'-Y_pred_Matern;

save('Y_pred_Matern.mat', 'Y_pred_Matern', 'Y_err_Matern')