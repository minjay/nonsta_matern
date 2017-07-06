addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

clear

load('data_EOF_regr_new.mat')
load('post_samples_real_reparam_nu4.mat')
load('Y_sim_need.mat')

beta = beta_hat(1:end-1);
tau = beta_hat(end);

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

theta_vec = theta(:);

% non-stationary variance function
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(length(theta_vec), 1) b_mat];

% get cov mat
cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A)+tau^2*eye(N);

% figure
T = 9;
L = chol(cov_mat, 'lower');
Z_sim = L\Y_sim_Gau_need';
% for t = 1:T
%     subplot(3, 3, t)
%     plot_pot(reshape(Y_sim_Gau_need(t, :), size(phi)), phi, theta, 1000, max(abs(Y_sim_Gau_need(t, :))));
% end

save('Z_sim_Gau_need.mat', 'Z_sim')
