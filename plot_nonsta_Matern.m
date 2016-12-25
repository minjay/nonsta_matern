load('data_EOF_regr_new.mat')
load('beta_hat_good_init.mat')

theta_vec = theta(:);
phi_vec = phi(:);
% stretching
[x, y, z] = trans_coord(theta_vec*4, phi_vec);

N = length(x);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance function
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(length(theta_vec), 1) b_mat];

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat_Matern = get_cov_nonsta_Matern(beta, r, b_mat)+tau^2*eye(N);

% figure
T = 9;
rng(1)
L = chol(cov_mat_Matern, 'lower');
Z_sim = randn(N, T)*1000;
Y_sim_Matern = (L*Z_sim)';

% for t = 1:T
%     subplot(3, 3, t)
%     plot_pot(reshape(Y_sim_Matern(t, :), size(phi)), phi, theta, 1000, max(abs(Y_sim_Matern(t, :))));
% end

save('Y_sim_Matern.mat', 'Y_sim_Matern')
