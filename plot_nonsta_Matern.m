load('data_EOF_regr_new.mat')
load('beta_hat.mat')

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

theta_vec = theta(:);
phi_vec = phi(:);
% stretching
[x, y, z] = trans_coord(theta_vec*4, phi_vec);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec*4);

b_mat(:, 1) = 1;

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat_Matern = get_cov_nonsta_Matern(beta, r, b_mat)+tau^2*eye(N);

save('cov_mat_Matern.mat', 'cov_mat_Matern', '-v7.3')

% figure
% T = 9;
% Y_sim_Matern = mvnrnd(zeros(T, N), cov_mat)*1000;
% rng(1)
% for t = 1:T
%     subplot(3, 3, t)
%     plot_pot(reshape(Y_sim_Matern(t, :), size(phi)), phi, theta, 1000, max(abs(Y_sim_Matern(t, :))));
% end
