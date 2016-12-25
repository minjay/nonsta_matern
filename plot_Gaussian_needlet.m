load('data_EOF_regr_new.mat')
load('post_samples_exp2.mat')

beta = beta_hat(1:end-1);
tau = beta_hat(end);

resid = resid_all(1, :);

rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
[pot_samples, index] = datasample(resid', 4000, 'Replace', false,...
    'Weights', w);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

% non-stationary variance function
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(length(theta_vec), 1) b_mat];
b_mat = b_mat(index, :);

B = 2;
j_min = 2;
j_max = 4;

[Npix, ~, A] = get_A_ss(B, j_min, j_max, theta_samples*4, phi_samples);

cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A);

figure
plot(theta_samples, sqrt(diag(cov_mat))*1e3, '.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Y_sim_Gau_need = mvnrnd(zeros(T, N), cov_mat)*1000;
% for t = 1:T
%     subplot(3, 3, t)
%     plot_pot(reshape(Y_sim_Gau_need(t, :), size(phi)), phi, theta, 1000, max(abs(Y_sim_Gau_need(t, :))));
% end

save('Y_sim_Gau_need.mat', 'Y_sim_Gau_need')
