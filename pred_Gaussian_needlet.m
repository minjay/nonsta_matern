clear

load('data_EOF_regr_new.mat')
load('post_samples_exp2.mat')

resid = resid_all(1, :);

beta = beta_hat(1:end-1);
tau = beta_hat(end);

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

% non-stationary variance function
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(N, 1) b_mat];

B = 2;
j_min = 2;
j_max = 4;

% get Npix
[Npix, ~, ~] = get_A_ss(B, j_min, j_max, 0, 0);

% get cov mat

cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A)+tau^2*eye(N);

%%% get pot_samples and index
rng(1)

% sampling
n = 4*1e3;
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
[pot_samples, index] = datasample(resid', n, 'Replace', false,...
    'Weights', w);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(pot_samples/1e3, n, 1);

SigmaP0 = cov_mat(:, index);
Y_pred_Gau_need = SigmaP0*tmp*1e3;

figure
plot_pot(reshape(Y_pred_Gau_need, size(phi)), phi, theta, 1000, max(abs(Y_pred_Gau_need)))
Y_err_Gau_need = resid'-Y_pred_Gau_need;
figure
plot_pot(reshape(Y_err_Gau_need, size(phi)), phi, theta, 1000, max(abs(Y_err_Gau_need)))

save('Y_pred_Gau_need.mat', 'Y_pred_Gau_need', 'Y_err_Gau_need')
