clear

load('data_EOF_regr_new.mat')
load('beta_hat_needlet.mat')

resid = resid_all(1, :);

beta = beta_hat(1:end-1);
tau = beta_hat(end);

theta_vec = theta(:);
phi_vec = phi(:);

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_vec*4);

b_mat(:, 1) = 1;

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
w(phi_vec<=pi/2) = 0;
[pot_samples, index] = datasample(resid', n, 'Replace', false,...
    'Weights', w);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(pot_samples/1e3, n, 1);

SigmaP0 = cov_mat(:, index);
Y_pred_Gau_need_far = SigmaP0*tmp*1e3;

figure
plot_pot(reshape(Y_pred_Gau_need_far, size(phi)), phi, theta, 1000, max(abs(Y_pred_Gau_need_far)))
Y_err_Gau_need_far = resid'-Y_pred_Gau_need_far;
figure
plot_pot(reshape(Y_err_Gau_need_far, size(phi)), phi, theta, 1000, max(abs(Y_err_Gau_need_far)))

save('Y_pred_Gau_need_far.mat', 'Y_pred_Gau_need_far', 'Y_err_Gau_need_far')
