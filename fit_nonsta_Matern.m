parpool(8)

load('data_EOF_regr_new.mat')
resid = resid_all(1, :);

load('post_samples_real_new_resid_j24_pen_B_spline_cubic_2knots_aurora_n4000_sigma_j_sq_0.01_0.0002_nu4_long_init0.mat')
eta_est = mean(post_samples.eta(:, 1501:end), 2);

rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
[pot_samples, index] = datasample(resid', 4000, 'Replace', false,...
    'Weights', w);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

% stretching
[x, y, z] = trans_coord(theta_samples*4, phi_samples);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance function
knots = [0 0 0 0 40/180 80/180 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_samples*4);

b_mat(:, 1) = 1;

m = size(b_mat, 2);

% rescale the observations
Y = pot_samples/1e3;

beta_init = [0 eta_est(2:end)' 2 10 0.1];
negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r, b_mat, Y);

lb = [-Inf eta_est(2:end)' 0 0 1e-3];
ub = [Inf eta_est(2:end)' 10 Inf Inf];

[beta_hat, f_min] = nonsta_Matern_fit(negloglik1, beta_init, lb, ub, true);

save('beta_hat.mat', 'beta_hat')

delete(gcp)
