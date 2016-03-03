parpool(8)

load('data_EOF_regr.mat')
resid = resid_all(1, :);

rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
index = rand_sampler_real(theta_vec*4);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);
pot_samples = resid(index)';

% stretching
[x, y, z] = trans_coord(theta_samples*4, phi_samples);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance function
knots = [0 0 0 0 0.25 0.5 0.75 1 1 1 1]*pi;
[b_mat, ~] = bspline_basismatrix(4, knots, theta_samples*4);

m = size(b_mat, 2);

% rescale the observations
Y = pot_samples/1e3;

beta_init = [zeros(1, m) 2 10 0.1];
negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r, b_mat, Y);

lb = [-Inf(1, m) 0 0 0];
ub = [Inf(1, m) 5 Inf Inf];

[beta_hat, f_min] = nonsta_Matern_fit(negloglik1, beta_init, lb, ub, true);

save('beta_hat.mat', 'beta_hat')

delete(gcp)
