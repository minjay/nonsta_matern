parpool(8)

load('data_regr.mat')

rng(1)

% sampling
theta_vec = theta(:);
phi_vec = phi(:);
index = rand_sampler(theta_vec*4, phi_vec);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);
pot_samples = resid(index)';

% stretching
[x, y, z] = trans_coord(theta_samples*4, phi_samples);

% get distance matrix
r = get_chordal_dist(x, y, z);

% non-stationary variance funcion
m = 4;
lambda_inv = 2.5;
b_mat = get_nonsta_var(m, lambda_inv, theta_samples*4);

% rescale the observations
Y = pot_samples/1e3;

beta_init = [zeros(1, m+1) 2 1 0.1];
negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r, b_mat, Y);

lb = [-Inf -Inf -Inf -Inf -Inf 0 0 0];
ub = [Inf Inf Inf Inf Inf 5 Inf Inf];

[beta_hat, f_min] = nonsta_Matern_fit(negloglik1, beta_init, lb, ub, true);

save('beta_hat.mat', 'beta_hat')

delete(gcp)
