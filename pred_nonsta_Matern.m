load('data_regr.mat')
load('beta_hat.mat')

rng(1)

% sampling
n = 1e3;
theta_vec = theta(:);
phi_vec = phi(:);
index = rand_sampler(theta_vec*4, phi_vec);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);
pot_samples = resid(index)';

% non-stationary variance funcion
m = 4;
lambda_inv = 2.5;
b_mat = get_nonsta_var(m, lambda_inv, theta_samples*4);

std_vec_est = exp(b_mat*reshape(beta_hat(1:m+1), m+1, 1));
plot(theta_samples, std_vec_est, '.')

n_r = 1e3;
r_vec = linspace(0, 2, n_r);
corr_vec = zeros(n_r, 1);
nu = beta_hat(m+2);
a = beta_hat(m+3);
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

% non-stationary variance funcion
b_mat = get_nonsta_var(m, lambda_inv, theta_vec*4);

beta = beta_hat(1:end-1);
tau = beta_hat(end);

% get cov mat
cov_mat = get_cov_nonsta_Matern(beta, r, b_mat)+tau^2*eye(N);

% kriging
Sigma00 = cov_mat(index, index);
tmp = Sigma00\reshape(pot_samples, n, 1);

SigmaP0 = cov_mat(:, index);
Y_pred = SigmaP0*tmp*1e3;

save('Y_pred.mat', 'Y_pred')