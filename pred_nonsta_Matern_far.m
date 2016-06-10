clear

load('data_EOF_regr_new.mat')
load('beta_hat.mat')
% load precomputed cov mat
load('cov_mat_Matern.mat')
resid = resid_all(1, :);

rng(1)

% sampling
n = 1000;
theta_vec = theta(:);
phi_vec = phi(:);
w = sin(theta_vec*4);
w(phi_vec<=pi/2 & theta_vec<=15/180*pi) = 0;
[pot_samples, index] = datasample(resid', n, 'Replace', false,...
    'Weights', w);
theta_samples = theta_vec(index);
phi_samples = phi_vec(index);

% kriging
Sigma00 = cov_mat_Matern(index, index);
tmp = Sigma00\reshape(pot_samples/1e3, n, 1);

SigmaP0 = cov_mat_Matern(:, index);
Y_pred_Matern_far = SigmaP0*tmp*1e3;

Y_err_Matern_far = resid'-Y_pred_Matern_far;

figure
subplot(1, 2, 1)
plot_pot(reshape(Y_pred_Matern_far, size(phi)), phi, theta, 1000, max(abs(Y_pred_Matern_far)))
colormap(jet)
subplot(1, 2, 2)
plot_pot_with_obs(reshape(Y_err_Matern_far, size(phi)), phi, theta, phi_samples, theta_samples, 1000)
colormap(jet)

save('Y_pred_Matern_far.mat', 'Y_pred_Matern_far', 'Y_err_Matern_far')