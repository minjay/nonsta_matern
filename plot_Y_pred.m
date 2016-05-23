load('data_EOF_regr_new.mat')
load('Y_pred.mat')
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

plot_pot(reshape(Y_pred, size(phi)), phi, theta, 1000, max(abs(Y_pred)))

Y_err = resid'-Y_pred;
plot_pot_with_obs(reshape(Y_err, size(phi)), phi, theta, phi_samples, theta_samples, 1000)
