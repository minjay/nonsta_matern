function f_value = negloglik_Gaussian_needlet(beta_all, b_mat, samples, B,...
    j_min, j_max, theta_samples, phi_samples)

n = length(samples);
beta = beta_all(1:end-1);
tau = beta_all(end);

cov_mat = get_cov_Gaussian_needlet(beta, b_mat, B, j_min, j_max,...
    theta_samples, phi_samples)+tau^2*eye(n);

[beta tau]

f_value = ecmnobj(reshape(samples, 1, n), zeros(n, 1), cov_mat);

end