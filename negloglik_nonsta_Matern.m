function f_value = negloglik_nonsta_Matern(beta_all, r, b_mat, samples)

n = length(samples);
beta = beta_all(1:end-1);
tau = beta_all(end);

cov_mat = get_cov_nonsta_Matern(beta, r, b_mat)+tau^2*eye(n);

[beta tau]

f_value = ecmnobj(reshape(samples, 1, n), zeros(n, 1), cov_mat);

end