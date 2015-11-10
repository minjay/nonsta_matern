function f_value = negloglik_needlet(beta_all, b_mat, samples, corr_mat)

n = length(samples);
eta = beta_all(1:end-1);
tau = beta_all(end);
std_vec = exp(b_mat*eta);

cov_mat = diag(std_vec)*corr_mat*diag(std_vec)+eye(n)*tau^2;

f_value = ecmnobj(samples', zeros(n, 1), cov_mat);

end