parpool(20)

load('Y_Gau_Matern_bootstrap.mat')
name = '1_50';
range = 1:50;
R = 200;
beta_fit_all = zeros(R, length(beta_hat));

parfor i = range
    negloglik1 = @(beta_all) negloglik_nonsta_Matern(beta_all, r_dist, b_mat, Y(i, :));
    [beta_fit, f_min] = nonsta_Matern_fit(negloglik1, beta_hat, lb, ub, false);
    beta_fit_all(i, :) = beta_fit;
end

save(['Gau_Matern_bootstrap_', name, '.mat'], 'beta_fit_all', 'beta_hat')

delete(gcp)
