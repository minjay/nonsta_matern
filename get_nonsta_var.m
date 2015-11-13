function b_mat = get_nonsta_var(m, lambda_inv, theta_samples)

n = length(theta_samples);
mu = pi/(m+1)*(1:m);
lambda = pi/(m+1)/lambda_inv;
b_mat = zeros(n, m+1);
b_mat(:, 1) = 1;
for i = 2:m+1
    b_mat(:, i) = exp(-(theta_samples-mu(i-1)).^2/2/lambda^2);
end

end