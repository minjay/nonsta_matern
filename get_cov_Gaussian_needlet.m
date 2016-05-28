function cov_mat = get_cov_Gaussian_needlet(beta, r, b_mat)

% number of basis functions
m = size(b_mat, 2)-1;
eta = beta(1:m+1);
nu = beta(m+2);
a = beta(m+3);
std_vec = exp(b_mat*reshape(eta, m+1, 1));
n = length(std_vec);

cov_mat = zeros(n);
for j = 1:n
    for i = 1:j-1
        value = nonsta_Matern(std_vec(i), std_vec(j), r(i, j), nu, a);
        cov_mat(i, j) = value;
        cov_mat(j, i) = value;
    end
    cov_mat(j, j) = nonsta_Matern(std_vec(j), std_vec(j), r(j, j), nu, a);
end

end