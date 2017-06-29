function cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A)

% number of basis functions
m = size(b_mat, 2)-1;
eta = [0 beta(1:m)];
std_vec = exp(b_mat*reshape(eta, m+1, 1)); 
sigma_sq = beta(m+1:end);

[N, M] = size(A);

DA = zeros(N, M);
for i = 1:N
    DA(i, :) = std_vec(i)*A(i, :);
end

len_j = length(Npix);
st = zeros(len_j, 1);
en = zeros(len_j, 1);
for j = 1:len_j
    st(j) = sum(Npix(1:j))-Npix(j)+1;
    en(j) = sum(Npix(1:j));
end

cov_mat = zeros(N);
for j = 1:len_j
    range = st(j):en(j);
    DA_sub = DA(:, range);
    cov_mat = cov_mat+sigma_sq(j)*(DA_sub*DA_sub');
end

end