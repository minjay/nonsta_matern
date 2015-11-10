load('data_regr.mat')

% beta_hat = [3.1530, 5.6248, -2.9758, 4.3551, -0.0394, 2.7891, 6.0510, 10.5671];
beta_hat = [0.8937, 6.2331, -3.5112, 5.1321, 0.0806, 4.9487, 19.9509, 11.2629];
std_vec_est = exp(b_mat*beta_hat(1:5)');
plot(theta_samples, std_vec_est, '.')

r_vec = linspace(0, 2, 1000);
value = zeros(length(r_vec), 1);
nu = beta_hat(6);
a = beta_hat(7);
tau = beta_hat(8);
for i = 1:length(r_vec)
    value(i) = Matern(r_vec(i), nu, a);
end
r_vec = asin(r_vec/2)*2;
plot(r_vec, value, 'LineWidth', 2)

theta_vec = theta(:);
phi_vec = phi(:);
[x, y, z] = trans_coord(theta_vec*4, phi_vec);

n = length(x);
r = zeros(n);
for j = 1:n
    for i = 1:j
        s = [x(i); y(i); z(i)];
        t = [x(j); y(j); z(j)];
        r(i ,j) = norm(s-t);
    end
end

% non-stationary variance funcion
m = 4;
mu = pi/(m+1)*(1:m);
lambda = pi/(m+1)*2.5/2;
b_mat = zeros(n, m+1);
b_mat(:, 1) = 1;
for i = 2:m+1
    b_mat(:, i) = exp(-(theta_vec*4-mu(i-1)).^2/2/lambda^2);
end

beta = beta_hat(1:end-1);
tau = beta_hat(end);
cov_mat = get_cov(beta, r, b_mat)+eye(n)*tau^2;

Sigma00 = cov_mat(index, index);
tmp = Sigma00\pot_samples';

SigmaP0 = cov_mat(:, index);
Y_pred = SigmaP0*tmp;