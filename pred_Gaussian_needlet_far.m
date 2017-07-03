clear

addpath(genpath('/home/minjay/NeedMat'))
addpath(genpath('/home/minjay/overcomplete'))
addpath(genpath('/home/minjay/div_curl'))
addpath(genpath('/home/minjay/model_output'))
addpath(genpath('/home/minjay/nonsta_matern'))
addpath(genpath('/home/minjay/bspline'))

load('data_EOF_regr_new.mat')
load('post_samples_real_exp3.mat')
load('post_samples_exp2.mat')

resid = resid_all(1, :)'/1e3;

beta = beta_hat(1:end-1);
tau = beta_hat(end);

load('mat_A.mat')
[N, M] = size(A);
A = A(361:N-360, :);
[N, M] = size(A);

% non-stationary variance function
load('ns.mat')
b_mat = kron(b_mat, ones(size(theta, 1), 1));
b_mat = [ones(N, 1) b_mat];

B = 2;
j_min = 2;
j_max = 4;

% get cov mat
cov_mat = get_cov_Gaussian_needlet(beta, b_mat, Npix, A)+tau^2*eye(N);

n = 1000;
theta_vec = theta(:);
phi_vec = phi(:);

% mimic the sampling design of SuperDARN real data
width = pi/2;
lat_low = 20/180*pi;
R = 100;

MSPE_Gau_needlet = zeros(R, 1);
MAE_Gau_needlet = zeros(R, 1);
CRPS_Gau_needlet = zeros(R, 1);
len_90_Gau_needlet = zeros(R, 1);
len_50_Gau_needlet = zeros(R, 1);
cp_90_Gau_needlet = zeros(R, 1);
cp_50_Gau_needlet = zeros(R, 1);

% set seed
rng(1)

for i = 1:R
    i
    
    % init weight vector w
    w = sin(theta_vec*4);
    % set the region of no data
    w(theta_vec>=lat_low) = 0;
    st = rand*2*pi;
    en = st+pi/2;
    % if part of the interval [st en] is outside of [0, 2*pi)
    if en>=2*pi
        w(phi_vec>=st) = 0;
        w(phi_vec<=en-2*pi) = 0;
    else
        w(phi_vec>=st & phi_vec<=en) = 0;
    end
        
    [pot_samples, index] = datasample(resid, n, 'Replace', false,...
        'Weights', w);
    theta_samples = theta_vec(index);
    phi_samples = phi_vec(index);

    % kriging
    Sigma00 = cov_mat(index, index);
    tmp = Sigma00\reshape(pot_samples, n, 1);

    index_pred = setdiff(1:N, index);
    SigmaP0 = cov_mat(:, index);
    Y_pred_Gau_needlet = SigmaP0*tmp;
    
    SigmaPP = cov_mat;
    Sigma0P = SigmaP0';
    Var_Y_pred_Gau_needlet = diag(SigmaPP-SigmaP0*(Sigma00\Sigma0P));

    Y_err_Gau_needlet = resid(index_pred)-Y_pred_Gau_needlet(index_pred);
    
    % get quantiles
    Y_lb_90 = Y_pred_Gau_needlet-norminv(0.95)*sqrt(Var_Y_pred_Gau_needlet);
    Y_ub_90 = Y_pred_Gau_needlet+norminv(0.95)*sqrt(Var_Y_pred_Gau_needlet);
    Y_lb_50 =  Y_pred_Gau_needlet-norminv(0.75)*sqrt(Var_Y_pred_Gau_needlet);
    Y_ub_50 = Y_pred_Gau_needlet+norminv(0.75)*sqrt(Var_Y_pred_Gau_needlet);
    
    MSPE_Gau_needlet(i) = mean(Y_err_Gau_needlet.^2);
    MAE_Gau_needlet(i) = mean(abs(Y_err_Gau_needlet));
    CRPS_Gau_needlet(i) = mean(CRPS(resid(index_pred), Y_pred_Gau_needlet(index_pred),...
        Var_Y_pred_Gau_needlet(index_pred)));
    
    % PI
    len_90_Gau_needlet(i) = mean(Y_ub_90(index_pred)-Y_lb_90(index_pred));
    len_50_Gau_needlet(i) = mean(Y_ub_50(index_pred)-Y_lb_50(index_pred));
    cp_90 = resid>=Y_lb_90 & resid<=Y_ub_90;
    cp_90_Gau_needlet(i) = mean(cp_90(index_pred));
    cp_50 = resid>=Y_lb_50 & resid<=Y_ub_50;
    cp_50_Gau_needlet(i) = mean(cp_50(index_pred));
end

save('pred_Gau_needlet_far.mat', 'MSPE_Gau_needlet', 'MAE_Gau_needlet', 'CRPS_Gau_needlet',...
    'len_90_Gau_needlet', 'len_50_Gau_needlet', 'cp_90_Gau_needlet', 'cp_50_Gau_needlet')