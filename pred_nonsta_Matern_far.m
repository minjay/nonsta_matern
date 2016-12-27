clear

load('data_EOF_regr_new.mat')
% load precomputed cov mat
load('cov_mat_Matern.mat')
% use unit kV
resid = resid_all(1, :)'/1e3;

n = 1000;
theta_vec = theta(:);
phi_vec = phi(:);
N = length(theta_vec);

% mimic the sampling design of SuperDARN real data
width = pi/2;
lat_low = 20/180*pi;
R = 20;

MSPE_Matern = zeros(R, 1);
MAE_Matern = zeros(R, 1);
CRPS_Matern = zeros(R, 1);
len_90_Matern = zeros(R, 1);
len_50_Matern = zeros(R, 1);
cp_90_Matern = zeros(R, 1);
cp_50_Matern = zeros(R, 1);

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
    Sigma00 = cov_mat_Matern(index, index);
    tmp = Sigma00\reshape(pot_samples, n, 1);

    index_pred = setdiff(1:N, index);
    SigmaP0 = cov_mat_Matern(:, index);
    Y_pred_Matern = SigmaP0*tmp;
    
    SigmaPP = cov_mat_Matern;
    Sigma0P = SigmaP0';
    Var_Y_pred_Matern = diag(SigmaPP-SigmaP0*(Sigma00\Sigma0P));

    Y_err_Matern = resid(index_pred)-Y_pred_Matern(index_pred);
    
    % get quantiles
    Y_lb_90 = Y_pred_Matern-norminv(0.95)*sqrt(Var_Y_pred_Matern);
    Y_ub_90 = Y_pred_Matern+norminv(0.95)*sqrt(Var_Y_pred_Matern);
    Y_lb_50 =  Y_pred_Matern-norminv(0.75)*sqrt(Var_Y_pred_Matern);
    Y_ub_50 = Y_pred_Matern+norminv(0.75)*sqrt(Var_Y_pred_Matern);
    
    MSPE_Matern(i) = mean(Y_err_Matern.^2);
    MAE_Matern(i) = mean(abs(Y_err_Matern));
    CRPS_Matern(i) = mean(CRPS(resid(index_pred), Y_pred_Matern(index_pred),...
        Var_Y_pred_Matern(index_pred)));
    
    % PI
    len_90_Matern(i) = mean(Y_ub_90(index_pred)-Y_lb_90(index_pred));
    len_50_Matern(i) = mean(Y_ub_50(index_pred)-Y_lb_50(index_pred));
    cp_90 = resid>=Y_lb_90 & resid<=Y_ub_90;
    cp_90_Matern(i) = mean(cp_90(index_pred));
    cp_50 = resid>=Y_lb_50 & resid<=Y_ub_50;
    cp_50_Matern(i) = mean(cp_50(index_pred));
end
