clear
clc

load('pred_Matern_far.mat')
load('pred_Gau_needlet_far.mat')

disp('MAE')
mean([MAE_Gau_needlet MAE_Matern])
std([MAE_Gau_needlet MAE_Matern])

disp('MSPE')
mean([MSPE_Gau_needlet MSPE_Matern])
std([MSPE_Gau_needlet MSPE_Matern])

disp('CRPS')
mean([CRPS_Gau_needlet CRPS_Matern])
std([CRPS_Gau_needlet CRPS_Matern])

disp('cp_50')
round(mean([cp_50_Gau_needlet cp_50_Matern])*100, 1)

disp('len_50')
round(mean([len_50_Gau_needlet len_50_Matern]), 2)

disp('cp_90')
round(mean([cp_90_Gau_needlet cp_90_Matern])*100, 1)

disp('len_90')
round(mean([len_90_Gau_needlet len_90_Matern]), 2)
