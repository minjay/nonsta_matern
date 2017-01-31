clear
clc

name = {'1_50', '51_100', '101_120', '121_140', '141_170', '171_200'};
for i = 1:length(name)
    load(['Gau_Matern_bootstrap_', name{i}])
    if i==1
        beta_fit_all_comb = beta_fit_all;
    else
        beta_fit_all_comb = beta_fit_all_comb + beta_fit_all;
    end
end

boxplot(beta_fit_all_comb)

est_std = std(beta_fit_all_comb, 0, 1);