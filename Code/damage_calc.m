function [ D, MISSIONS, TH, Nb, Sc ] = damage_calc(profile, mat, a0, a1, a2, logsp_avg, blocks, CV)

Sa_eq_profile = profile.Sa./(1 - profile.Sm/mat.UTS);
nbins_profile = round((max(Sa_eq_profile)-min(Sa_eq_profile))/blocks);
[ Nb, Sc ] = hist( Sa_eq_profile, nbins_profile );
alpha_N = max(profile.N);
Nb = ceil(alpha_N*Nb); 

Sc_3CV = (1 + 3*CV) * Sc; 
Nf = 10.^(a0 + a1*(log10(Sc_3CV) - logsp_avg) - 3*a2);

D = sum(Nb./Nf);
MISSIONS = 1/D;

TH = zeros(length(Sc),3);
TH(:,1) = Nb;
TH(:,2) = Sc;
TH(:,3) = -Sc;

end