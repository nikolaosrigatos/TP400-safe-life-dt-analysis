function [Pf_series, life_profile, hpf] = prob_damage_tol( n_missions, crack_geom )

POE = @(x) 5.42e-6 .* exp(-61.546.*(x)); % x in inches

aa = exp(linspace(log(0.001*25.4),log(0.1*25.4),20)); % Crack size in mm
aa_inch = aa/25.4;
area = pi*crack_geom.d * crack_geom.t;
area_inch = area/(25.4.^2);

total_number_holes = 32;

POE_calc = POE(aa_inch); 
hpf = POE_calc * area_inch; % Probability of failure for a single hole
Pf_series = 1 - (1 - hpf).^(total_number_holes);
Pf_goal = 2e-5;

life_profile = interp1(Pf_series,n_missions,Pf_goal);

end