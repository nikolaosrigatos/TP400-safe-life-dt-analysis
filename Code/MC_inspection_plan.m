function [ Pf_noinsp, Pf_lf, Pf_ld, Pf_eddy, N_cycles ] = MC_inspection_plan(cdf_lf, cdf_ld, cdf_eddy, N_sim, n_missions,num_insp)

aa = exp(linspace(log(0.001*25.4),log(0.1*25.4),20)); 

u_a0 = unifrnd(0, 1, 1, N_sim);
a0 = expinv(u_a0, 1/61.546);
a0 = a0.*25.4;
N_a0 = interp1(aa, n_missions, a0, 'spline', 'extrap');

%%

failures_noinsp = zeros(N_sim,2);
inspection_results_lf = [];
inspection_results_ld = [];
inspection_results_eddy = [];
failures_lf = [];
failures_ld = [];
failures_eddy = [];

k = 1;
l = 1;
m = 1;

    for i=1:N_sim

         if num_insp == 0
             if N_a0(i) < 3000
                  failures_noinsp(i,1) = 1;
                  failures_noinsp(i,2) = N_a0(i);
             else
                  failures_noinsp(i,1) = 0;
                  failures_noinsp(i,2) = N_a0(i);
             end

         else

            avsN_a0 = interpolate_curves(aa,a0(i));
            det_results_lf = PND_MC( avsN_a0, cdf_lf, num_insp, N_a0(i) );
            det_results_ld = PND_MC( avsN_a0, cdf_ld, num_insp, N_a0(i) );
            det_results_eddy = PND_MC( avsN_a0, cdf_eddy, num_insp, N_a0(i) );

            inspection_results_lf = vertcat(inspection_results_lf, det_results_lf);
            inspection_results_ld = vertcat(inspection_results_ld, det_results_ld);
            inspection_results_eddy = vertcat(inspection_results_eddy, det_results_eddy);

            if det_results_lf(end,1) == 0
                failures_lf(k,:) = det_results_lf(end,:);   
                k = k+1;
            end
            if det_results_ld(end,1) == 0
                failures_ld(l,:) = det_results_ld(end,:);   
                l = l+1;
            end
            if det_results_eddy(end,1) == 0
                failures_eddy(m,:) = det_results_eddy(end,:);   
                m = m+1;
            end

         end

    end

Pf_noinsp = length(failures_noinsp(failures_noinsp(:,1)==1))/N_sim;

if num_insp ~= 0
N_cycles = linspace(0,max(N_a0),3000);

failed_before_N_cycles_lf = zeros(length(N_cycles),1);
failed_before_N_cycles_ld = zeros(length(N_cycles),1);
failed_before_N_cycles_eddy = zeros(length(N_cycles),1);
for i=1:length(N_cycles)
failed_before_N_cycles_lf(i) = length(failures_lf(failures_lf(:, 2) <= N_cycles(i)));
failed_before_N_cycles_ld(i) = length(failures_ld(failures_ld(:, 2) <= N_cycles(i)));
failed_before_N_cycles_eddy(i) = length(failures_eddy(failures_eddy(:, 2) <= N_cycles(i)));
end

% Calculate the conditional probability
Pf_lf = failed_before_N_cycles_lf ./ N_sim;
Pf_ld = failed_before_N_cycles_ld ./ N_sim;
Pf_eddy = failed_before_N_cycles_eddy ./ N_sim;

% Pf_series_lf = 1 - (1 - Pf_lf).^(32);
% Pf_series_ld = 1 - (1 - Pf_ld).^(32);
% Pf_series_eddy = 1 - (1 - Pf_eddy).^(32);

end

end






