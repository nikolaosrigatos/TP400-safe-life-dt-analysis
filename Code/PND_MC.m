function det_results = PND_MC( crack_growth_profile, cdf_datapoints, n_inspections, N_a0  )

detection_result = zeros(n_inspections,1);
inspection_interval = crack_growth_profile(end,1)/(n_inspections + 1);
inspections_det = inspection_interval * (1:n_inspections);

inspections_prob = zeros(n_inspections,1);
for i = 1:n_inspections

    u = unifrnd(0, 1, 1, 1);
    inspections_prob(i) = norminv(u, inspections_det(i), inspections_det(i)*0.033);

end


for i = 1:n_inspections

    crack_size_calc = interp1(crack_growth_profile(:, 1), crack_growth_profile(:, 2), inspections_prob(i));
    POD_calc = interp1(cdf_datapoints(:, 1), cdf_datapoints(:, 2), crack_size_calc);

    detection_result = rand() <= POD_calc;

end

det_results = zeros(n_inspections,2);
det_results(:,1) = detection_result;
det_results(:,2) = inspections_prob;
det_results(end,2) = N_a0;

end
