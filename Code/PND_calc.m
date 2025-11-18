function [ n_inspections, inspections, PND ] = PND_calc( crack_growth_profile, cdf_datapoints )

target_POF = 2e-5;
PND = 1e5;
n_inspections = 0;


while PND > target_POF
    
    n_inspections = n_inspections + 1;
    inspection_interval = crack_growth_profile(end,1)/(n_inspections + 1);
    inspections = inspection_interval * (1:n_inspections);
    
    crack_size_calc = interp1(crack_growth_profile(:,1), crack_growth_profile(:,2), inspections,'spline');
    POD_calc = interp1(cdf_datapoints(:, 1), cdf_datapoints(:, 2), crack_size_calc,'spline') ;

    POND = 1 - POD_calc;
    PND(n_inspections)  = prod(POND);

end
    
end



