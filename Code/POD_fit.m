function cdf_data = POD_fit( POD_data, max_crack_size )

mu = interp1(POD_data(:,2), log10(POD_data(:,1)),0.5);
mu_2sigma = interp1(POD_data(:,2), log10(POD_data(:,1)),0.975);
sigma = (mu_2sigma - mu)/2;
crack_size = linspace(0,max_crack_size,1000);

cdf_normal = cdf('normal', log10(crack_size), mu, sigma);
cdf_data(:,1) = crack_size; 
cdf_data(:,2) = cdf_normal;

end