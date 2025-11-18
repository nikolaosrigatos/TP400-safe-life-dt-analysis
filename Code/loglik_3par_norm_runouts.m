function [loglik] = loglik_3par_norm_runouts(xx)
% Function to compute - log-likelihood function for normal distribution.
% Log-likelihood is a function of the three parameters of the distribution,
% contianed in the input vector xx, that has to be minimized
% GN* 2023

global fatigue_data RO % recall the global varible defined in the main script

x1=xx(1);
x2=xx(2);
x3=xx(3);

logsp_avg = mean(log10(fatigue_data(:,1)));

for jj=1:length(fatigue_data)
    
    S(jj)=fatigue_data(jj,1); % stress level for the jj data point
    Nf(jj)=fatigue_data(jj,2); % Nf for the jj data point
    y(jj)=log10(Nf(jj)); % log of the data (normal distribution)
    
    mu=x1+x2*( log10(S(jj))-logsp_avg ); % mu for the given S
    sigma=x3; % sigma constant
    
    z(jj)=(y(jj)-mu)/sigma; % standardized variable
    
    if Nf(jj) >= RO
        % If the data point is a runout (censored), use the survival function
        L(jj) = log(1-normcdf(z(jj))); % contribution of censored data to ML
    elseif Nf(jj) < RO
        L(jj)=-0.5*log(2*pi)-log(sigma)-0.5*z(jj).^2; % contribution of failure to ML
    end
        
end

loglik=-sum(L);

end

