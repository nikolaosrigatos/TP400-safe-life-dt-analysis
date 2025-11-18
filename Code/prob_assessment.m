function Pf = prob_assessment( Nb, Sc, num_simulations, num_missions, a0, a1, a2, logsp_avg )


u_c = unifrnd(0, 1, 1, num_simulations);
c_par = 10.^norminv(u_c, a0, a2);
m_par = -a1;

Nf_mission = zeros(length(Nb), num_simulations);
damage_mission = zeros(1, num_simulations);
accumulated_damage = zeros(1, num_simulations);
Pf = zeros(num_missions, 1);

for i = 1:num_simulations

        Nf_mission(:, i) = (c_par(i) ./ (Sc / 10^logsp_avg).^m_par); 
        damage_mission(i) = sum( Nb' ./ Nf_mission(:, i) );

end

C = 10^a0;

result = 0;
for i=1:length(Nb)
    result = Nb(i)*(Sc(i)./10^logsp_avg).^m_par + result;
end

dmg = 0:1e-2:1.0; 
Pf_an = zeros(num_missions,1);

for i = 1:num_missions
    
    % Monte Carlo
    accumulated_damage = accumulated_damage + damage_mission;
    Nf = sum ( accumulated_damage >= 1 );
    Pf(i) = Nf / num_simulations;

    % Analytical Solution
    loading = i*result;
    mean = loading/C;
    Damage = 1 - logncdf(dmg,log10(mean),a2);
    Pf_an(i) = Damage(end);

end

figure()
hold on
box on
grid on
plot(1:num_missions, Pf);
% plot(1:num_missions,Pf1)
plot(1:num_missions, Pf_an);
xlabel('Number of missions')
ylabel('Probability of failure')
title('Probabilistic assessment of damage')
legend('MC','Analytical','Location','Southeast')
hold off

end