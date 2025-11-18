function a_crit = FAD_config1(crack_geom, mat, a_max)

KIC = mat.KIC * sqrt(1000); 
J = pi*(crack_geom.De^4 - crack_geom.Di^4)/64;
S = crack_geom.Mmax/J * crack_geom.De/2;
Rm = crack_geom.Dm/2;

alpha = @(a) a/Rm;
S_ref = @(alpha) pi*S/(pi - alpha - 2 * (sin(alpha))^2/(pi - alpha) - sin(2*alpha)/2);
Lr = @(S_ref) S_ref/mat.SY;
lambda_tab = @(a) (12*(1-mat.nu^2))^0.25 * a /sqrt(Rm*crack_geom.t);
K = @(a) S*sqrt(pi*a);

Lr_max = 0.5*(mat.SY + mat.UTS)/mat.SY;


aa = linspace(0,a_max,1000);
Lr_calc = zeros(size(aa));
KJ_calc = zeros(size(aa));
f_Lr = zeros(size(aa));

switch mat.plateau

    case 'yes'
        Deps = 0.0375 * (1 - mat.SY/1000);
        lambda = 1 + (mat.E*Deps)/mat.SY;
        N = 0.3* ( 1 - mat.SY/mat.UTS);

        for ii=1:length(aa)

            a = aa(ii);
            alpha_calc = alpha(a);
            S_ref_calc = S_ref(alpha_calc);
            Lr_calc(ii) = Lr(S_ref_calc);

            lambda_tab_calc = lambda_tab(a);
            M3_calc = interp1(crack_geom.tab(:,1), crack_geom.tab(:,2), lambda_tab_calc, 'spline', 'extrap');
            M4_calc = interp1(crack_geom.tab(:,1), crack_geom.tab(:,3), lambda_tab_calc, 'spline', 'extrap');
            Mb = M3_calc + M4_calc;
            K_calc = Mb * K(a);

            if Lr_calc(ii) <=1
                f_Lr(ii)= 1./sqrt(1+0.5*Lr_calc(ii).^2);
            elseif Lr_calc(ii) <= Lr_max
                f_Lr(ii)= 1/sqrt(lambda+1/(2*lambda)).*Lr_calc(ii).^((N-1)/(2*N));
            else
                f_Lr(ii)= 1/sqrt(lambda+1/(2*lambda)).*Lr_max.^((N-1)/(2*N));
            end
    
            KJ_calc(ii) = K_calc /(f_Lr(ii));

        end

    case {'no'}
        mu = min([0.0001*mat.E/mat.SY, 0.6]);
        N = 0.3* ( 1 - mat.SY/mat.UTS); 

        for ii=1:length(aa)

            a = aa(ii);
            alpha_calc = alpha(a);
            S_ref_calc = S_ref(alpha_calc);
            Lr_calc(ii) = Lr(S_ref_calc);

            lambda_tab_calc = lambda_tab(a);
            M3_calc = interp1(crack_geom.tab(:,1), crack_geom.tab(:,2), lambda_tab_calc, 'spline', 'extrap');
            M4_calc = interp1(crack_geom.tab(:,1), crack_geom.tab(:,3), lambda_tab_calc, 'spline', 'extrap');
            Mb = M3_calc + M4_calc;
            K_calc = Mb * K(a);
             
            if Lr_calc(ii) >=0 && Lr_calc(ii) <=1
                f_Lr(ii) = (1 + 0.5*Lr_calc(ii)^2)^(-0.5) * (0.3 + 0.7*exp(-mu*Lr_calc(ii)^6));
            elseif Lr_calc(ii) > 1 && Lr_calc(ii) <= Lr_max
                f_Lr(ii) = (1 + 0.5)^(-0.5) * (0.3 + 0.7*exp(-mu))*Lr_calc(ii)^((N-1)/(2*N));
            else
                f_Lr(ii)= 1/sqrt(1.5).*(0.3+0.7*exp(-mu)).*Lr_max.^((N-1)/(2*N)); 
            end
    
            KJ_calc(ii) = K_calc/(f_Lr(ii));

        end

end

[~, index1] = min(abs(KIC - KJ_calc));
a_crit1 = aa(index1);

[~, index2] = min(abs(Lr_max - Lr_calc));
a_crit2 = aa(index2);

a_crit = min([a_crit1, a_crit2]);

figure()
grid on
box on 
hold on
yyaxis right
ylim([0.3, 1.3]);
yticks(0.3:0.1:1.3)
ylabel('L_r');
ax = gca;
ax.YColor = 'k';
plot(aa,Lr_calc,'k')
yline(Lr_max,'--k')
xlim([0 a_max])
yyaxis left
ylim([0, 6000]);
yticks(0:1000:6000)
ylabel('K_J [MPa*sqrt(mm)]');
ax = gca;
ax.YColor = 'b';
plot(aa,KJ_calc)
yline(KIC,'--b')
xline(a_crit,'--k')
xlabel('Crack length a [mm]');
title('Configuration 1');
hold off

end