function a_crit = FAD_config2(crack_geom, mat, a_max)

KIC = mat.KIC * sqrt(1000); 
J = pi*(crack_geom.De^4 - crack_geom.Di^4)/64;
S = crack_geom.Mmax/J * crack_geom.De/2;
Y  = @(a) (sec(pi*a/crack_geom.W))^0.5;
K  = @(Y,a,S) Y * S * sqrt(pi*a);
Lr = @(a,S) S/(1-2*a/crack_geom.W)/mat.SY;

Lr_max = 0.5*(mat.SY + mat.UTS)/mat.SY;

aa = linspace(0,a_max,1000);
Lr_calc = zeros(size(aa));
KJ_calc = zeros(size(aa));
f_Lr = zeros(size(aa));

switch mat.plateau

    case 'yes'
        delta_e = 0.0375 * (1 - mat.SY/1000);
        lambda = 1 + (mat.E*delta_e)/mat.SY;
        N = 0.3* ( 1 - mat.SY/mat.UTS);

        for ii=1:length(aa)

            a = aa(ii);
            Lr_calc(ii) = Lr(a,S);
            Y_calc = Y(a);
            K_calc = K(Y_calc,a,S);
             
            if Lr_calc(ii) >=0 && Lr_calc(ii) <=1
                f_Lr(ii) = (1 + 0.5*Lr_calc(ii)^2)^(-0.5);
            elseif Lr_calc(ii) == 1
                f_Lr(ii) = (lambda + 0.5*lambda)^(-0.5);
            elseif Lr_calc(ii) > 1 && Lr_calc(ii) <= Lr_max
                f_Lr(ii) = (lambda + 0.5*lambda)^(-0.5) *Lr_calc(ii)^((N-1)/(2*N));
            end
    
            KJ_calc(ii) = K_calc /(f_Lr(ii));

        end

    case {'no'}
        mu = min([0.001*mat.E/mat.SY, 0.6]);
        N = 0.3* ( 1 - mat.SY/mat.UTS); 

        for ii=1:length(aa)

            a = aa(ii);
            Lr_calc(ii) = Lr(a,S);
            Y_calc = Y(a);
            K_calc = K(Y_calc,a,S);
             
            if Lr_calc(ii) >=0 && Lr_calc(ii) <=1
                f_Lr(ii) = (1 + 0.5*Lr_calc(ii)^2)^(-0.5) * (0.3 + 0.7*exp(-mu*Lr_calc(ii)^6));
            elseif Lr_calc(ii) > 1 && Lr_calc(ii) <= Lr_max
                f_Lr(ii) = (1 + 0.5)^(-0.5) * (0.3 + 0.7*exp(-mu))*Lr_calc(ii)^((N-1)/(2*N));
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
ylim([0.3, 1.1]);
yticks(0:0.1:1.1)
ylabel('L_r');
ax = gca;
ax.YColor = 'k';
plot(aa,Lr_calc,'k')
yline(Lr_max,'--k')
yyaxis left
xlim([0 a_max])
ylim([0, 6000]);
yticks(0:1000:6000)
ylabel('K_J [MPa*sqrt(mm)]');
ax = gca;
ax.YColor = 'b';
plot(aa,KJ_calc)
yline(KIC,'--b')
xline(a_crit,'--k')
xlabel('Crack length a [mm]');
title('Configuration 2');
hold off

end