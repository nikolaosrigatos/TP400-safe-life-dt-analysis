%% SAFE LIFE AND DAMAGE TOLERANCE ANALYSIS OF THE TP400 PROPELLER SHAFT

close all
clear all                                                                                                     %#ok<CLALL>
clc

%% DATA

global fatigue_data RO                                                                                        %#ok<GVMIS>
addpath( 'Fatigue Data', 'Mission Profiles', 'SF_NASGRO', 'POD Curves', 'NASGRO', 'Monte Carlo' )

% 1. FATIGUE DATA

% -- S-N curve                                                                                                              
fatigue_data = xlsread("fatigue_data1.xlsx");                                                                 %#ok<*XLSRD>
S_fatigue_data = fatigue_data(:,1); 
N_fatigue_data = fatigue_data(:,2);
RO = max(N_fatigue_data); % Runout cycles

% 2. MATERIAL DATA

mat = struct();
mat.E = 210e3;   % Young's modulus [MPa]
mat.nu = 0.3;    % Poisson ratio [-]
mat.SY = 972.2;  % Yielding stress [MPa]
mat.UTS = 1103;  % Ultimate tensile strength [MPa]
mat.KIC = 112.1; % Fracture strength [MPa m^0.5]

% 3. ROTATIONAL SPEED (constant) 

rot_speed = 900/60; % Shaft rotational speed [rps]

% 4. MISSION PROFILES DATA

profiles = struct();

% MISSION PROFILE 1
profiles.profile1 = struct();
profiles.profile1.t =  xlsread('profile1.xlsx','A:A');  % [s]
profiles.profile1.Sa = xlsread('profile1.xlsx','B:B');  % [MPa]
profiles.profile1.Sm = xlsread('profile1.xlsx','C:C');  % [MPa]

profiles.profile1.N = zeros(length(profiles.profile1.t),1);
for i=1:length(profiles.profile1.t)-1
    profiles.profile1.N(i) = rot_speed*(profiles.profile1.t(i+1) - profiles.profile1.t(i));
end

% MISSION PROFILE 2
profiles.profile2 = struct();
profiles.profile2.t =  xlsread('profile2.xlsx','A:A');  % [s]
profiles.profile2.Sa = xlsread('profile2.xlsx','B:B');  % [MPa]
profiles.profile2.Sm = xlsread('profile2.xlsx','C:C');  % [MPa]

profiles.profile2.N = zeros(length(profiles.profile2.t),1);
for i=1:length(profiles.profile2.t)-1
    profiles.profile2.N(i) = rot_speed*(profiles.profile2.t(i+1) - profiles.profile2.t(i));
end

% MISSION PROFILE 3
profiles.profile3 = struct();
profiles.profile3.t =  xlsread('profile3.xlsx','A:A');  % [s]
profiles.profile3.Sa = xlsread('profile3.xlsx','B:B');  % [MPa]
profiles.profile3.Sm = xlsread('profile3.xlsx','C:C');  % [MPa]

profiles.profile3.N = zeros(length(profiles.profile3.t),1);
for i=1:length(profiles.profile3.t)-1
    profiles.profile3.N(i) = rot_speed*(profiles.profile3.t(i+1) - profiles.profile3.t(i));
end

% 5. CRACK GEOMETRY
crack_geom = struct();
crack_geom.d = 6;
crack_geom.t = 8;
crack_geom.Di = 180;
crack_geom.De = crack_geom.Di + 2*crack_geom.t;
crack_geom.Dm = (crack_geom.Di + crack_geom.De)/2;
crack_geom.W = pi*crack_geom.Dm/2; 
crack_geom.Mmax = 70e6;
mat.plateau = 'yes';
crack_geom.tab = zeros(11,3);
crack_geom.tab(:,1) = [0;0.251;0.502;1.505;2.257;3.261;4.515;5.518;6.772;7.776;9.032];
crack_geom.tab(:,2) = [0;0.021;0.028;0.054;0.063;0.069;0.074;0.079;0.088;0.1;0.119];
crack_geom.tab(:,3) = [1;0.828;0.733;0.544;0.45;0.364;0.299;0.264;0.23;0.205;0.179];

%% PART 1

% Safe - Life assessment
% 1.1) Fit the SN diagram with a 3 parameters model, considering also the runouts.
% 1.2) Calculate the damage associated to a mission composed only by profile 1. Assume constant scatter for 
% the SN diagram, use the μ-3σ curve. Lump the stress history in blocks of approx 15MPa. Assume a CV of 5%
% for the loads and use the μ+3σ percentile.
% 1.3) Calculate the damage associated to a standard mission (profile 1 + profile 2)
% 1.4) Which is the maximum number of hours of flight for standard mission?
% 1.5) Probabilistic question: Plot a curve of failure probability Pf against life.

%% PART 1 - TASK 1.1

%S-N diagram fitting with a 3 parameters model (also considering the runouts)
x0=[5 -2 std(log(N_fatigue_data))];
x1=fminsearch('loglik_3par_norm_runouts',x0); % Find the parameters
[par,~,~,~,~,H]=fminunc('loglik_3par_norm_runouts',x1); % Find the parameters
a0 = par(1);
a1 = par(2);
a2 = par(3);

logsp_avg = mean(log10(S_fatigue_data));

% Calculate the mean curve
sp50 = linspace(440,1000,200);
p50 = (a0 + a1*(log10(sp50)-logsp_avg));
s50 = 10.^p50;

% Calculate the μ - 3σ curve
sp3s = linspace(440,1000,200);
p3s = p50 - 3*a2;
s3s = 10.^p3s;


fprintf('\n --- PART 1 ---\n\n');
fprintf('\nPART 1 - TASK 1.1\nThe parameters of the model are: a0 = %0.4f, a1 = %0.4f, a2 = %0.4f\n', a0, a1, a2);

% Plot the SN diagram
figure()
P1 = loglog(fatigue_data(:,2),fatigue_data(:,1),'ko');
hold on
grid on
P2 = loglog(s50,sp50,'k','LineWidth',2);
P3 = loglog(s3s,sp3s,'r','LineWidth',2);
xlabel('N [cycles]');
ylabel('S [MPa]');
title('S-N curve');
set(gca,'FontSize',12);
set(gca, 'XScale', 'log', 'YScale', 'log');
legend('Fatigue Data','μ','μ-3σ')
xlim([10^2 1e7])
ylim([420 1080])
hold off

%% PART 1 - TASK 1.2

blocks = 15; 
CV = 5/100;
[ D_profile1, missions_profile1, TH_profile1, Nb1, Sc1 ] = damage_calc(profiles.profile1, mat, a0, a1, a2, logsp_avg, blocks, CV);

fprintf('\nPART 1 - TASK 1.2\nThe damage associated to Profile 1 is: %e\n', D_profile1);

%% PART 1 - TASK 1.3

[ D_profile2, missions_profile2, TH_profile2, Nb2, Sc2 ] = damage_calc(profiles.profile2, mat, a0, a1, a2, logsp_avg, blocks, CV);
D_std_mission = D_profile1 + D_profile2;

fprintf('\nPART 1 - TASK 1.3\nThe damage associated to a standard mission is: %e\n', D_std_mission);

%% PART 1 - TASK 1.4

missions_std_mission = 1/D_std_mission;
t_max_std_missions = missions_std_mission * 4;

fprintf('\nPART 1 - TASK 1.4\nThe maximum hours of flight for standard missions is: %0.4f hours\n', t_max_std_missions);

%% PART 1 - TASK 1.5

Nb_comb = horzcat(Nb1, Nb2);
Sc_comb = horzcat(Sc1, Sc2);

num_simulations = 1e4;

% Profile 1
num_missions1 = 3e5;
Pf_profile1 = prob_assessment( Nb1, Sc1, num_simulations, num_missions1, a0, a1, a2, logsp_avg );

% Standard Mission
num_missions2 = 3.5e3;
Pf_std = prob_assessment( Nb_comb, Sc_comb, num_simulations, num_missions2, a0, a1, a2, logsp_avg );

%% PART 2

% Deterministic damage tolerance analysis
% 2.1) Determine the critical crack size using the EPFM theory for configuration 1 and configuration 2
% 2.2) Find the maximum number of missions for Profile 1. Lump the signal and randomize the sequence
% 2.3) Find the maximum number of standard missions (Profile 1 and Profile 3). Lump the signal and randomize
% the sequence.
% 2.4) If the target life of 3000h is not met, find the number of inspections ensuring a probability compliant
% with regulations (2e-5). Use the POD curves in the Advisory Circular AC 33-70.2 (Appendix 4).

%% PART 2 - TASK 2.1

% CONFIGURATION 1
a_max = 100;
a_crit_config1 = FAD_config1(crack_geom, mat, a_max);

% CONFIGURATION 2 
a_max = 60;
a_crit_config2 = FAD_config2(crack_geom, mat, a_max);

% Maximum admissible crack to use for NASGRO
max_adm_crack = min([a_crit_config1, a_crit_config2]);

fprintf('\n --- PART 2 ---\n\n');
fprintf('\nPART 2 - TASK 2.1\nThe critical crack size for configuration 1 is: %0.4f mm\n', a_crit_config1);
fprintf('The critical crack size for configuration 2 is: %0.4f mm\n', a_crit_config2);
fprintf('The maximum admissible crack is: %0.4f mm\n', max_adm_crack);

%% PART 2 - TASK 2.2

fprintf('\nPART 2 - TASK 2.2\nFile generation: ');
toggle = 'no';
num_permutations = 10;
folder_path = 'D:\Academic Career\Politecnico di Milano\Studies\2nd Year\3rd Semester\Structural Reliability of Aerospace Components\Projects\Turboprop shaft\SF_NASGRO\PROFILE1\';
write_SF_files(TH_profile1, toggle, num_permutations, folder_path)

missions_NASGRO_profile1 = 406 ;
fprintf('\nThe maximum number of missions for Profile 1 is: %i missions\n', missions_NASGRO_profile1);

%% PART 2 - TASK 2.3

fprintf('\nPART 2 - TASK 2.3\nFile generation: ');
TH_profile_comb = vertcat( TH_profile1, TH_profile2 );

toggle = 'no';
folder_path = 'D:\Academic Career\Politecnico di Milano\Studies\2nd Year\3rd Semester\Structural Reliability of Aerospace Components\Projects\Turboprop shaft\SF_NASGRO\PROFILE_COMB\';
write_SF_files(TH_profile_comb, toggle, num_permutations, folder_path)

missions_NASGRO_profile_comb = 150; 
fprintf('\nThe maximum number of missions for a standard mission is: %i missions\n', missions_NASGRO_profile_comb);

%% PART 2 - TASK 2.4 

fprintf('\nPART 2 - TASK 2.4\n');

POD_lf = load('POD_liquid_full.txt');
POD_ld = load('POD_liquid_directed.txt');
POD_eddy = load('POD_eddy.txt');

max_crack_size = max_adm_crack;
cdf_lf = POD_fit( POD_lf, max_crack_size );
cdf_ld = POD_fit( POD_ld, max_crack_size );
cdf_eddy = POD_fit( POD_eddy, max_crack_size );

% PROFILE 1 

crack_growth_profile1 = load('aN_profile1.txt');
crack_growth_profile1(:,1) = crack_growth_profile1(:,1) * 3;

[ n_inspections_lf1, inspections_lf1, PND_lf1 ] = PND_calc( crack_growth_profile1, cdf_lf );
[ n_inspections_ld1, inspections_ld1, PND_ld1 ] = PND_calc( crack_growth_profile1, cdf_ld );
[ n_inspections_eddy1, inspections_eddy1, PND_eddy1 ] = PND_calc( crack_growth_profile1, cdf_eddy );

fprintf('\nPROFILE 1\n');
fprintf(['The number of inspections needed to reach the target PND for Penetrant liquid - full-field is %i ' ...
    'and the inspection interval is %0.4f hours\n'], n_inspections_lf1, inspections_lf1(2)-inspections_lf1(1));
fprintf(['The number of inspections needed to reach the target PND for Penetrant liquid - directed is %i ' ...
    'and the inspection interval is %0.4f hours\n'], n_inspections_ld1, inspections_ld1(2)-inspections_ld1(1));
fprintf(['The number of inspections needed to reach the target PND for Eddy current is %i ' ...
    'and the inspection interval is %0.4f hours\n'], n_inspections_eddy1, inspections_eddy1(2)-inspections_eddy1(1));

max_inspections1 = max( [n_inspections_eddy1, n_inspections_ld1, n_inspections_lf1] );

figure()
box on
hold on
plot(crack_growth_profile1(:,1),crack_growth_profile1(:,2), '-k')
xline(inspections_lf1,'-b')
legend('Crack growth', 'Inspections LF');
title('Profile 1 - Penetrant liquid - full-field')
xlabel('Time (hours)')
ylabel('Crack dimension (mm)')
hold off

figure()
box on
hold on
plot(crack_growth_profile1(:,1),crack_growth_profile1(:,2), '-k')
xline(inspections_ld1,'-r')
legend('Crack growth', 'Inspections LD');
title('Profile 1 - Penetrant liquid - directed')
xlabel('Time (hours)')
ylabel('Crack dimension (mm)')
hold off

figure()
box on
hold on
plot(crack_growth_profile1(:,1),crack_growth_profile1(:,2), '-k')
xline(inspections_eddy1,'-g')
legend('Crack growth', 'Inspections Eddy');
title('Profile 1 - Eddy current')
xlabel('Time (hours)')
ylabel('Crack dimension (mm)')
hold off

figure()
box on
grid on
hold on
plot(1:n_inspections_lf1, PND_lf1, 'o--b')
plot(1:n_inspections_ld1, PND_ld1, 'o--r')
plot(1:n_inspections_eddy1, PND_eddy1, 'o--k')
yline(2e-5)
xticks(1:max_inspections1);
set(gca,'XLim',[0.8 (max_inspections1+0.2)])
set(gca, 'YScale', 'log')
set(gca,'YLim',[10^(-15) 10^1],'YTick',[10^(-15) 10^(-10) 10^(-5) 10^0])
title('Profile 1')
xlabel('Number of inspections')
ylabel('Failure probability (PND)')
legend('Penetrant liquid - full field', 'Penetrant liquid - directed', 'Eddy current', 'Target','Location','Southwest')
hold off

% STANDARD MISSION

crack_growth_std_mission = load('aN_stdmission.txt');
crack_growth_std_mission(:,1) = crack_growth_std_mission(:,1) * 4;

[ n_inspections_lf2, inspections_lf2, PND_lf2 ] = PND_calc( crack_growth_std_mission, cdf_lf );
[ n_inspections_ld2, inspections_ld2, PND_ld2 ] = PND_calc( crack_growth_std_mission, cdf_ld );
[ n_inspections_eddy2, inspections_eddy2, PND_eddy2 ] = PND_calc( crack_growth_std_mission, cdf_eddy );

fprintf('\nSTANDARD MISSION\n');
fprintf(['The number of inspections needed to reach the target PND for Penetrant liquid - full-field is %i ' ...
    'and the inspection interval is %0.4f hours\n'], n_inspections_lf2, inspections_lf2(2)-inspections_lf2(1));
fprintf(['The number of inspections needed to reach the target PND for Penetrant liquid - directed is %i ' ...
    'and the inspection interval is %0.4f hours\n'], n_inspections_ld2, inspections_ld2(2)-inspections_ld2(1));
fprintf(['The number of inspections needed to reach the target PND for Eddy current is %i ' ...
    'and the inspection interval is %0.4f hours\n'], n_inspections_eddy2, inspections_eddy2(2)-inspections_eddy2(1));

max_inspections2 = max( [n_inspections_eddy2, n_inspections_ld2, n_inspections_lf2] );

figure()
box on
hold on
plot(crack_growth_std_mission(:,1),crack_growth_std_mission(:,2), '-k')
xline(inspections_lf2,'-b')
legend('Crack growth', 'Inspections LF');
title('Standard mission - Penetrant liquid - full-field')
xlabel('Time (hours)')
ylabel('Crack dimension (mm)')
hold off

figure()
box on
hold on
plot(crack_growth_std_mission(:,1),crack_growth_std_mission(:,2), '-k')
xline(inspections_ld2,'-r')
legend('Crack growth', 'Inspections LD');
title('Standard mission - Penetrant liquid - directed')
xlabel('Time (hours)')
ylabel('Crack dimension (mm)')
hold off

figure()
box on
hold on
plot(crack_growth_std_mission(:,1),crack_growth_std_mission(:,2), '-k')
xline(inspections_eddy2,'-g')
legend('Crack growth', 'Inspections Eddy');
title('Standard mission - Eddy current')
xlabel('Time (hours)')
ylabel('Crack dimension (mm)')
hold off

figure()
box on
grid on
hold on
plot(1:n_inspections_lf2, PND_lf2, 'o--b')
plot(1:n_inspections_ld2, PND_ld2, 'o--r')
plot(1:n_inspections_eddy2, PND_eddy2, 'o--k')
yline(2e-5)
xticks(1:max_inspections2);
set(gca,'XLim',[0.8 (max_inspections2+0.2)])
set(gca, 'YScale', 'log')
set(gca,'YLim',[10^(-25) 10^1],'YTick',[10^(-20) 10^(-10) 10^(0)])
title('Combined mission')
xlabel('Number of inspections')
ylabel('Failure probability (PND)')
legend('Penetrant liquid - full field', 'Penetrant liquid - directed', 'Eddy current', 'Target','Location','Southwest')
hold off

%% PART 3 

% Probabilistic damage tolerance assesment 
% 3.1) Starting from the anomalies distrubution for circular holes in the AC 33-70.2, perform a 
% probabilistic damage toleranse assessment, considering that we have 8 holes on each of the 4 shafts 
% (weakest link). Which is the life with a probability of being exceeded of 2e-5 considering only Profile 1?
% And for a standard mission?
% 3.2) If the target failure probability is not ensured, apply an inspection plan in a probabilistic 
% framework. Use a Monte Carlo simulation according to the scheme proposed in the AC 33-70.2 (Optional)

%% PART 3 - TASK 3.1

fprintf('\n --- PART 3 ---\n\n');

% Profile 1
n_missions_profile1 = [7811 4854 3089 2047 1451 1092 864 709 599 518 459 414 373 337 304 273 242 210 177 141]*3;
[Pf_WL_profile1, life_profile1, hpf1] = prob_damage_tol( n_missions_profile1, crack_geom);

fprintf('\nPART 3 - TASK 3.1\n');
fprintf('\nPROFILE 1\n');
fprintf('The life with a probability of being exceeded of 2e-5 is: %0.4f hours\n', life_profile1); 

figure()
box on
plot(n_missions_profile1, Pf_WL_profile1)
set(gca, 'YScale', 'log')
hold on
plot(n_missions_profile1, hpf1)
xline(life_profile1,'--','Life@Target');
yline(2e-5,'--','Target');
title('Probabilistic damage tolerance - Profile 1')
xlabel('Time [hours]')
ylabel('Probability of failure')
legend('8 holes - 4 parts - WL', 'Single Hole','Location','Southeast')
grid on
hold off


% Standard Mission
n_missions_std = [483 427 379 336 300 268 240 217 197 179 164 152 140 130 120 111 101 91 79 67]*4;

[Pf_WL_std, life_std, hpf2] = prob_damage_tol( n_missions_std, crack_geom);

fprintf('\nSTANDARD MISSION\n');
fprintf('The life with a probability of being exceeded of 2e-5 is: %0.4f hours\n', life_std);

figure()
box on
plot(n_missions_std, Pf_WL_std)
set(gca, 'YScale', 'log')
hold on
plot(n_missions_std, hpf2)
xline(life_std,'--','Life@Target');
yline(2e-5,'--','Target');
title('Probabilistic damage tolerance - Standard Mission')
xlabel('Time [hours]')
ylabel('Probability of failure')
legend('8 holes - 4 parts - WL', 'Single Hole','Location','Southeast')
grid on
hold off

%% PART 3 - TASK 3.2

% N_sim = 1e5;
% 
% % Profile 1
% max_num_insp1 = 9;
% MC_profile1 = struct('Pf_lf', cell(1, max_num_insp1), 'Pf_ld', cell(1, max_num_insp1), 'Pf_eddy', cell(1, max_num_insp1));
% 
% for i=1:max_num_insp1
% 
% [ Pf_noinsp, Pf_lf1, Pf_ld1, Pf_eddy1, N_cycles ] = MC_inspection_plan( cdf_lf, cdf_ld, cdf_eddy, N_sim, n_missions_profile1, i );
% MC_profile1(i).Pf_lf = Pf_lf1; 
% MC_profile1(i).Pf_ld = Pf_ld1; 
% MC_profile1(i).Pf_eddy = Pf_eddy1; 
% 
% end

% % Standard Mission
% max_num_insp2 = 4;
% Pf_lf2= zeros(max_num_insp1,1);
% Pf_ld2 = zeros(max_num_insp1,1);
% Pf_eddy2 = zeros(max_num_insp1,1);
% 
% for i=1:max_num_insp1
% 
% [ Pf_noinsp, Pf_lf2, Pf_ld2, Pf_eddy2, N_cycles ] = MC_inspection_plan( cdf_lf, cdf_ld, cdf_eddy, N_sim, n_missions_profile_std, i );
% 
% end


