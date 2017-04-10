% author: manaswi, keren, hemma
% description: force displacement simulation of arytenoid cartilage
%              based on material from the book
%              Myoelastic Aerodynamic Throry of Phonation
%              by Ingo R. Titze

% NOTE: All units in mks system. angles in radians unless specified

% Tabula rasa
clear all;
close all;
clc;

%% DEBUGGING commands
% dbstop in computeStrainRateDependentComponentOfActiveStress at 6 if epsilon_dot>epsilon_dot_m

%% initialize variables

% assume no strain and only passive stress in cadaveric position.
sigma_IA = computePassiveStressInMuscle(0,C.muscleConstants_IA);
sigma_LCA = computePassiveStressInMuscle(0,C.muscleConstants_LCA);
sigma_PCA = computePassiveStressInMuscle(0,C.muscleConstants_PCA);
sigma_TA = computePassiveStressInMuscle(0,C.muscleConstants_TA);
sigma_CT = computePassiveStressInMuscle(0,C.muscleConstants_CT);

F_IA = sigma_IA * C.muscleConstants_IA(1,2);
F_LCA = sigma_LCA * C.muscleConstants_LCA(1,2);
F_PCA = sigma_PCA * C.muscleConstants_PCA(1,2);
F_TA = sigma_TA * C.muscleConstants_TA(1,2);
F_CT = sigma_CT * C.muscleConstants_CT(1,2);
F_lig = 0;
F_muc = 0;

epsilon_r = 0;
epsilon_dot_r = 0;
epsilon_t = 0;
epsilon_dot_t = 0;
xi_a = 0;
xi_dot_a = 0;
psi_a = 0;
psi_dot_a = 0;
theta_a = 0;
theta_dot_a = 0;

X_state = [F_IA, F_LCA, F_PCA, F_TA, F_CT, F_lig, F_muc, ...
           sigma_IA, sigma_LCA, sigma_PCA, sigma_TA, sigma_CT, ...
           epsilon_r, epsilon_dot_r, epsilon_t, epsilon_dot_t, ...
           xi_a, xi_dot_a, psi_a, psi_dot_a, theta_a, theta_dot_a];

%% define constants
% constants have been moved to the C.m file
% anything that is of the form C.Name is a constant
% import C

%% define time and timestep
tMin = 0;
tMax = 1;
dt = 1e-4;
nSteps = round((tMax - tMin) / dt); % must be a natural number

store = zeros(nSteps,5);
storeState = zeros(nSteps, 23);
storeDebug1 = zeros(nSteps, 64);
storeDebug2 = zeros(nSteps, 64);
storeDebug3 = zeros(nSteps, 64);
storeDebug4 = zeros(nSteps, 64);

drawTableOfMuscleConstants();

f1 = figure(1);
f1.Position = [10 570 560 420];
plot(C.x_CAJ, C.y_CAJ, 'mo');
hold on
plot(C.x_V,C.y_V,'b*')

xlim([-0.02, 0.02]);
ylim([-0.02, 0.02]);
line([-0.05,0.05],[0, 0]);
line([0, 0],[-0.05,0.05]);

% WriteLog('open');

% loop
for count = 1:1:nSteps
    
    time = tMin+count*dt;
    
    % muscle activation
%     U = zeros(1,5);
%     U = [0.9, 0.5, 0, 0.2, 0.2];
    U = [0 0 1 0 0];
%     U = [0 1 0 0 0];
    
    % update state
    [X_state, debugData] = simStep(X_state, U, dt);
    
    % we are interested in 
    xi_V  = C.x_CAJ + cos(X_state(1,21)) * (C.x_V - C.x_CAJ) - sin(X_state(1,21)) * (C.y_V - C.y_CAJ) + X_state(1,17);
    psi_V = C.y_CAJ + sin(X_state(1,21)) * (C.x_V - C.x_CAJ) + cos(X_state(1,21)) * (C.y_V - C.y_CAJ) + X_state(1,19) + C.Lo*(X_state(1,13) + X_state(1,15));
    
    store(count,:) = [time, xi_V, psi_V, X_state(1,17)+C.x_CAJ, X_state(1,19)+C.y_CAJ];
    storeState(count,:) = [time, X_state];
    storeDebug1(count,:) = [time, debugData(1,:)];
    storeDebug2(count,:) = [time, debugData(2,:)];
    storeDebug3(count,:) = [time, debugData(3,:)];
    storeDebug4(count,:) = [time, debugData(4,:)];
    
%     figure(1)
%     plot(xi_V, psi_V, 'b.')
%     hold on
%     plot(X_state(1,17)+C.x_CAJ, X_state(1,19)+C.y_CAJ, 'r.')
end

figure(1) % x y scatter plot of arytenoid and vocal process
plot(store(:,2), store(:,3), 'b.')
hold on
plot(store(:,4), store(:,5), 'r.')
    
f2 = figure(2); % x y coordinates of vocal process with time
f2.Position = [590 570 560 420];
plot(store(:,1),store(:,2))
hold on
plot(store(:,1),store(:,3))

f3 = figure(3); % forces in muscles
f3.Position = [1170 570 560 420];
plot(storeState(:,1), storeState(:,2))
hold on
plot(storeState(:,1), storeState(:,3))
plot(storeState(:,1), storeState(:,4))
plot(storeState(:,1), storeState(:,5))
plot(storeState(:,1), storeState(:,6))
plot(storeState(:,1), storeState(:,7))
plot(storeState(:,1), storeState(:,8))
legend('F_{IA}','F_{LCA}','F_{PCA}','F_{TA}','F_{CT}','F_{lig}','F_{muc}')
title('forces')

f4 = figure(4); % stress in muscles
f4.Position = [10 50 560 420];
plot(storeState(:,1), storeState(:,9))
hold on
plot(storeState(:,1), storeState(:,10))
plot(storeState(:,1), storeState(:,11))
plot(storeState(:,1), storeState(:,12))
plot(storeState(:,1), storeState(:,13))
plot(storeState(:,1), storeState(:,14))
plot(storeState(:,1), storeState(:,15))
legend('\sigma_{IA}','\sigma_{LCA}','\sigma_{PCA}','\sigma_{TA}','\sigma_{CT}','\sigma_{lig}','\sigma_{muc}')
title('stress')

f5 = figure(5); % strain in muscles
f5.Position = [590 50 560 420];
plot(storeDebug1(:,1), storeDebug1(:,2))
hold on
plot(storeDebug1(:,1), storeDebug1(:,3))
plot(storeDebug1(:,1), storeDebug1(:,4))
plot(storeDebug1(:,1), storeDebug1(:,5))
plot(storeDebug1(:,1), storeDebug1(:,6))
plot(storeDebug1(:,1), storeDebug1(:,7))
plot(storeDebug1(:,1), storeDebug1(:,8))
legend('\epsilon_{IA}','\epsilon_{LCA}','\epsilon_{PCA}','\epsilon_{TA}','\epsilon_{CT}','\epsilon_{lig}','\epsilon_{muc}')
title('strain')

f6 = figure(6); % passive stress
f6.Position = [1170 50 560 420];
plot(storeDebug1(:,1), storeDebug1(:,16))
hold on
plot(storeDebug1(:,1), storeDebug1(:,17))
plot(storeDebug1(:,1), storeDebug1(:,18))
plot(storeDebug1(:,1), storeDebug1(:,19))
plot(storeDebug1(:,1), storeDebug1(:,20))
plot(storeDebug1(:,1), storeDebug1(:,21))
plot(storeDebug1(:,1), storeDebug1(:,22))
legend('\sigma_p_{IA}','\sigma_p_{LCA}','\sigma_p_{PCA}','\sigma_p_{TA}','\sigma_p_{CT}','\sigma_p_{lig}','\sigma_p_{muc}')
title('passive stress')

f7 = figure(7); % stress-strain
f7.Position = [-800 50 560 420];
plot(storeDebug1(:,2).*100, storeDebug1(:,16))
hold on
plot(storeDebug1(:,3).*100, storeDebug1(:,17))
plot(storeDebug1(:,4).*100, storeDebug1(:,18))
plot(storeDebug1(:,5).*100, storeDebug1(:,19))
plot(storeDebug1(:,6).*100, storeDebug1(:,20))
plot(storeDebug1(:,7).*100, storeDebug1(:,21))
plot(storeDebug1(:,8).*100, storeDebug1(:,22))
legend('IA','LCA','PCA','TA','CT','lig','muc')
title('passive stress - strain %')