% unpack state variables into individual components
function [F_IA, F_LCA, F_PCA, F_TA, F_CT, F_lig, F_muc, ...
     sigma_IA, sigma_LCA, sigma_PCA, sigma_TA, sigma_CT, ...
     epsilon_r, epsilon_dot_r, epsilon_t, epsilon_dot_t, ...
     eta_a, eta_dot_a, psi_a, psi_dot_a, theta_a, theta_dot_a] = unpackStateVariables(X)
 
    F_IA = X(1,1);
    F_LCA = X(1,2);
    F_PCA = X(1,3);
    F_TA = X(1,4);
    F_CT = X(1,5);
    F_lig = X(1,6);
    F_muc = X(1,7);
    sigma_IA = X(1,8);
    sigma_LCA = X(1,9);
    sigma_PCA = X(1,10);
    sigma_TA = X(1,11);
    sigma_CT = X(1,12);
    epsilon_r = X(1,13);
    epsilon_dot_r = X(1,14);
    epsilon_t = X(1,15);
    epsilon_dot_t = X(1,16);
    eta_a = X(1,17);
    eta_dot_a = X(1,18);
    psi_a = X(1,19);
    psi_dot_a = X(1,20);
    theta_a = X(1,21);
    theta_dot_a = X(1,22);
end