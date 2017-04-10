% define constants
classdef C
   properties (Constant)
      

        % Direction cosines and directional moment arms for muscles of CA joint
        % Data from canines
        alpha_IA    = -0.697;
        alpha_LCA   = -0.198;
        alpha_PCA   = -0.1;     
        alpha_TA    =  0.015;
        beta_IA     = -0.644;
        beta_LCA    =  0.886;
        beta_PCA    = -0.8;     
        beta_TA     =  0.990;
        gamma_IA    = -3.30  * 1e-3;
        gamma_LCA   =  3.915 * 1e-3;
        gamma_PCA   = -5.49  * 1e-3;
        gamma_TA    =  0.8   * 1e-3;

        % damping coefficients of CA joint
        dx = 1;%0.001;
        dy = 1;%0.001;
        dr = 1;%0.1;

        % damping coefficents in vocal cords
%         tr = 0.1;
%         tt = 0.1;

        % Constants for CT joint mechanics
        w = 11.1 * 1e-3;   % m    human
        h = 16.1 * 1e-3;   % m    human
        cos_phi = 0.76;    % phi = 41 deg human

        % muscle constants
        % order:
        % L, Ac, M, rho, sigma_0, sigma_2, epsilon_1, epsilon_2, B, sigma_m, epsilon_m, b, epsilon_dot_m, ti, tp, ts
        muscleConstants_IA  = [9.3*1e-3, 24.5*1e-6, 0.121*1e-3, 1.04*1e3, 2000, 30000, -0.5, 0, 3.5, 96000, 0.4, 1.25, 4, 0.01, 0.1, 0.09];
        muscleConstants_LCA = [14.4*1e-3, 11.9*1e-6, 0.314*1e-3, 1.04*1e3, 3000, 59000, -0.5, 0.05, 4, 96000, 0.4, 2.37, 4, 0.01, 0.1, 0.09];
        muscleConstants_PCA = [15*1e-3, 48.1*1e-6, 0.544*1e-3, 1.04*1e3, 5000, 55000, -0.5, 0.1, 5.3, 96000, 0.4, 1.86, 4, 0.01, 0.1, 0.09]; % **** epsilon m = 0.0
        muscleConstants_TA  = [18.3*1e-3, 40.9*1e-6, 0.8232*1e-3, 1.04*1e3, 1000, 1500, -0.5, -0.05, 6.5, 105000, 0.2, 1.07, 6, 0.01, 0.1, 0.06];
        muscleConstants_CT  = [13.8*1e-3, 41.9*1e-6, 0.9423*1e-3, 1.04*1e3, 2200, 5000, -0.5, -0.06, 7, 87000, 0.2, 2.4, 2.2, 0.02, 0.1, 0.08];
        muscleConstants_lig = [16*1e-3, 5*1e-6, 0, 0, 1000, 9000, -0.5, -0.35, 4.4, 0, 0, 0, 0, 0, 0.01, 0.08];
        muscleConstants_muc = [16*1e-3, 6.1*1e-6, 0, 0, 1000, 1400, -0.5, 0, 17, 0, 0, 0, 0, 0, 0.1, 0.08];

        % cadaveric lengths
        Lo =  C.muscleConstants_lig(1,1); % vocal fold cadaveric length
        Lc =  C.muscleConstants_CT(1,1); % CT muscle cadaveric length
        
        % positions
        x_CAJ    =   7.1 * 1e-3; % m, x position of CAJ from origin
        y_CAJ    =  -7.1 * 1e-3; % m, y position of CAJ from origin
        x_V =   4.0 * 1e-3; % m, x position of vocal process from origin in cadaveric position
        y_V =   0.0 * 1e-3; % m, y position of vocal process from origin in cadaveric position
        
        % mass and moment of inertia
        Ir = 113261;        % kg / m2
        Mt = 1.7524 * 1e-3; % kg
        Ia = 682;           % kg / m2
        Ma = 0.0861 * 1e-3; % kg

        % more constants
        xi_bar_02  = C.x_V;
        psi_bar_02 = C.y_V;
   end
end