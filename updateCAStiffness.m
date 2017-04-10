% calculate non linear stiffness of CA joint
function [kx, ky, Kappa] = updateCAStiffness(xi_a, psi_a, theta_a)
    
    if abs(xi_a) <= 0.003 % m
        kx = 60 * (1 + 2*1e5*xi_a^2);  % N/m
    else kx = 60 * (1 + 2*1e5*xi_a^2 + 2*1e13*(abs(xi_a) - 0.003)^4);
    end

    if abs(psi_a) <= 0.003 % m
        ky = 200 * (1 + 2*1e5*psi_a^2);    % N/m
    else ky = 200 * (1 + 2*1e5*psi_a^2 + 2*1e13*(abs(psi_a) - 0.003)^4);
    end

    if abs(theta_a) < 0.3   % rad
        Kappa = 0.005 * (1 + 10 * theta_a^2);     % Nm/rad
    else Kappa = 0.005 * (1 + 10 * theta_a^2 + 10^4*(abs(theta_a) - 0.3)^4);
    end

end