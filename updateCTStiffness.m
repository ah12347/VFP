% calculate non linear stiffness for CT joint
function [kr, kt] = updateCTStiffness(theta, Lo, epsilon_t, epsilon_r,h)
    
    if abs(theta) <= 0.2    %rad
        kr = 0.05 * (1 + 1500*abs(epsilon_t*Lo)) * (1 + 40*theta^2);     % Nm/rad
    else kr = 0.05 * (1 + 1500*abs(epsilon_t*Lo)) * (1 + 40*theta^2 + 2*10^5*(abs(theta) - 0.2)^4);
    end
    
    if abs(Lo*epsilon_t) <= 1.5 * 1e-3 % m
        kt = 500 * (1 + 10^6*Lo^2*epsilon_t^2) * (1 + 3*Lo*abs(epsilon_r)/(pi*h));    % N/m
    else kt = 500 * (1 + 10^6*Lo^2*epsilon_t^2 + 2*10^14*(abs(Lo*epsilon_t) - 0.015)^4)*(1 + 3*Lo*abs(epsilon_r)/(pi*h));
    end
end