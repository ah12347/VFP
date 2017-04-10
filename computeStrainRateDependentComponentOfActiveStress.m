% strain rate dependent component of time independent active stress
function g = computeStrainRateDependentComponentOfActiveStress(epsilon_dot, muscleConstants)
    % unpack muscle constants
    epsilon_dot_m = muscleConstants(1,13);
    
    if(epsilon_dot > epsilon_dot_m)
        disp('strain rate exceeds max');
        epsilon_dot = epsilon_dot_m;
    end
    
    if epsilon_dot <= 0
        g = max(0, ((epsilon_dot/epsilon_dot_m) + 1) / (1 - (3*epsilon_dot/epsilon_dot_m)));
    else g = ((epsilon_dot/epsilon_dot_m) + (1/9)) / (((5/9)*epsilon_dot/epsilon_dot_m) + (1/9));
    end
end