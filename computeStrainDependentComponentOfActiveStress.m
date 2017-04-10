% strain dependent component of time independent active stress
function f = computeStrainDependentComponentOfActiveStress(epsilon, muscleConstants)
    % unpack muscle constants
    epsilon_m = muscleConstants(1,11);
    b = muscleConstants(1,12);
    
    if (epsilon > epsilon_m), disp('error?'); end
    
    % strain dependent component of time independent active stress
    f = max(0, 1-b*(epsilon-epsilon_m)^2);
end