% calculate passive stress in muscle
function sigma_passive = computePassiveStressInMuscle(epsilon, muscleConstants)

    % unpack muscle constants
    epsilon_1 = muscleConstants(1,7);
    epsilon_2 = muscleConstants(1,8);
    sigma_0 = muscleConstants(1,5);
    sigma_2 = muscleConstants(1,6);
    B = muscleConstants(1,9);
    
    % non linear passive stress
    if epsilon > epsilon_2
        sigma_passive = (-sigma_0/epsilon_1) * (epsilon - epsilon_1) + sigma_2*(exp(B*(epsilon - epsilon_2)) - 1 - B*(epsilon - epsilon_2));
    else sigma_passive = (-sigma_0/epsilon_1) * (epsilon - epsilon_1);
    end
end