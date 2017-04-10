% calculate non linear elasticity modulus of muscle
function E = computeElasticityModulus(epsilon, muscleConstants)
    
    % unpack muscle constants
    epsilon_1 = muscleConstants(1,7);
    epsilon_2 = muscleConstants(1,8);
    sigma_0 = muscleConstants(1,5);
    sigma_2 = muscleConstants(1,6);
    B = muscleConstants(1,9);
    
    % modulus of elasticity
    if epsilon > epsilon_2
        E = -sigma_0/epsilon_1 + B*sigma_2*(exp(B*(epsilon - epsilon_2)) - 1);
    else E = -sigma_0/epsilon_1;
    end
end