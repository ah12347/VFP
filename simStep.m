function [X_next, debugData] = simStep(X, U, dt)
    % simple 4th order Runge-Kutta method <-- check this
    [k1,d1] = continuousDynamics(X,U);
    [k2,d2] = continuousDynamics(X+k1*dt/2,U);
    [k3,d3] = continuousDynamics(X+k2*dt/2,U);
    [k4,d4] = continuousDynamics(X+k3*dt,U);
    
    X_next = X + (k1 + 2*k2 + 2*k3 + k4)*dt/6; 
    debugData = [d1;d2;d3;d4];
end