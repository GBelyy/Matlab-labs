function rho = airDencity(y)
    T = 300;
    R = 8.31;
    M = 0.029;
    g = 9.8;
    rho = 1;
    rho = rho*exp((-M*g*y)/(R*T));
end

