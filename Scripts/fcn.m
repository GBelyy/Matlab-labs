function dudt = fcn(t, u)

    y = u(1);
    v = u(2);

    m = 0.2;
    V = 1;
    g = 9.8;
    k = 0.02;
    
    dudt = zeros(size(u));
 
    F_arch = g * V * airDencity(y);
    F_resist = -k * airDencity(y)*v;
    F_grav = -m * g;

    dudt(1) = v;
    dudt(2) = (F_arch + F_resist + F_grav)/m;
    
end

