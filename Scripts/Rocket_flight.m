% Проверки функций физической модели
test_resistanceForce = resistanceForce([0; 1], 1, linspace(-pi, pi));
test_dissipationMomentum = dissipationMomentum(linspace(-1, 1), 1, 1);

test_y = [0 1 10 100]; % Высота
test_airDensity = airDensity(test_y);

test_m = 1; % Масса
test_g = 9.8; % Ускорение свободного падения
test_gravityForce = gravityForce(test_m, test_g);

test_mu = -20; % Расход массы ракеты
test_vjet = [0; -2e4]; % Вектор скорости потока продуктов горения
test_jetForce = jetForce(test_mu, test_vjet);

test_jetForce_norm = 1; % Длина вектора реактивной силы
test_L = 70; % Длина ракеты
test_beta = linspace(-0.01, 0.01); % Углы отклонения вектора тяги
%test_jetMomentum = jetMomentum(test_jetForce_norm, test_L / 2, test_beta);

% Проверки численной модели в точке старта
test_fcn_u0 = fcn(0, [20e3; 0; 0; 0; 0; 0; 0]); % Проверка сил в покое
test_fcn_u1 = fcn(0, [20e3; 0; 0; 0; 100; 0; 0]); % Проверка сил при движении вверх
test_fcn_u2 = fcn(0, [20e3; 0; 0; 0; 0; 0; 0.1]); % Проверка моментов сил вращении

%% Решение ode45
options = odeset('MaxStep', 1);
m0 = 20000;
u0 = zeros(7,1);
u0(1) = m0;

[t, u] = ode45(@fcn, [0 600], u0, options);

m = u(:,1);
x = u(:,2);
y = u(:,3);
vx = u(:,4);
vy = u(:,5);
alpha = u(:,6);
omega = u(:,7);
v1 = sqrt(vx.^2 + vy.^2);
r = u(:,2:3);
v = u(:,4:5);

%% Визуализация решения plot

subplot(2,3,1)
plot(t,x,t,y);
title('Координаты');
xlabel('t, c');
ylabel('Расстояния x и y, км');
grid on;

subplot(2,3,2)
plot(t,vx,'b-',t,vy,'r-', t, v1, 'k');
title('Компонены скоростей и ее модуль');
xlabel('t, c');
ylabel('Скорости vx,vy,|v|, м/с');
grid on;

subplot(2,3,3)
plot(x,y,'-');
title('Траектория');
xlabel('x, м');
ylabel('y, м');
grid on;

subplot(2,3,4)
plot(t,alpha,'-');
title('Угол отклонения от вертикали');
xlabel('t, c');
ylabel('alpha, rad');
grid on;

subplot(2,3,5)
plot(t,omega,'-');
title('Угловая скорость');
xlabel('t, c');
ylabel('omega, c^-1');
grid on;

subplot(2,3,6)
plot(t,m,'-');
title('Масса');
xlabel('t, c');
ylabel('m, кг');
axis([0 600 7500 20000]);
grid on;

%% Функция физической модели fcn
function dudt = fcn(t,u)
    %вытаскиваем параметры из u
    m = u(1);
    x = u(2);
    y = u(3);
    vx = u(4);
    vy = u(5);
    alpha = u(6);
    omega = u(7);
    %получаем плотность
    rho = airDensity(y);
    %задаем некоторые константы
    L = 70;
    v_jet = -20000; %в м/с
    g = -9.8;
    %устанавливаем зависимости расхода и угла отклонения двигателя
    massFlow(t < 500) = -20;
    massFlow(t >= 500) = 0;
    beta(t<80) = 0;
    beta(t>=80 & t<120) = -0.0001;
    beta(t>=120 & t<200) = 0;
    beta(t>=200 & t<240) = 0.00007;
    beta(t>=240) = 0;
    %собираем скорость в вектор
    v = [vx; vy];
    %обеспечиваем поворот вектора выброса газов
    Fr0 = [0; v_jet*massFlow];
    A = [cos(beta) sin(beta); -sin(beta) cos(beta)];
    S = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
    Fr = A*S*Fr0;
    %угол между вектором скорости и ракетой
    gamma = acos((vx*sin(alpha)+vy*cos(alpha))/sqrt(vx.^2 + vy.^2));
    %формируем вектор результирующих сил
    F = Fr + resistanceForce(v, rho, gamma) + [0; m*g];
    %формируем моменты
    M = -dissipationMomentum(omega, rho, L) + norm(Fr0)*sin(-beta)*L/2;
    I = 1/12*m*L^2;
    dudt = zeros(size(u));
    %задаем производные
    dudt(1) = massFlow;
    dudt(2) = vx;
    dudt(3) = vy;
    dudt(4) = F(1)/m;
    dudt(5) = F(2)/m;
    dudt(6) = omega;
    dudt(7) = M/I;
end

%% Вспомогательные функции
function rho = airDensity(y)
    T = 300;
    R = 8.31;
    M = 0.029;
    g = 9.8;
    rho0 = 1;
    rho = rho0 * exp(-(M*g*y)/(R*T));
end

function F_grav = gravityForce(m, g)
    F_grav = [0;-m * g];
end

function F = resistanceForce(v, rho, gamma)
    gamma(gamma > pi) = 2 * pi - gamma(gamma > pi);

    k = zeros(size(gamma));
    k(gamma <= pi/4) = 1;
    k(gamma > pi/4 & gamma <= 3*pi/4) = 10;
    k(gamma > 3*pi/4) = 2;

    F = -v * norm(v) * k * rho * 10;
    end

    function F_jet = jetForce(mu, vjet)
    F_jet = vjet * mu;
end

function M = dissipationMomentum(omega, rho, L)
    k = 500;
    M = k * omega * rho * L/2;
end