%% Полет аэростата ode45
clc
clear variables;
close all;

u0 =[0;0];
[t,u] = ode45(@fcn,[0 600],u0);
y = u(:,1);
v = u(:,2);

subplot(1,3,1)
plot(t,y,'.')
grid on;
set(gca,'XTick',[0:100:600]);
title('Высота aэростата');
xlabel('Время,с');
ylabel('Высота,м');

subplot(1,3,2)
plot(t,v,'.')
grid on;
set(gca,'XTick',[0:100:600]);
title('Cкорость aэростата');
xlabel('Время,с');
ylabel('Скорость,м/с');

subplot(1,3,3)
plot(t,'.')
grid on;
title('Расчетная сетка');
xlabel('Номер шага моделирования');
ylabel('Шаг по времени,с');

