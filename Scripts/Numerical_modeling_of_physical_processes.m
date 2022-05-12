%% Свободное падение тела (кинематика)
clc
clear variables;
close all;
dt = 0.1;
t = 0:dt:10;
y0 = 25;
v0 = 30;
g = 9.8;
y = y0 + v0*t - g*t.^2/2;
y_max = max(y)
plot(t,y,'.')
grid on;
set(gca,'XTick',[0:1:10],'YTick',[-200:50:100]);
title({'Зависимость высоты свободно падающего тела от времени',...
        'Максимальная высота y_{max} = 70.911 м'});
xlabel('Время,с');
ylabel('Высота,м');
%% Свободное падение тела (динамика)
clc
clear variables;
close all;
dt = 0.1;
t = 0:dt:10;
v0 = 30;
y0 = 25;
g = 9.8;

v = v0 + ...
cumsum(repmat(-g,size(t))*dt);

y = y0 + ...
cumsum(v * dt);

subplot(1,2,1)
plot(t,v,'.')
grid on;
set(gca,'XTick',[0:1:10],'YTick',[-70:10:30]);
title('Cкорость свободно падающего тела');
xlabel('Время,с');
ylabel('Скорость,м/с');

subplot(1,2,2)
plot(t,y,'.')
grid on;
set(gca,'XTick',[0:1:10],'YTick',[-200:50:100]);
title('Высота свободно падающего тела');
xlabel('Время,с');
ylabel('Высота,м');
%% Полет аэростата
clc
clear variables;
close all;
m = 0.2;
V = 1;
g = 9.8;
k = 0.01;
dt = 1;
t = 0:dt:600;
v = zeros(size(t));
y = zeros(size(t));

for i=1:length(t)-1
    F_arch = g * V * airDencity(y(i));
    F_resist = -k * v(i);
    F_grav = -m * g;
    v(i + 1) = v(i) + ...
        (F_arch + F_grav + F_resist) / m * dt;
    y(i + 1) = y(i) + v(i + 1) * dt;
end

subplot(1,2,1)
plot(t,y,'.')
grid on;
set(gca,'XTick',[0:100:600]);
title('Высота aэростата');
xlabel('Время,с');
ylabel('Высота,м');

subplot(1,2,2)
plot(t,v,'.')
grid on;
set(gca,'XTick',[0:100:600]);
title('Cкорость aэростата');
xlabel('Время,с');
ylabel('Скорость,м/с');