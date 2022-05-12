%% Параметры моделируемого сигнала
clear variables;
close all;
T = 3;			% желаемая длительность генерируемого сигнала
fd = 44100;		% частота дискретизации
f = 5000;		% частота генерируемого сигнала

tt = linspace(0, T, T * fd);	        % создание вектора времени
N = 2 ^ (nextpow2(length(tt)) - 1);     % оптимизация количества элементов (выравнивание до ближайшей степени двойки)
tt = tt(1:N);
ff = linspace(0, fd, N);				% расчёт частот гармоник

%% Моделирование сигнала

signal = sin(2*pi*f*tt);	% создание вектора сигнала
%% Моделирование шума
rng('default');

noise = randn(size(tt));
plot(tt,signal,tt,noise)
%% Моделирование суперпозиции сигнала и шума

mix = 0.1 * signal + noise;		% сложение 0.1 * сигнал + шум
%% разложение сигнала в спектр

MIX = fft(mix);
plot(ff,abs(MIX))
%% подавление шума

MIX_filtered = MIX;
MIX_filtered(abs(MIX_filtered)<2000) = 0;
plot(ff,abs(MIX_filtered))
%% Восстановление сигнала

mix_filtered = ifft(MIX_filtered);
plot(tt,mix_filtered)
%% Визуализация
subplot(6,1,1)
plot(tt,signal)
xlabel('Время,сек')
title('Исходный сигнал')

subplot(6,1,2)
plot(tt,noise)
xlabel('Время,сек')
title('Шум приёмника')

subplot(6,1,3)
plot(tt,mix)
xlabel('Время,сек')
title('Зашумленный сигнал')

subplot(6,1,4)
plot(ff,abs(MIX))
xlabel('Частота,Гц')
title('Спектр зашумленного сигнала')

subplot(6,1,5)
plot(ff,abs(MIX_filtered))
xlabel('Частота,Гц')
title('Спектр сиггнала после фильтрации')

subplot(6,1,6)
plot(tt,mix_filtered)
xlabel('Время,сек')
title('Восстановленный после фильтрациии сигнал')