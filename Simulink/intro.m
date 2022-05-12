m = 1;
l = 1;
g = 9.8;
b = 0.1;
A = 0.1;
w = 73;

startDegree = 170;
sim('Kapitsa_Pendulum');
x1 = x;
y1 = y;
plot(x1,y1);

hold on;

startDegree = 10;
sim('Kapitsa_Pendulum');
x2 = x;
y2 = y;
plot(x2,y2);
