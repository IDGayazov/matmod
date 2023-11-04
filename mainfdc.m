clc; clear vars; close all;

bvp1;

n = 11;
x = linspace(a, b, n)';

tic

y = solvebvp(p, r, q, f, bc, x);

time = toc;% секунды

disp(time);

% точное решение
ya = u(x);

figure;
plot(x, ya, "-red", x, y, "-blue*");
xlabel("x");
ylabel("solution");
legend("real solution", "my solution");
title("Решение задачи");
