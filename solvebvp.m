function y = solvebvp(p, r, q, f, bc, x)

x = x(:);

n = numel(x);

% i = 2:n;
% h = [0, x(i) - x(i - 1)];

h = [0, diff(x')];

hp = 0.5*[h(2), (x(3: n - 1)' - x(2:n-2)'), h(n)];

xp = [h(2) / 2, 0.5*(x(2:n)' + x(1:n - 1)')];

% disp(xp);

px = p(xp);
qx = q(xp);
rx = r(xp);

% disp(qx);

a = zeros(n, 1);
b = zeros(n, 1);
c = zeros(n, 1);
F = zeros(n, 1);

if bc(1, 1) == 1
    a(1) = 1;
    c(2) = 0;
    F(1) = bc(1, 3);
else
    a(1) = (px(1) - rx(1) * hp(1))/h(2) + hp(1) * qx(1) / 2 + bc(1, 2);
    c(1) = (-px(1) + rx(1) * hp(1))/h(2) + hp(1) * qx(1)/2;
    F(1) = hp(1) * f(xp(1)) + bc(1, 3);
end
for i = 2:n-1
    a(i) = (px(i + 1) - rx(i + 1) * hp(i)) / h(i+1) ...
        + (px(i) + rx(i) * hp(i)) / h(i+1) ...
        + hp(i) * (qx(i) + qx(i + 1)) / 2;
    b(i - 1) = (-px(i) - rx(i) * hp(i)) / h(i) + hp(i) * qx(i) / 2;
    c(i + 1) = (-px(i+1) + rx(i+1) * hp(i)) / h(i+1) + hp(i) * qx(i+1) / 2;
    F(i) = hp(i) * (f(i+1) + f(i)) / 2;
end
if bc(2, 1) == 1
    a(n) = 1;
    b(n) = 0;
    F(n) = bc(2, 3);
else
    a(n) = (px(n) + rx(n) * hp(n))/h(n) + hp(n) * qx(n) / 2 + bc(1, 2);
    b(n) = (-px(n) - rx(n) * hp(n))/h(n) + hp(n) * qx(n)/2;
    F(n) = hp(n) * f(xp(n)) + bc(1, 3);
end

% disp(F);
% disp(b);
% disp(a);
% disp(c);

A = spdiags([b a c], -1:1, n, n);
% full(A)
y = A \ F;
end