a = -1; b = 1;

eps = 1;
p = @(x) ones(size(x));
r = @(x) -4*x./(eps+x.^2);
q = @(x) -2./(eps+x.^2);
f = @(x) 0 * ones(size(x));

u = @(x) 1./(eps+x.^2);

bc = [1 0 u(a)
      1 0 u(b)];
