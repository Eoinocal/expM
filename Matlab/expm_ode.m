
A = rand(5);

options = odeset('RelTol',1e-11,'AbsTol',1e-11, 'stats','on');

f = @(t, x) A*x;

[w, h] = size(A);
eA = zeros(size(A));

for j = 1:w
    x0 = zeros(h,1);
    x0(j) = 1;
    [t, y] = ode15s(f, [0 1], x0, options);
    eA(:,j) = y(end,:);
end

eA
expm(A)

norm(eA - expm(A), inf)

