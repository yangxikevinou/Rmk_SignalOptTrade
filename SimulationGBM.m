clear;
hold on;

%%%
%%% Parameters
%%%
kappa = 0.01;
T     = 1;
x     = 100;
c     = 50;
phi   = 1;
rho   = 1;
beta  = sqrt(phi/kappa);

%%%
%%% Grid
%%%
dp = 0.01;
p  =0:c/499:c;

figure(1);
hold on;
r = zeros(size(p));
for i = 1:length(p)
    R1 = (2*kappa*beta*cosh(beta*T) + 2*rho*sinh(beta*T))^(-1);
    R2 = 2*x*(rho*beta*cosh(beta*T) + phi*sinh(beta*T));
    fun = @(u) (beta.*cosh(beta.*(T-u)) + rho.*sinh(beta.*(T-u))./kappa).*(1./sqrt(u).*normpdf(0.5.*sqrt(u)-1./sqrt(u).*log(c./p(i)),0,1) + 0.5.*normcdf(0.5.*sqrt(u)-1./sqrt(u).*log(c./p(i)),0,1));
    R3 = integral(fun,0,T);

    r(i) = R1*(R2+p(i).*R3);
end
plot(p,r);

title([ ...
    'T = ', num2str(T), '; ' ...
    'x = ', num2str(x), '; ' ...
    'c = ', num2str(c), '; ' ...
    'kappa = ', num2str(kappa), '; ' ...
    'phi = ', num2str(phi), '; ' ...
    'rho = ', num2str(rho) ...
]);
axis([min(p) max(p) 0 1.5*max(r)]);
