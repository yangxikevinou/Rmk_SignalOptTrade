%% Dimensional Analysis
% There are three base quantities in finance, i.e. money, share, and time,
% resembling those three in mechanics, i.e. mass, length, and time. We
% shall perform dimensional analysis on the variables in the optimal
% trading problem, and identify the dimensionless variables. We shall use
% $``[\dot]''$ to denote the dimension of a variable.
% 
% $[T]=[t]=time$
% 
% $[x]=[X]=share$
% 
% $[V]=money$
% 
% $[c]=[P]=[S]=money \cdot share^{-1}$
% 
% $[r]=share \cdot time^{-1}$
% 
% $[\kappa]=money \cdot time \cdot share^{-2}$
% 
% $[\phi]=moeny \cdot time^{-1} \cdot share^{-2}$
% 
% $[\rho]=money \cdot share^{-2}$
% 
% $[\sigma]=time^{-1/2}$
% 
% In this model, we shall use $T,x$ and $P_0=S_0$ as the basis variables and
% express the rest in terms of dimensionless quantities.
% 
% $\tilde{V}=\frac{V}{xP_0}$
% 
% $\tilde{c}=\frac{c}{P_0}, \tilde{P}=\frac{P}{P_0}, \tilde{S}=\frac{S}{S_0}$
% 
% $\tilde{r}=\frac{r}{x/T}$
% 
% $\tilde{\kappa}=\frac{\kappa}{P_0 T/x}$
% 
% $\tilde{\phi}=\frac{\phi}{P_0 /(xT)}$
% 
% $\tilde{\rho}=\frac{\rho}{P_0 / x}$
% 
% $\tilde{\sigma}=\sigma T^{1/2}$
% 
% In the arithmetic Brownian motion model, we introduce the dimensionless 
% ``moneyness'' variable $m=\frac{c-P_t}{S_0}$, while in the geometric
% Brownian motion model, we instead use the dimensionless ``log-moneyness''
% variable $l=\log(1+\frac{c-P_t}{S_t})$. Both variables are bounded below
% by zero, and they are exactly zero iff the asset price hits the cap.


%% Scripts

%%%
%%% Inputs
%%%
% base quantities
T     = 1; % longest trading horizon
x     = 1; % position
S0    = 1; % initial price
% other variables
% c     = 100; % price cap
sigma = 0.5; % volatility parameter
kappa = 0.1; % temporary price impact
phi   = 1; % running inventory penalty
rho   = 1; % terminal inventory penalty

% theta of lookback call
% arithmetic Brownian motion
% theta=@(m,u) S0*sigma^2*normpdf(m,0,sigma.*sqrt(u));
% geometric Brownian motion
theta=@(m,u) S0*sigma^2*(normpdf(m,sigma.^2.*u/2,sigma.*sqrt(u))+normcdf(m,sigma.^2.*u/2,sigma.*sqrt(u),'upper'));


%%%
%%% Intermediate variables
%%%
beta  = sqrt(phi/kappa);
G=@(u) beta.*cosh(beta.*u)+rho./kappa.*sinh(beta.*u);
dG=@(u) phi./kappa.*sinh(beta.*u)+beta.*rho./kappa.*cosh(beta.*u);

%%%
%%% Grid
%%%
t=linspace(0.01,T,100); % time remaining, y-axis
m=linspace(0,3,50); % moneyness, x-axis

% difference between optimal rate and AC rate
Rdiff = zeros(numel(t),numel(m));
for i = 1:numel(t)
    for j = 1:numel(m)
        fun = @(u) G(t(i)-u).*theta(m(j),u);
        Rdiff(i,j) = integral(fun,0,t(i));
    end
end
Rdiff=0.5/kappa*Rdiff;
Rmart=x*dG(t);

figure;
surf(m,t,bsxfun(@rdivide,Rdiff,Rmart'));
title([ ...
    '$x = $', num2str(x), '; ' ...
    '$P_0 = $', num2str(S0), '; ' ...
    '$\tilde{\kappa} = $', num2str(kappa*x/S0/T), '; ' ...
    '$\tilde{\phi} = $', num2str(phi*x*T/S0), '; ' ...
    '$\tilde{\rho} = $', num2str(rho*x/S0), '; ' ...
    '$\tilde{\sigma} = $', num2str(sigma*sqrt(T))
],'Interpreter','latex');
xlabel('(log-)moneyness');
ylabel('time remaining');
zlabel('(optimal rate - AC rate) / AC rate');

% optimal rate
Ropt=bsxfun(@plus,Rdiff,Rmart');
Ropt=bsxfun(@rdivide,Ropt,G(t'));
figure;
surf(m,t,bsxfun(@rdivide,Ropt,x./t'));
title([ ...
    '$x = $', num2str(x), '; ' ...
    '$P_0 = $', num2str(S0), '; ' ...
    '$\tilde{\kappa} = $', num2str(kappa*x/S0/T), '; ' ...
    '$\tilde{\phi} = $', num2str(phi*x*T/S0), '; ' ...
    '$\tilde{\rho} = $', num2str(rho*x/S0), '; ' ...
    '$\tilde{\sigma} = $', num2str(sigma*sqrt(T))
],'Interpreter','latex');
xlabel('(log-)moneyness');
ylabel('time remaining');
zlabel('optimal rate / linear rate');