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

%
% Parameter sets
%
phi  = 1e-2; % running inventory penalty
plot = 1;    % enable/disable plot

%%%
%%% Inputs
%%%
% base quantities
T     = 5;    % longest trading horizon
x     = 1; % position
S0    = 1; % initial price
% other variables
kappa = 1e-1; % temporary price impact
rho   = phi;  % terminal inventory penalty
sigma = 0.5;  % volatility parameter

% theta of lookback call
% arithmetic Brownian motion
thetaABM=@(m,u) S0*sigma^2*normpdf(m,0,sigma.*sqrt(u));
% geometric Brownian motion
thetaGBM=@(m,u) S0*sigma^2*(normpdf(m,sigma.^2.*u/2,sigma.*sqrt(u))+normcdf(m,sigma.^2.*u/2,sigma.*sqrt(u),'upper'));


%%%
%%% Intermediate variables
%%%
beta  = sqrt(phi/kappa);
G=@(u) beta.*cosh(beta.*u)+rho./kappa.*sinh(beta.*u);
dG=@(u) phi./kappa.*sinh(beta.*u)+beta.*rho./kappa.*cosh(beta.*u);

%%%
%%% Grid
%%%
t=linspace(0.01,T,50); % time remaining, y-axis
m=linspace(0,3,25); % moneyness, x-axis

% difference between optimal rate and AC rate
RdiffABM = zeros(numel(t),numel(m));
RdiffGBM = zeros(numel(t),numel(m));
for i = 1:numel(t)
    for j = 1:numel(m)
        funABM = @(u) G(t(i)-u).*thetaABM(m(j),u);
        RdiffABM(i,j) = integral(funABM,0,t(i));
        funGBM = @(u) G(t(i)-u).*thetaGBM(m(j),u);
        RdiffGBM(i,j) = integral(funGBM,0,t(i));
    end
end
RdiffABM=0.5/kappa*RdiffABM;
RdiffGBM=0.5/kappa*RdiffGBM;
Rmart=x*dG(t);

figure(1);
surf(m,t,bsxfun(@rdivide,RdiffABM,Rmart'));
%title([ ...
%    'Arithmetic Brownian Motion: ', ...
%    '$x = $', num2str(x), '; ' ...
%    '$P_0 = $', num2str(S0), '; ' ...
%    '$\tilde{\kappa} = $', num2str(kappa*x/S0/T), '; ' ...
%    '$\tilde{\phi} = $', num2str(phi*x*T/S0), '; ' ...
%    '$\tilde{\rho} = $', num2str(rho*x/S0), '; ' ...
%    '$\tilde{\sigma} = $', num2str(sigma*sqrt(T))
%],'Interpreter','latex');
title([...
    '$\gamma = ', num2str(phi*x*T/S0), '$ and ' ...
    '$\Gamma = ', num2str(rho*x/S0), '$'
],'interpreter','latex');
xlabel('Moneyness');
ylabel('Remaining Time');
zlabel('Relative Increase');
view([50 40]);
colorbar;

if plot == 1
    cleanfigure;
    matlab2tikz(['figDiffABM_', num2str(phi), '.tikz'], 'height', '\figheight', 'width', '\figwidth', 'extraTikzpictureOptions', 'font=\footnotesize');
end


figure(2);
surf(m,t,bsxfun(@plus,RdiffABM,Rmart'));
%title([ ...
%    'Arithmetic Brownian Motion: ', ...
%    '$x = $', num2str(x), '; ' ...
%    '$P_0 = $', num2str(S0), '; ' ...
%    '$\tilde{\kappa} = $', num2str(kappa*x/S0/T), '; ' ...
%    '$\tilde{\phi} = $', num2str(phi*x*T/S0), '; ' ...
%    '$\tilde{\rho} = $', num2str(rho*x/S0), '; ' ...
%    '$\tilde{\sigma} = $', num2str(sigma*sqrt(T))
%],'Interpreter','latex');
title([...
    '$\gamma = ', num2str(phi*x*T/S0), '$ and ' ...
    '$\Gamma = ', num2str(rho*x/S0), '$'
],'interpreter','latex');
xlabel('Moneyness');
ylabel('Remaining Time');
zlabel('Optimal Rate');
view([50 40]);
colorbar;

if plot == 1
    cleanfigure;
    matlab2tikz(['figABM_', num2str(phi), '.tikz'], 'height', '\figheight', 'width', '\figwidth', 'extraTikzpictureOptions', 'font=\footnotesize');
end

figure(3);
surf(m,t,bsxfun(@rdivide,RdiffGBM,Rmart'));
%title([ ...
%    'Geometric Brownian Motion: ', ...
%    '$x = $', num2str(x), '; ' ...
%    '$P_0 = $', num2str(S0), '; ' ...
%    '$\tilde{\kappa} = $', num2str(kappa*x/S0/T), '; ' ...
%    '$\tilde{\phi} = $', num2str(phi*x*T/S0), '; ' ...
%    '$\tilde{\rho} = $', num2str(rho*x/S0), '; ' ...
%    '$\tilde{\sigma} = $', num2str(sigma*sqrt(T))
%],'Interpreter','latex');
title([...
    '$\gamma = ', num2str(phi*x*T/S0), '$ and ' ...
    '$\Gamma = ', num2str(rho*x/S0), '$'
],'interpreter','latex');
xlabel('Log-Moneyness');
ylabel('Remaining Time');
zlabel('Relative Increase');
view([50 40]);
colorbar;

if plot == 1
    cleanfigure;
    matlab2tikz(['figDiffGBM_', num2str(phi), '.tikz'], 'height', '\figheight', 'width', '\figwidth', 'extraTikzpictureOptions', 'font=\footnotesize');
end


figure(4);
surf(m,t,bsxfun(@plus,RdiffGBM,Rmart'));
%title([...
%    'Geometric Brownian Motion: ' ...
%    '$x = $', num2str(x), '; ' ...
%    '$P_0 = $', num2str(S0), '; ' ...
%    '$\lambda = $', num2str(kappa*x/S0/T), '; ' ...
%    '$\gamma = $', num2str(phi*x*T/S0), '; ' ...
%    '$\Gamma = $', num2str(rho*x/S0), '; ' ...
%    '$\sigma = $', num2str(sigma*sqrt(T))
%],'interpreter','latex');
title([...
    '$\gamma = ', num2str(phi*x*T/S0), '$ and ' ...
    '$\Gamma = ', num2str(rho*x/S0), '$'
],'interpreter','latex');
xlabel('Log-Moneyness','interpreter','latex');
ylabel('Remaining Time','interpreter','latex');
zlabel('Optimal Rate','interpreter','latex');
view([50 40]);
colorbar;

if plot == 1
    cleanfigure;
    matlab2tikz(['figGBM_', num2str(phi), '.tikz'], 'height', '\figheight', 'width', '\figwidth', 'extraTikzpictureOptions', 'font=\footnotesize');
end

% optimal rate
%Ropt=bsxfun(@plus,Rdiff,Rmart');
%Ropt=bsxfun(@rdivide,Ropt,G(t'));
%figure;
%surf(m,t,bsxfun(@rdivide,Ropt,x./t'));
%title([ ...
%    '$x = $', num2str(x), '; ' ...
%    '$P_0 = $', num2str(S0), '; ' ...
%    '$\tilde{\kappa} = $', num2str(kappa*x/S0/T), '; ' ...
%    '$\tilde{\phi} = $', num2str(phi*x*T/S0), '; ' ...
%    '$\tilde{\rho} = $', num2str(rho*x/S0), '; ' ...
%    '$\tilde{\sigma} = $', num2str(sigma*sqrt(T))
%],'Interpreter','latex');
%xlabel('(log-)moneyness');
%ylabel('time remaining');
%zlabel('optimal rate / linear rate');