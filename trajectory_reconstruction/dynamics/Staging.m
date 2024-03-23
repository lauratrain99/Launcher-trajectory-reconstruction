function [m0, m_subR, m_stage, m_str, m_prop, DV_req, DV_stage, ...
    lambda, MR] = Staging(N, Isp, eps, m_pl, h_orbit)                       % Provides all the masses for the different stages

global mu R_earth g0

a = h_orbit+R_earth;                                                        % Semimajor axis [m] 
DV_loss = 2000;                                                             % Delta_V for losses [m/s]
DV_req = sqrt(mu/a)+DV_loss;                                                % Delta_V required [m/s]
C = Isp*g0;                                                                 % Exhaust velocity [m/s]  
eps_mean = mean(eps);                                                       % Mean structural coefficient 
C_mean = mean(C);                                                           % Mean characteristic velocity [m/s]

% Initial Lagrange Multiplier 
u = DV_req/C_mean;                                                          % Guess value
lambda_tot = ((exp(-u/N)-eps_mean)/(1-eps_mean))^N;                         % Total payload ratio 
lambda_i = lambda_tot^(1/N);                                                % Stage payload ratio 
MRi = 1/(eps_mean*(1-lambda_i)+lambda_i);                                   % Stage Mass ratio 

P0 = 1/(MRi*C_mean*eps_mean-C_mean);                                        % Lagrange multiplier first guess


% Finding Lagrange Multiplier 
fun = @(x)-DV_req;   

for j = 1:N
    fun = @(x)fun(x)+C(j)*log((1 + x*C(j))/(C(j)*x*eps(j)));
end

P = fzero(fun,P0);                                                          % Exact Lagrange Multiplier

for i = N:-1:1                                                              % Mass and Payload Ratio for Stages
    MR(i) = (1+P*C(i))/(P*C(i)*eps(i));                                                          
    lambda(i) = (1-eps(i)*MR(i))/((1-eps(i))*MR(i));                         
end


lambda_tot = prod(lambda);                                                  % Total payload ratio
DV_stage = -C.*log(eps.*(1-lambda)+lambda);                                 % Stage DeltaV [m/s]

%% Errors

if abs(sum(DV_stage)-DV_req) > 1E-5                                         % Check final Delta_V
    error('Computational error in staging!')
end

dev2 = -(1+P*C)./(MR.^2)+(eps./(1-eps.*MR)).^2;                       
if any(dev2 < 0)
    error('Not a MIN')    
end

%% Outputs 

m0 = m_pl/lambda_tot;                                                       % Initial mass [kg]

m_subR = zeros(1,N+1);
m_subR(end) = m_pl;

for i = N:-1:1
    m_subR(i) = m_subR(i+1)/lambda(i);                                      % Subrocket mass [kg]
    m_stage(i) = m_subR(i)-m_subR(i+1);                                     % Stage mass [kg]
    m_str(i) = m_stage(i)*eps(i);                                           % Structural mass [kg] 
    m_prop(i) = m_stage(i)-m_str(i);                                        % Fuel or propellant mass [kg] 
end

end