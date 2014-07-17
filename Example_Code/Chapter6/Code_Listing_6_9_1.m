%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates a solution to the nonlinear 1D Serre-Green Nagdhi eqn  %-%
%-% given by the equations:                                             %-%
%-% n_t + ((d+n)u)_x = 0, (1)                                           %-%
%-% q_t + (qu -0.5u^2+gn-0.5(d+n)^2u_x^2)_x = 0, (2)                    %-%
%-% q - u + (1/3)(d+n)^2u_xx + (d+n)n_xu_x = 0. (3)                     %-%
%-% The spatial domain is [-100,100] with zero flux boundary conditions.%-%
%-% The approximation is done by a Gaussian radial basis function pseudo%-%
%-% spectral method with shape parameter 1. Needs gau.m and RHS.m.      %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_9_1()
    clear all;  close all;  clc;
    dt = 1e-1;  Tfin = 3;
    tspan = 0:dt:Tfin;

    n = 400;
    L  = 50;

    x = ( linspace ( -L, L, n ) )';
    Nx = length(x);
    cx = (2)*ones(Nx,1);  %blanket shape parameter

    [Ax,D1x,D2x] = deal(zeros(Nx));
    for j = 1 : Nx
        [Ax(:,j),D1x(:,j),D2x(:,j)] = gau(x,x(j),cx(j));
    end
    D1x = D1x / ( Ax );  D2x = D2x / ( Ax );
    D1x(1,:) = zeros(size(D1x(1,:)));  D1x(end,:) = zeros(size(D1x(1,:)));
    D2x(1,:) = zeros(size(D1x(1,:)));  D2x(end,:) = zeros(size(D1x(1,:)));

    % Various physical paramters
    d = 0.5;
    g = (1/(0.45*sqrt(d)))^2;
    a = 0.025;
    BETA = 1.0 / 3.0;
    c = sqrt(g*(d+a));
    kappa = sqrt(3*a)/(d*sqrt(a+d));

    % Exact solutions
    eta  = @(x,t) a * sech( 0.5*kappa*(x - c * t)).^2;
    u    = @(x,t) c*eta(x,t) ./ (d + eta(x,t));

    %Initial conditions
    eta0  = eta(x,0.0);
    u_x   = (D1x*(u(x,0)));
    u_xx  = (D2x*(u(x,0))); 
    eta_x = D1x*(eta0);
    q0    = u(x,0) - (d+eta0).*((BETA)*(d+eta0).*u_xx + eta_x.*u_x);
    Q = q0;  ETA = eta0;
    init = [ETA; Q];

    %for higher accuracy set RelTol to 2.3e-14
    options = odeset('RelTol',2.3e-11,'AbsTol',eps);
    [t,w] = ode113(@(t,q) RHS(t,q,0,D1x,D2x,g,d,BETA), tspan,init,options);
    W  = w(end,:);  ETA = W(1:n)';

    error3 = norm( ETA - eta(x,Tfin) , inf) 
    error4 = norm( ETA - eta(x,Tfin) , inf)/norm( eta(x,Tfin) , inf)

    figure(1),  plot(x,eta0,x,eta(x,Tfin),'go',x,ETA,'k')
    figure(2),  plot(x,abs(eta(x,Tfin) - ETA))
end

function [phi,phi1,phi2,phi3,phi4] = gau(x,xc,c)
    % 1-D guassian radial basis function
    f = @(r,c) exp(-(c*r).^2);
    r = x - xc;
    phi = f(r,c);
    if nargout > 1
    % 1-st derivative    
        phi1 = -2*r*c^2.*exp(-(c*r).^2);
        if nargout > 2
        % 2-nd derivative
            phi2 = 2*c^2*exp(-c^2*r.^2).*(2*c^2*r.^2 - 1);
            if nargout > 3
            % 3-rd derivative    
                phi3 = -4*c^4*r.*exp(-c^2*r.^2).*(2*c^2*r.^2 - 3);
                if nargout > 4
                % 4-th derivative        
                    phi4 = 4*c^4*exp(-c^2*r.^2).*(4*c^4*r.^4 - 12*c^2*r.^2 + 3);
                end
            end
        end
    end
end

function f = RHS(t,q,dummy,D1x,D2x,g,d,BETA) 

    n  = length(D1x);  I = eye(n);  
    ETA  = q(1:n);  ETA_x = D1x*ETA;
    Q = q(n+1 : 2*n);
    
    %Solve discretized elliptic equation at time level t_n
    L   = ( BETA*diag((d+ETA).^2)*D2x + diag((d+ETA).*ETA_x)*D1x - I );
    U   = L \ (-Q);

    rhs1 = -D1x*((d+ETA) .*U );
    rhs2 = -D1x*(Q.*U - 0.5*(U).^2 + g*ETA - 0.5*(d+ETA).^2.*((D1x*U).^2));
    f = [rhs1;rhs2];
end