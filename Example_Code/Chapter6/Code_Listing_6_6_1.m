%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates the solution to the BBM equation which is given by     %-%
%-% u_t + u_x + uu_x - u_{xxt} = 0. The spatial domain is [-100,100]    %-%
%-% with no flux boundary conditions. The approximation is done by a    %-%
%-% Gaussian radial basis function spectral method with shape parameter %-%
%-% 0.8. ODE113 is used.  A uniform grid is used for the RBF centers.   %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_6_1()
    clear all; close all; clc;
    n = 300; L1 = 100;
    dt = 1e-1; Tfin = 10;
    tspan = 0 : dt : Tfin; x = linspace ( -L1, L1, n )';
    shape = 0.8; %Shape parameter
    cx = (shape)*ones(n,1); [Ax,D1x,D2x] = deal(zeros(n));
    for j = 1 : n
        [Ax(:,j),D1x(:,j),D2x(:,j)] = gau(x,x(j),cx(j));
    end
    D1x = D1x / ( Ax ); D2x = D2x / ( Ax );
    %Zero flux
    D1x(1,:) = zeros(size(D1x(1,:))); D1x(end,:) = zeros(size(D1x(1,:)));
    D2x(1,:) = zeros(size(D1x(1,:))); D2x(end,:) = zeros(size(D1x(1,:)));
    c_s = (5/3)*(1/10) + 1; %Profile speed
    U = @(x,t) 3*(c_s-1)*(sech( 0.5*sqrt((c_s-1)/c_s)*(x - c_s*t) )).^2;
    RHS_u = @(t,u) - 0.5*D1x*(u.^2) - D1x*u; %BBM right hand side
    I = eye(n);
    LL = I - D2x; %Mass matrix (I-D2x)*U =-D1x*U -0.5*D1x*(U^2)
    init = U(x,0.0); %Initial conditon
    options = odeset('RelTol',2.3e-14,'AbsTol',1e-16,'Mass',LL);
    [t,w] = ode113(@(t,u) RHS_u(t,u),tspan,init,options);
    W1 = w(end,:); W1 = W1';
    Error = norm(W1 - U(x,Tfin),inf)
    plot(x,W1,'r.',x,U(x,Tfin))
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