%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates the solution to the KdV eqn (u_t + uu_x + u_xxx = 0).  %-%
%-% The spatial domain is [-100,100] with zero flux boundary conditions.%-%
%-% The approximation is done by a Gaussian radial basis function       %-%
%-% spectral method with shape parameter 1. A uniform grid is used for  %-%
%-% the RBF centers.                                                    %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_4_1()
    clear all; close all; clc;
    n = 350;
    L1 = 100; 
    x = linspace ( -L1, L1, n )';
    dt = 1e-1; 
    Tfin = 10;  
    tspan = 0 : dt: Tfin;
    shape = 1;
    cx = (shape)*ones(n,1);
    %Set up RBF interpolation and differentiation matrices.
    [Ax,D1x,D3x] = deal(zeros(n));
    for j=1:n
        [Ax(:,j),D1x(:,j),~,D3x(:,j)] = gau(x,x(j),cx(j));
    end
    D1x = D1x /( Ax ); D3x = D3x /( Ax );
    %Zero flux boundaries
    D1x(1,:) = zeros(size(D1x(1,:))); D1x(end,:) = zeros(size(D1x(1,:)));
    D3x(1,:) = zeros(size(D1x(1,:))); D3x(end,:) = zeros(size(D1x(1,:)));
    %Exact solution
    A = 1/sqrt(6); L = 1; x0 = 0;
    U = @(x,t) 3*A^2*sech(A*L*(x - x0/L)/2 - A^3*t/2).^2;
    %Right hand side of u_t=-u_x-u_{xxx}
    RHS_u = @(t,u) -u.*(D1x*u) - D3x*u;
    init = U(x,0.0);
    options = odeset('RelTol',2.3e-14,'AbsTol',1e-16);
    [t,w] = ode113(@(t,u) RHS_u(t,u),tspan,init,options);
    W1 = w(end,:); W1 = W1';
    Error = norm(W1 - U(x,Tfin),inf)
    plot(x,W1,'r.',x,U(x,Tfin),x,init,'go')
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