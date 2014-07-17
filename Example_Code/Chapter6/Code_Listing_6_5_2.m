%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates a solution to the KdV equation on the domain [-100,100 %-%
%-% with zero flux boundary conditions. A local Gaussian RBF spectral   %-%
%-% method is implemented with a shape parameter of 1. The stencil size %-%
%-% must be odd.  Needs the files gau.m and constructW.m.               %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_5_2()
    clear all; close all; clc;
    dt = 1e-1; Tfin = 10;
    L1 = 100; tspan = 0 : dt: Tfin;
    n = 350; x = linspace ( -L1, L1, n )';
    sten = 253; x1 = x(1:sten); shape = 1;
    cx1 = (shape)*ones(sten,1); [A,D1,D3] = deal(zeros(sten));
    for j=1:sten
        [A(:,j),D1(:,j),~,D3(:,j)] = gau(x1,x1(j),cx1(j));
    end
    %Construct local RBF differentiation matrices
    [W1] = constructW(A,D1,sten,n); [W2] = constructW(A,D3,sten,n);
    W1x = W1; W3x = W2;
    %Zero flux boundary conditions
    W1x(1,:) = zeros(size(W1x(1,:))); W1x(end,:) = zeros(size(W1x(1,:)));
    W3x(1,:) = zeros(size(W1x(1,:))); W3x(end,:) = zeros(size(W1x(1,:)));
    A = sqrt(1)/sqrt(6); L = 1; x0 = 0;
    U = @(x,t) 3*A^2*sech(A*L*(x - x0/L)/2 - A^3*t/2).^2; %Exact solution
    RHS_u = @(t,u) -u.*(W1x*u) - W3x*u; %Right hand side, u_t=-uu_x-u_{xxx}
    init = U(x,0.0);
    options = odeset('RelTol',2.3e-14,'AbsTol',1e-16);
    [t,w] = ode113(@(t,u) RHS_u(t,u),tspan,init,options);
    w1 = w(end,:); w1 = w1';
    Error = norm(w1 - U(x,Tfin),inf)
    plot(x,w1,'r.',x,U(x,Tfin),'-')
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