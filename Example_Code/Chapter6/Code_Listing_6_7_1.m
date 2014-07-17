%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates the solution to a Boussinesq system given by           %-%
%-% n_t + u_x + (un)_x + a*u_{xxx} - b*n_{xxt} = 0 (1)                  %-%
%-% u_t + n_x + (uu_x) + c*n_{xxx} - d*u_{xxt} = 0, (2)                 %-%
%-% where n = n(x,t), u = u(x,t), a = b = c = 0, and d = 1/3.           %-%
%-% See Chen 1997 for information about soliton solutions to the system.%-%
%-% The spatial domain is [-40,80] with zero flux boundary conditions.  %-%
%-% The approximation is done by a multiquadratic radial basis function %-%
%-% spectral method with shape parameter 0.4. Matlab's ODE113 is used   %-%
%-% to evolve the equations in time. A uniform grid is used for the RBF %-%
%-% centers.                                                            %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_7_1()
    clear all; close all; clc;
    N = 350; dt = 1e-1; Tfin = 1; tspan = 0 : dt : Tfin;
    L = 80; d = 1/3; C_s = (21/3)*(1/10) + 1; rho = 0.5; x0 = 0;
    u = @(x,t) (1-d*rho)*C_s+3*d*C_s*rho*sech(0.5*sqrt(rho)*(x+x0-C_s*t)).^2;
    eta= @(x,t) -1*ones(size(x));
    x = linspace(-L/2,L,N)'; Nx = length(x);
    cx = (0.4)*ones(Nx,1); [Ax,D1x,D2x] = deal(zeros(Nx));
    for j=1:Nx
        [Ax(:,j),D1x(:,j),D2x(:,j)] = mq(x,x(j),cx(j));
    end
    %Set up interpolation and differentiation matrices
    D1x = D1x / (Ax); D2x = D2x / (Ax);
    D1x(1,:) = zeros(size(D1x(1,:))); D1x(end,:) = zeros(size(D1x(1,:)));
    D2x(1,:) = zeros(size(D1x(1,:))); D2x(end,:) = zeros(size(D1x(1,:)));
    I = eye(N); ETA = eta(x,0); U = u(x,0); init = [ETA; U];
    %Right hand side function for ODE solver
    RHS_u = @(t,q) [-D1x*q(N+1:end)-D1x*(q(N+1:end).*q(1:N));...
    -D1x*q(1:N) - 0.5*D1x*(q(N+1:end).^2)];
    L2 = I - d*D2x; L2 = blkdiag(I,L2);%need a mass matrix for u_{xxt}
    options = odeset('RelTol',2.3e-14,'AbsTol',1e-16,'Mass',L2);
    [t,w] = ode113(@(t,q) RHS_u(t,q),tspan,init,options);
    W = w(end,:)'; U = W(N+1:end);
    Err = norm(U - u(x,Tfin),inf)
    figure(1), plot(x, abs(U-u(x,Tfin) ),'r.-')
end

function [phi,phi1,phi2,phi3,phi4] = mq(x,xc,c)
    % 1-D multiquadric radial basis function
    f = @(r,c) sqrt((c*r).^2 + 1);
    %f = @(r,c) exp(-(c*r).^2);
    %f = @(r,c) 1 ./ ((c*r).^2 + 1);
    r = x - xc;
    phi = f(r,c);

    if nargout > 1
    % 1-st derivative    
     phi1 = (c^2)*r./phi;
     %phi1 = -2*r*c^2.*exp(-(c*r).^2);
     %phi1 = -2*(c.^2*r) ./ ((c*r).^2 + 1).^2;
        if nargout > 2
        % 2-nd derivative
        phi2 = (c^2)./(phi.^3);
            if nargout > 3
            % 3-rd derivative    
            phi3 = -3*(c^4)*r./(phi.^5);
                if nargout > 4
                % 4-th derivative        
                phi4 = 12*(c^4)*((c*r).^2-0.25)./(phi.^7);
                end
            end
        end
    end
end