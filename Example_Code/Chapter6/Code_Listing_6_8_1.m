%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates a solution to the 1D SW equations using MQ RBFs.       %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_8_1()

    clear all; close all; clc;
    n = 2^(8); dt = 1e-2; Tfin = 15; nmax = floor(Tfin/dt); g = 9.80;
    L = 200; x = ( linspace ( 0, L, n ) )';
    shape = 5; %shape parameter
    cx = (shape)*ones(n,1);
    [Ax,D1x] = deal(zeros(n));
    for j = 1 : n
        [Ax(:,j),D1x(:,j)] = mq(x,x(j),cx(j));
    end
    D1x = D1x / (Ax);

    %Initial conditions
    h = 2.0 + sin ( (2*pi/L) * x ); v = zeros(size(h));
    u = v; uh = u.*h;
    %Right hand side functions h_t = ..., (uh)_t = ...
    RHS_h = @(h,uh) -D1x*(uh);
    RHS_uh = @(h,uh) -D1x*( (uh).^2 ./ h + 0.5*g*h.^2 );
    for jj = 1 : nmax
        k1 = dt*RHS_h(h ,uh); % RK4 step for h
        k2 = dt*RHS_h(h + k1/2,uh);
        k3 = dt*RHS_h(h + k2/2,uh);
        k4 = dt*RHS_h(h + k3 ,uh);
        h = h + (1/6)*( k1 + 2*(k2 + k3) + k4 );
        k1 = dt*RHS_uh(h,uh); % RK4 step for uh
        k2 = dt*RHS_uh(h,uh + k1/2);
        k3 = dt*RHS_uh(h,uh + k2/2);
        k4 = dt*RHS_uh(h,uh + k3);
        uh = uh + (1/6)*( k1 + 2*(k2 + k3) + k4 );
        h(1) = h(end-1); h(end) = h(2); %Periodic boundary conditions
        uh(1) = uh(end-1); uh(end) = uh(2);
        H(:,jj) = h;
    end
    [X,T] = meshgrid(x,linspace(0,Tfin,size(H,2)));
    figure(1), surf(X,T,H'), shading interp
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