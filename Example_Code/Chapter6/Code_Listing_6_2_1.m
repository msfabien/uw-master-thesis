%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates the advection equation u_t=-u_x on [-1,1]. The RBF is  %-%
%-% a Gaussian RBF with a shape parameter equal to 15. Needs the file   %-%
%-% gau.m.                                                              %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_2_1()
    clear all; close all; clc;
    %Exact
    f = @(x,t) exp(-50*(x+0.5-t).^2);
    a = -1; b = -a; %plotting domain
    n = 100; x = linspace(a,b,n)';
    shape = 15; %shape parameter
    cx = shape*ones(n,1);
    [A,D1] = deal(zeros(n));
    for j = 1 : n
        [A(:,j),D1(:,j)] = gau(x,x(j),cx(j)); %construct RBF matrices
    end
    D1x = D1 / ( A );
    Tfin = 1; dt = 1/300; nmax = floor(Tfin/dt);
    D1x = D1x(2:end,2:end); %this handles the BC u(x,0)=0
    x = x(2:end); %BC
    u = f(x,0); I = eye(n-1);
    for jj = 1 : nmax
        u = (I + 0.5*dt*D1x) \ (u - 0.5*dt*(D1x*u)); %Crank-Nicolson
    end
    norm(u-f(x,Tfin),inf), plot(x,f(x,0),'g',x,f(x,Tfin),'b',x,u,'r.')
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