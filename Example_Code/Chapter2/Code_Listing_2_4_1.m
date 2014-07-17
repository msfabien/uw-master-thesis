%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates a solution to the KdV equation u_t + uu_x + u_xxx = 0  %-%
%-% on the domain [-pi,pi] with periodic boundary conditions via a pse- %-%
%-% udo-spectral method with RK4 time-stepping scheme.                  %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_2_4_1()
    clear all; close all; clc;
    N = 2^7; % Mesh refinement
    L = 2*pi; % Domain length
    dt = 0.05*(L/N)^3; % Time step
    x = (2*pi/N)*(-N/2:N/2-1)';
    x0 = -2; % Soliton initial position
    A = 16; % Soliton amplitude
    U = @(x,t) 3*A^2*sech(A*(x - x0)/2 - A^3*t/2).^2; %Exact solution
    u = U(x,0.0);  % Initial condition
    tf = 0.006; % Final time
    nt = round(tf/dt); % Total timesteps
    for index = 1:nt
        a = -dt*RHS_KdV(u); % March forward du/dt = -RHS_KdV(u)
        b = -dt*RHS_KdV(u + 0.5*a); % using the RK4 scheme
        c = -dt*RHS_KdV(u + 0.5*b);
        d = -dt*RHS_KdV(u + c);
        u = u + (a + 2*b + 2*c + d)/6; % Update solution u
    end
    plot(x,u,x,U(x,0.0),x,U(x,tf),'.')
    norm(U(x,tf) - u,inf)
end

function U = RHS_KdV(u)
    % This function pseudospectrally calculates the spatial derivatives in the
    % KdV equation u_t = - ( uu_x + u_xxx ). We can deal with the nonlinear
    % term uu_x easily this way.
    N = length(u);
    uhat = fft(u);
    k = [0:N/2-1 0 -N/2+1:-1]';
    ux = real(ifft(1i*k.*uhat));
    uxxx = real(ifft(-1i*k.^3.*uhat));
    U = u.*ux + uxxx;
end