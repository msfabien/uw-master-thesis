%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-%  Approximate the solution to the KdV equation u_t + uu_x + u_xxx = 0%-%
%-%  on the domain [-pi,pi] with periodic boundary conditions via a sp- %-%
%-%  ectral integrating factor method with a RK4 (IFR4) time-stepping   %-%
%-%  scheme.                                                            %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_2_4_3()
    clear all; close all; clc;
    % Set up grid and single-soliton initial data:
    N = 2^7; % Mesh refinement
    L = 2*pi; % Domain length
    dt = 0.4/N^2; % Time step
    x = (L/N)*(-N/2:N/2-1)';
    x0 = -2; % Soliton initial position
    A = 5; % Soliton amplitude
    tf = 0.1; % Final time
    U = @(x,t) 3*A^2*sech(A*(x - x0)/2 - A^3*t/2).^2;%Exact solution
    u = U(x,0); % Initial condition
    UF = U(x,tf); % Final solution
    v = fft(u);
    k = [0:N/2-1 0 -N/2+1:-1]'; %wave numbers
    ik3 = 1i*k.^3;
    nmax = round(tf/dt);
    for index = 1:nmax
        g = -.5*1i*dt*k;
        E = exp(dt*ik3/2); % Integrating Factor
        E2 = E.^2;
        a = g.*fft(real( ifft( v ) ).^2);
        b = g.*fft(real( ifft(E.*(v+a/2)) ).^2); % 4th-order Runge-Kutta
        c = g.*fft(real( ifft(E.*v + b/2) ).^2); % with integrating factor
        d = g.*fft(real( ifft(E2.*v+E.*c) ).^2);
        v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6; % Update approximation
        u = real(ifft(v));
    end
    plot(x,u,x,UF,'r.')
    norm(u - UF,inf)
end