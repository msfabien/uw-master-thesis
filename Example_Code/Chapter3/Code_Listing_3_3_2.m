%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximate the solution to the BBM equation, which is given by     %-%
%-% u_t + u_x + uu_x - u_{xxt} = 0. The spatial domain is [-100,100]    %-%
%-% and periodic boundary conditions are assumed. The approximation is  %-%
%-% done by a integrating factor spectral method with RK4 timestepping. %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_3_3_2()
    clear all; close all; clc;
    %Set spatial and temporal information
    N = 2^(8);
    L = 200; %spatial domain length
    b = L/(2*(pi));
    x = (L/N)*(-N/2:N/2-1)';
    T = 15; %Final time value
    dt = 1e-1; %Large time step (compare with Code_Listing_3_3_1)
    nmax = round(T/dt);
    %Parameters for exact solution to the BBM equation
    c_s = 2; %Profile speed
    U = @(x,t) 3*(c_s-1)*(sech( 0.5*sqrt((c_s-1)/c_s)*(x - c_s*t) )).^2;
    k = [0:N/2-1 0 -N/2+1:-1]';
    ik = 1i*b* k./(b^2 + k.^2);
    g = -0.5*dt*ik;
    E = exp(-dt*ik/2); E2 = E.^2; %Integrating factor exponentials
    u = U(x,0.0); v = fft(u); %Initial condition
    for n = 1:nmax
        a = g.*fft(real( ifft( v ) ).^2);
        b = g.*fft(real( ifft(E.*(v+a/2)) ).^2); % IFRK4
        c = g.*fft(real( ifft(E.*v + b/2) ).^2);
        d = g.*fft(real( ifft(E2.*v+E.*c) ).^2);
        v = E2.*v + (E2.*a + 2*E.*(b +c) + d)/6;
    end
    u = real(ifft(v));
    Err = norm(u - U(x,T),inf)
    plot(x,u,'r.',x,U(x,T),x,U(x,0),'g')
end