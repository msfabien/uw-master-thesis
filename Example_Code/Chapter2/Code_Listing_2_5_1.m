%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximate the solution to the Small Dispersion Limit KdV equation %-%
%-% u_t + 6uu_x + ep^2u_xxx = 0 on the domain [-5,5] with periodic bou- %-%
%-% ndary conditions via a IFRK4 scheme. Initial condition of -sec^2(x).%-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_2_5_1()
    clear all; close all; clc;
    N = 2^(9); 
    tmax = 0.4; 
    dt = 1e-3;
    ep = 10^(-1.5);
    L = 10; b = L/(2*(pi));
    B = b; % Accounting for domain
    x = (L/N)*(-N/2:N/2-1)'; % length change
    U = @(x)  1./cosh(x).^2;
    u = U(x);
    v = fft(u);
    k = [0:N/2-1 0 -N/2+1:-1]'; 
    ik3 = 1i*k.^3/b^3;
    for n = 1:round(tmax/dt)
        g = -3*1i*dt*k/B;       % Accounting for domain length change
        E = exp(ep^2*dt*ik3/2); 
        E2 = E.^2;
        a = g.*fft(real( ifft( v ) ).^2);
        b = g.*fft(real( ifft(E.*(v+a/2)) ).^2); % 4th-order
        c = g.*fft(real( ifft(E.*v + b/2) ).^2); % Runge-Kutta
        d = g.*fft(real( ifft(E2.*v+E.*c) ).^2);
        v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
        u = real(ifft(v));
    end
    plot(x,u,'r',x,U(x))
end