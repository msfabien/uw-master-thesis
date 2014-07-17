%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates the solution to the small dispersion CH equation,given %-%
%-% u_t +2*Ku_x+ 3uu_x - ep^2(u_{xxt}+2u_xu_{xx}+uu_{xxx})= 0. The      %-%
%-% spatial domain is [-5,5] and periodic boundary conditions are       %-%
%-% assumed.  The approximation is done by a Fourier spectral method    %-%
%-% with RK4 timestepping. (K = 1.2).                                   %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_4_4_1()
    clear all; close all; clc
    N = 2^(9);         %number of grid points
    T = 1.0;           %final time simulation
    dt = 1e-3;         %time step
    nmax = round(T/dt);
    L = 10;            %domain length
    ep = 10^(-1.5);    %small dispersion parameter
    kappa = 1.2;       %CH parameter kappa
    x = (L/N)*(-N/2:N/2-1);
    M = [0:N/2-1 0 -N/2+1:-1];
    k = 2*pi*M/L; %Domain length incorporated into wavenumbers
    V = @(x,t) -sech(x).^2;
    c_1 = (-1i*k/2).*(ep^2+ep^2*k.^2)./(1 + ep^2*k.^2);
    c_2 = -1i*k./(1 + ep^2*k.^2);
    RHS1 = @(v) c_1.*fft(real(ifft(v)).^2);
    RHS2 = @(v) c_2.*fft(real(ifft(v)).^2 + ep^2*0.5*real(ifft(1i*k.*v)).^2 );
    RHS3 = @(v) 2*kappa*c_2.*v;
    RHS = @(v) RHS1(v) + RHS2(v) + RHS3(v);
    u0 = V(x,0.0);
    v = fft(u0);
    for n = 1:nmax
        a = dt*RHS(v);
        b = dt*RHS(v + a/2);
        c = dt*RHS(v + b/2);
        d = dt*RHS(v + c);
        v = v + (a + 2*b + 2*c + d)/6;
    end
    u = real(ifft(v));
    plot(x,u0,x,u,'r')
end