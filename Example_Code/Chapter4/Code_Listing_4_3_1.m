%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates the solution to the CH equation, which is given by     %-%
%-% u_t - u_{xxt} + 3uu_x = 2u_xu_{xx}+uu_{xxx}. The spatial domain is  %-%
%-% [0,50] and periodic boundary conditions are assumed. The approxim-  %-%
%-% ation is done by a Fourier spectral method with RK4 time-stepping.  %-%
%-% When kappa=0 the CH equation loses its physical significance.       %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_4_3_1()
    clear all; close all; clc;
    N = 2^11;   %number of grid points
    L = 50;     %domian length
    T = 3.2;    %final time
    c = 1;      %peakon speed
    dt = 5e-2;  %time step
    x = L*(0:(N-1))/N;
    nmax = round(T/dt);
    M = [0:(N/2-1) (-N/2):(-1)];
    k = 2*pi*M/L; %Domain length incorporated into wavenumbers
    V = @(x,t) c*exp(-abs(x-5-c*t)); %Exact solution
    c_1 = -1i*k/2;
    c_2 = -1i*k./(1+k.^2);
    RHS1 = @(v) c_1.*fft(real(ifft(v)).^2);
    RHS2 = @(v) c_2.*fft(real(ifft(v)).^2 + 0.5*real(ifft(1i*k.*v)).^2 );
    RHS = @(v) RHS1(v) + RHS2(v); %Right hand side for v_t = RHS(v)
    u = V(x,0.0);
    v = fft(u); %Initial condition
    for n = 1:nmax %Time stepping loop
        for jj = 1 : length(k) % 2/3's dealiasing rule
            if ( abs(M(jj)) > (2*pi/L)*2*N/3 )
                v(jj) = 0;
            end
        end
        a = dt*RHS(v); %RK4 update
        b = dt*RHS(v + a/2);
        c = dt*RHS(v + b/2);
        d = dt*RHS(v + c);
        v = v + (a + 2*b + 2*c + d)/6;
    end
    u = real(ifft(v));
    norm(u - V(x,T),inf)
    plot(x,u,'r.',x,V(x,T),x,V(x,0.0),'g')
end