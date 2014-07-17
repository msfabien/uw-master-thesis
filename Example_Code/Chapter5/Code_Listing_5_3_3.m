%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates a solution to the 1D SW equations given by the system  %-%
%-% (h)_t + (uh)_x = 0 (1)                                              %-%
%-% (uh)_t + (u^2h + 0.5gh^2)_x = 0, (2)                                %-%
%-% dam-break I.C. h(x >= L/2) = 1; h(x <= L/2) = 3; and initial        %-%
%-% velocity of u(x,0) = 0. Periodic boundary conditions on [0,100] are %-%
%-% assumed, and the approximation is done by a Fourier pseudospectral  %-%
%-% method, with a RK4 time stepping. Hyperviscosity is used.           %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_5_3_3()
    clear all; close all; clc;
    n = 2^(8); 
    g = 1.0;
    L = 100.0; 
    B = L/(2*pi); %Domain rescale factor
    Tfin = 8.0;
    dt = 1e-1; 
    nmax = floor(Tfin/dt);
    x = linspace( 0, L, n )';
    h(x >= L/2) = 1; h(x <= L/2) = 3; h = h'; %Initial depth
    h0 = h;
    u = zeros(size(h)); %Initial velocity
    uh = u.*h; %Initial discharge
    k = [0:n/2-1 0 -n/2+1:-1]';
    %Pseudospectrally evaluate spatial derivatives.
    RHS_h = @(h,uh) -real(ifft(1i*k.*fft(uh))) - ... %RHS for (1)
    n^(-0.65)*real(ifft(k.^2.*fft(h))); %(hyperviscosity term)
    RHS_uh = @(h,uh) -real(ifft(1i*k.*fft(uh.^2./h+0.5*g*h.^2))); %RHS for (2)
    for jj = 1 : nmax
        k1 = (dt/B)*RHS_h(h ,uh); %Update h via RK4
        k2 = (dt/B)*RHS_h(h + k1/2,uh);
        k3 = (dt/B)*RHS_h(h + k2/2,uh);
        k4 = (dt/B)*RHS_h(h + k3 ,uh);
        h = h + (1/6)*( k1 + 2*(k2 + k3) + k4 );
        k1 = (dt/B)*RHS_uh(h,uh); %Update uh via RK4
        k2 = (dt/B)*RHS_uh(h,uh + k1/2);
        k3 = (dt/B)*RHS_uh(h,uh + k2/2);
        k4 = (dt/B)*RHS_uh(h,uh + k3);
        uh = uh + (1/6)*( k1 + 2*(k2 + k3) + k4 );
    end
    plot(x,h0,'r',x,h)
end