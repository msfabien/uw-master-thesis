%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates a solution to the 1D SW equations given by the system  %-%
%-% (h)_t + (uh)_x = 0 (1)                                              %-%
%-% (uh)_t + (u^2h + 0.5gh^2)_x = 0, (2)                                %-%
%-% sine wave initial depth h(x,0) = 2 + sin((pi/100)x) initial         %-%
%-% velocity of u(x,0) = 0. Periodic boundary conditions on [0,200] are %-%
%-% assumed, andthe approximation is done by a Fourier pseudospectral   %-%
%-% method, with a RK4 time stepping. No hyperviscosity is used         %-%
%-% (solution stays smooth).                                            %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Code_Listing_5_3_4()
    clear all; close all; clc;
    N = 2^(7); g = 9.81;
    L = 200.0; B = L/(2*pi); %Domain rescale factor
    Tfin = 15.0; dt = 1e-1; nmax = floor(Tfin/dt);
    x = linspace( 0, L, N )';
    h = 2.0 + sin((2*pi/L) * x); %Initial depth
    u = zeros(size(h)); %Initial velocity
    uh = u.*h; %Initial discharge
    k = [0:N/2-1 0 -N/2+1:-1]';  qq = 1;
    %Pseudospectrally evaluate spatial derivatives.
    RHS_h = @(h,uh)  -real(ifft(1i*k.*fft(uh))); %RHS for (1)
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
        if (mod(jj,5) == 0)
            H(:,qq) = h;    %depth
            UH(:,qq) = uh;  %momentum
            qq = qq + 1;
        end
    end
    [X,T] = meshgrid(x,linspace(0,Tfin,size(H,2)));
    figure(1), surf(X,T,H'), shading interp
    figure(2), plot(x,h)
end