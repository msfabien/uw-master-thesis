%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates a solution to the 1D SW equations given by the system  %-%
%-% (h)_t + (uh)_x = 0 (1)                                              %-%
%-% (uh)_t + (u^2h + 0.5gh^2)_x = -ghB_x, (2)                           %-%
%-% sine wave initial depth h(x,0) = 1 + sin((pi/10)x) initial velocity %-%
%-% of u(x,0) = 0. Periodic boundary conditions on [0,10] are assumed,  %-%
%-% and the approximation is done by a Fourier pseudospectral method,   %-%
%-% with a RK4 time stepping. Hyperviscosity is used, and nonzero       %-%
%-% topography is imposed.                                              %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Code_Listing_5_3_5()
clear all; close all; clc;
N = 2^(7); g = 9.81; L = 10; b = L/(2*pi);
Tfin = 1.5; dt = 1e-3; nmax = floor(Tfin/dt);
x = linspace(0,L,N)'; h = 1 + sin((1*pi/L)*x);
u = zeros(size(h)); uh = u.*h;
k = [0:N/2-1 0 -N/2+1:-1]';
% Linear topography 1-x/10
p = polyfit([0 10],[1 0.0],1);
B = @(x,t) (-1/10)*x + 1;
B_x = @(x,t) (-1/10);
RHS_h = @(h,uh,t) -real(ifft(1i*k.*fft(uh))) - ...
(5e-1*N^(-0.5))*real(ifft(-(1i*k).^2.*fft(h))); % Hyperviscosity
RHS_uh = @(h,uh,t) -real(ifft(1i*k.*fft(uh.^2 ./ h + 0.5*g*h.^2))) - ...
g*h.*B_x(x,t); % Linear topography
for jj = 1 : nmax
t = jj*dt;
k1 = (dt/b)*RHS_h(h ,uh,t); % RK4
k2 = (dt/b)*RHS_h(h + k1/2,uh,t);
k3 = (dt/b)*RHS_h(h + k2/2,uh,t);
k4 = (dt/b)*RHS_h(h + k3 ,uh,t);
h = h + (1/6)*( k1 + 2*(k2 + k3) + k4 );
k1 = (dt/b)*RHS_uh(h,uh,t);
k2 = (dt/b)*RHS_uh(h,uh + k1/2,t);
k3 = (dt/b)*RHS_uh(h,uh + k2/2,t);
k4 = (dt/b)*RHS_uh(h,uh + k3,t);
uh = uh + (1/6)*( k1 + 2*(k2 + k3) + k4 );
end
plot(x,h)
end