%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% A 2D SWE simulation. Water droplet in a rectangular paralelipiped.  %-%
%-% The boundary conditions are periodic. Grid: [-100,100]x[-100,100].  %-%                                    %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_5_6_1()
    clear all; close all; clc;
    N = 2^(6); L = 200; x = linspace(-L/2,L/2,N)'; y = x; g = 9.8;
    Tfin = 45; dt = 5e-2; nmax = floor(Tfin/dt);
    B = L/(2*pi); %Domain scaling parameter
    k = [0:N/2-1 0 -N/2+1:-1]'; %1D wave numbers
    [kx,ky] = meshgrid(k); %2D wave numbers
    [xx, yy] = meshgrid(x,y);
    %Initial condition
    h = 1.0 + exp(-(xx.^2+yy.^2)/(L/2)); uh = zeros(size(h)); vh = uh;
    RHS_h = @(h,uh,vh) -real(ifft2(1i*kx.*fft2(uh))) - ...
    real(ifft2(1i*ky.*fft2(vh)));
    RHS_uh = @(h,uh,vh) -real(ifft2(1i*kx.*fft2(uh.^2 ./ h + 0.5*g*h.^2))) -...
    -real(ifft2(1i*ky.*fft2(((vh)./h).*uh)));
    RHS_vh = @(h,uh,vh) -real(ifft2(1i*ky.*fft2(vh.^2 ./ h + 0.5*g*h.^2))) -...
    -real(ifft2(1i*kx.*fft2(((uh)./h).*vh)));
    grid = surf(h); %plot simulation
    axis([1 N 1 N 1 3]); view(-20,55); hold all; %plot simulation
    for jj = 1 : nmax
        set(grid ,'zdata', h); drawnow %plot simulation
        h = h + (dt/B)*RHS_h(h ,uh,vh); %naive forward euler
        uh = uh + (dt/B)*RHS_uh(h,uh,vh);
        vh = vh + (dt/B)*RHS_vh(h,uh,vh);
    end
end