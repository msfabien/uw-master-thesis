%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximate the solution to the KdV equation u_t + uu_x + u_xxx = 0 %-%
%-% on the domain [-pi,pi] with periodic boundary conditions via a spe- %-%
%-% ctral modified exponential time differencing method with a (mETDR4) %-%
%-% time-stepping scheme.                                               %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_2_4_4()
    clear all; close all; clc;
    N = 2^7; % Mesh refinement
    L = 2*pi; % Domain length
    dt = 0.4/N^2; h = dt; % Time step
    x = (L/N)*(-N/2:N/2-1)';
    x0 = -2; % Soliton initial position
    A = 5; % Soliton amplitude
    tf = 0.1; % Final time
    U = @(x,t) 3*A^2*sech(A*(x - x0)/2 - A^3*t/2).^2; %Exact solution
    u = U(x,0); % Initial condition
    UF = U(x,tf); % Final solution
    % Precompute mETDRK4 scalar quantities (due to Kassam-Trefethen):
    k = [0:N/2-1 0 -N/2+1:-1]'; % Wave numbers
    L = 1i*k.^3; % Fourier multipliers
    E = exp(h*L); E2 = exp(h*L/2);
    M = 64; % Number of points for complex means
    r = exp(2i*pi*((1:M)-0.5)/M); % Roots of unity
    LR = h*L(:,ones(M,1))+r(ones(N,1),:);
    Q = h*mean( (exp(LR/2)-1)./LR ,2);
    f1 = h*mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3,2);
    f2 = h*mean( (4+2*LR+exp(LR).*(-4+2*LR))./LR.^3,2);
    f3 = h*mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3,2);
    g = -.5i*k;
    % March forward in time via ETDRK4 formula (due to Cox-Matthews):
    v = fft(u);
    nmax = round(tf/dt); % Total timesteps
    for index = 1:nmax
        Nv = g.*fft(real(ifft(v)).^2);
        a = E2.*v+Q.*Nv; Na = g.*fft(real(ifft(a)).^2);
        b = E2.*v+Q.*Na; Nb = g.*fft(real(ifft(b)).^2);
        c = E2.*a+Q.*(2*Nb-Nv); Nc = g.*fft(real(ifft(c)).^2);
        v = E.*v+(Nv.*f1+(Na+Nb).*f2+Nc.*f3);
        u = real(ifft(v));
    end
    plot(x,u,x,UF,'r.'), norm(u - UF,inf)
end