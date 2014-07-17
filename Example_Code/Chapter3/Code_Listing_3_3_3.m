%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximate the solution to the BBM equation, which is given by     %-%
%-% u_t + u_x + uu_x - u_{xxt} = 0. The spatial domain is [-100,100]    %-%
%-% and periodic boundary conditions are assumed. The approximation is  %-%
%-% done by a exponential time differencing spectral method with RK4    %-%
%-% time-stepping.                                                      %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_3_3_3()
    clear all; close all; clc;
    %Set spatial and temporal information
    N = 2^(8);
    L = 200; b = L/(2*(pi)); %spatial domain length and rescaling
    x = (L/N)*(-N/2:N/2-1)';
    T = 15; %Final time value
    dt = 1e-1; h = dt; %Large time step
    nmax = round(T/dt);
    %Parameters for exact solution to the BBM equation
    c_s = 2; %Profile speed
    U = @(x,t) 3*(c_s-1)*(sech( 0.5*sqrt((c_s-1)/c_s)*(x - c_s*t) )).^2;
    % Precompute ETDRK4 scalar quantities (due to Kassam-Trefethen):
    k = [0:N/2-1 0 -N/2+1:-1]'; % Wave numbers
    L = -1i*b*k./(b^2+k.^2); % Fourier multipliers
    E = exp(h*L); E2 = exp(h*L/2);
    M = 2^6; % Number of points for complex means
    r = exp(2i*pi*((1:M)-0.5)/M); % Roots of unity
    LR = h*L(:,ones(M,1))+r(ones(N,1),:);
    Q = h*mean( (exp(LR/2)-1)./LR ,2);
    f1 = h*mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3,2);
    f2 = h*mean( (4+2*LR+exp(LR).*(-4+2*LR))./LR.^3,2);
    f3 = h*mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3,2);
    g = 0.5*L;
    u = U(x,0.0); v = fft(u); %Initial conditon
    
    % March forward in time via ETDRK4 formula (due to Cox-Matthews):
    for index = 1:nmax
        Nv = g.*fft(real(ifft(v)).^2); a = E2.*v+Q.*Nv;
        Na = g.*fft(real(ifft(a)).^2); b = E2.*v+Q.*Na;
        Nb = g.*fft(real(ifft(b)).^2); c = E2.*a+Q.*(2*Nb-Nv);
        Nc = g.*fft(real(ifft(c)).^2);
        v = E.*v+(Nv.*f1+(Na+Nb).*f2+Nc.*f3);
    end
    u = real(ifft(v));
    Err = norm(u - U(x,T),inf)
    plot(x,u,'r.',x,U(x,T),x,U(x,0),'g')
end