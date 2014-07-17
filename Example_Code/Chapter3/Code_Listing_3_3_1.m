%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximate the solution to the BBM equation, which is given by     %-%
%-% u_t + u_x + uu_x - u_{xxt} = 0. The spatial domain is [-100,100]    %-%
%-% and periodic boundary conditions are assumed. The approximation is  %-%
%-% done by a direct-spectral method with RK4 time-stepping.            %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_3_3_1()
    clear all; close all; clc
    %Set spatial and temporal information
    N = 2^(8);
    L = 200; %spatial domain length
    b = L/(2*(pi)); %Rescaling paramter
    x = (L/N)*(-N/2:N/2-1)';
    T = 15; %Final time value
    dt = 1e-1; %Large time step
    nmax = round(T/dt);
    
    %Parameters for exact solution to the BBM equation
    c_s = 2; %Profile speed
    U = @(x,t) 3*(c_s-1)*(sech( 0.5*sqrt((c_s-1)/c_s)*(x - c_s*t) )).^2;
    
    %Set up right hand side function for the problem dw/dt = RHS(w)
    k = [0:N/2-1 0 -N/2+1:-1]'; %Wave numbers
    g1 = -0.5*1i* b*k./(b^2 + k.^2); g2 = 2*g1;
    RHS = @(v) g1.*fft(real(ifft(v)).^2) + g2.*fft(real(ifft(v)));
    u = U(x,0.0); v = fft(u); %Initial condition
    for n = 1:nmax
        a = RHS(v);
        b = RHS(v + a*(dt/2));
        c = RHS(v + b*(dt/2));
        d = RHS(v + c*(dt));
        v = v + (dt/6)*(a + 2*(b + c) + d);
    end
    u = real(ifft(v));
    Err = norm(u - U(x,T),inf)
    plot(x,u,'r.',x,U(x,T),x,U(x,0),'g')
end