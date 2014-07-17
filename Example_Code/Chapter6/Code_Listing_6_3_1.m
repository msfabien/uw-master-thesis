%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% Approximates the first derivative of the Runge function 1/(1+25x^2) %-%
%-% on the domain [-1,1] via a local Multiquadratic RBF with a shape    %-%
%-% parameter of 2. The stencil size can be any odd number less than or %-%
%-% equal to n.                                                         %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_3_1()
    clear all; close all; clc;
    f = @(x) 1./(1+25*x.^2); %Exact Runge function
    fp = @(x) (-50*x).*(f(x)).^2; %Exact Runge function derivative
    sten = 7; %Stencil size must be odd, and less than or equal to n.
    n = 50; x = linspace(-1,1,n)'; Hh = zeros(n,sten);
    cx1 = (2)*ones(sten,1); x1 = x(1:sten); [A,D1] = deal(zeros(sten));
    
    for j=1:sten  %RBF interpolant pattern repeats for uniform grids
        [A(:,j),D1(:,j)] = mq(x1,x1(j),cx1(j)); 
    end
    
    for ii = 1 : n %construct stencil weights at center ii
        if ( ii <= 0.5*(sten-3)+1 )
            D1para = ii; h = D1(D1para,:); Hh(ii,:) = h;
        elseif ( ii >= n - (0.5*(sten-3)) )
            D1para = ii - n + sten; h = D1(D1para,:); Hh(ii,:) = h;
        else
            D1para = (sten+1)/2; h = D1(D1para,:); Hh(ii,:) = h;
        end
    end
    WW = Hh / (A); %construct stencil weights at every center
    W = deal(zeros(n)); index1 = 1; index2 = sten;
    for ii = 1 : n %construct sparse differentiation matrix
        if ( ii <= 0.5*(sten-3)+1 )
            aa = 1; bb = sten;
            W(ii,aa:bb) = WW(ii,:);
        elseif ( ii >= n - (0.5*(sten-3)) )
            aa = n - sten+1; bb = n;
            W(ii,aa:bb) = WW(ii,:);
        else
            aa = index1; bb = index2; index1 = index1 + 1; index2 = index2 + 1;
            W(ii,aa:bb) = WW(ii,:);
        end
    end
    norm(W*f(x)-fp(x),inf)
    
    subplot(2,1,1), plot(x,W*f(x),'r.',x,fp(x))
    subplot(2,1,2), spy(W)
end

function [phi,phi1,phi2,phi3,phi4] = mq(x,xc,c)
    % 1-D multiquadric radial basis function
    f = @(r,c) sqrt((c*r).^2 + 1);
    r = x - xc;
    phi = f(r,c);

    if nargout > 1
    % 1-st derivative    
     phi1 = (c^2)*r./phi;
        if nargout > 2
        % 2-nd derivative
        phi2 = (c^2)./(phi.^3);
            if nargout > 3
            % 3-rd derivative    
            phi3 = -3*(c^4)*r./(phi.^5);
                if nargout > 4
                % 4-th derivative        
                phi4 = 12*(c^4)*((c*r).^2-0.25)./(phi.^7);
                end
            end
        end
    end
end