% Approximates the Runge function 1/(1+25x^2) on [-1,1] divided into three
% equally spaced points. The RBF is a inverse multiquadratic with a shape
% parameter equal to 5. Needs the file imq.m.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-%  Approximates the Runge function 1/(1+25x^2) on [-1,1] using three  %-%
%-%  different RBF centers. The RBF chosen is a multiquadratic with a   %-%
%-%  shape parameter of 1.25. Needs the function file mq.m.             %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_1_2()
    clear all; close all; clc;
    f = @(x) 1./(1+25*x.^2);                %Exact
    fp = @(x) (-50*x).*(f(x)).^2;           %Exact
    fpp = @(x) 50*(75*x.^2-1) .* (f(x)).^3; %Exact
    a = -1;  b = -a;  xx = linspace(a,b,100); %plotting domain
    n = 31;
    x = linspace(a,b,n)';  shape = 5; %shape parameter
    cx = shape*ones(n,1);  [A,D1,D2] = deal(zeros(n));
    for j = 1 : n
        [A(:,j),D1(:,j),D2(:,j)] = imq(x,x(j),cx(j)); %construct RBF matrices
    end
    D1x = D1 / ( A ); D2x = D2 / ( A ); %solve for RBF derivative matrices
    lam = A \ f(x);
    Err1 = norm(A*lam - f(x),inf) %Interpolation error
    Err2 = norm(D1x*f(x) - fp(x),inf) %1st derivative error
    Err3 = norm(D2x*f(x) - fpp(x),inf) %2nd derivative error
    plot(xx,f(xx)),hold on,  plot(x,A*lam,'ro'),hold off
end

function [phi,phi1,phi2,phi3,phi4] = imq(x,xc,c)
    % 1-D inverse multiquadric radial basis function
    f = @(r,c) 1./((c*r).^2 + 1);
    r = x - xc;
    phi = f(r,c);

    if nargout > 1
    % 1-st derivative    
    phi1 = -2*r*c^2./((c*r).^2 + 1).^2;
        if nargout > 2
        % 2-nd derivative
            phi2 = (6*c^4*r.^2 - 2*c^2)./(c^2*r.^2 + 1).^3;
            if nargout > 3
            % 3-rd derivative    
            phi3 = -(24*c^4*r.*(c^2*r.^2 - 1))./(c^2*r.^2 + 1).^4;
                if nargout > 4
                % 4-th derivative        
                phi4 = 24*(5*c^8*r.^4 - 10*c^6*r.^2 + c^4)./(c^2*r.^2 + 1).^5;
                end
            end
        end
    end
end