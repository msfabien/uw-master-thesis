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


function Code_Listing_6_1_1()
    clear all; close all; clc;
    f = @(x) 1./(1+25*x.^2);%Runge function
    a = -1; b = 1;
    MM = 100;
    xp = linspace(a,b,MM)'; %fine grid to smooth out RBF approximations
    INDEX = [3,5,30]; %number of RBF centers
    shape = 1.25; %RBF blanket shape parameter
    %Variables for plots
    IdX = [0,-0.1,-.2]; ccc = {'r','b','k'}; cc2 = {'ro','bo','ko'};
    hold on,
    for ii = 1 : length(INDEX)
        nmm = INDEX(ii);
        uu = linspace(a,b,nmm)'; cx = shape*ones(nmm,1);
        n = INDEX(ii);
        [Ap] = deal(zeros(MM,n));
        for j = 1 : nmm
            [Ap(:,j)] = mq(xp,uu(j),cx(j));
        end
        alpha = Ap \ f(xp);
        y4 = Ap*alpha;
        plot(xp,y4,ccc{ii}),
        plot(uu,ones(size(uu))*IdX(ii),cc2{ii})
    end
    plot(xp,f(xp),'go')
    hold off
    legend('3 RBF Approx','3 RBF cenetrs','5 RBF Approx','5 RBF cenetrs',...
    '30 RBF Approx','30 RBF cenetrs','Exact')
    title('RBF approximation'), ylim([-0.3 1.2])
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
