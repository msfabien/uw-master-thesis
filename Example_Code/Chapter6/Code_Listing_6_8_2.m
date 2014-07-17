%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-% ------------------------------------------------------------------- %-%
%-% A 2D SWE simulation. Water droplet in a rectangular paralelipiped.  %-%
%-% ------------------------------------------------------------------- %-%
%-% Author: Maurice S. Fabien, University of Washington (Jan-Jun 2014)  %-%
%-%                          , Rice University          (2014-    )     %-%
%-% Email : fabien@rice.edu                                             %-%
%-% GitHub: https://github.com/msfabien/                                %-%
%-% ------------------------------------------------------------------- %-%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Code_Listing_6_8_2()

    clear all; close all; clc;
    n = 80;  dt = 1e-2;  g = 9.8;
    L1 = 100; x = ( linspace ( -L1, L1, n ) )'; y = x;
    [XX,YY] = meshgrid(x,y); xx = XX(:); yy = YY(:);
    shape = 0.2393; %shape parameter
    cx = shape*ones(n,1);
    %Set up RBF matrices
    [Ax,D1x] = deal(zeros(n));
    for j = 1 : n
        [Ax(:,j),D1x(:,j)] = mq(x,x(j),cx(j));
    end
    D1x = D1x / (Ax); I = eye(n); %Compute 1D derivatives
    Lx = sparse(kron(I,D1x)); Ly = sparse(kron(D1x,I)); %Compute partial derivatives
    %Initial conditions
    h = 1.0 + exp( -(xx.^2+yy.^2)/(2*L1) ); H = zeros(n, n); U = H; V = H;
    H = H + reshape(h,n , n); vh = zeros(size(h)); uh = vh;
    grid = surf(H); axis([0 n 0 n 1 3]); hold all; %plot simulation
    RHS_h = @(h,uh,vh) -Lx*(uh) - Ly*(vh);
    RHS_uh = @(h,uh,vh) -Lx*( (uh).^2 ./ h + 0.5*g*h.^2 ) - Ly*((uh./h).*vh);
    RHS_vh = @(h,uh,vh) -Ly*( (vh).^2 ./ h + 0.5*g*h.^2 ) - Lx*((uh./h).*vh);
    while ( 1 == 1 ) %for jj = 1 : nmax
        set(grid ,'zdata', H); drawnow %plot simulation
        h = h + (dt)* RHS_h(h,uh,vh);
        uh = uh + (dt)*RHS_uh(h,uh,vh);
        vh = vh + (dt)*RHS_vh(h,uh,vh);
        u = uh ./ h; U(:,:) = reshape(u,n, n); %Reflecting boundary conditions
        U(:,1) = 0.0; U(:,end) = 0.0;
        U(1,:) = -U(2,:); U(end,:) = -U(end-1,:);
        uh = reshape(U,n^2, 1) .* h;
        v = vh ./ h; V(:,:) = reshape(v,n, n); %Reflecting boundary conditions
        V(1,:) = 0.0; V(end,:) = 0.0;
        V(:,1) = -V(:,2); V(:,end) = -V(:,end-1);
        vh = reshape(V,n^2, 1) .* h;
        H(:,:) = reshape(h,n, n); %Consistency conditions
        H(:,1) = H(:,2); H(:,n) = H(:,n-1); %(see Cleve Moler reference)
        H(:,n) = H(:,n-1);
        H(1,:) = H(2,:); H(n,:) = H(n-1,:);
    end
end

function [phi,phi1,phi2,phi3,phi4] = mq(x,xc,c)
    % 1-D multiquadric radial basis function
    f = @(r,c) sqrt((c*r).^2 + 1);
    %f = @(r,c) exp(-(c*r).^2);
    %f = @(r,c) 1 ./ ((c*r).^2 + 1);
    r = x - xc;
    phi = f(r,c);

    if nargout > 1
    % 1-st derivative    
     phi1 = (c^2)*r./phi;
     %phi1 = -2*r*c^2.*exp(-(c*r).^2);
     %phi1 = -2*(c.^2*r) ./ ((c*r).^2 + 1).^2;
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