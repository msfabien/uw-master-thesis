% This function constructs the local RBF differentiation matrices.
% INPUT VARIABLES:
% **A , RBF interpolation matrix of size (sten x sten)
% **L , RBF differentiation matrix of size (sten x sten)
% **sten , Stencil size, must be odd and less than or equal to n
% **n , number of RBF centers
function [W] = constructW(A,L,sten,n)
Hh = zeros(n,sten);
for ii = 1 : n %Construct RBF differentiation weights at center ii
if ( ii <= 0.5*(sten-3)+1 )
Lpara = ii; h = L(Lpara,:); Hh(ii,:) = h;
elseif ( ii >= n - (0.5*(sten-3)) )
Lpara = ii - n + sten; h = L(Lpara,:); Hh(ii,:) = h;
else
Lpara = (sten+1)/2; h = L(Lpara,:); Hh(ii,:) = h;
end
end
WW = Hh / (A); %Construct all RBF differentiation weights
W = deal(zeros(n)); index1 = 1; index2 = sten;
for ii = 1 : n %Construct RBF differentiation matrix W
if ( ii <= 0.5*(sten-3)+1 )
aa = 1; bb = sten; W(ii,aa:bb) = WW(ii,:);
elseif ( ii >= n - (0.5*(sten-3)) )
aa = n - sten+1; bb = n; W(ii,aa:bb) = WW(ii,:);
else
aa = index1; bb = index2; index1 = index1 + 1; index2 = index2 + 1;
W(ii,aa:bb) = WW(ii,:);
end
end