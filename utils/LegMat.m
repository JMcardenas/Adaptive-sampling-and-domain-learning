%-------------------------------------------------------------------------%
% Filename: LegMat.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description: Computes the 1D Legendre matrix of size length(pts) x k, 
% where pts is the set of data points
%
% Inputs:
% pts - vector of points 
% k - dimension of vector of points
%
% Output:
% A - 1D Legendre matrix of size pts x k
%-------------------------------------------------------------------------%

function[A] = LegMat(pts,k)

A      = zeros(length(pts),k);
A(:,1) = 1;
A(:,2) = pts*sqrt(3);

for i=2:k-1
    A(:,i+1) = (pts.* (2*i - 1).*A(:,i)./sqrt(i-1/2)- (i - 1).*A(:,i-1)./sqrt(i-3/2)).*sqrt(i+1/2)/i;
end

end

