%-------------------------------------------------------------------------%
% Filename: g_dim.m
% Authors: Juan M. Cardenas  
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description:  computes g function for function f=f_3 in paper. 
%
% Inputs:
% d - dimension 
%
% Output:
% output - g(d)
%-------------------------------------------------------------------------%

function[output] = g_dim(d)

output = 1 - ((d-2)/100).*(8 + (d-3).*(d-7));

end