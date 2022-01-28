%-------------------------------------------------------------------------%
% Filename: const_dom.m
% Authors: Juan M. Cardenas  
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description:  This function indicates if f(z) lies in the domain. 
%
% Inputs:
% f - function
% domain - domain of interest
% tol - tolerance
% a,b - parameter for domain
%
% Output:
% output - 1 (inside domain) or 0 (outside domain)
%-------------------------------------------------------------------------%


function[output] = const_dom(f,domain,tol,a,b)
 
    if domain == 1
        if f >= tol
            output = 1;
        else
            output = 0;
        end
    else 
        % domain 2
        if (a <= f) && (f <= b)
            output = 1;
        else
            output = 0;
        end
    end
end
        