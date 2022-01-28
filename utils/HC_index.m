%-------------------------------------------------------------------------%
% Filename: HC_index.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description: Computes a hyperbolic cross index set in d dimensions with 
% prod (i_j+1)<= K+1. 
%
% Inputs:
% K - degree of polynomial
% d - dimension 
%
% Output:
% I - the multi-indices in a matrix I where each column is a multi-index
%-------------------------------------------------------------------------%

function I=HC_index(K,d)

% initialize 1D matrix I
I=0:K;

if d>=2
    for t=2:d
        % for each new index i, test the jth index of I to see if in the HC
        J=[];
        for i=0:K
            
            [a,b]=size(I);
            for j=1:b
                
                % get the jth index of I
                z=I(:,j);
                
                % test the product and add the new vector to the new index set
                if (i+1)*prod(z+ones(size(z)))<=K+1
                    
                    z=[I(:,j) ; i ];
                    J=[J z];
                    
                end
            end
        end
        I=J;
    end
    
end

end
