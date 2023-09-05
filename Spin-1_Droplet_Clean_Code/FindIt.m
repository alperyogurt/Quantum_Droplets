% This function gives the result of the I(t) for different t values. If
% t<10000, it finds the closest data point in the look-up table. If t value
% is greater then 10000 (t>10000), then uses an analytic expression of
% sqrt(t)


%Inputs: 

% t:                 q/n_1 c_1        It might also be a vector of t values. 

%tVectorExtended:    the range of t values covered in look up .mat files

%phi1VectorExtended: the range of I(t) values covered in look up .mat files


%Output: 

% It:              I(t) value or vector if the input t is a vector. 
function [It] = FindIt(t,tVectorExtended,phi1VectorExtended) 

    It = zeros(size(t));
    
    %Find the indexes with t values greater then 10000
    indGreat = find(ge(t,10000));
    
    %Find the indexes with t values less then 10000
    indLower = find(t<10000);
    
    % Calculate I(t) for t>10000 as sqrt(t)
    It(indGreat) = (t(indGreat).^(0.5));
    
    % Calculate I(t) for t<10000 by finding the nearest point from the
    % look-up file
    t(indLower) = round(t(indLower)*10)/10;
    It(indLower) = phi1VectorExtended(t(indLower)*10+1);
end