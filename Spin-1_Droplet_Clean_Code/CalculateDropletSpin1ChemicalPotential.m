%%This function returns to the muVector, mu value and the variance of the
%%muVector for any wavefunction provided as input. 

% Inputs:
% rStep:     (1x1) The step size of radial coordinate

% r:               The radius vector

% phiGround:       Input wavefunction

% ItDummyVector:   I(t) values for the different regions of the input
% wavefunction

% ItDerDummyVector: I'(t) values for the different regions of the input
% wavefunction

% t1:              q/(n_1 c_1)

% Outputs: 

% muVector:         mu vector for each radial points

% mu:               calculated mu value. The mean of the muVector.

% varMu:            variance of the muVector

%muPercentageQuality: the ratio of the varMu and mu
%%
function [muVector, mu, varMu, muPercentageQuality] = CalculateDropletSpin1ChemicalPotential(r,rStep,phiGround,ItDummyVector,ItDerDummyVector,t1)
    
    % Calculate first and second derivatives of the wavefunction
    derivative = zeros(size(phiGround));
    
    for ii= 2:length(r)-1
        derivative(ii) = (phiGround(ii+1)-phiGround(ii))/rStep;
    end
    derivative(length(r)) = derivative(length(r)-1);
    
    derivative2 = zeros(size(derivative));
    
    for ii= 2:length(r)-1
        derivative2(ii) = (derivative(ii+1)-derivative(ii))/rStep;
    end
    
    %Handle the end point
    derivative2(length(r)) = derivative2(length(r)-1);
    
    %Resize the ItDummyVectors for proceeding calculations
    ItDummyVector = [ItDummyVector ItDummyVector(length(ItDummyVector))];
    ItDerDummyVector = [ItDerDummyVector ItDerDummyVector(length(ItDerDummyVector))];
    
    
    
    
    %Calculate Hamiltionian*phiGround whose eigenvalues are chemical
    %potentials
    A = -derivative(2:length(r))./r(2:length(r)) - derivative2(2:length(r))/2 -3*phiGround(2:length(r)).^3 +...
        (5/2*ItDummyVector.*phiGround(2:length(r)).^3-t1*ItDerDummyVector.*phiGround(2:length(r))).*phiGround(2:length(r));
    
    %Create the muVector for each space coordinate r
    muVector = A./phiGround(2:length(r));
    
    % Handle the last part of the muVector by hand in order to terminate
    % faster
    muVector(length(r)-1) = muVector(length(r)-6);
    
    %Calculate the output variables
    varMu = var(muVector);
    mu = mean(muVector);
    muPercentageQuality = varMu/abs(mu)*100;
end
