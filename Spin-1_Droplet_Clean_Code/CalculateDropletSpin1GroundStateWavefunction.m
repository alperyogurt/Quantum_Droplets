%% This function returns to the wavefunction of the spin-1 droplet with the
% input arguments: 

% t1:        (1x1) q/(n_1 c_1) determines the strength of the second order
% Zeeman effect q coefficient

% rStep:     (1x1) The step size of radial coordinate

% rInterval: (2x1) The range of the radial coordinates to be examined

% tStep:     (1x1) The step size of the 'infinitesimal' time iterations
% within imaginary time evolution

% tIteration (1x1) Once decided that the calculated state is not
% sufficiently good to be the ground state, the difference equation of the
% imaginary time method is executed 'tIteration' amount more. 

% N:         (1x1) The total number of particles

% expectedPercentageQuality: (1x1) A termination criteria for the obtained
% wavefunction to be the ground state of the Hamiltonian. Calculates the
% percentage difference between the very last value of the chemical potential mu and
% the previous value of the mu, and pass to the next iteration if the
% difference is greater than this variable. 

% expectedMuQuality: (1x1) Calculates the variance of the mu vector and
% divides it to the mean value. If this 'quality factor' is greater than
% the input, continues to the next iteration. 


% Outputs: 
% phiGround: (1xlength(r)) The calculated ground state wavefunction

% r:          The radius vector
% mu:         Calculated chemical potential

%varMu:       Calculated variance of the chemical potential ('Variance'
%means the spatial noise of the chemical potential due to the division of
%H\psi / \psi for each radial point. 

%muPercentageQuality: The difference between the very last chemical
%potential and the current chemical potential values. We require the
%difference between the latest iterations to be lowest possible. 
%%
function [phiGround,r,mu,varMu,muPercentageQuality] = CalculateDropletSpin1GroundStateWavefunction(t1,rStep,rInterval,tStep,tIteration,N,expectedPercentageQuality,expectedMuQuality) 
  
  %Load the 'look up table' that contains the I(t) values up to t=10000. 
  load('ItLookUpTable.mat');
  load('ItDerivativeLookUpTable.mat');

  % Create the radius variable
  r = [rInterval(1):rStep:rInterval(2)];

  %Initialize the ground state wavefunction as a Gaussian (non-normalized
  %for the moment) 
  phiGround = 0.25*exp(-r.^2);
  
  
  
  %Keep counter for each loop of iterations
  Counter=0;
  
  %Assign some 'big' value to quality factors 
  muPercentageChange = 5;
  muQuality = 5;
  
  %Imaginary time derivation iterations to be executed as long as both of  
  %the quality criterias are satisfied
  while (muPercentageChange>expectedPercentageQuality || muQuality>expectedMuQuality)
      %Skip this part if it is the very first iteration
      if Counter ~= 0
        muOld = muNew;
        varMuOld = varMuNew;
      end
      
      %Imaginary Time Iterations on the time domain. Number of iterations
      %is tIteration
      for ii=1:tIteration
          
        %Normalize the phiGround at each iteration
        A = trapz(r,(r.^2).*(phiGround.^2));
        phiGround = phiGround*sqrt(N/(4*pi*A));
        
        % Calculate I(t) and I'(t) for different regions of the
        % wavefunction
        ItDummyVector = FindIt(t1./(phiGround(2:length(r)-1).^2),tVectorExtended,phi1VectorExtended);
        ItDerDummyVector = FindItDerivative(t1./(phiGround(2:length(r)-1).^2),tVectorExtended,phi1DerVectorExtended);

        %Calculate the infinitesimal change 'phiSpaceChange' and iterate the phiGround
        phiSpaceChange = (phiGround(3:length(r))-phiGround(2:length(r)-1))./(r(2:length(r)-1)*rStep) + ...
            (1/2)*(phiGround(3:length(r))-2*phiGround(2:length(r)-1)+phiGround(1:length(r)-2))./(rStep^2) + ...
            3*phiGround(2:length(r)-1).^3 - ...
            (5/2*ItDummyVector.*phiGround(2:length(r)-1).^3-t1*ItDerDummyVector.*phiGround(2:length(r)-1)).*phiGround(2:length(r)-1);

        %Complete the iteration
        phiGround(2:length(r)-1) = phiGround(2:length(r)-1) + tStep.*phiSpaceChange;
        
        %Handle the very first and very last points of the wavefunction
        %'phiGround'
        phiGround(length(r)) = phiGround(length(r)-1);
        phiGround(1) = phiGround(2);
  end
   
    %Calculate the muVector, mean value of the muVector and the variance
    [muVector,mu, varMu, muPercentageQuality] = CalculateDropletSpin1ChemicalPotential(r,rStep,phiGround,ItDummyVector,ItDerDummyVector,t1);
    
    %Print some values to follow the variables to be swept at the moment
    muNew =mu
    N
    varMuNew = varMu;
    varMu
    t1
    
    
    %Complete the quality controls and renew the quality variables in order
    %to be checked on the next 'while' turn. 
    
    if Counter ~= 0
        muPercentageChange = (muOld-muNew)/muNew*100;
        muQuality = varMuNew/muNew;
    end
    
    %If the imaginary time evolution kept the execution until the absolute 
    %magnitude of the chemical potential value is less then 0.01, stop the
    %evolution and exit the while loop. Interpret this wavefunction as the
    %'expanded gas'
   
    if abs(muNew) <0.01
        muQuality = 0;
        muPercentageChange = 0;
    end 
    
    Counter = Counter + 1;
  end
 end