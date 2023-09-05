function [E, MF, LHY, KE] = dropletSpin1TotalEnergy(t1,r,rStep,phiGround,tVectorExtended,phi1VectorExtended) 
    derivative = zeros(size(phiGround));
    
    for ii= 2:length(r)-1
        derivative(ii) = (phiGround(ii+1)-phiGround(ii))/rStep;
    end
    ItDummyVector = FindIt(t1./(phiGround.^2),tVectorExtended,phi1VectorExtended);
    
    KEDensity = (1/2*(derivative.^2)).';
    MFDensity = (-3/2*phiGround.^4).';
    LHYDensity = (ItDummyVector.*(phiGround.^5)).';
    EDensity = KEDensity + MFDensity + LHYDensity;
    MF = 4*pi*trapz(r,r.^2.*MFDensity);
    KE = 4*pi*trapz(r,r.^2.*KEDensity);
    LHY = 4*pi*trapz(r,r.^2.*LHYDensity);
    E = 4*pi*trapz(r,r.^2.*EDensity);
end