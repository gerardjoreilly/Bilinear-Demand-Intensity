function [lam_bl, lam1, lam2]=maf_bilinear_EDPbased(C,a1,a2,b1,b2,slim,beta_DR,beta_DU,beta_CR,beta_CU,k0,k1,k2)
% [lam_bl, lam1, lam2]=maf_bilinear_EDPbased(theta_c,a1,a2,b1,b2,slim,beta_DR,beta_DU,beta_CR,beta_CU,k0,k1,k2)
% This function computes the MAF of exceedance of a limit state given an
% EDP theta_c with the 2nd order hazard curve as per Vamvatsikos [2013] and
% the bilinear demand-intensity model proposed by O'Reilly
% 
% Inputs:
% C:            median capacity
% a1:           EDP-IM relation term (Part 1)
% b1:           EDP-IM relation term (Part 1)
% a2:           EDP-IM relation term (Part 2)
% b2:           EDP-IM relation term (Part 2)
% slim:         Limit state intensity
% beta_DR:      aleatory uncertainty associated with demand
% beta_DU:      epistemic uncertainty associated with demand
% beta_CR:      aleatory uncertainty associated with capacity
% beta_CU:      epistemic uncertainty associated with capacity
% k0:           second order hazard model parameter
% k1:           second order hazard model parameter
% k2:           second order hazard model parameter
%
% Output
% lam_bl:       MAF of exceedance using bilinear model
% lam1:         MAF of exceedance using part 1 linear model
% lam2:         MAF of exceedance using part 2 linear model

% Compute the total uncertainty
beta_T=sqrt(beta_DR^2+beta_DU^2+beta_CR^2+beta_CU^2);

% Check the continuity of the bilinear model
Clim1=a1*slim^b1;
Clim2=a2*slim^b2;
if abs(Clim1-Clim2)>0.005
    error('Lack of continuity in blinear demand-intensity model, refit the coefficients a and b');
end

% Compute the median intensity
smed1=(C/a1)^(1/b1);
smed2=(C/a2)^(1/b2);
H1=k0*exp(-k1*log(smed1)-k2*log(smed1)^2);
H2=k0*exp(-k1*log(smed2)-k2*log(smed2)^2);

% Compute the phi' term (Equation 30)
phi1=1/(1+2*k2*beta_T^2/(b1^2));
phi2=1/(1+2*k2*beta_T^2/(b2^2));

% Compute the lognormal distribution parameters and 'weights'
mu1=phi1*((log(C)-log(a1))/b1-k1*beta_T^2/(b1^2)); % Equation 33
mu2=phi2*((log(C)-log(a2))/b2-k1*beta_T^2/(b2^2));
sigma1=beta_T^2*sqrt(phi1)/b1; % Equation 32
sigma2=beta_T^2*sqrt(phi2)/b2;
F1=logncdf(slim,mu1,sigma1); % Equation 34
F2=logncdf(slim,mu2,sigma2);

G1=sqrt(phi1)*k0^(1-phi1)*H1^phi1*exp(k1^2*phi1*beta_T^2/(2*b1^2)); %Equation 29
G2=sqrt(phi2)*k0^(1-phi2)*H2^phi2*exp(k1^2*phi2*beta_T^2/(2*b2^2));

% Compute the lambda terms
lam_bl=F1*G1+(1-F2)*G2;
lam1=G1;
lam2=G2;
