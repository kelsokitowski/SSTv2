[EX1_1,EX2_1,Q_1,f1_1,f2_1,f3_1]= ETDRKcoefficients(-2.0*kVals.^2*nu,kVals,h);
[EX1_2,EX2_2,Q_2,f1_2,f2_2,f3_2]= ETDRKcoefficients(-2.0*kVals.^2*nu,kVals,h);
[EX1_3,EX2_3,Q_3,f1_3,f2_3,f3_3]= ETDRKcoefficients(-2.0*kVals.^2*nu,kVals,h);
[EX1_4,EX2_4,Q_4,f1_4,f2_4,f3_4]= ETDRKcoefficients(-2.0*kVals.^2*D,kVals,h);
[EX1_5,EX2_5,Q_5,f1_5,f2_5,f3_5]= ETDRKcoefficients(-2.0*kVals.^2*D,kVals,h);
[EX1_6,EX2_6,Q_6,f1_6,f2_6,f3_6]= ETDRKcoefficients(-1.0*kVals.^2*(D+nu),kVals,h);
function [EX1,EX2,Q,f1,f2,f3]= ETDRKcoefficients(linearOperator,kVals,h)
L = linearOperator;
N = length(kVals); %number of kVals
%L = -2.0*kVals.^2*nu; %linear operator (excluding time derivative)
EX1 = exp(h*L); EX2 = exp(h*L/2); %surprisingly these are the only two exponentials
%%...they will be part of a messy algebraic system, but still

M = 16; %number of points to use for the periodic trapezoidal rule for poles.
r = exp(1i*pi*((1:M)-0.5)/M)'; %16 points on the unit circle in the complex plane "roots of unity"
%there will be a pole at z = L. Use cauchy integral formula

for kj = 1:N
    for mj = 1:M
    LR(mj,kj)= h*L(kj)+r(mj); %generates h*L+r for each k, doing 16 r values to be used in periodic trapezoidal rule
    end
end


%ETDRK4 coefficients
Q = h*real(mean( (exp(LR/2)-1.0)./LR ,1));
f1 = h*real(mean( (-4.0-LR+exp(LR).*(4.0-3.0*LR+LR.^2.0))./LR.^3.0 ,1));
f2 = h*real(mean( (2.0+LR+exp(LR).*(-2.0+LR))./LR.^3.0 ,1));
f3 = h*real(mean( (-4.0-3.0*LR-LR.^2.0+exp(LR).*(4.0-LR))./LR.^3.0 ,1));
end
