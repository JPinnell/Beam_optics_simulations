function U = HG(X,Y,M,N,weights,w0)
% This function computes a superposition of HG mode at the plane z = 0
% X,Y are matrices, e.g. x = -1:0.1:1, y = -1:0.1:1, [X,Y] = meshgrid(x,y)
% "weights" is a weight vector for the coefficients in the superposition
% E.g. generate mode U = 5*HG_0^-5 + I*HG_0^3 - HG_1^4 + 2*HG_3^4 then,
% p = [0,0,1,3], l = [-5,3,4,4], weights = [5,I,-1,2]

U = zeros(size(X)); %initialise electric field
 
for i = 1:length(weights)
    U = U + weights(i).*(1/w0).*(sqrt(2.^(1-N(i)-M(i))./(pi.*factorial(N(i)).*factorial(M(i)))).*Hermite(M(i),sqrt(2).*X./w0).*Hermite(N(i),sqrt(2).*Y./w0).*exp(-(X.^2+Y.^2)./w0^2));
end

end

function y = Hermite(n,x)
% Computes Hermite polynomial (faster than hermiteh)
y = zeros(size(x));
for m = 0:floor(n/2)
    y = y + (-1)^m/(factorial(m)*factorial(n-2*m)).*(2.*x).^(n-2*m);
end
y = y.*factorial(n);

end