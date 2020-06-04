function U = LG(R,Phi,p,l,weights,w0)
% This function computes a superposition of LG_p^l modes at the plane z=0
% R,Phi are matrices, e.g. x = -1:0.1:1; y = -1:0.1:1; [X,Y] = meshgrid(x,y); [Phi,R] = cart2pol(X,Y);
% "weights" is a weight vector for the coefficients in the superposition
% E.g. generate mode U = 5*LG_0^-5 + I*LG_0^3 - LG_1^4 + 2*LG_3^4 then,
% p = [0,0,1,3], l = [-5,3,4,4], weights = [5,I,-1,2]

U = zeros(size(R)); % initialise electric field

for i = 1:length(weights)
    Norm = sqrt(2*factorial(p(i))/(pi*factorial(p(i)+abs(l(i)))))./w0; % normalisation factor
    Radial = (sqrt(2).*R./w0).^(abs(l(i))).*Laguerre(p(i),abs(l(i)),2.*R.^2./(w0.^2)); % radial term
    Gaussian = exp(-R.^2./w0^2); % Gaussian envelope
    Azimuthal = exp(1i.*l(i).*Phi); % Azimuthal phase
    U = U + weights(i).*Norm.*Radial.*Gaussian.*Azimuthal;
end
end

function y = Laguerre(p,l,x)
% Computes the associated Laguerre polynomials (faster than laguerrel)
y = zeros(size(x));
if p == 0
    y = ones(size(x));
elseif p == 1
    y = 1 + l - x;
else
    for i = 0:p
        y = y + (-1)^i*factorial(p+l)/(factorial(p-i)*factorial(l+i)*factorial(i)).*x.^i;
    end
end
end
