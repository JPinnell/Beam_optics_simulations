function Z = Zernike(n,m,rho,phi)
% Computes the single Zernike phase map Z_n^m (see wikipedia "Zernike polynomials")
% We must have n>=0 and n>=|m|. Also m = -n,-n+2,...,n-2,n otherwise Z = 0
% rho,phi are meshgrid polar coordinates

[V,H] = size(rho); % (row,col) dimensions 

% normalise rho to unit disk (take largest possible circle)
if V > H % if #rows > #cols
    rho = rho./rho(V/2+1,H); 
else
    rho = rho./rho(V,H/2+1);
end

% Compute radial polynomials R_n^m
R = zeros(V,H); % initialise
m1 = abs(m);
if mod(n-m1,2) == 0 % n - m must be even, else R is 0
    for k = 0:((n-m1)/2)
        R = R + (-1)^k*nchoosek(n-k,k)*nchoosek(n-2*k,(n-m1)/2-k).*rho.^(n-2*k); 
    end
end

if m >= 0
    Z = sqrt(2*n+2).*R.*cos(m1.*phi);
else
    Z = sqrt(2*n+2).*R.*sin(m1.*phi);
end
Z(rho>1) = 0;

end
