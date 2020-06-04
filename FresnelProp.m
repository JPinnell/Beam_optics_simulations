function Out = FresnelProp(Input, dx, lambda, z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out = FresnelProp(Input, dx, lambda, z)
% v1 J.Pinnell 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Descrition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs Fresnel propagation of an input field.
% Specifically, given a field U(x,y,0), it finds the field at U(x,y,z).
% It actually uses the Angular Spectrum method of propagation under the 
% Fresnel approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input  - input field (must be a square matrix)
% dx     - input sample size (pixel size)
% lambda - wavelength of light
% z      - propagation distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out - field at the output plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make spatial frequency coordinates (sample at Nyquist frequency)
[N,~] = size(Input);                   % Input must be square matrix
L = N*dx;                              % side length of Input
kx = -1/(2*dx) : 1/L : 1/(2*dx) - 1/L; % frequency coordinates sampled at Nyquist frequency
[KX, KY] = meshgrid(kx,-kx);           % for vectorisation

% Propagate using Angular spectrum method
G = fft2(fftshift(Input));                           % compute angular spectrum of input field
psi = fftshift(exp(-1i*pi*lambda*z.*(KX.^2+KY.^2))); % transfer function of plane waves
U = G.*psi;                                          % propagate plane waves in frequency space
Out = ifftshift(ifft2(U));                           % revert back to spatial coordinates
end