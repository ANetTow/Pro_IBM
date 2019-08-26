function r = calc_resp_sine(mu_cell, daylength, delta_t)
% CALC_RESP_SINE calculates respiration rate of Prochlorococcus with given
% cellular growth rate and daylength under LD conditions with sinusoidal
% light and discrete timesteps. Based on Zinser et al (2009) measurement
% that R = (1/3)*G.
%
% INPUT :
%   mu_cell =   cellular growth rate/per cell photosynthesis rate (d^{-1}) 
%   daylength = length of light period (h)
%   delta_t =   time step of model (d)
%
% OUTPUT :
%   r = respiration rate (d^-1)
%
% Usage:
%   r = calc_resp_sine(mu_cell, daylength, delta_t)
%
% Started:  20/Jun/2013 Annette Hynes, UGA

a = 2*pi*24./(2*daylength);                         % convert to days
N = floor((daylength/24)*(1/delta_t));              % This should be 50 for 12-h day
xi = sin(0.5*(N-1).*a*delta_t).*sin(0.5*(N).*a*delta_t)./sin(0.5*a*delta_t); % Expansion of sum of sines

r = -log(2/3 + (1/3)*exp(-mu_cell*delta_t*xi));