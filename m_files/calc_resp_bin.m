function r = calc_resp_bin(mu_cell, daylength)
% CALC_RESP_BIN alculates respiration rate of Prochlorococcus with given
% cellular growth rate and daylength under LD conditions with binary light.
% Based on Zinser et al (2009) measurement that R = 1/3 G
%
% r = calc_resp(mu_cell, daylength)
%
% INPUT:
%   mu_cell = cellular growth rate/ per cell photosynthesis rate (d^{-1})
%   daylength = length of light period (h)
%
% OUTPUT:
%   r = respiration rate (d^-1)
% 
% Started:  20/Jun/2013 Annette Hynes, UGA
% Modified: 

d = daylength/24;       % convert to days

%r = log((1/3)*exp(-mu_cell*d) + 2/3)./log(exp(-d) + exp(-d - 1) - exp(-2*d));
r = -log(2/3 + (1/3)*exp(-mu_cell*d));