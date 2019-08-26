function light = light_sine(time, daylength)
% LIGHT_SINE gives light output as a sine wave, depending on the time of
% day. 
%
% INPUT :
%   time =  current time (d); 0 = dawn, fraction of a day daylength =
%           length of light period (h) 
%
% OUTPUT :
%   light = 0 for night
%                 1 for high noon
%
% Usage:
%   light = light_sine(time, daylength)
%
% Started:  15/Jul/2010 Annette Hynes, UGA

lighttime = time - fix(time);
n = length(lighttime);
light = zeros(n, 1);
c = 1/(2*daylength/24);         % Coefficient to get correct daylength

% dawn at 0, noon at daylength/2, dusk at daylength

light(lighttime <= daylength/24) = sin(2*pi*c*lighttime(lighttime <= daylength/24));