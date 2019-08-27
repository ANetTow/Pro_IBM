function cell_next = exp_growth(cell_now, mu, time, light)
% EXP_GROWTH increases the size of cells using a smooth exponential
% function and according to current size, growth rates, and light levels.  
%
% INPUT:
%   cell_now =  the cell size before growing (fg C; scalar or vector)
%   mu =        exponention growth rate (Malthusian parameter; d^-1; scalar or vector)
%   time =      amount of time passed (d)
%   light =     light level from 0 to 1
%
% OUTPUT:
%   cell_next = the cell size after growing (same dimensions as cell_now)
% 
% Usage:
%   cell_next = exp_growth(cell_now, mu, time, light)
%
% Started: 09/Jun/2010 Annette Hynes, UGA

cell_next = cell_now.*(exp((mu.* time)*light)); % exponential growth scaled by light level