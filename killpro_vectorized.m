function [data_store_out, npro_out, num_cells_grazed] = killpro_vectorized(current_time, data_store, npro, t_step, g)
% KILLPRO_VECTORIZED takes a structure of Prochlorococcus and "grazes" them
% by randomly eliminating a fixed proportion of cells. Based on Brad
% Blythe's original killpro.m.
%
% INPUT:
%   current_time =  current index
%   data_store =    structure containing information on the cells
%   npro =          the number of cells
%   t_step =        the time step size (d)
%   g =             the grazing rate (d^-1)
%
% OUTPUT:
%   pro_out =       the structure of cells after grazing
%
% Usage:
%   [data_store_out, npro_out, num_cells_grazed] = killpro_vectorized(current_time, data_store, npro, t_step, g)
%
% Started:  12/May/2010 Annette Hynes, UGA

% Select cells for grazing
indx = randperm(npro);                          % First x cells in ramdon list are grazed
grazing_rate = g*t_step;
num_cells_grazed = floor(npro*grazing_rate);

if num_cells_grazed > 0
    indx_cells_grazed = indx(1:num_cells_grazed);
    data_store(current_time).size(indx_cells_grazed) = [];
    data_store(current_time).dna(indx_cells_grazed) = [];
    data_store(current_time).time_2(indx_cells_grazed) = [];
    data_store(current_time).mu(indx_cells_grazed) = [];
    data_store_out = data_store;
    npro_out = length(data_store_out(current_time).size);
else
    data_store_out = data_store;
    npro_out = npro;
end