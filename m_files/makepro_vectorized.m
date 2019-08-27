function data_store = makepro_vectorized(data_store, index_div, current_time)
% MAKEPRO_VECTORIZED takes cells that are eligible to divide, performs
% binary fission, and passes the cells back. Based on Brad Blythe's
% original makepro.m.
%
% INPUT:
%   data_store =    structure containing cellular information 
%   index_div =     indices of the cells eligible to divide 
%   current_time =  index of current timepoint
%
% OUTPUT:
%   data_store =    updated structure of cellular information
%
% Usage:
%   data_store = makepro_vectorized (data_store, index_div, current_time)
%
% Started:  12/May/2010 Annette Hynes, UGA

current_size = data_store(current_time).size;                                   % Size before dividing
n = length(current_size);                                                       % Number of cells before dividing
n_div = length(index_div);                                                      % Number of cells dividing
new_cells = [index_div; (n + 1:n + n_div)'];                                    % indices of new cells, second daughter cell tacked on to the end

data_store(current_time).size(index_div) = current_size(index_div)*0.5;                 % Daughter cell #1
data_store(current_time).size(end + 1: end + n_div) = current_size(index_div) * 0.5;    % Daughter cell #2 added on the end

data_store(current_time).dna(new_cells) = 1;                                    % One copy of DNA
data_store(current_time).time_2(new_cells) = 0;                                 % Reset the dividing clock
data_store(current_time).mu(1:n + n_div) = data_store(current_time - 1).mu(1);  % Constant growth rate
