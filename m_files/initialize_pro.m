function [t_span, t_step, pro_structure] = initialize_pro(growthrate)
% INITIALIZE_PRO initializes the Prochlorococcus model. Based on Brad
% Blythe's paramters.m and makepro.m. 
% 
% INPUT:
%   growthrate =    cellular exponential growth rate (Malthusian parameter; d^-1)
%
% OUTPUT:
%   t_span =        Run time for model (d; scalar)
%   t_step =        Time step for model (d)
%   pro_structure = Structure to store cellular information
%
% Usage:
%   [t_span, t_step, pro_structure] = initialize_pro(growthrate)
%
% Started:  14/May/2010 Annette Hynes, UGA

% TIME
t_span = 31;                                                            % Run time for model (d)
t_step = 0.01;                                                          % Time step (d; = 14.4 min) 
n_steps = t_span/t_step;                                                % # timesteps

% PROCHLOROCOCCUS
pro_structure = struct('size', {}, 'dna', {}, 'mu', {}, 'time_2', {});  % Structure holding cell data
pro_structure = repmat(pro_structure, n_steps + 1, 1);                  % Structure to store info at each time step

npro = 10000;                                                           % Initial number of cells
pro_structure(1).mu = growthrate*ones(npro, 1);

av_size = 50;                                                           % fg C per cell
spread = 20;
pro_structure(1).size = abs(av_size*ones(npro, 1) + spread*randn(npro, 1));
pro_structure(1).dna = 1*ones(npro, 1);                                 % All in G1 (genome copies)
pro_structure(1).time_2 = zeros(npro, 1);                               % Time when cell reaches dna = 2, used to determine time lag for division after replication (d) 
