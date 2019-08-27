function [data_store, av_cell_size, av_cell_dna, tm, mu_pop, cells, S, G2, index] = Pro_IBM(mu_cell, T_S,...
    T_G2, daylength, light_regime, Cg, Ps_width, Ps_zero )
% PRO_IBM runs the Prochlorococcus diel cycling model .
% 
% Based on Brad Blythe?s driver_July_2_2008.m in model_V_0.3.
%
% INPUT:
%   mu_cell =       Cellular growth rate (d^(?1))
%   T_S =           Duration of S phase (d)
%   T_G2 =          Duration G2 phase (d)
%   daylength =     Length of daylight , sunrise to sunset (h)
%   light_regime =  String denoting whether light is binary ('binary'), sinusoidal ('sine'), or constant ('constant') 
%   Cg =            Circadian gate (d)
%   Ps_width =      Parameter that controls the "width" of the probability function for cells entering S phase based on their size (fg C). Typically 85.
%   Ps_zero =       Parameter that states where the probability for cells entering S is zero. Typically 45.
%
% OUTPUT:
%   data_store =    Structure containing cellular information
%   av_cell_size =  Average cell size (fg C per cell )
%   av_cell_dna =   Average DNA per cell (genome copies per cell )
%   tm =            Time vector (d)
%   mu_pop =        Population growth rate calculated from final days when grazing is turned off 
%   cells =         number of cells
%   S =             fraction of cells in S phase (1 < dnapercell < 2) 
%   G2 =            fraction of cells in G2 phase (dna = 2)
%   index =         indices for the last three days
%
% Other files used: 
%   initialize_pro.m
%   makepro_vectorized.m 
%   killpro_vectorized.m 
%   light_sine.m 
%   exp_growth.m 
%   calc_resp_bin.m 
%   calc_resp_sine.m 
%   calc_resp_cont.m
%
% Usage:
%   [data_store, av_cell_size, av_cell_dna, tm, mu_pop, cells, S, G2, index] = Pro_IBM(mu_cell, T_S, T_G2, daylength, light_regime, Cg, Ps_width, Ps_zero);
%
% Started: 24/Jun/2013 Annette Hynes, UGA

tic                                                         % Start the clock

% Initialize data structure
if strcmp(light_regime, 'binary') || strcmp(light_regime, 'sine')  % LD
    [t_span, t_step, data_store] = initialize_pro(mu_cell);
elseif strcmp(light_regime, 'constant')                             % LL
    [t_span, t_step, data_store] = initialize_pro_constant(mu_cell);
end

% Initialize storage arrays
cnt = 1;                                                    % Timestep counter
n_steps = t_span/t_step;

npro = length(data_store(1).size);                          % Number of cells          

G2_cells = zeros(n_steps + 1, 1);                           % Have 2+ copies of DNA
S_cells = zeros(n_steps + 1, 1);                            % Duplicating: 1 < DNA > 2

av_cell_size = zeros(n_steps + 1,1); 
av_cell_size(1) = mean([data_store(1).size]); 
av_cell_dna = zeros(n_steps + 1,1); 
av_cell_dna(1) = mean([data_store(1).dna]);

cells = zeros(n_steps + 1, 1);                              % number of cells
cells(1) = npro;
lt = zeros(n_steps + 1, 1);                                 % Light
tm = zeros(n_steps + 1, 1);                                 % Time
grazed = zeros ( n_steps + 1, 1);   

% Time stepping through the model
for t = t_step:t_step:t_span                                % t in days
    cnt = cnt + 1;
    
    % Initialize the next step of the data structure
    data_store(cnt).size = data_store(cnt - 1).size; 
    data_store(cnt).dna = data_store(cnt - 1).dna; 
    data_store(cnt).time_2 = data_store(cnt - 1).time_2; 
    data_store(cnt).mu = data_store(cnt - 1).mu;

    % Light
    lighttime = single(t - fix(t));                         % Time since dawn (d)
    if strcmp(light_regime, 'binary')                       % Binary light
        light = light_binary(t, daylength);
        p1 = 0.002866*daylength - 0.01245;                  % Parameters relating daylength and light field to growth rate
        p2 = 0.0299*daylength - 0.01857;
        t_zero = single(Cg);
        r = calc_resp_bin(mu_cell , daylength);             % Respiration rate
    elseif strcmp(light_regime, 'sine')                     % Sinusoidal light level
        light = light_sine(t ,daylength);
        p1 = 0.0008077*daylength + 0.0002555;               % Parameters relating daylength and light field to growth rate

        p2 = 0.0191*daylength - 0.01667;
        t_zero = single(Cg);
        r = calc_resp_sine(mu_cell, daylength, t_step);     % Respiration rate
    elseif strcmp(light_regime , 'constant')                % Constant light
        light = light_constant(t, daylength); 
        p1 = 0.05377;                                       % Parameters relating daylength and light field to growth rate
        p2 = 0.02665;
        t_zero = 0;                                         % No circadian gate! 
        r = calc_resp_const(mu_cell);                       % Respiration rate
    end
    
    mu_pop_est = p1*mu_cell^2 + p2*mu_cell; 
    npro_now = npro;
    
    % Phase fraction data
    this_dna = data_store(cnt - 1).dna;
    G2_cells(cnt - 1) = length(this_dna(this_dna >= 2));
    S_cells(cnt - 1) = length(this_dna(this_dna > 1 & this_dna < 2));
    
    % Growth
    data_store(cnt).size = exp_growth(data_store(cnt - 1).size, data_store(cnt - 1).mu, t_step , light);    % exponential

    % Respiration
    data_store(cnt).size = data_store(cnt).size - r.*t_step.*data_store(cnt).size; 
    
    % Determine which cells will replicate DNA
    chosen = [];

    size = data_store(cnt).size(:);
    p_enterS = exp((size - Ps_zero )/( Ps_width/log(2))) - 1;   % probability of entering S increases with size .
    p_enterS(p_enterS > 1) = 1;                             % can?t have probability > 1
 
    if lighttime >= t_zero                                  % Circadian gate
        enterS_dice_roll = rand(npro , 1);                  % Roll the S dice, based on a uniform distribution 
        ind_chosen = find (p_enterS >= enterS_dice_roll);   % Determine which cells are chosen to enter
        not_dbl = find(data_store(cnt - 1).dna(:) < 1.001 & data_store(cnt - 1).dna(:) > 0.999);    % Only 1 copy of DNA
        chosen = ind_chosen(ismembc(ind_chosen, not_dbl));  % cells actually chosen to initiate DNA replication
    end

    dna_plus = find((data_store(cnt - 1).dna(:) - 1) > 1e-6);   % Already started replicating
    indx_no_rep = setxor(1:npro_now, [chosen; dna_plus]);       % Won?t replicate even though they qualify

    % Increase DNA
    DNA_inc = t_step /T_S; 
    if ~isempty(chosen)                                     % fraction increase in DNA per timestep given T_S
        data_store(cnt).dna(chosen) = data_store(cnt - 1).dna(chosen) + (DNA_inc);  % Initialize DNA increase
    end
    
    data_store(cnt).dna(dna_plus) = data_store(cnt - 1).dna(dna_plus) + (DNA_inc);  % Continue DNA increase
    data_store(cnt).dna(indx_no_rep) = data_store(cnt - 1).dna(indx_no_rep);    % DNA remains the same
     
    % Cell division
    two = find(data_store(cnt).dna >= 2);                   % Finished replicating DNA
    index_time2 = find(data_store(cnt).time_2 == 0);        % Haven?t previously entered G2
    index_uptime = two(ismembc(two, index_time2));          % Fit both criteria
    data_store(cnt).time_2(index_uptime) = t;               % Time they reached G2
    data_store(cnt).dna(two) = 2;                           % Not greater than two copies of DNA

    lag_test = t - data_store(cnt).time_2(two);             % Only test those cells with more than 2 copies of DNA
    lag_ok = lag_test >= T_G2;                              % Find cells that have completed G2 phase
    dividing = two(lag_ok);                                 % Index of dividing cells

    if ~isempty(dividing)
        data_store = makepro_vectorized(data_store, dividing, cnt); % Divide by binary fission
    end
    
    npro = npro + length(dividing);                         % Update cell number to include newly divided daughters
    
    non_zero_size = nnz(data_store(cnt).size);
    non_zero_dna = nnz(data_store(cnt).dna);
    
    av_cell_size(cnt) = sum(data_store(cnt).size)/non_zero_size;
    av_cell_dna(cnt) = sum(data_store(cnt).dna)/non_zero_dna;   
    
    % Grazing
    grazing_rate = mu_pop_est; 
    
    if tm <= t_span - 6
        [data_store, npro, n_grazed] = killpro_vectorized(cnt, data_store, npro, t_step, grazing_rate);
    else                                                    % No grazing during the last few days
        n_grazed = 0;
    end
    
    grazed ( cnt ) = n_grazed ;

    % Update arrays
    lt(cnt) = light; 
    tm(cnt) = t; 
    cells(cnt) = npro;
    
    % Save memory
    if cnt > 3
    data_store(cnt - 3).size = []; 
    data_store(cnt - 3).dna = []; 
    data_store(cnt - 3).mu = []; 
    data_store(cnt - 3).time_2 = [];
    end
end

% Percent cells in each cell cycle phase
G2 = G2_cells./cells; 
S = S_cells./cells;

index = (tm >= t_span - 4 & tm <= t_span - 1);              % Last few days

% Calculate mu_pop
t1 = tm(end) - 4 + Cg - 0.05; 
dif_t1 = tm - t1;
N1 = cells(abs(dif_t1) == min(abs(dif_t1)));

t2 = tm(end) - 1 + Cg - 0.05; 
dif_t2 = tm - t2;
N2 = cells(abs(dif_t2) == min(abs(dif_t2))); 
mu_pop = log(N2/N1)/(t2 - t1);

toc % Stop the timer