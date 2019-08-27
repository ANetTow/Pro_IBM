# Pro_IBM
An individual-based model simulating the cell cycle of the marine cyanobacterium *Prochlorococcus*.

## Introduction

*Prochlorococcus* is the smallest but mightiest of the phytoplankton, responsible for about 20% of global primary production. Their daily cycle of growth and division can be used to calculate division rates of field populations, and this individual-based model (IBM) simulates *Prochlorococcus* diel patterns and estimates cell cycle parameters including division rate (Hynes et al, 2015a; Hynes et al, 2015b).

<img src = https://github.com/ANetTow/Pro_IBM/blob/master/Pro_IBM_flowchart.png title="Pro IBM Flowchart" align="left" style="float" width="300">*Flowchart of model events.  From Hynes et al (2015a), Fig. 1.*

The IBM is fully described in Hynes et al (2015a). The flowchart shows the processes and decisions in the model.  All cells grow and respire.  Cells that are large anough and have been randomly selected can begin DNA replication, entering S phase.  Cells that have completed DNA replication and have endured the duration of G2 can divide. All cells can be randomly selected to be grazed, and the grazing rate is chosen to balance population growth so the model does not become numerically overwhelmed.  

The model is run for 30 days.  During the last 3 days, grazing is turned off for analysis.

<br clear="left"/>

## Functions

The model is written in MATLAB (devloped in R2012b, still runs in R2018b).  Download the files in the [m_files](https://github.com/ANetTow/Pro_IBM/tree/master/m_files) directory and either add them to your MATLAB search path OR put them in your working directory.  Variables and usage for each function can be found by calling `help function_name` at the prompt.

* `Pro_IBM` runs the model with user-defined input parameters.

* `initialize_pro` initializes the structure that stores cell data and the starting population.

* `makepro_vectorized` takes cells ready to divide and turns them into two daughter cells.

* `killpro_vectorized` randomly grazes cells at a fixed proportion.

* `light_sine` simulates light as a truncated sine wave given the time of day and the day length.

* `exp_growth` increases the size of cells according to current size, cellular growth rates, and light levels.

* `calc_resp_bin` calculates cellular respiration rate according to cellular growth rate and daylength under LD conditions with binary light.

* `calc_resp_sine` calculates cellular respiration rate according to cellular growth rate and daylength under LD conditions with sinusoidal light.

* `calc_resp_cont` calculates cellular respiration rate under continuous light conditions.

## Variables

Cell cycle parameters and other variables are defined by the user:  

`[data_store, av_cell_size, av_cell_dna, tm, mu_pop, cells, S, G2, index] = Pro_IBM(mu_cell, T_S, T_G2, daylength, light_regime, Cg, Ps_width, Ps_zero)`

Input/output variables and reasonable ranges for *Prochlorococcus* are summarized in the table below:

Input:

|Name|    Description|                  Values |
|---|     ---|                          ---|
|mu_cell| Maximum cellular growth rate| 1 - 3.5 d^(-1)|
|T_S|     Duration of S phase|          0.1 - 0.33 d|
|T_G2|    Duration G2 phase|            0.05 - 0.30 d|
|daylength|Length of daylight, sunrise to sunset|10 - 14 h|
|light_regime|String denoting whether light is binary, sinusoidal, or constant| 'binary', 'sine', or 'constant'|
|Cg|      Circadian gate|               0.1 - 0.35 d|
|Ps_width|Parameter that controls the "width" of the probability function for cells entering S phase based on their size| 85 fg C|
|Ps_zero|Parameter that states where the probability for cells entering S is zero|45|

Output:

  - `data_store` =    Structure containing cellular information
  - `av_cell_size` =  Average cell size (fg C per cell )
  - `av_cell_dna` =   Average DNA per cell (genome copies per cell )
  - `tm` =            Time vector (d)
  - `mu_pop` =        Population growth rate calculated from final days when grazing is turned off 
  - `cells` =         number of cells
  - `S` =             fraction of cells in S phase (1 < dnapercell < 2) 
  - `G2` =            fraction of cells in G2 phase (dna = 2)
  - `index` =         indices for the last three days

## References

- Hynes A. M., B. J. Blythe, and B. J. Binder (2015a).  "An individual-based model for the analysis of *Prochlorococcus* diel cycle behavior," *Ecol. Model.* **301**:1 - 15, [doi: 10.1016/j.ecolmodel.2015.01.011](https://doi.org/10.1016/j.ecolmodel.2015.01.011).

- Hynes A. M., K. L. Rhodes, and B. J. Binder (2015b).  "Assessing cell cycle-based methods of measuring *Prochlorococcus* division rates using an individual-based model," *Limnol. Oceanogr. Methods* **13**:640 - 650, [doi: 10.1002/lom3.10054]( https://doi.org/10.1002/lom3.10054).
