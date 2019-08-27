# Pro_IBM
An individual-based model simulating the cell cycle of the marine cyanobacterium *Prochlorococcus*.

## Introduction

*Prochlorococcus* is the smallest but mightiest of the phytoplankton, responsible for about 20% of global primary production. Their daily cycle of growth and division can be used to calculate division rates of field populations, and this individual-based model (IBM) simulates *Prochlorococcus* diel patterns and estimates cell cycle parameters including division rate (Hynes et al, 2015a; Hynes et al, 2015b).

<img src = https://github.com/ANetTow/Pro_IBM/blob/master/Pro_IBM_flowchart.png title="Pro IBM Flowchart" align="left" style="float" width="300">*Flowchart of model events.  From Hynes et al (2015a), Fig. 1.*

The IBM is fully described in Hynes et al (2015a). The flowchart shows the processes and decisions in the model.  All cells grow and respire.  Cells that are large anough and have been randomly selected can begin DNA replication, entering S phase.  Cells that have completed DNA replication and have endured the duration of G2 can divide. All cells can be randomly selected to be grazed, and the grazing rate is chosen to balance population growth so the model does not become numerically overwhelmed.  

The model is run for 30 days.  During the last 3 days, grazing is turned off for analysis.

## Functions

## Variables

## References

-Hynes A. M., B. J. Blythe, and B. J. Binder (2015a).  "An individual-based model for the analysis of *Prochlorococcus* diel cycle behavior," *Ecol. Model.* **301**:1 - 15, doi:10.1016/j.ecolmodel.2015.01.011.

-Hynes A. M., K. L. Rhodes, and B. J. Binder (2015b).  "Assessing cell cycle-based methods of measuring *Prochlorococcus* division rates using an individual-based model," *Limnol. Oceanogr. Methods* **13**:640 - 650, doi:10.1002/lom3.10054.
