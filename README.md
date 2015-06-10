# RNAseq-on-PacBio

Real-Time, Single Molecule RNA-seq 

This is a model for cDNA sequencing based on public-domain information 
DNA polymerase kinetics and PacBio's SMRT architecture using ZMW 
confinement.  The basic modeling elements of the SMRT architecture are

* Context-averaged DNA polymerase incorporation kinetics with the 
standard deviation of the base incorporation time equal to the mean.
* Constant probability DNA photodamage model
* Circular template with sense/antisense strands and hairpin adapters.
The RNA-seq application with such a system involves a number of sample 
preparation and instrument protocol steps, including 
** reverse-transcriptase conversion to cDNA, 
** repair, A-tailing, and adapter ligation
** polymerase binding, 
** complex loading into single-molecule measurement wells, 
** addition of reagents and chelates
** lights on, acquisition on, and streaming data compression
** postprocessing and analysis
The information structure of the RNA transcripts is a critical aspect of 
this analysis, as is the time delay between steps (6) and (7).  Modeling
of (a)-(c) determines the stochastic dynamics of the situation.
