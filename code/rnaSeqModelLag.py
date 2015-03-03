# -*- coding: utf-8 -*-
"""
Real-Time, Single Molecule RNA-seq 

This is a model for cDNA sequencing based on public-domain information 
DNA polymerase kinetics and PacBio's SMRT architecture using ZMW 
confinement.  The basic modeling elements of the SMRT architecture are

(a) Context-averaged DNA polymerase incorporation kinetics with the 
standard deviation of the base incorporation time equal to the mean.

(b) Constant probability DNA photodamage model

(c) Circular template with sense/antisense strands and hairpin adapters.

The RNA-seq application with such a system involves a number of sample 
preparation and instrument protocol steps, including 

(1) reverse-transcriptase conversion to cDNA, 
(2) repair, A-tailing, and adapter ligation
(4) polymerase binding, 
(5) complex loading into single-molecule measurement wells, 
(6) addition of reagents and chelates
(7) lights on, acquisition on, and streaming data compression
(8) postprocessing and analysis

The information structure of the RNA transcripts is a critical aspect of 
this analysis, as is the time delay between steps (6) and (7).  Modeling
of (a)-(c) determines the stochastic dynamics of the situation.
 
Created on Sun Feb 12 10:03:23 2012

@author: akhalak
"""

import numpy as np
import matplotlib.pylab as plt
#import pdb as debug

#import StringIO, Image
#import win32com.client

import samplePrepModelLib as akS

#
# run parameter flags
#
savePlots = True

#
# define grid of readlength tau and of transcript size
#
transcriptVals = np.r_[500:5000:250];
readlengthVals = np.r_[1000:10000:1000];
transcriptVals = [4000]
lagTimeVals = np.array([10, 60, 5*60, 10*60, 15*60])

#lagTimeVals = np.array([15*60])
N = len(transcriptVals)
M = len(readlengthVals)
K = len(lagTimeVals)
fullUtrFracAll = np.zeros((M,N,K))


# Algorithm Summary
#
# 1) The start position is based on a Gaussian with a mean at the
# base rate times the time between reaction start and acquisition start
# The variance on the Gaussian is also that time, based on a random walk
# model.
#
# 2) Convert the distribution from (1) to the required unrolled RL to
# cover both UTRs.  This is based on the insert size, from which the remainder
# required to cover both UTRs follows the following formula.
#
# 3) Compute the unrolled Readlength distribution.  This is based on the
# readlength tau.
#
# 4) Integrate the product of the CDF of required unrolled RL and the PDF
# of the readlength distribution to get the fraction of complete
# transcripts of that size.
#


#
# loop over grid of transcripts lengths and readlengths
#
for n in np.r_[0:N]:
    for m in np.r_[0:M]:
        for k in np.r_[0:K]:
            
            #
            # assign transcript size and readlength, and compute
            # total template length
            #
            transcriptLen = transcriptVals[n]
            photodamageTau = readlengthVals[m]
            halfCycleLen = akS.compHalfCycleLen(transcriptLen)
            RLgrid = np.r_[0:(2*halfCycleLen)]
    
            #
            # compute PDFs for the starting base in normal and hotstart modes.
            # Need an extension (probably using an MCMC model) to have
            # different speeds with and without strand displacement.
            #
            startBase = akS.startBaseModel(RLgrid, transcriptLen, lagTime = lagTimeVals[k]);
    
            #
            # calculate the required number of bases for each member of the RLgrid
            #
            basesToGo = akS.RnaCoverageModel(RLgrid, transcriptLen)
    
            #
            # compute basesToGo PDF on RLgrid
            #
            basesToGoPDF = np.zeros(len(RLgrid))
            for i in np.r_[0:len(RLgrid)]:
                basesToGoPDF[basesToGo[i]] = basesToGoPDF[basesToGo[i]] + startBase[i]
    
            #
            #  build up readlength CDF over same grid.  Since the
            #  ReadlengthModel is expecting a grid that spans the support of
            #  the distribution, it should be longer than what is actually
            #  needed.
            #
            RL_perc = akS.photodamageModel(photodamageTau, RLgrid=np.r_[0:np.max((2*halfCycleLen,photodamageTau*3))])
            RLcdf   = 1 - RL_perc[RLgrid]
    
            #
            # compute hotstart and nonhotstart fractions
            #
            fullUtrFracAll[m,n,k] = np.sum(RLcdf * basesToGoPDF)
            

#
# Non-hotstart summary plots of transcript size vs readlength
#
fullUtrFrac    = np.zeros((M,N))
for k in np.r_[0:K]:

    fullUtrFrac = fullUtrFracAll[:,:,k]
    plt.close(k+1)
    plt.figure(k+1)
    plt.clf()
    plt.contourf(transcriptVals, readlengthVals, np.abs(fullUtrFrac), levels=np.r_[0:1:0.1])
    plt.xlabel('Transcript Coding Region Length')
    plt.ylabel('Readlength Tau')
    plotlabel = 'Fraction of Transcripts Read, Lag = %.0fsec' % lagTimeVals[k]
    plt.title(plotlabel)
    plt.colorbar()
    if savePlots:
        fig = plt.gcf()    
        filename = 'Contour_lag%.0f.png' % lagTimeVals[k] 
        fig.savefig(filename, format='png')


plt.show()
