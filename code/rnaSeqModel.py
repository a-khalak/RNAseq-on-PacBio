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
import pdb as debug

import samplePrepModelLib as akS

#
# parameter Flags
#
showPlots = True
savePlots = True

#
# define grid of readlength tau and of transcript size
#
transcriptVals = np.r_[500:5000:500]
readlengthVals = np.r_[1000:10000:1000]
transcriptVals = [3500]
#readlengthVals = [2000]
N = len(transcriptVals)
M = len(readlengthVals)
fullUtrFrac = np.zeros((M,N))
fullUtrFracHot = np.zeros((M,N))

hsm = akS.hsm
hsm['nonHotstartLagTime'] = 900
hsm['HotstartLagTime'] = 10


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

#if showPlots:
#    print "Click Upper Left to Quit, Upper Right to Advance, Lower Plot to Save"

#
# loop over grid of transcripts lengths and readlengths
#
for n in np.r_[0:N]:
    for m in np.r_[0:M]:
        
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
        startBase = akS.startBaseModel(RLgrid, transcriptLen, lagTime = hsm['nonHotstartLagTime']);
        hotstartBase = akS.startBaseModel(RLgrid, transcriptLen, lagTime = hsm['hotstartLagTime']);

        #
        # calculate the required number of bases for each member of the RLgrid
        #
        basesToGo = akS.RnaCoverageModel(RLgrid, transcriptLen)

        #
        # compute basesToGo PDF on RLgrid
        #
        basesToGoPDF = np.zeros(len(RLgrid))
        basesToGoPDFHot = np.zeros(len(RLgrid))
        for i in np.r_[0:len(RLgrid)]:
            basesToGoPDF[basesToGo[i]] = basesToGoPDF[basesToGo[i]] + startBase[i]
            basesToGoPDFHot[basesToGo[i]] = basesToGoPDFHot[basesToGo[i]] + hotstartBase[i]

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
        fullUtrFrac[m,n] = np.sum(RLcdf * basesToGoPDF)
        fullUtrFracHot[m,n] = np.sum(RLcdf * basesToGoPDFHot)

        if showPlots:
            plotlabel = "cDNA/RNA-seq Model, CodingLen=%.0f, ReadlengthTau=%.0f " % (transcriptLen, photodamageTau)
    
            plt.figure('Single experiment')
            plt.clf()
            plt.subplot(211);
            plt.plot(RLgrid, startBase, 'b')
            ax1 = plt.gca()
            plt.plot (RLgrid, hotstartBase, 'r')
            plt.xlabel('start position')
            plt.ylabel('probability')
            plt.ylim((0, np.max(startBase)*2.5))
            ax2 = plt.twinx(ax1)
            plt.plot(RLgrid, basesToGo[RLgrid], 'g')
            plt.ylabel('Sequencing Bases Required for both UTRs')
    
            # plot some indicator lines for primer locations
            plt.axes(ax1)
            forwPrSite = max(RLgrid) * np.array([1, 1])
            revPrSite = halfCycleLen * np.array([1, 1])
            probrange = np.array([0, np.max(hotstartBase)])
            plt.plot(np.array([0, 0]), probrange, 'k--')
            plt.plot(forwPrSite, probrange, 'k--')
            plt.plot(revPrSite, probrange, 'k--')
            plt.title(plotlabel)
    
            plt.subplot(212)
            plt.plot(RLgrid, basesToGoPDF)
            ax1 = plt.gca()
            plt.plot(RLgrid, basesToGoPDFHot, 'r')
            plt.xlabel('Min. Readlength Required to Hit Both UTRs')
            plt.ylabel('Probability Density')
            plt.ylim((0, np.max(basesToGoPDF)*2.5))
            ax2 = plt.twinx(ax1)
            plt.plot(RLgrid, RLcdf)
            plt.ylabel('Readlength CDF')
            plt.ylim((0, 1))
            plt.show()

            if savePlots:
                filename = 'Detail_Read%.0f_Trans%.0f.png' %  (photodamageTau, transcriptLen)
                fig = plt.gcf()
                fig.savefig(filename, format='png')

    
            #
            # click the top plot to advance and bottom plot to save
            #
            if 0:
                debug.set_trace()
                x = plt.ginput()            
                print (x[0])
                if (x[0][1] < 1 and savePlots):
                    fig = plt.gcf()    
                    filename = 'Detail_Read%.0f_Trans%.0f.png' %  (photodamageTau, transcriptLen)
                    fig.savefig(filename, format='png')
                elif ((x[0][1]>1) & (x[0][0]<1000)):
                    raise Exception ("quitting")
            elif 0:
                x = raw_input("(q)uit, (a)dvance, or (s)ave: ")
                if (x == "s" or x == "S"):
                    fig = plt.gcf()
                    filename = 'Detail_Read%.0f_Trans%.0f.png' %  (photodamageTau, transcriptLen)
                    fig.savefig(filename, format='png')
                elif (x == "q" or x == "Q"):
                    raise Exception  ("quitting")



