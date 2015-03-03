# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 01:44:20 2012

@author: akhalak
"""

import numpy as np
import scipy
import scipy.stats.distributions as dist
import matplotlib.mlab as mlab

#
# define polymerase properties in a single parameter structure
#
hsm = dict()
hsm['polymerase speed'] = 2.2
hsm['nonHotstartLagTime'] = 60*15    # fifteen minutes of lag
hsm['hotstartLagTime'] = 10
hsm['movieTime'] = 50*60             # 50 minutes of movie
hsm['fivePrUtrlen'] = 130
hsm['threePrUtrlen'] = 250
hsm['hairpin'] = 44
hsm['minMapBases'] = 80


#
# photodamage model is an exponential that is movie limited
#
def photodamageModel(tau, acqLength=hsm['movieTime'], 
                     baseRate=hsm['polymerase speed'], RLgrid=np.r_[0:10000]):
    readLen = dist.expon.pdf(RLgrid, scale=tau)
    acqReadMax = baseRate * acqLength
    readLen[RLgrid > acqReadMax] = 0
    normalization = scipy.integrate.trapz(readLen, RLgrid)
    readLen = readLen / normalization
    readLen = np.cumsum(readLen)
    return readLen

#
# develop functions that map unrolled start position to length required
# to hit both UTRs.
#
def forw_catchUtr(x,l,fivePrUtrlen=hsm['fivePrUtrlen'], minMapBases=hsm['minMapBases']):
    return fivePrUtrlen - x + l + minMapBases
        
def forw_missUtr(x,l, fivePrUtrlen=hsm['fivePrUtrlen'], 
                 threePrUtrlen=hsm['threePrUtrlen'], minMapBases=hsm['minMapBases'],
                 hairpin=hsm['hairpin']):
    return fivePrUtrlen - x + minMapBases + 2*l + 2*threePrUtrlen + hairpin

def rev_catchUtr(x,l, threePrUtrlen = hsm['threePrUtrlen'], minMapBases = hsm['minMapBases']):
    return threePrUtrlen - x + l + minMapBases
    
def rev_missUtr(x,l, fivePrUtrlen=hsm['fivePrUtrlen'], 
                 threePrUtrlen=hsm['threePrUtrlen'], minMapBases=hsm['minMapBases'],
                 hairpin=hsm['hairpin']):
    return threePrUtrlen - x + minMapBases + 2*l + 2*fivePrUtrlen + hairpin



def startBaseModel(RLgrid, transcriptLen, fivePrUtrlen=hsm['fivePrUtrlen'],
                   threePrUtrlen=hsm['threePrUtrlen'], hairpin=hsm['hairpin'],
                   lagTime=hsm['nonHotstartLagTime'],
                   polymerase_speed=hsm['polymerase speed']):

    #
    # these define the distribution of unrolled start position for alignments
    #
    totalAdvanceMean = lagTime * polymerase_speed
    totalAdvanceStd = np.sqrt(totalAdvanceMean)
    halfCycleLen = fivePrUtrlen + threePrUtrlen + hairpin + transcriptLen
    
    #
    # allow the non-hotstart 'walk-in' to go up to 3 full laps.  This
    # needs to be generated in a more sophisticated manner to handle
    # different average speeds for strand displacement vs. non-strand
    # displacement.  The easiest way would be an MCMC approach.
    #    
    startBase = 0.5 * dist.norm.pdf(RLgrid, loc=totalAdvanceMean, scale=totalAdvanceStd) +  \
    0.5 * dist.norm.pdf(RLgrid, loc=(totalAdvanceMean + halfCycleLen), scale=totalAdvanceStd) 

    startBase = startBase + 0.5 * dist.norm.pdf(RLgrid+2*halfCycleLen, loc=totalAdvanceMean, scale=totalAdvanceStd) + \
    0.5 * dist.norm.pdf(RLgrid+2*halfCycleLen, loc=(totalAdvanceMean + halfCycleLen), scale=totalAdvanceStd)        

    startBase = startBase + 0.5 * dist.norm.pdf(RLgrid+4*halfCycleLen, loc=totalAdvanceMean, scale=totalAdvanceStd) + \
    0.5 * dist.norm.pdf(RLgrid+4*halfCycleLen, loc=(totalAdvanceMean + halfCycleLen), scale=totalAdvanceStd)        

    startBase = startBase / np.sum(startBase)
    
    return startBase

def compHalfCycleLen(transcriptLen, fivePrUtrlen=hsm['fivePrUtrlen'], 
                     threePrUtrlen=hsm['threePrUtrlen'], 
                     hairpin=hsm['hairpin']):
    halfCycleLen = fivePrUtrlen + threePrUtrlen + hairpin + transcriptLen
    return halfCycleLen


def RnaCoverageModel(RLgrid, transcriptLen, fivePrUtrlen=hsm['fivePrUtrlen'], 
                     threePrUtrlen = hsm['threePrUtrlen'], 
                     minMapBases = hsm['minMapBases'],
                     hairpin = hsm['hairpin']):
    # compute expected value of bases to go for each startLoc
    # which follows a separate logic for forward and reverse (average
    # these).  This follows four cases: two for each primer side,
    # by another two for whether a sufficient number of bases to
    # recognize the nearest UTR have been sequenced.
    #
    basesToGo = np.zeros(len(RLgrid))
    halfCycleLen = compHalfCycleLen(transcriptLen, fivePrUtrlen=fivePrUtrlen,
                                    threePrUtrlen=threePrUtrlen, hairpin=hairpin)
    fcycleSize = np.remainder(RLgrid, halfCycleLen*2)
    hcycleSize = np.remainder(RLgrid, halfCycleLen)

    #
    # case 1: you catch the 5' UTR on the forward strand
    #
    indCatchNorm = mlab.find(fcycleSize < (fivePrUtrlen - minMapBases))
    basesToGo[indCatchNorm] = forw_catchUtr(hcycleSize[indCatchNorm], transcriptLen)

    #
    # case 2: you miss the 5' UTR on the forward strand
    #
    indCatchNorm = mlab.find((fcycleSize >= fivePrUtrlen - minMapBases) & (fcycleSize < halfCycleLen))
    basesToGo[indCatchNorm] = basesToGo[indCatchNorm] + forw_missUtr(hcycleSize[indCatchNorm], transcriptLen)

    #
    # case 3: you catch the 3' UTR on the complementary strand
    #
    indCatchNorm = mlab.find ((fcycleSize >= halfCycleLen) &
        (fcycleSize < (halfCycleLen + threePrUtrlen - minMapBases)))
    basesToGo[indCatchNorm] = basesToGo[indCatchNorm] + rev_catchUtr(hcycleSize[indCatchNorm], transcriptLen)

    #
    # case 4: you miss the 3' UTR on the complementary strand
    #
    indCatchNorm = mlab.find(fcycleSize >= (halfCycleLen + threePrUtrlen - minMapBases))
    basesToGo[indCatchNorm] = basesToGo[indCatchNorm] + rev_missUtr(hcycleSize[indCatchNorm], transcriptLen)

    return basesToGo
