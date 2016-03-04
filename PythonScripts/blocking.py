import numpy as np
import matplotlib.pyplot as plt
import sys

def extract(filename):

    infile = open(filename, 'r')
			
    # initialize list
    energies = []

    # skip comment line
    infile.readline()

    # extract energies from file
    for line in infile:
        energies.append(float(line))

    infile.close()

    return energies


def blocking(energies, numberOfBlocks, blockSize):

    # find mean of each block for given block size
    blockMeans = np.zeros(numberOfBlocks)
    for j in xrange(numberOfBlocks):
        blockMeans[j] = mean(energies[j*blockSize:(j+1)*blockSize])

    # find mean and variance of block means
    finalMean = sum(blockMeans)/numberOfBlocks
    finalMean2 = sum(blockMeans**2)/numberOfBlocks
    finalVariance = finalMean2 - finalMean**2
    #finalVariance = sum((blockMeans-finalMean)**2)/numberOfBlocks
    return finalMean, finalVariance


def mean(blockEnergies):
    return sum(blockEnergies)/len(blockEnergies)
    
        



def main():

    # extract data and convert to array
    energies = np.array(extract('energy.dat'))

    # block parameters
    numberOfEnergies = len(energies)
    blockSizeStepLength = 50 #(maximumBlockSize - minimumBlockSize) / numberOfBlockSizes
    minimumBlockSize = 1
    maximumBlockSize = 10000
    numberOfBlockSizes = (maximumBlockSize - minimumBlockSize) / blockSizeStepLength
	

    print 'No of energies = ', numberOfEnergies
    print 'Min size = '      , minimumBlockSize
    print 'Max size = '      , maximumBlockSize
    print 'Step length = '   , blockSizeStepLength
    print 'Number of block sizes = ', numberOfBlockSizes

    # loop over block sizes
    meanOfGivenBlockSize 		= np.zeros(numberOfBlockSizes)
    varianceOfGivenBlockSize 	= np.zeros(numberOfBlockSizes)
    allBlockSizes 				= np.zeros(numberOfBlockSizes)
    allNumberOfBlocks 			= np.zeros(numberOfBlockSizes)
    for i in xrange(numberOfBlockSizes):
        blockSize = minimumBlockSize + i*blockSizeStepLength
        numberOfBlocks = numberOfEnergies / blockSize
        allBlockSizes[i] = blockSize
        allNumberOfBlocks[i] = numberOfBlocks
        finalMean, finalVariance = blocking(energies, numberOfBlocks, blockSize)
        meanOfGivenBlockSize[i] = finalMean
        varianceOfGivenBlockSize[i] = finalVariance

        # show progress
        sys.stdout.write("\r{0}".format((float(i)/numberOfBlockSizes)*100))
        sys.stdout.flush()

    return meanOfGivenBlockSize, varianceOfGivenBlockSize, allBlockSizes, allNumberOfBlocks



# main
meanOfGivenBlockSize, varianceOfGivenBlockSize, allBlockSizes, allNumberOfBlocks = main()
plt.plot(allBlockSizes, varianceOfGivenBlockSize)
plt.show()
print






