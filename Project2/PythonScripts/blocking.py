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


def extractBinary(path):

    energies = np.fromfile(path)
    print "Done loading file"

    """
    mydt = np.dtype([('x', np.double),('y', np.double)])
    positions = np.fromfile(path, mydt)
    print "Done loading file"
    r = np.zeros(len(positions))
    for i in xrange(len(positions)):
        r[i] = np.sqrt([positions[i][0]*positions[i][0] + positions[i][1]*positions[i][1]])
    """

    return energies;



def blocking(energies, numberOfBlocks, blockSize):

    # find mean of each block for given block size
    blockMeans = np.zeros(numberOfBlocks)
    for j in xrange(numberOfBlocks):
        blockMeans[j] = mean(energies[j*blockSize:(j+1)*blockSize])

    # find mean and variance of block means
    finalMean = sum(blockMeans)/numberOfBlocks
    finalMean2 = sum(blockMeans**2)/numberOfBlocks
    finalVariance = (finalMean2 - finalMean**2)/numberOfBlocks
    return finalVariance



def mean(blockEnergies):
    return sum(blockEnergies)/len(blockEnergies)
   
       

def main(path, binary):

    if binary:
        energies = extractBinary(path)

    else:
        energies = np.array(extract(path))
    
    # block parameters
    numberOfEnergies = len(energies)
    blockSizeStepLength = 50
    minimumBlockSize = 1
    #maximumBlockSize = len(energies) / 250
    maximumBlockSize = 3000
    numberOfBlockSizes = (maximumBlockSize - minimumBlockSize) / blockSizeStepLength
    
    print 'No of energies = ', numberOfEnergies
    print 'Min size = '      , minimumBlockSize
    print 'Max size = '      , maximumBlockSize
    print 'Step length = '   , blockSizeStepLength
    print 'Number of block sizes = ', numberOfBlockSizes
    
    varianceOfGivenBlockSize 	= np.zeros(numberOfBlockSizes)
    allBlockSizes 				= np.zeros(numberOfBlockSizes)

    # loop over block sizes
    for i in xrange(numberOfBlockSizes):
        blockSize = minimumBlockSize + i*blockSizeStepLength
        numberOfBlocks = numberOfEnergies / blockSize
        allBlockSizes[i] = blockSize
        finalVariance = blocking(energies, numberOfBlocks, blockSize)
        varianceOfGivenBlockSize[i] = finalVariance

        # show progress
        sys.stdout.write("\r%2d %% complete" % ((float(i)/numberOfBlockSizes)*100))
        sys.stdout.flush()

    print
    return varianceOfGivenBlockSize, allBlockSizes



# ----- main -----
path = "../dataFiles/energiesN2w100Ms6.bin"
binary = True
varianceOfGivenBlockSize, allBlockSizes = main(path, binary)

plt.plot(allBlockSizes, varianceOfGivenBlockSize)
plt.xlabel(r'$n_b$', fontsize=30)
plt.ylabel(r'$\sigma$', fontsize=30)
plt.title('N = 100')
plt.show()






