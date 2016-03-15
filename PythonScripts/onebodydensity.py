import numpy as np
import matplotlib.pyplot as plt
import sys

def extract(filename):

    infile = open(filename, 'r')
			
    # initialize list
    distances = []

    # extract energies from file
    for line in infile:
        words = line.split()
        x = float(words[0])
        y = float(words[1])
        z = float(words[2])
        distances.append(np.sqrt(x**2 + y**2 + z**2))

    infile.close()

    return np.array(distances)


def density():

    distances = extract('positionsInteraction6.dat')
    
    plt.hist(distances, bins=100)
    plt.show()

    
    

# ----- main -----
density()


