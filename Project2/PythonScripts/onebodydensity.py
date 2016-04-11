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
    mean_distance = sum(distances)/len(distances)
    
    figure = plt.hist(distances, bins=200, color='b', normed=1)
    plt.axvline(mean_distance, color='r', linestyle='dashed', label='Average distance')
    plt.legend(fontsize=20)
    plt.xlabel(r'$r/a_{ho}$', fontsize=20)
    plt.ylabel('Normalized number of particles', fontsize=20)
    plt.title('System 2', fontsize=20)
    plt.show()

    
    

# ----- main -----
density()


