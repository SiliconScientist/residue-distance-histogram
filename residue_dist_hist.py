import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import squareform, pdist

def remove_ipynb_checkpoints(path):
    flag = False
    
    for i in range(len(os.listdir(path))):
        if (os.listdir(path)[i] == '.ipynb_checkpoints'):
            flag = True
    if flag:
        os.rmdir(f'{path}/.ipynb_checkpoints')

def pdb_to_coords(path, index):
    # Generate data from from a chosen dirrectory
    ppdb = PandasPdb().read_pdb(f'{path}/{os.listdir(path)[index]}')


    protein_df = ppdb.df['ATOM']
    df_new = protein_df[(protein_df['atom_name'] == 'CA') & (protein_df['residue_name'] == amino_acid) & (protein_df['chain_id'] == chain_id)]

    # Pull the 3D coordinates out of the dataframe
    # 11:14 are the columns with the x, y, z coordinates.
    monomer_coords = df_new.values[:, 11:14]
    return monomer_coords

# Generate symmetry mate
def generate_dimer_coords(coords, matrix, ajustment):
    
    # List comprehension
    anti_coords = [np.dot(matrix, coord) for coord in coords]

    # # Old method
    # anti_coords = []
    # for coord in coords:
    #     anti_coords.append(np.dot(matrix, coord))

    # Convert to numpy array
    anti_coords = np.array(anti_coords)

    mirror_coords = anti_coords + ajustment

    dimer_coords = np.append(coords, mirror_coords, axis=0)
    
    return dimer_coords

def calc_dist(coords, cutoff):

    # dist is a square 2D numpy array (i.e. a matrix with the same # of rows and columns) 
    # dist has a diagonal of zeros and is symmetric (has duplicate values on either side of the diagonal)
    dist = squareform(pdist(coords, 'euclidean'))

    # dist looks like:
    # [0, 1, 2]
    # [1, 0, 3]
    # [2, 3, 0]

    # Remove the duplicate values by extracting the upper triangle of the dist matrix
    upper_triangle = np.triu(dist)

    # upper_triangle looks like:
    # [0, 1, 2]
    # [0, 0, 3]
    # [0, 0, 0]

    # Remove zeros
    distances = upper_triangle[upper_triangle != 0]

    # distances looks like: [1, 2, 3]
    
    distances = distances[distances < cutoff]

    return distances

def display_hist(distances, histo_bins, hist_color, title):
    plt.hist(distances, bins = histo_bins, density = True, color = hist_color)
    plt.xlabel("Distance in angstroms")
    plt.ylabel("Number of occurrences")
    plt.title(title)
    plt.show()




# Parameters
amino_acid = 'ASN'
taxa = 'psychro'

# You can find the following two parameter on a pdb code at REMARK 350
# Rotational symmetry operator
matrix_1ab4 = np.array([[-1, 0, 0],
                        [0, -1, 0],
                        [0, 0, 1]])

# Translational symmetry operator
ajustment_1ab4 = [119.63, 119.63, 0]

# defaults
chain_id = 'A'
atomic_dist_cutoff = 20
histo_bins = 20
hist_color = 'turquoise'


distances_combined = []

pdb_directory = f'/Users/averyhill/MyDocuments/schoeffler_research_summer_2021/pdbs/{taxa}philes/{taxa}_gyra_pdb_folder'

pdb_count = len(os.listdir(pdb_directory))

# Jupyter notebook will often put checkipoints into folders, this functions removes them
remove_ipynb_checkpoints(pdb_directory)

for i in range(pdb_count):
    # Get monomer coordinates from pdb_to_coords()
    monomer_coords = pdb_to_coords(pdb_directory, i)

    # Dimerize the monomer coordinates using generate_dimer()
    dimer_coords = generate_dimer_coords(monomer_coords, matrix_1ab4, ajustment_1ab4)

    # Convert numpy array to dataframe because calc_dist only accepts dataframes
    dimer_coord_array = pd.DataFrame(dimer_coords, columns = ['x_coord','y_coord','z_coord'])

    # Calculate the distance between the coordinates using calc_dist()
    distances = calc_dist(dimer_coord_array, atomic_dist_cutoff)

    # Sum the distances
    distances_combined = np.append(distances_combined, distances)

# Show the histogram of the summed distances using display_hist
display_hist(distances_combined, histo_bins, hist_color, taxa)