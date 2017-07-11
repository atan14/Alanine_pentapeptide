import sys

from ANN_simulation import *

option = sys.argv[1]
folder = sys.argv[2]
autoencoder_file = sys.argv[3]
folder_to_store_files = sys.argv[4]
num_of_dihedrals = int(sys.argv[5])
dimensionality = int(sys.argv[6])
starting_index_last_few_frames = int(sys.argv[7])
ending_index = int(sys.argv[8])
random_flag = int(sys.argv[9])

a=pickle.load(open(autoencoder_file, 'rb'))

if num_of_dihedrals == 3:
    dihedral_angle_range = [0, 1, 2]
elif num_of_dihedrals == 2:
    dihedral_angle_range = [1, 2]

a.generate_mat_file_for_WHAM_reweighting(folder, mode=option,folder_to_store_files=folder_to_store_files,
    input_data_type='Cartesian', dihedral_angle_range=dihedral_angle_range, dimensionality=dimensionality,
     starting_index_of_last_few_frames=starting_index_last_few_frames,
     ending_index_of_frames=ending_index, random_dataset=random_flag)
