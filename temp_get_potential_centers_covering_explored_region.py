import matplotlib
matplotlib.use("agg")

from ANN_simulation import *

import sys

autoencoder_pkl_file = sys.argv[1]
autoencoder_coef_file = sys.argv[2]
folder_containing_output = sys.argv[3]
force_constant = sys.argv[4]
temperature = sys.argv[5]
interval = sys.argv[6]
total_steps = sys.argv[7]
solvent_type = sys.argv[8]
ensemble_type = sys.argv[9]

a = pickle.load(open(autoencoder_pkl_file, 'rb'))
dimensionality = len(a.get_PCs()[0])
# print "dimensionality = %d" % dimensionality
list_of_points = []
method = 1

if method == 1:
    num_of_points_each_dim = int(sys.argv[10])
    if dimensionality == 2:
        for xx in np.linspace(-1, 1, num_of_points_each_dim):
            for yy in np.linspace(-1, 1, num_of_points_each_dim):
                list_of_points += [[xx,yy]]
    elif dimensionality == 3:
        for xx in np.linspace(-1, 1, num_of_points_each_dim):
            for yy in np.linspace(-1, 1, num_of_points_each_dim):
                for zz in np.linspace(-1, 1, num_of_points_each_dim):
                    list_of_points += [[xx, yy, zz]]
    else:
        raise Exception('error')
    res = a.get_proper_potential_centers_for_WHAM(list_of_points, 0.1, 3)
elif method == 2:
    total_num_of_potential_centers = int(sys.argv[10])  # for method 2
    res = a.get_proper_potential_centers_for_WHAM_2(total_num_of_potential_centers)   # for method 2
else:
    raise Exception('error')

res = np.array(res)
try:
    fig, ax = plt.subplots()
    temp_PCs = a.get_PCs()
    ax.scatter(temp_PCs.T[0], temp_PCs.T[1], s=4)
    ax.scatter(res.T[0], res.T[1], s=20, marker="X", c='red')
    fig.savefig('temp_potential_center.png')
except:
    pass


for device_index, item in enumerate(res):
    if isinstance(molecule_type, Trp_cage):
        pc_string = 'pc_' + ','.join([str(round(_1, 2)) for _1 in item])
        print "python ../src/biased_simulation_Trp_cage.py %s %s %s %s %s %s %s %s --platform CUDA --temperature %s --equilibration_steps 10000 --minimize_energy 1 --data_type_in_input_layer 1 --device %d --fast_equilibration 1 --starting_checkpoint auto" % \
                                (interval, total_steps, force_constant, folder_containing_output, autoencoder_coef_file,
                                    pc_string, solvent_type, ensemble_type, 
                                    temperature, device_index % 2)
    elif isinstance(molecule_type, Alanine_dipeptide):
        pc_string = 'pc_' + ','.join([str(round(_1, 2)) for _1 in item])
        print "python ../src/biased_simulation.py %s %s %s %s %s %s  --temperature %s --data_type_in_input_layer 1 --platform CPU" % \
              (interval, total_steps, force_constant, folder_containing_output, autoencoder_coef_file,
               pc_string, temperature)
    elif isinstance(molecule_type, Alanine_pentapeptide):
        pc_string = 'pc_' + ','.join([str(round(_1, 2)) for _1 in item])
        print "python ../src/biased_simulation_general.py Alanine_pentapeptide %s %s %s %s %s %s %s %s --platform CUDA --temperature %s --equilibration_steps 10000 --minimize_energy 1 --data_type_in_input_layer 1 --fast_equilibration 1 --starting_checkpoint none" % \
              (interval, total_steps, force_constant, folder_containing_output, autoencoder_coef_file,
               pc_string, solvent_type, ensemble_type, temperature)
