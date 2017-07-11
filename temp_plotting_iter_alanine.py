# fix autoencoder, compare different datasets
import matplotlib
matplotlib.use('agg')

from ANN_simulation import *
import sys

iter_index = int(sys.argv[1])
include_new_biased_data = int(sys.argv[2])

pkl_file = '../resources/Alanine_pentapeptide/network_%d.pkl' % iter_index
a = Sutils.load_object_from_pkl_file(pkl_file)
_1 = coordinates_data_files_list(['../target/Alanine_pentapeptide/unbiased/']
                         + ['../target/Alanine_pentapeptide/network_%d' % item for item in range(1, iter_index + include_new_biased_data)])
_1 = _1.create_sub_coor_data_files_list_using_filter_conditional(lambda x: not 'aligned' in x)

pdb_file_list = _1.get_list_of_corresponding_pdb_files()
scaling_factor = 10.0
coor_data = _1.get_coor_data(scaling_factor)
coor_data = Sutils.remove_translation(coor_data)
PCs = np.array(a.get_PCs(coor_data))
dihedrals = Alanine_pentapeptide.get_many_dihedrals_from_cossin(
    Alanine_pentapeptide.get_many_cossin_from_coordinates(coor_data))

phi = [item[5] for item in dihedrals]
psi = [item[4] for item in dihedrals]
dis_between_two_ends = np.sqrt(np.sum((coor_data[:,:3] - coor_data[:, -3:]) ** 2, axis=1))
temp_RMSD = Alanine_pentapeptide.metric_RMSD_of_atoms(pdb_file_list, '../resources/Alanine_pentapeptide.pdb', 'name CA')
colors = [phi, temp_RMSD, dis_between_two_ends]
name_of_colors = ['phi', 'temp_RMSD', 'dis_between_two_ends']

fig, axes = plt.subplots(2, 3)
for item in range(3):
    ax = axes[0][item]
    im = ax.scatter(PCs.T[0], PCs.T[1], c=colors[item], s=4, cmap='gist_rainbow')
    fig.colorbar(im, ax=ax)
    ax.set_title(name_of_colors[item])
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    try:
        b = plotting(a)
        b.plotting_potential_centers(fig, ax, coordinates_data_files_list(
                                         ['../target/Alanine_pentapeptide/network_%d' % (iter_index)])._list_of_coor_data_files)
    except Exception as e:
        print e
        pass

    ax = axes[1][item]
    im = ax.scatter(phi, psi, c=colors[item], s=4, cmap='gist_rainbow')
    fig.colorbar(im, ax=ax)
    ax.set_title(name_of_colors[item])
    ax.set_xlabel("phi"); ax.set_ylabel("psi")
    ax.set_xlim([-4,4])
    ax.set_ylim([-4,4])

fig.set_size_inches((15, 10))

if include_new_biased_data:
    include_new_biased_data_string = "_new_data_included"
else:
    include_new_biased_data_string = ""

# fig.savefig('iter_%02d.png' % iter_index)
fig.savefig('iter_%02d%s.png' % (iter_index, include_new_biased_data_string), bbox_inches = 'tight')
