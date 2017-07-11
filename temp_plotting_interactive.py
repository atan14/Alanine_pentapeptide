import sys

from ANN_simulation import *

autoencoder_pkl = sys.argv[1]

molecule_type = Alanine_pentapeptide()

step_interval = 1

fig, ax = plt.subplots()
_1 = coordinates_data_files_list(['../target/Alanine_pentapeptide/'])
_1 = _1.create_sub_coor_data_files_list_using_filter_conditional(lambda x: not 'aligned' in x)

my_file_list = _1.get_list_of_coor_data_files()
pdb_file_list = _1.get_list_of_corresponding_pdb_files()
input_data = _1.get_coor_data(CONFIG_49)
input_data = Sutils.remove_translation(input_data)

coloring = molecule_type.metric_RMSD_of_atoms(pdb_file_list, step_interval=step_interval)

a = Sutils.load_object_from_pkl_file(autoencoder_pkl)

b = plotting(a)

temp_fig, _, _ = b.plotting_with_coloring_option("PC", fig, ax,
                                                 input_data_for_plotting=input_data,
                                                 color_option='other',
                                                 contain_colorbar=True,
                                                 other_coloring=coloring,
                                                 smoothing_using_RNR=False,
                                                 enable_mousing_clicking_event=True,
                                                 related_coor_list_obj=_1,
                                                 )
# b.plotting_potential_centers(fig, ax,
#                              coordinates_data_files_list(
#                                  ['../target/Trp_cage/network_%d' % (i)])._list_of_coor_data_files)

plt.show()


