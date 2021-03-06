from config import *
from molecule_spec_sutils import *  # import molecule specific unitity code
from coordinates_data_files_list import *
from sklearn.cluster import KMeans
from keras.models import Sequential
from keras.optimizers import *
from keras.layers import Dense, Activation, Lambda, Reshape
from keras.regularizers import l2
from keras.callbacks import EarlyStopping
import random

##################    set types of molecules  ############################

molecule_type = Sutils.create_subclass_instance_using_name(CONFIG_30)

##########################################################################


class autoencoder(object):
    """the neural network for simulation
    this class includes abstract methods, that must be implemented by subclasses
    """
    def __init__(self,
                 index,  # the index of the current network
                 data_set_for_training,
                 output_data_set = None,  # output data may not be the same with the input data
                 autoencoder_info_file=None,  # this might be expressions, or coefficients
                 training_data_interval=CONFIG_2,
                 in_layer_type=LinearLayer,
                 hidden_layers_types=CONFIG_17,
                 out_layer_type=LinearLayer,  # different layers
                 node_num=CONFIG_3,  # the structure of ANN
                 max_num_of_training=CONFIG_5,
                 filename_to_save_network=CONFIG_6,
                 hierarchical=CONFIG_44,
                 network_verbose=CONFIG_46,
                 output_as_circular=CONFIG_47,
                 *args, **kwargs           # for extra init functions for subclasses
                 ):

        self._index = index
        self._data_set = data_set_for_training
        self._output_data_set = output_data_set
        self._training_data_interval = training_data_interval
        if autoencoder_info_file is None:
            self._autoencoder_info_file = "../resources/%s/autoencoder_info_%d.txt" % (CONFIG_30, index)
        else:
            self._autoencoder_info_file = autoencoder_info_file

        if not in_layer_type is None: self._in_layer_type = in_layer_type
        if not hidden_layers_types is None: self._hidden_layers_type = hidden_layers_types
        if not out_layer_type is None: self._out_layer_type = out_layer_type

        self._node_num = node_num
        self._max_num_of_training = max_num_of_training
        if filename_to_save_network is None:
            self._filename_to_save_network = "../resources/%s/network_%s.pkl" % (
            CONFIG_30, str(self._index))  # by default naming with its index
        else:
            self._filename_to_save_network = filename_to_save_network

        self._hierarchical = hierarchical
        self._network_verbose = network_verbose
        num_of_PC_nodes_for_each_PC = 2 if self._hidden_layers_type[1] == CircularLayer else 1
        self._num_of_PCs = self._node_num[2] / num_of_PC_nodes_for_each_PC
        self._connection_between_layers_coeffs = None
        self._connection_with_bias_layers_coeffs = None
        self._output_as_circular = output_as_circular
        self._init_extra(*args, **kwargs)
        return

    @abc.abstractmethod
    def _init_extra(self):
        """must be implemented by subclasses"""
        pass

    def save_into_file(self, filename=CONFIG_6, fraction_of_data_to_be_saved = 1.0):
        if filename is None:
            filename = self._filename_to_save_network

        if fraction_of_data_to_be_saved != 1.0:
            number_of_data_points_to_be_saved = int(self._data_set.shape[0] * fraction_of_data_to_be_saved)
            print ("Warning: only %f of data (%d out of %d) are saved into pkl file" % (fraction_of_data_to_be_saved,
                                                                                        number_of_data_points_to_be_saved,
                                                                                        self._data_set.shape[0]))
            self._data_set = self._data_set[:number_of_data_points_to_be_saved]
            if not self._output_data_set is None:        # for backward compatibility
                self._output_data_set = self._output_data_set[:number_of_data_points_to_be_saved]

        if os.path.isfile(filename):  # backup file if previous one exists
            os.rename(filename, filename.split('.pkl')[0] + "_bak_" + datetime.datetime.now().strftime(
                "%Y_%m_%d_%H_%M_%S") + '.pkl')

        with open(filename, 'wb') as my_file:
            pickle.dump(self, my_file, pickle.HIGHEST_PROTOCOL)
        return

    def get_expression_of_network(self):
        """
        this function generates expression of PCs in terms of inputs
        """
        # FIXME: the expression no longer works, since I made input list for input layer of autoencoder consistent
        # for both alanine dipeptide and trp-cage, always [cos, sin, cos, sin ....],
        # which is consistent with ANN_Force, instead of [cos, cos, cos, cos, sin, sin, sin, sin]
        type_of_middle_hidden_layer = self._hidden_layers_type[1]

        node_num = self._node_num
        expression = ""

        # 1st part: network
        for i in range(2):
            expression = '\n' + expression
            mul_coef = self._connection_between_layers_coeffs[i].reshape(node_num[i + 1], node_num[i])
            bias_coef = self._connection_with_bias_layers_coeffs[i]

            for j in range(np.size(mul_coef, 0)):
                temp_expression = 'in_layer_%d_unit_%d = ' % (i + 1, j)

                for k in range(np.size(mul_coef, 1)):
                    temp_expression += ' %f * out_layer_%d_unit_%d +' % (mul_coef[j, k], i, k)

                temp_expression += ' %f;\n' % (bias_coef[j])
                expression = temp_expression + expression  # order of expressions matter in OpenMM

            if i == 1 and type_of_middle_hidden_layer == CircularLayer:
                for j in range(np.size(mul_coef, 0) / 2):
                    temp_expression = 'out_layer_%d_unit_%d = ( in_layer_%d_unit_%d ) / radius_of_circular_pair_%d;\n' % \
                                      (i + 1, 2 * j, i + 1, 2 * j, j)
                    temp_expression += 'out_layer_%d_unit_%d = ( in_layer_%d_unit_%d ) / radius_of_circular_pair_%d;\n' % \
                                       (i + 1, 2 * j + 1, i + 1, 2 * j + 1, j)
                    temp_expression += 'radius_of_circular_pair_%d = sqrt( in_layer_%d_unit_%d * in_layer_%d_unit_%d + in_layer_%d_unit_%d * in_layer_%d_unit_%d );\n' \
                                       % (j, i + 1, 2 * j, i + 1, 2 * j, i + 1, 2 * j + 1, i + 1, 2 * j + 1)

                    expression = temp_expression + expression
            else:
                for j in range(np.size(mul_coef, 0)):
                    temp_expression = 'out_layer_%d_unit_%d = tanh( in_layer_%d_unit_%d );\n' % (i + 1, j, i + 1, j)
                    expression = temp_expression + expression

        # 2nd part: relate PCs to network
        if type_of_middle_hidden_layer == CircularLayer:
            temp_expression = 'PC0 = acos( out_layer_2_unit_0 ) * ( step( out_layer_2_unit_1 ) - 0.5) * 2;\n'
            temp_expression += 'PC1 = acos( out_layer_2_unit_2 ) * ( step( out_layer_2_unit_3 ) - 0.5) * 2;\n'
            expression = temp_expression + expression
        elif type_of_middle_hidden_layer == TanhLayer:
            temp_expression = 'PC0 = out_layer_2_unit_0;\nPC1 = out_layer_2_unit_1;\n'
            expression = temp_expression + expression

        # 3rd part: definition of inputs
        expression += molecule_type.get_expression_for_input_of_this_molecule()

        return expression

    def write_expression_into_file(self, out_file=None):
        if out_file is None: out_file = self._autoencoder_info_file

        expression = self.get_expression_of_network()
        with open(out_file, 'w') as f_out:
            f_out.write(expression)
        return

    def write_coefficients_of_connections_into_file(self, out_file=None):
        if out_file is None: out_file = self._autoencoder_info_file

        with open(out_file, 'w') as f_out:
            for item in [0, 1]:
                f_out.write(str(list(self._connection_between_layers_coeffs[item])))
                f_out.write(',\n')

            for item in [0, 1]:
                f_out.write(str(list(self._connection_with_bias_layers_coeffs[item])))
                f_out.write(',\n')
        return

    @abc.abstractmethod
    def get_PCs(self, input_data=None):
        """must be implemented by subclasses"""
        pass

    @abc.abstractmethod
    def train(self):
        """must be implemented by subclasses"""
        pass

    @abc.abstractmethod
    def get_output_data(self, input_data=None, num_of_PCs=None):
        """must be implemented by subclasses"""
        pass

    @abc.abstractmethod
    def get_mid_result(self, input_data=None):
        """must be implemented by subclasses"""
        pass

    def get_training_error(self, num_of_PCs=None):
        """
        :param num_of_PCs: this option only works for hierarchical case, indicate you would like to get error with
        a specific number of PCs (instead of all PCs)
        """
        input_data = np.array(self._data_set)
        actual_output_data = self.get_output_data(num_of_PCs)
        if hasattr(self, '_output_data_set') and not self._output_data_set is None:
            expected_output_data = self._output_data_set
        else:
            expected_output_data = input_data

        return np.linalg.norm(expected_output_data - actual_output_data) / sqrt(self._node_num[0] * len(input_data))

    def get_fraction_of_variance_explained(self, num_of_PCs=None):
        """ here num_of_PCs is the same with that in get_training_error() """
        input_data = np.array(self._data_set)
        actual_output_data = self.get_output_data(num_of_PCs)
        if hasattr(self, '_output_data_set') and not self._output_data_set is None:
            expected_output_data = self._output_data_set
        else:
            expected_output_data = input_data

        var_of_input = sum(np.var(expected_output_data, axis=0))
        var_of_err = sum(np.var(actual_output_data - expected_output_data, axis=0))
        return 1 - var_of_err / var_of_input

    def get_commands_for_further_biased_simulations(self, list_of_potential_center=None,
                                                    num_of_simulation_steps=None,
                                                    autoencoder_info_file=None,
                                                    force_constant_for_biased=None,
                                                    bias_method=CONFIG_65
                                                    ):
        """this function creates a list of commands for further biased simulations that should be done later,
        either in local machines or on the cluster
        """
        if num_of_simulation_steps is None:
            num_of_simulation_steps = CONFIG_8
        if autoencoder_info_file is None:
            autoencoder_info_file = self._autoencoder_info_file
        if bias_method == "US":
            PCs_of_network = self.get_PCs()
            if self._hidden_layers_type[1] == CircularLayer:
                assert (len(PCs_of_network[0]) == self._node_num[2] / 2)
            else:
                assert (len(PCs_of_network[0]) == self._node_num[2])
            if list_of_potential_center is None:
                list_of_potential_center = molecule_type.get_boundary_points(list_of_points=PCs_of_network)
            if force_constant_for_biased is None:
                if isinstance(molecule_type, Trp_cage):
                    folder_state_coor_file = '../resources/1l2y_coordinates.txt'
                elif isinstance(molecule_type, Alanine_dipeptide):
                    folder_state_coor_file = '../resources/alanine_dipeptide_coordinates.txt'
                elif isinstance(molecule_type, Alanine_pentapeptide):
                    folder_state_coor_file = '../resources/Alanine_pentapeptide_coordinates.txt'
                else:
                    raise Exception('molecule type error')

                if CONFIG_53 == "fixed":
                    force_constant_for_biased = [CONFIG_9 for _ in list_of_potential_center]
                elif CONFIG_53 == "flexible":
                    input_folded_state = np.loadtxt(folder_state_coor_file) / CONFIG_49
                    PC_folded_state = self.get_PCs(Sutils.remove_translation(input_folded_state))[0]
                    print("PC_folded_state = %s" % str(PC_folded_state))
                    force_constant_for_biased = [2 * CONFIG_54 / np.linalg.norm(np.array(item) - PC_folded_state) ** 2
                                                 for item in list_of_potential_center]
                elif CONFIG_53 == "truncated":
                    input_folded_state = np.loadtxt(folder_state_coor_file) / CONFIG_49
                    PC_folded_state = self.get_PCs(Sutils.remove_translation(input_folded_state))[0]
                    print("PC_folded_state = %s" % str(PC_folded_state))
                    force_constant_for_biased = [min(2 * CONFIG_54 / np.linalg.norm(np.array(item) - PC_folded_state) ** 2,
                                                     CONFIG_9) for item in list_of_potential_center]
                else:
                    raise Exception("error")

            todo_list_of_commands_for_simulations = []
            if CONFIG_48 == 'Cartesian':
                input_data_type = 1
            elif CONFIG_48 == 'cossin':
                input_data_type = 0
            else:
                raise Exception("error input data type")

            for index, potential_center in enumerate(list_of_potential_center):
                if isinstance(molecule_type, Alanine_dipeptide):
                    parameter_list = (str(CONFIG_16), str(num_of_simulation_steps), str(force_constant_for_biased[index]),
                                      '../target/Alanine_dipeptide/network_%d' % self._index,
                                      autoencoder_info_file,
                                      'pc_' + str(potential_center).replace(' ', '')[1:-1],
                                      input_data_type
                                      # need to remove white space, otherwise parsing error
                                      )
                    command = "python ../src/biased_simulation.py %s %s %s %s %s %s --data_type_in_input_layer %d" % parameter_list
                    if CONFIG_42:  # whether the force constant adjustable mode is enabled
                        command = command + ' --fc_adjustable --autoencoder_file %s --remove_previous ' % (
                            '../resources/Alanine_dipeptide/network_%d.pkl' % self._index)

                elif isinstance(molecule_type, Trp_cage):
                    fast_equilibration_flag = CONFIG_72
                    parameter_list = (str(CONFIG_16), str(num_of_simulation_steps), str(force_constant_for_biased[index]),
                                      '../target/Trp_cage/network_%d/' % self._index,
                                      autoencoder_info_file,
                                      'pc_' + str(potential_center).replace(' ', '')[1:-1],
                                      CONFIG_40, CONFIG_51, input_data_type, index % 2, fast_equilibration_flag)
                    command = "python ../src/biased_simulation_general.py Trp_cage %s %s %s %s %s %s %s %s --data_type_in_input_layer %d --device %d --fast_equilibration %d" % parameter_list
                    if CONFIG_42:
                        command = command + ' --fc_adjustable --autoencoder_file %s --remove_previous' % (
                            '../resources/Trp_cage/network_%d.pkl' % self._index)

                elif isinstance(molecule_type, Alanine_pentapeptide):
                    fast_equilibration_flag = CONFIG_72
                    parameter_list = (
                    str(CONFIG_16), str(num_of_simulation_steps), str(force_constant_for_biased[index]),
                    '../target/Alanine_pentapeptide/network_%d/' % self._index,
                    autoencoder_info_file,
                    'pc_' + str(potential_center).replace(' ', '')[1:-1],
                    CONFIG_40, CONFIG_51, input_data_type, index % 2, fast_equilibration_flag)
                    command = "python ../src/biased_simulation_general.py Alanine_pentapeptide %s %s %s %s %s %s %s %s --data_type_in_input_layer %d --device %d --fast_equilibration %d" % parameter_list
                    if CONFIG_42:
                        command = command + ' --fc_adjustable --autoencoder_file %s --remove_previous' % (
                            '../resources/Alanine_pentapeptide/network_%d.pkl' % self._index)

                else:
                    raise Exception("molecule type not defined")

                todo_list_of_commands_for_simulations += [command]
        elif bias_method == "MTD":
            todo_list_of_commands_for_simulations = []
            self.write_expression_script_for_plumed()
            dimensionality = CONFIG_36
            pc_string = 'pc_' + ','.join(['0' for _ in range(dimensionality)])
            if isinstance(molecule_type, Alanine_dipeptide):
                for mtd_sim_index in range(5):
                    parameter_list = (str(CONFIG_16), str(num_of_simulation_steps), str(mtd_sim_index),
                                      '../target/Alanine_dipeptide/network_%d/' % self._index,
                                      self._autoencoder_info_file, pc_string)
                    command = "python ../src/biased_simulation.py %s %s %s %s %s %s --data_type_in_input_layer 1 --bias_method MTD" % parameter_list
                    todo_list_of_commands_for_simulations += [command]
            elif isinstance(molecule_type, Trp_cage):
                for mtd_sim_index in range(6):
                    parameter_list = (str(CONFIG_16), str(num_of_simulation_steps), str(mtd_sim_index),
                                      '../target/Trp_cage/network_%d/' % self._index, self._autoencoder_info_file,
                                      pc_string, CONFIG_40, CONFIG_51, mtd_sim_index % 2)
                    command = "python ../src/biased_simulation_general.py Trp_cage %s %s %s %s %s %s %s %s --data_type_in_input_layer 1 --bias_method MTD --device %d" % parameter_list
                    todo_list_of_commands_for_simulations += [command]
            elif isinstance(molecule_type, Alanine_pentapeptide):
                for mtd_sim_index in range(6):
                    parameter_list = (str(CONFIG_16), str(num_of_simulation_steps), str(mtd_sim_index),
                                      '../target/Alanine_pentapeptide/network_%d/' % self._index, self._autoencoder_info_file,
                                      pc_string, CONFIG_40, CONFIG_51, mtd_sim_index % 2)
                    command = "python ../src/biased_simulation_general.py Alanine_pentapeptide %s %s %s %s %s %s %s %s --data_type_in_input_layer 1 --bias_method MTD --device %d" % parameter_list
                    todo_list_of_commands_for_simulations += [command]
            else:
                raise Exception("molecule type not defined")
        else:
            raise Exception("bias method not found")

        return todo_list_of_commands_for_simulations

    def get_proper_potential_centers_for_WHAM(self, list_of_points, threshold_radius, min_num_of_neighbors):
        """
        This function selects some 'proper' potential centers within the domain from list_of_points, by "proper"
        we mean there are at least min_num_of_neighbors data points that are located within the radius of threshold_radius
        of the specific potential center.
        Typically list_of_points could be evenly distributed grid points in PC space
        """
        data_points = np.array(self.get_PCs())
        list_of_points = np.array(list_of_points)
        assert (data_points.shape[1] == list_of_points.shape[1])
        proper_potential_centers = []

        for item in list_of_points:
            neighbors_num = sum([np.dot(item - x, item - x) < threshold_radius * threshold_radius for x in data_points])

            if neighbors_num >= min_num_of_neighbors:
                proper_potential_centers += [item]

        return proper_potential_centers

    def get_proper_potential_centers_for_WHAM_2(self, total_number_of_potential_centers):
        data_points = np.array(self.get_PCs())
        kmeans = KMeans(init='k-means++', n_clusters=total_number_of_potential_centers, n_init=10)
        kmeans.fit(data_points)
        return kmeans.cluster_centers_

    def generate_mat_file_for_WHAM_reweighting(self,
                                               directory_containing_coor_files,
                                               mode="Bayes",  # mode = "standard" or "Bayes"
                                               folder_to_store_files='./standard_WHAM/', dimensionality=2,
                                               input_data_type='cossin',  # input_data_type could be 'cossin' or 'Cartesian'
                                               scaling_factor=CONFIG_49,  # only works for 'Cartesian'
                                               dihedral_angle_range=[1,2],  # only used for alanine dipeptide
                                               starting_index_of_last_few_frames=0,  # number of last few frames used in calculation, 0 means to use all frames
                                               ending_index_of_frames = 0,  # end index, for FES convergence check
                                               random_dataset = False,  # pick random dataset to estimate variance
                                               num_of_bins = 20
                                               ):
        """
        note: 
        dihedral_angle_range, starting_index_of_last_few_frames, ending_index_of_frames, random_dataset 
        may not work for Bayes mode
        num_of_bins only works for Bayes mode
        """
        if folder_to_store_files[-1] != '/':
            folder_to_store_files += '/'
        if not os.path.exists(folder_to_store_files):
            subprocess.check_output(['mkdir', folder_to_store_files])

        if mode == "Bayes":
            for item in ['bias', 'hist', 'traj', 'traj_proj']:
                directory = folder_to_store_files + item
                subprocess.check_output(['mkdir', '-p', directory])
                assert (os.path.exists(directory))
        else: pass

        list_of_coor_data_files = coordinates_data_files_list(
            [directory_containing_coor_files])._list_of_coor_data_files
        force_constants = []
        harmonic_centers = []
        window_counts = []
        coords = []
        umbOP = []
        for item in list_of_coor_data_files:
            # print('processing %s' %item)
            temp_force_constant = float(item.split('output_fc_')[1].split('_pc_')[0])
            force_constants += [[temp_force_constant] * dimensionality  ]
            temp_harmonic_center_string = item.split('_pc_[')[1].split(']')[0]
            harmonic_centers += [[float(item_1) for item_1 in temp_harmonic_center_string.split(',')]]
            if input_data_type == 'cossin':
                temp_coor = self.get_PCs(molecule_type.get_many_cossin_from_coordinates_in_list_of_files([item]))
            elif input_data_type == 'Cartesian':
                temp_coor = self.get_PCs(Sutils.remove_translation(np.loadtxt(item) / scaling_factor))
            else:
                raise Exception('error input_data_type')

            if random_dataset:
                # data_index_list = random.sample(range(temp_coor.shape[0]), int(0.5 * temp_coor.shape[0]))  # nonrepeated
                data_index_list = [random.choice(range(temp_coor.shape[0])) for _ in range(temp_coor.shape[0])]  # allow repeated data
                # print "random data_index_list"
            else:
                data_index_list = np.arange(temp_coor.shape[0])
                data_index_list = data_index_list[starting_index_of_last_few_frames:]
                if ending_index_of_frames != 0: data_index_list = data_index_list[:ending_index_of_frames]

            temp_coor = temp_coor[data_index_list]
            assert len(temp_coor) == len(data_index_list)
            temp_window_count = temp_coor.shape[0]
            window_counts += [float(temp_window_count)]   # there exists problems if using int

            coords += list(temp_coor)
            if isinstance(molecule_type, Alanine_dipeptide):
                temp_angles = np.array(molecule_type.get_many_dihedrals_from_coordinates_in_file([item]))[data_index_list]
                temp_umbOP = [[a[temp_dihedral_index] for temp_dihedral_index in dihedral_angle_range] for a in temp_angles]
                assert (temp_window_count == len(temp_umbOP)), (temp_window_count, len(temp_umbOP))
                assert (len(dihedral_angle_range) == len(temp_umbOP[0]))
                umbOP += temp_umbOP

        if mode == "standard":
            max_of_coor = map(lambda x: round(x, 1) + 0.1, map(max, list(zip(*coords))))
            min_of_coor = map(lambda x: round(x, 1) - 0.1, map(min, list(zip(*coords))))
            interval = 0.1

            window_counts = np.array(window_counts)
            sciio.savemat(folder_to_store_files + 'WHAM_nD__preprocessor.mat', {'window_counts': window_counts,
                                                                                'force_constants': force_constants,
                                                                                'harmonic_centers': harmonic_centers,
                                                                                'coords': coords, 'dim': dimensionality,
                                                                                'temperature': 300.0,
                                                                                'periodicity': [[0.0] * dimensionality],
                                                                                'dF_tol': 0.001,
                                                                                'min_gap_max_ORIG': [
                                                                                    [min_of_coor[item_2], interval,
                                                                                     max_of_coor[item_2]] for item_2 in range(dimensionality)]
                                                                                })
            sciio.savemat(folder_to_store_files + 'umbrella_OP.mat',
                          {'umbOP': umbOP
                           })

        elif mode == "Bayes":
            # write info into files
            # 1st: bias potential info
            with open(folder_to_store_files + 'bias/harmonic_biases.txt', 'w') as f_out:
                for item in range(len(force_constants)):
                    f_out.write('%d\t' % (item + 1))
                    for write_item in harmonic_centers[item]:
                        f_out.write('%f\t' % write_item)
                    for write_item in force_constants[item]:
                        f_out.write('%f\t' % write_item)
                    f_out.write("\n")

            # 2nd: trajectory, and projection trajectory in phi-psi space (for reweighting), and histogram
            num_of_bins = num_of_bins
            binEdges = np.linspace(-1, 1, num_of_bins)
            with open(folder_to_store_files + 'hist/hist_binEdges.txt', 'w') as f_out:
                for _ in range(dimensionality):
                    for item in binEdges:
                        f_out.write('%f\t' % item)
                    f_out.write('\n')

            num_of_bins_proj = 40
            binEdges_proj = np.array([np.linspace(-np.pi, np.pi, num_of_bins_proj),
                                      np.linspace(-np.pi, np.pi, num_of_bins_proj)])
            with open(folder_to_store_files + 'hist/hist_binEdges_proj.txt', 'w') as f_out:
                for row in binEdges_proj:
                    for item in row:
                        f_out.write('%f\t' % item)
                    f_out.write('\n')

            end_index = 0
            for item, count in enumerate(window_counts):
                start_index = int(end_index)
                end_index = int(start_index + count)
                with open(folder_to_store_files + 'traj/traj_%d.txt' % (item + 1), 'w') as f_out_1, \
                        open(folder_to_store_files + 'traj_proj/traj_%d.txt' % (item + 1), 'w') as f_out_2, \
                        open(folder_to_store_files + 'hist/hist_%d.txt' % (item + 1), 'w') as f_out_3:
                    for line in coords[start_index:end_index]:
                        for item_1 in line:
                            f_out_1.write('%f\t' % item_1)

                        f_out_1.write("\n")

                    for line in umbOP[start_index:end_index]:
                        for item_1 in line:
                            f_out_2.write('%f\t' % item_1)

                        f_out_2.write("\n")

                    temp_hist, _ = np.histogramdd(np.array(coords[start_index:end_index]),
                                                     bins=[binEdges.tolist()] * dimensionality)
                    for _1 in temp_hist.flatten():
                        f_out_3.write('%d\t' % _1)
        else:
            raise Exception("error mode!")

        return


class neural_network_for_simulation(autoencoder):
    """
    A better class name would be 'autoencoder_pybrain', to indicate that it is constructed using Pybrain, but for
    **backward compatibility** reason I keep this one.
    """
    def _init_extra(self,
                    trainer=None,
                    network_parameters=CONFIG_4,  # includes [learningrate,momentum, weightdecay, lrdecay]
                    ):
        self._network_parameters = network_parameters
        self._trainer = trainer  # save the trainer so that we could train this network step by step later
        self._molecule_net = None
        return

    def get_mid_result(self, input_data=None):
        if input_data is None: input_data = self._data_set
        mid_result = []
        for item in input_data:
            self._molecule_net.activate(item)
            mid_result.append([list(layer.outputbuffer[0]) for layer in
                               self._molecule_net.modulesSorted[5:]])  # exclude bias nodes and input layer
        return mid_result

    def get_PCs(self, input_data=None):
        """
        write an independent function for getting PCs, since it is different for TanhLayer, and CircularLayer
        """
        if input_data is None: input_data = self._data_set
        num_of_hidden_layers = len(self._hidden_layers_type)
        index_of_bottleneck_hidden_layer = (
                                           num_of_hidden_layers - 1) / 2  # it works for both 3-layer and 5-layer structure
        type_of_middle_hidden_layer = self._hidden_layers_type[index_of_bottleneck_hidden_layer]
        temp_mid_result = self.get_mid_result(input_data=input_data)
        mid_result_1 = [item[index_of_bottleneck_hidden_layer] for item in temp_mid_result]

        if type_of_middle_hidden_layer == CircularLayer:
            PCs = [[acos(item[2 * _1]) * np.sign(item[2 * _1 + 1]) for _1 in range(self._node_num[2] / 2)] for item in
                   mid_result_1]
            assert (len(PCs[0]) == self._node_num[2] / 2), (len(PCs[0]), self._node_num[2] / 2)
        else:
            PCs = mid_result_1
            assert (len(PCs[0]) == self._node_num[2])

        return PCs

    def get_output_data(self, input_data=None, num_of_PCs=None):
        if input_data is None: input_data = self._data_set
        output_data = np.array([self._molecule_net.activate(item) for item in input_data])
        dim_of_output = self._node_num[-1]
        if (not self._hierarchical) or num_of_PCs is None:
            output_data = [item[- dim_of_output:] for item in output_data]
        else:
            output_data = [item[(num_of_PCs - 1) * dim_of_output: num_of_PCs * dim_of_output] for item in output_data]
        return output_data

    def train(self):
        ####################### set up autoencoder begin #######################
        node_num = self._node_num
        num_of_PC_nodes_for_each_PC = 2 if self._hidden_layers_type[1] == CircularLayer else 1
        num_of_PCs = node_num[2] / num_of_PC_nodes_for_each_PC

        in_layer = self._in_layer_type(node_num[0], "IL")
        num_of_hidden_layers = len(self._hidden_layers_type)

        if self._hierarchical:
            if num_of_hidden_layers == 3:  # 5-layer autoencoder
                hidden_layers = [(self._hidden_layers_type[0])(node_num[1], "HL0"),
                                 (self._hidden_layers_type[1])(node_num[2], "PC"),
                                 (self._hidden_layers_type[2])(node_num[3] * num_of_PCs, "HL2")]
                bias_layers = [BiasUnit("B0"), BiasUnit("B1"), BiasUnit("B2"), BiasUnit("B3")]
            else:
                raise Exception("this num of hidden layers is not implemented")

            out_layer = self._out_layer_type(node_num[num_of_hidden_layers + 1] * num_of_PCs, "OL")

            parts_of_PC_layer = [ModuleSlice(hidden_layers[1], outSliceFrom=item * num_of_PC_nodes_for_each_PC,
                                             outSliceTo=(item + 1) * num_of_PC_nodes_for_each_PC)
                                 for item in range(num_of_PCs)]
            parts_hidden_2_layer = [ModuleSlice(hidden_layers[2],
                                                inSliceFrom=item * node_num[3], inSliceTo=(item + 1) * node_num[3],
                                                outSliceFrom=item * node_num[3], outSliceTo=(item + 1) * node_num[3])
                                    for item in range(num_of_PCs)]
            parts_output_layer = [ModuleSlice(out_layer, inSliceFrom=item * node_num[4],
                                              inSliceTo=(item + 1) * node_num[4])
                                  for item in range(num_of_PCs)]
            self._parts_of_PC_layer = parts_of_PC_layer
            self._parts_hidden_2_layer = parts_hidden_2_layer
            self._parts_output_layer = parts_output_layer

            layers_list = [in_layer] + hidden_layers + [out_layer]

            molecule_net = FeedForwardNetwork()

            molecule_net.addInputModule(in_layer)
            for item in (hidden_layers + bias_layers):
                molecule_net.addModule(item)

            molecule_net.addOutputModule(out_layer)

            connection_between_layers = list(range(num_of_hidden_layers + 1))
            connection_with_bias_layers = list(range(num_of_hidden_layers + 1))

            # set up full connections
            for i in range(2):
                connection_between_layers[i] = FullConnection(layers_list[i], layers_list[i + 1])
                connection_with_bias_layers[i] = FullConnection(bias_layers[i], layers_list[i + 1])
                molecule_net.addConnection(connection_between_layers[i])  # connect two neighbor layers
                molecule_net.addConnection(connection_with_bias_layers[i])

            # set up shared connections
            connection_with_bias_layers[2] = MotherConnection(node_num[3])
            connection_with_bias_layers[3] = MotherConnection(node_num[4])
            connection_between_layers[3] = MotherConnection(node_num[3] * node_num[4])
            connection_between_layers[2] = [MotherConnection(node_num[3] * num_of_PC_nodes_for_each_PC) for _ in
                                            range(num_of_PCs)]

            for _1 in range(num_of_PCs):
                molecule_net.addConnection(SharedFullConnection(connection_with_bias_layers[2],
                                                                bias_layers[2], parts_hidden_2_layer[_1]))
                molecule_net.addConnection(SharedFullConnection(connection_with_bias_layers[3],
                                                                bias_layers[3], parts_output_layer[_1]))

            for _1 in range(num_of_PCs):
                molecule_net.addConnection(SharedFullConnection(connection_between_layers[3], parts_hidden_2_layer[_1],
                                                                parts_output_layer[_1]))

            for _1 in range(num_of_PCs):
                for _2 in range(_1, num_of_PCs):
                    molecule_net.addConnection(SharedFullConnection(connection_between_layers[2][_1],
                                                                    parts_of_PC_layer[_1], parts_hidden_2_layer[_2]))
        else:
            if num_of_hidden_layers == 3:  # 5-layer autoencoder
                hidden_layers = [(self._hidden_layers_type[0])(node_num[1], "HL1"),
                                 (self._hidden_layers_type[1])(node_num[2], "HL2"),
                                 (self._hidden_layers_type[2])(node_num[3], "HL3")]
                bias_layers = [BiasUnit("B1"), BiasUnit("B2"), BiasUnit("B3"), BiasUnit("B4")]
            elif num_of_hidden_layers == 1:
                hidden_layers = [(self._hidden_layers_type[0])(node_num[1], "HL1")]
                bias_layers = [BiasUnit("B1"), BiasUnit("B2")]
            else:
                raise Exception("this num of hidden layers is not implemented")

            out_layer = self._out_layer_type(node_num[num_of_hidden_layers + 1], "OL")

            layers_list = [in_layer] + hidden_layers + [out_layer]

            molecule_net = FeedForwardNetwork()

            molecule_net.addInputModule(in_layer)
            for item in (hidden_layers + bias_layers):
                molecule_net.addModule(item)

            molecule_net.addOutputModule(out_layer)

            connection_between_layers = list(range(num_of_hidden_layers + 1))
            connection_with_bias_layers = list(range(num_of_hidden_layers + 1))

            for i in range(num_of_hidden_layers + 1):
                connection_between_layers[i] = FullConnection(layers_list[i], layers_list[i + 1])
                connection_with_bias_layers[i] = FullConnection(bias_layers[i], layers_list[i + 1])
                molecule_net.addConnection(connection_between_layers[i])  # connect two neighbor layers
                molecule_net.addConnection(connection_with_bias_layers[i])

        molecule_net.sortModules()  # this is some internal initialization process to make this module usable

        ####################### set up autoencoder end #######################

        trainer = BackpropTrainer(molecule_net, learningrate=self._network_parameters[0],
                                  momentum=self._network_parameters[1],
                                  weightdecay=self._network_parameters[2],
                                  lrdecay=self._network_parameters[3],
                                  verbose=self._network_verbose)

        sincos = self._data_set[::self._training_data_interval]  # pick some of the data to train
        data_as_input_to_network = sincos

        if self._hierarchical:
            data_set = SupervisedDataSet(node_num[0], num_of_PCs * node_num[num_of_hidden_layers + 1])
            for item in data_as_input_to_network:
                data_set.addSample(item, list(item) * num_of_PCs)
        else:
            data_set = SupervisedDataSet(node_num[0], node_num[num_of_hidden_layers + 1])
            for item in data_as_input_to_network:
                data_set.addSample(item, item)

        training_print_info = '''training network with index = %d, training maxEpochs = %d, structure = %s, layers = %s, num of data = %d,
parameter = [learning rate: %f, momentum: %f, weightdecay: %f, lrdecay: %f]\n''' % \
                              (self._index, self._max_num_of_training, str(self._node_num),
                               str(self._hidden_layers_type).replace("class 'pybrain.structure.modules.", ''),
                               len(data_as_input_to_network),
                               self._network_parameters[0], self._network_parameters[1], self._network_parameters[2],
                               self._network_parameters[3],)

        print("Start " + training_print_info)

        trainer.trainUntilConvergence(data_set, maxEpochs=self._max_num_of_training)

        self._connection_between_layers_coeffs = [item.params for item in connection_between_layers]
        self._connection_with_bias_layers_coeffs = [item.params for item in connection_with_bias_layers]

        print('Done ' + training_print_info)

        self._trainer = trainer
        self._molecule_net = molecule_net
        return self


class autoencoder_Keras(autoencoder):
    def _init_extra(self,
                    network_parameters = CONFIG_4,
                    batch_size = 100,
                    enable_early_stopping=True
                    ):
        self._network_parameters = network_parameters
        self._batch_size = batch_size
        self._enable_early_stopping = enable_early_stopping
        self._molecule_net_layers = None              # why don't I save molecule_net (Keras model) instead? since it it not picklable:
                                                      # https://github.com/luispedro/jug/issues/30
                                                      # https://keras.io/getting-started/faq/#how-can-i-save-a-keras-model
        return

    def get_output_data(self, input_data=None, num_of_PCs = None):
        if input_data is None: input_data = self._data_set
        temp_model = Sequential()
        for item in self._molecule_net_layers:
            temp_model.add(item)

        return temp_model.predict(input_data)

    def get_mid_result(self, input_data=None):
        """The out format of this function is different from that in Pybrain implementation"""
        if input_data is None: input_data = self._data_set
        temp_model = Sequential()
        temp_model_bak = temp_model
        result = []
        for item in self._molecule_net_layers[:-2]:
            temp_model = temp_model_bak
            temp_model.add(item)
            temp_model_bak = copy.deepcopy(temp_model)   # this backup is required to get the correct results, no idea why
            result.append(temp_model.predict(input_data))

        return result

    def get_PCs(self, input_data=None):
        if input_data is None: input_data = self._data_set
        temp_model = Sequential()
        if hasattr(self, '_output_as_circular') and self._output_as_circular:  # use hasattr for backward compatibility
            for item in self._molecule_net_layers[:-5]:
                temp_model.add(item)
        else:
            for item in self._molecule_net_layers[:-2]:
                temp_model.add(item)

        if self._hidden_layers_type[1] == CircularLayer:
            PCs = [[acos(item[2 * _1]) * np.sign(item[2 * _1 + 1]) for _1 in range(len(item) / 2)]
                   for item in temp_model.predict(input_data)]
            assert (len(PCs[0]) == self._node_num[2] / 2), (len(PCs[0]), self._node_num[2] / 2)
        elif self._hidden_layers_type[1] == TanhLayer:
            print len(input_data)
            PCs = temp_model.predict(input_data)
            assert (len(PCs[0]) == self._node_num[2])
        else:
            raise Exception("PC layer type error")

        return PCs

    def train(self):    
        node_num = self._node_num
        num_of_PC_nodes_for_each_PC = 2 if self._hidden_layers_type[1] == CircularLayer else 1
        num_of_PCs = node_num[2] / num_of_PC_nodes_for_each_PC
        data = self._data_set
        if hasattr(self, '_output_data_set') and not self._output_data_set is None:
            print ("outputs different from inputs")
            output_data_set = self._output_data_set
        else:
            output_data_set = data

        num_of_hidden_layers = len(self._hidden_layers_type)
        if self._hierarchical:
            raise Exception('hierarchical version not implemented')
        elif num_of_hidden_layers != 3:
            raise Exception('not implemented for this case')
        else:
            molecule_net = Sequential()
            molecule_net.add(Dense(input_dim=node_num[0], output_dim=node_num[1], activation='tanh',W_regularizer=l2(self._network_parameters[4][0])))   # input layer
            if self._hidden_layers_type[1] == CircularLayer:
                molecule_net.add(Dense(input_dim=node_num[1], output_dim=node_num[2], activation='linear',W_regularizer=l2(self._network_parameters[4][1])))
                molecule_net.add(Reshape((node_num[2] / 2, 2), input_shape=(node_num[2],)))
                molecule_net.add(Lambda(temp_lambda_func_for_circular_for_Keras))  # circular layer
                molecule_net.add(Reshape((node_num[2],)))
                molecule_net.add(Dense(input_dim=node_num[2], output_dim=node_num[3], activation='tanh',W_regularizer=l2(self._network_parameters[4][2])))
                molecule_net.add(Dense(input_dim=node_num[3], output_dim=node_num[4], activation='linear',W_regularizer=l2(self._network_parameters[4][3])))

            elif self._hidden_layers_type[1] == TanhLayer:
                molecule_net.add(Dense(input_dim=node_num[1], output_dim=node_num[2], activation='tanh',W_regularizer=l2(self._network_parameters[4][1])))
                molecule_net.add(Dense(input_dim=node_num[2], output_dim=node_num[3], activation='tanh',W_regularizer=l2(self._network_parameters[4][2])))
                molecule_net.add(Dense(input_dim=node_num[3], output_dim=node_num[4], activation='linear',W_regularizer=l2(self._network_parameters[4][3])))
            else:
                raise Exception ('this type of hidden layer not implemented')

            if hasattr(self, '_output_as_circular') and self._output_as_circular:
                molecule_net.add(Reshape((node_num[4] / 2, 2), input_shape=(node_num[4],)))
                molecule_net.add(Lambda(temp_lambda_func_for_circular_for_Keras))  # circular layer
                molecule_net.add(Reshape((node_num[4],)))

            molecule_net.compile(loss='mean_squared_error', metrics=['accuracy'],
                                 optimizer=SGD(lr=self._network_parameters[0],
                                               momentum=self._network_parameters[1],
                                               decay= self._network_parameters[2],
                                               nesterov=self._network_parameters[3])
                                 )

            training_print_info = '''training network with index = %d, training maxEpochs = %d, structure = %s, layers = %s, num of data = %d,
parameter = [learning rate: %f, momentum: %f, lrdecay: %f, regularization coeff: %s], output as circular = %s\n''' % \
                                  (self._index, self._max_num_of_training, str(self._node_num),
                                   str(self._hidden_layers_type).replace("class 'pybrain.structure.modules.", ''),
                                   len(data),
                                   self._network_parameters[0], self._network_parameters[1],
                                   self._network_parameters[2], str(self._network_parameters[4]), str(self._output_as_circular))

            print("Start " + training_print_info + str(datetime.datetime.now()))
            call_back_list = []
            earlyStopping = EarlyStopping(monitor='val_loss', patience=30, verbose=0, mode='min')
            if self._enable_early_stopping:
                call_back_list += [earlyStopping]


            molecule_net.fit(data, output_data_set, nb_epoch=self._max_num_of_training, batch_size=self._batch_size,
                             verbose=int(self._network_verbose), validation_split=0.2, callbacks=call_back_list)

            dense_layers = [item for item in molecule_net.layers if isinstance(item, Dense)]
            for _1 in range(len(dense_layers)):
                assert(dense_layers[_1].get_weights()[0].shape[0] == node_num[_1]), (dense_layers[_1].get_weights()[0].shape[1], node_num[_1])   # check shapes of weights

            self._connection_between_layers_coeffs = [item.get_weights()[0].T.flatten() for item in molecule_net.layers if isinstance(item, Dense)]  # transpose the weights for consistency
            self._connection_with_bias_layers_coeffs = [item.get_weights()[1] for item in molecule_net.layers if isinstance(item, Dense)]

            print('Done ' + training_print_info + str(datetime.datetime.now()))

            self._molecule_net_layers = molecule_net.layers

        return self

def temp_lambda_func_for_circular_for_Keras(x):
    """This has to be defined at the module level here, otherwise the pickle will not work
    """
    return x / ((x ** 2).sum(axis=2, keepdims=True).sqrt())
