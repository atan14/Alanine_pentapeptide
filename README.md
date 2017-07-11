# Alanine_pentapeptide

python main_work.py

python auto_qsub.py "bash temp_plotting_iter_alanine_all.sh" --gpu 0 --submit

python temp_get_potential_centers_covering_explored_region.py ../resources.Alanine_pentapeptide/network_10.pkl ../resources/Alanine_pentapeptide/autoencoder_info_10.txt FES $fc $temperature $interval $total_steps $solvent_type $ensemble_type $num_of_points_each_dim > temp_simulation.sh

python auto_qsub.py "bash temp_simulation.sh" --gpu 0 --submit

python generate_coordinates.py Alanine_pentapeptide --path FES

python temp_generating_files_for_WHAM.py Bayes FES ../resources/Alanine_pentapeptide/network_10.pkl WHAM $num_of_dihedrals $dimensionality 0 0 0

bash build.sh WHAM 1

python BayesWHAM_plotter.py 
OR use BayesWHAM_plotter.m
