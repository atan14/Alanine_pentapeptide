#!/bin/bash

wham_dir=$1
whether_reweight=$2
dim=3 
periodicity=[0,0,0] 
T=300
harmonicBiasesFile="${wham_dir}/bias/harmonic_biases.txt"
histBinEdgesFile="${wham_dir}/hist/hist_binEdges.txt"
histDir="${wham_dir}/hist"
tol_WHAM=1E-8 
maxIter_WHAM=1E6 
steps_MH=1E4 
saveMod_MH=1E3 
printMod_MH=1E3 
maxDpStep_MH=5E-6 
seed_MH=200184 
prior="none" 
alpha=1.0 

python BayesWHAM.py $dim $periodicity $T $harmonicBiasesFile $histBinEdgesFile $histDir $tol_WHAM $maxIter_WHAM $steps_MH $saveMod_MH $printMod_MH $maxDpStep_MH $seed_MH $prior $alpha


if [[ $2 == "1" ]]; then
	T=300
	dim_UMB=3
	periodicity_UMB=[0,0,0]
	harmonicBiasesFile_UMB=${wham_dir}/bias/harmonic_biases.txt
	trajDir_UMB="${wham_dir}/traj"
	histBinEdgesFile_UMB="${wham_dir}/hist/hist_binEdges.txt"
	histDir_UMB="${wham_dir}/hist"
	fMAPFile_UMB="f_MAP.txt"
	fMHFile_UMB="f_MH.txt"
	trajDir_PROJ="${wham_dir}/traj_proj"
	histBinEdgesFile_PROJ="${wham_dir}/hist/hist_binEdges_proj.txt"

	python BayesReweight.py $T $dim_UMB $periodicity_UMB $harmonicBiasesFile_UMB $trajDir_UMB $histBinEdgesFile_UMB $histDir_UMB $fMAPFile_UMB $fMHFile_UMB $trajDir_PROJ $histBinEdgesFile_PROJ

fi
