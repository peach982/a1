
source activate /Lustre01/tangqianzi/anaconda3forSeurat4/envs/R4.0.2
export HDF5_USE_FILE_LOCKING='FALSE'

cd /Lustre03/data/tangqianzi/Gut_unroll/outputs/monocle23_time_course/nutrient_MEs/
Rscript 08.plot_pseudotime_and_ME.R sugar
Rscript 08.plot_pseudotime_and_ME.R amino_acid
Rscript 08.plot_pseudotime_and_ME.R bile_salt
Rscript 08.plot_pseudotime_and_ME.R inorganic_solutes
Rscript 08.plot_pseudotime_and_ME.R lipid
Rscript 08.plot_pseudotime_and_ME.R metal_ion
Rscript 08.plot_pseudotime_and_ME.R Nucleotide
Rscript 08.plot_pseudotime_and_ME.R organic_solutes
Rscript 08.plot_pseudotime_and_ME.R vitamin_absorption
Rscript 08.plot_pseudotime_and_ME.R water
