make
#ls ~/tools/SKIPPERtools/correlation_simulator/files/sim_vs_sim_halo0_hits_MC_N0_0_DC457_R*.root | xargs -L1 -P8 ./runKS_SENSEI.sh
for i in {0..99}
do
	#ls ~/tools/SKIPPERtools/correlation_simulator/files/sim_vs_sim_halo60_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR*.root | xargs -L1 -P8 ./runKS_SENSEI.sh
	ls ~/tools/SKIPPERtools/correlation_simulator/files/sim_vs_sim_halo0_hits_MC_N0_0_DC457_R18*.root | xargs -L1 -P1 ./runKS_SENSEI.sh
	hadd ./files/pvalues__${i}.root ./files/pvalues_run0_sim_vs_sim_halo0_hits_MC_N0_0_DC457_R*t
	sleep 10
done

