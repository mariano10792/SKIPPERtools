#for j in {10,20,40,60,80}
#for j in 60
#do
#	for i in ~/tools/SKIPPERtools/correlation_simulator/sim_halo${j}_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR600_*.root 0; 
#		do ./KS_SENSEI.exe $i; 
#	done
#done


#for r in {0..99}
#for r in 0
#do
#	for j in 60
#	do
#		#for i in ~/tools/SKIPPERtools/correlation_simulator/sim_vs_sim_halo${j}_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run12_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR600_1_189.root ;
#		for i in ~/tools/SKIPPERtools/correlation_simulator/files/sim_vs_sim_halo60_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run14_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR600_1_101.root;
#			do ./KS_SENSEI.exe $i $r; 
#		done
#	done
#done

make
#for i in ~/tools/SKIPPERtools/correlation_simulator/files/sim_vs_sim_halo60_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_1_*t;
for i in ~/tools/SKIPPERtools/correlation_simulator/files/sim_halo60_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_1_*t; 
#for i in ~/tools/SKIPPERtools/correlation_simulator/files/sim_halo60_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_2_*t; 
#for i in ~/tools/SKIPPERtools/correlation_simulator/files/sim_halo0_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_2_*t; 
#for i in ~/tools/SKIPPERtools/correlation_simulator/files/sim_halo0_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_1_*t; 
	do ./KS_SENSEI.exe $i 0; 
done
