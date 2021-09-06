# Script used to run CCD_simulation_correlation.exe
# First argument is module
# Second argument is run number (9,12,14,15) -- (9,12,14) for 0hs
# The for is across the halos we want to analyse.

#for j in {10,20,40,60,80}
#for j in {100,110,120,130,140,150}

#for j in 60
#do 
#	make 
#	for i in /mnt/mariano/MINOSnewfirmware/commissioning_data/1pix1e/hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run$2*EXPOSURE72000*_$1_*t 
#		do ./CCD_simulation_correlation.exe $i $j
#	done
#done


make


for j in 0
do 
#	for i in /mnt/mariano/MINOSnewfirmware/commissioning_data/1pix1e/hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run$2*EXPOSURE72000*_$1_*t
	for i in /mnt/mariano/MINOSnewfirmware/commissioning_data/minos2/hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR600_$1_*t
		do ./CCD_simulation_correlation.exe $i $j
	done
#	for i in /mnt/mariano/MINOSnewfirmware/commissioning_data/1pix1e/hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run$2*EXPOSURE72000*_$1_*t
	for i in /mnt/mariano/MINOSnewfirmware/commissioning_data/minos3/hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR1800_$1_*t
		do ./CCD_simulation_correlation.exe $i $j
	done
done


