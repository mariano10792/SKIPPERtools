while true;
do
	#make && ./milo.exe hits_skp_moduleC40_41-ssc16_17-lta20_60_test10_NROW50_NBINROW1_NCOL200_NBINCOL16_EXPOSURE0_CLEAR0_2_6.fz.root
	#make && ./milo.exe hits_skp_module4-lta2-mistica_TEMP135K-run1_NROW3200_NBINROW1_NCOL470_NBINCOL1_EXPOSURE0_CLEAR600_212.root
	make && ./milo.exe hits_skp_module2-lta5-ssc4_TEMP135K-run2_NROW500_NBINROW2_NCOL240_NBINCOL2_EXPOSURE0_CLEAR600_322.root
	sleep 5
	
done
