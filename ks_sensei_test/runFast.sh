make
#for i in {0..9999}
#do
#	#ls $0 | xargs -L1 -P1 ./runKS_SENSEI.sh $i
#	./KS_SENSEI_test.exe $i
#done

cat 10000.txt | xargs -L1 -P8 ./runKS_SENSEI.sh

