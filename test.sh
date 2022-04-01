BinDir="/annoroad/data/PD/pipeline_PM/CisTrans"
SPDir="/annoroad/data/PD/results_PM/QB/QB_0238_20200621"
SP="QB20CT0144-1-A2"
perl ${BinDir}/CisTransOfAnnoFile_v2.pl \
	-var ${SPDir}/${SP}/searchDB/chemoAndTarget/targetResults/PanCancer_All/PanCancer_All_QB20CT0144-1-A2_snv.TargetMut.xls \
	-var ${SPDir}/${SP}/searchDB/chemoAndTarget/targetResults/PanCancer_All/PanCancer_All_QB20CT0144-1-A2_sindel.TargetMut.xls \
	-bam ${SPDir}/${SP}/alignment/${SP}.final.bam \
	-dir ${BinDir}/test
