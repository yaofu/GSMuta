#!/bin/sh 

function usage
{
	echo "
	Usage : bash driver_detection.sh -i xx.maf -o output_dir -m method['Poisson'] -bed coverage_file[''] -frac [0.6] -core [10] -rc [2] -bin [1] -interval [3] -pathway pathway_file['data/kegg.tab'] 
	
		-i	[Required] Input MAF file of somatic mutations
		-o	[Required] Output directory
		-m	[Optional] 'Poisson' or 'NB' (negative binomial), default: Poisson
		-bed	[Optional] Directory for sample specific coverage information, default: NA
		-frac	[Optional] Fraction of genes included for non-silent mutations to construct background, default: 0.6
		-core	[Optional] Number of cores used for parallelization, default: 10
		
		-rc	[Optional] For oncogene analysis, the lowest number mutations on gene to be included, default: 3
		-bin	[Optional] For oncogene analysis, the binning size, default: 1
		-interval	[Optional] For oncogene analysis, the distance to merger clusters, default: 3
		-pathway	[Optional] Input pathway, default: KEGG
"
}


OUTPUT_DIR=res
METHOD=Poisson
BEDS=NA
FRAC=0.6
CORE=10
OG_NUM=3
BIN=1
INTERVAL=3
PATHWAY=data/kegg.tab


if [ -e $1 ];then
usage
exit
fi

while [ "$1" != "" ]; do
    case $1 in
        -i | --file )           shift
                                INPUT_MAF=$1
				;;
	-o)			shift
				OUTPUT_DIR=$1
				;;
	-m)			shift 
				METHOD=$1
				;;
	-bed)			shift
				BEDS=$1
				;;
	-frac)			shift
				FRAC=$1;
				;;
	-core)			shift
				CORE=$1;
				;;
	-rc)			shift
				OG_NUM=$1;
				;;
	-bin)			shift
				BIN=$1;
				;;
	-interval)		shift
				INTERVAL=$1;
				;;
	-pathway)		shift
				PATHWAY=$1;
				;;
	-h | --help )           usage
                                exit

        esac
        shift
done

if [[ $INPUT_MAF == "" ]];then
usage
exit
fi


[ -d $OUTPUT_DIR ] || mkdir $OUTPUT_DIR

tag=`echo $INPUT_MAF | awk -F "[./]" '{print $(NF-1)}'`

MAF_cleanup_1="$OUTPUT_DIR/$tag-clean1.maf"
MAF_cleanup_2="$OUTPUT_DIR/$tag-clean2.maf"


echo "Parameters: -i $INPUT_MAF -o $OUTPUT_DIR -m $METHOD -bed $BEDS -frac $FRAC -core $CORE -rc $OG_NUM -bin $BIN -interval $INTERVAL -pathway $PATHWAY"
echo "" 

Rscript --no-save --no-restore  Code/MAF_cleanup.R $INPUT_MAF $MAF_cleanup_1
perl Code/Sequence_Retrieve.pl $MAF_cleanup_1 $MAF_cleanup_2
Rscript --no-save --no-restore Code/GSMuta.R $MAF_cleanup_2 $FRAC $OUTPUT_DIR $METHOD $CORE $BEDS $OG_NUM $BIN $INTERVAL $PATHWAY

