# Args:
#   [1] : FASTQ - FASTQ file
#   [2] : OUT_HTML - File where .html report will be written
#   [3] : OUT_ZIP - File where .zip report files will be written

RSCRIPT=$1

RLIBS=$2

RUN_HMM=$3

TUM_FILE=$4

GC_FILE=$5

MAP_FILE=$6

OUT_DIR=$7

OUT_BASENAME=$8

MAP_CUTOFF=$9

NUM_STATES=${10}

MU=${11}

M=${12}

KAPPA=${13}

E=${14}

S=${15}


export R_LIBS=$RLIBS:$R_LIBS


$RSCRIPT $RUN_HMM --tumour_file=$TUM_FILE --gc_file=$GC_FILE --map_file=$MAP_FILE --out_dir=$OUT_DIR --out_basename=$OUT_BASENAME  --map_cutoff=$MAP_CUTOFF --num_states=$NUM_STATES --param_mu=$MU --param_m=$M --param_k=$KAPPA --param_e=$E --param_s=$S
