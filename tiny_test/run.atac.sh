BASEDIR=/scratch/cwl_working_dir/base
OUTDIR=./test_run/
WORKDIR=/scratch/cwl_working_dir/work
TMPDIR=/scratch/cwl_working_dir/tmp
TMPOUTDIR=/scratch/cwl_working_dir/tmp_out
CWLSCRIPT=../CWL/workflows/ATACseq_pipeline.cwl
INPUT=./test_main.atac.yml

SINGULARITY_LOCALCACHEDIR=/scratch/.singularity/cache
SINGULARITY_PULLFOLDER=/scratch/.singularity/pull
SINGULARITY_TMPDIR=/scratch/.singularity/tmp
export SINGULARITY_LOCALCACHEDIR
export SINGULARITY_PULLFOLDER
export SINGULARITY_TMPDIR

if [ ! -d "$OUTDIR" ]
then
	mkdir "$OUTDIR"
fi

cwltool --debug --singularity \
	--tmp-outdir-prefix $TMPOUTDIR --tmpdir-prefix $TMPDIR \
	--basedir $BASEDIR --outdir $OUTDIR \
	$CWLSCRIPT $INPUT

#cwltoil --singularity --logDebug --clean always --cleanWorkDir always \
#	--tmp-outdir-prefix $TMPOUTDIR --tmpdir-prefix $TMPDIR \
#	 --workDir $WORKDIR --basedir $BASEDIR --outdir $OUTDIR \
#	$CWLSCRIPT $INPUT
