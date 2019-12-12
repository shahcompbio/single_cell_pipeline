

```
samtools view -b -L tests/split_bams/snv_regions.bed ../remixt/HCC1395_chr15.bam > tests/split_bams/HCC1395_chr15.bam
samtools view -b -L tests/split_bams/snv_regions.bed ../remixt/HCC1395BL_chr15.bam > tests/split_bams/HCC1395BL_chr15.bam
samtools index tests/split_bams/HCC1395_chr15.bam
samtools index tests/split_bams/HCC1395BL_chr15.bam
```

```
docker run -w `pwd` -v `pwd`:`pwd` -v /var/run/docker.sock:/var/run/docker.sock -v /usr/bin/docker:/usr/bin/docker --env PYTHONPATH=`pwd`:`pwd`/pypeliner singlecellpipeline/single_cell_pipeline:v0.5.4 single_cell split_wgs_bam --input_yaml tests/split_bams/inputs_normal.yaml --out_dir tests/variant_calling/inputs/normal/ --submit local --loglevel DEBUG --config_file tests/split_bams/config.yaml --context_config context_config.yaml --tmpdir tmp
```

```
docker run -w `pwd` -v `pwd`:`pwd` -v /var/run/docker.sock:/var/run/docker.sock -v /usr/bin/docker:/usr/bin/docker --env PYTHONPATH=`pwd`:`pwd`/pypeliner singlecellpipeline/single_cell_pipeline:v0.5.4 single_cell split_wgs_bam --input_yaml tests/split_bams/inputs_tumour.yaml --out_dir tests/variant_calling/inputs/tumour/ --submit local --loglevel DEBUG --config_file tests/split_bams/config.yaml --context_config context_config.yaml --tmpdir tmp
```
