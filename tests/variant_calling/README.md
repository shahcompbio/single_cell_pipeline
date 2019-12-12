

```docker run -w `pwd` -v `pwd`:`pwd` -v /var/run/docker.sock:/var/run/docker.sock -v /usr/bin/docker:/usr/bin/docker --env PYTHONPATH=`pwd`:`pwd`/pypeliner:`pwd`/biowrappers singlecellpipeline/single_cell_pipeline:v0.5.4 single_cell variant_calling --input_yaml tests/variant_calling/input.yaml --out_dir tests/variant_calling/outputs --submit local --loglevel DEBUG --config_file tests/variant_calling/config.yaml --context_config context_config.yaml --tmpdir tmp --nocleanup --rerun
```
