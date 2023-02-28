# Annotation Stats pipeline

Generate statistics for a (Rapid Release) Genebuild annotation.


## setup environment

create Python 3.10 or 3.9 virtual environment
```
# e.g.
pyenv virtualenv 3.10.9 annotation_stats
```

install Python dependencies
```
poetry install
```


## run the pipeline

```
module load nextflow-22.10.1-gcc-11.2.0-ju5saqw

nextflow annotation_stats.nf
```
