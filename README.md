# parse_genexus_archive
1. Сбилдить

```
docker build -t parse_genexus_archive ~/path/to/code
```

2. Запустить

```
docker run -v /host/path/to/data:/data -t -i parse_genexus_archive python3 main.py --genexus_archive /data/docker/path/to/zip  --reference /data/docker/path/to/fasta --output /data/out.fa
```

Префикс перед первым символом '_' из имени файла --genexus_archive идет в заголовок fasta файла

```
usage: main.py [-h] --genexus_archive GENEXUS_ARCHIVE --reference REFERENCE [--cores_num CORES_NUM] [--min_coverage MIN_COVERAGE] [--quality_threshold QUALITY_THRESHOLD] --output OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  --genexus_archive GENEXUS_ARCHIVE
                        Genexus output zip file
  --reference REFERENCE
                        Reference genome fasta file
  --cores_num CORES_NUM
                        Cores num for bowtie2, default = 16
  --min_coverage MIN_COVERAGE
                        Minimum coverage per base, default = 1
  --quality_threshold QUALITY_THRESHOLD
                        Minimum quality score (in phred) a base has to reach to be counted, default=15
  --output OUTPUT       Output fasta file name
```
