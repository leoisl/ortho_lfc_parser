# ortho_lfc_parser

A tool that parses an orthogroup file and a set of differential expression files and outputs positive, negative,
both positive and negative, and non significant events.

# Dependencies

* `python 3.7+`

# Running on paper data

```
python ortho_lfc_parser/main.py --orthogroups data/Orthogroups.csv data/Emu_LarvaXAdult_DifExp.csv data/Hmi_LarvaXAdult_DifExp.csv data/Mco_LarvaXAdult_DifExp.csv
```

Output can be found in dir `out`. These will be `tsv` files containing the positive, negative, both positive and negative,
non significant and all events.

# Usage

```
usage: main.py [-h] --orthogroups ORTHOGROUPS [--outdir OUTDIR] DE_files [DE_files ...]

Orthogroup LFC parser.

positional arguments:
  DE_files              Differential expression files

optional arguments:
  -h, --help            show this help message and exit
  --orthogroups ORTHOGROUPS
                        The orthogroups file
  --outdir OUTDIR       Output directory

```
