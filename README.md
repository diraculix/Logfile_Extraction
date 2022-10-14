# Logfile_Extraction 
> _Proton therapy irradiation record file data extraction, conversion, analysis and treatment dose reconstruction._

1. [logfile_extractor.py](##logfile_extractor.py)
    1. [init( )](###init())
    1. [prepare_dataframe( )](###prepare_dataframe())
1. [logzip_converter.py](##logzip_converter.py)
  

## logfile_extractor.py

Main data extraction script for pre-processed (unpacked) .csv log-files, following a hierarchy:

```
Logfiles
└───Fraction-ID(20210911)
    └───Beam(01)
    └───Beam(02)
    └───Beam(03)
        └───..
            └───beam.*.csv
            └───beam_config.*.csv
            └───events.*.csv
            └───map_record_XXX*
                └───*_part_XX*.csv
                └───*_tuning_XX.csv
            └───map_specif_XXX*
                └───*_part_XX*.csv
                └───*_tuning_XX.csv
```

File | Function
:--- | :---
`beam.csv` | World
`beam_config.csv` | World
`map_record*.csv` | World
`map_specif*.csv` | World

### init( )

Init is init.

### prepare_dataframe( )

hey


## logzip_converter.py
hola
