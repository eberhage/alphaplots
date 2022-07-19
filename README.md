# alphaplots.py
This script will scan a folder and use the Pickle (.pkl) files of an AlphaFold output to generate the MSA, the per-AA plDDT distribution and (for each model) a predicted alignment error (PAE) plot.

## Usage
```
python3 alphaplots.py -i /path/to/alphafold/result/output/folder
```
## Parameters
### Required
```
-i <input_dir>
--input_dir <input_dir>
```
### Optional
To set a destination for the generated plots
```
-o <output_dir>
--output_dir <output_dir>
(default = <input_dir>)
```
To add a prefix to the output file names
```
-n <prefix>
--name <prefix>
(default = None)
```
To limit the plot generation to the first (alphabetical) <number_of_models> models. Otherwise all pkl files will be plotted.
```
-m <number_of_models>
--models <number_of_models>
(default = 0 = all)
```
To generate a json file containing the information of pLDDT and PAE
```
--jsondump
```
To read a json file in the input directory instead of the pkl files (usually faster)
```
--jsonload <file>
```
Version
```
-v
--version
```
Help
```
-h
--help
```
