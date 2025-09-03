# alphaplots.py
This script will scan a folder and use the Pickle (.pkl) files of an AlphaFold 2 output to generate the MSA, the per-AA plDDT distribution and (for each model) a predicted alignment error (PAE) plot.

See [alphaplots3](https://github.com/eberhage/alphaplots3) for the AlphaFold 3 version.

## Example Usage
Make plots from Pickle files. That's it.
```
python3 alphaplots.py -i /path/to/alphafold/result/output/folder
```
Make plots from Pickle files. Make JSON from Pickle files. Delete Pickle files.
```
python3 alphaplots.py -i /path/to/alphafold/result/output/folder --jsondump --rmpkl
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
To read a json file in the input directory (or relative to it) instead of the pkl files (usually faster)
```
--jsonload <file>
```
To auto-delete all model pkl files (not features.pkl) in the input directory. Should only be used after or in combination with jsondump
```
--rmpkl
```
To skip the plotting alltogether. Only makes sense with jsondump.
```
--noplot
```
To auto-answer every question with "yes". Use with caution. For experienced users only.
```
--yes
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
