__version_info__ = (1,1,1)
__version__ = '.'.join(map(str, __version_info__))
__author__ = 'Jan Eberhage, Institute for Biophysical Chemistry, Hannover Medical School (eberhage.jan@mh-hannover.de)'

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import argparse
import pickle
import json

def convert(x):
  if hasattr(x, "tolist"):  # numpy arrays have this
    return x.tolist()
  raise TypeError(x)

def find_pkl_models(input_dir, model_num=0):
  model_names = []
  for path, currentDirectory, files in os.walk(input_dir):
    for file in files:
      if file.startswith('result') and file.endswith('.pkl'):
        model_names.append(os.path.join(path, file))
  model_names.sort()
  print(f'Found {str(len(model_names))} models')
  if not model_names:
    if os.path.exists(os.path.join(input_dir, 'pae_plddt.json')):
      print(f'JSON file found. Try "python3 '+' '.join(sys.argv)+' --jsonload pae_plddt.json"')
    sys.exit('Exiting')
  if model_num:
    print(f'Unpickling the first {model_num} models. Skipping the rest. Please wait.')
    model_names = model_names[:model_num]
  else:
    print('Unpickling all models. Please wait.')
  return model_names

def get_pae_plddt_from_pkl(model_names, input_dir):
  out = {}
  for i,name in enumerate(model_names):
    shortname = name.replace(os.path.join(input_dir, 'result_'),'').replace('multimer_v2_pred_','').replace('pred_','').replace('.pkl','')
    print(f'Loading {name}')
    d = pickle.load(open(name,'rb'))
    out[name] = {'short_name': shortname, 'plddt': d['plddt']}
    if 'predicted_aligned_error' in d:
      out[name]['pae'] = d['predicted_aligned_error']
  return out

def get_pae_plddt_from_json(json_path):
  print(f'Reading {json_path}')
  with open(json_path) as handle:
    content = json.loads(handle.read())
  print('Found '+str(len(content))+' models in the provided JSON file')
  return content

def generate_json_dump(pae_plddt_per_model, out_dir):
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)
  print('Generating pae_plddt.json in the output directory for further usage. Remember to also keep features.pkl.')
  with open(os.path.join(out_dir, 'pae_plddt.json'), 'w') as f:
    json.dump(pae_plddt_per_model, f, indent = 2, default=convert)  

def generate_output_images(feature_dict, out_dir, name, pae_plddt_per_model):
  print('Generating plots in '+out_dir)
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)
  msa = feature_dict['msa']
  indexes = feature_dict['residue_index']
  chain_starts = [i for i,x in enumerate(indexes) if x == 0]
  seqid = (np.array(msa[0] == msa).mean(-1))
  seqid_sort = seqid.argsort()
  non_gaps = (msa != 21).astype(float)
  non_gaps[non_gaps == 0] = np.nan
  final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

  ##################################################################
  plt.figure(figsize=(14, 4), dpi=100)
  ##################################################################
  plt.subplot(1, 2, 1)
  plt.title("Sequence coverage")
  plt.imshow(final,
             interpolation='nearest', aspect='auto',
             cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
  plt.plot((msa != 21).sum(0), color='black')
  plt.xlim(-0.5, msa.shape[1] - 0.5)
  plt.ylim(-0.5, msa.shape[0] - 0.5)
  if len(chain_starts) > 1:
    for chain_break in chain_starts[1:]:
      plt.plot([chain_break,chain_break],[-0.5, msa.shape[0] - 0.5],color="black",linewidth=1)
  plt.colorbar(label="Sequence identity to query", )
  plt.xlabel("Positions")
  plt.ylabel("Sequences")	
  ##################################################################
  plt.subplot(1, 2, 2)
  plt.title("Predicted LDDT per position")
  for model_name, value in pae_plddt_per_model.items():
    plt.plot(value["plddt"], label=value["short_name"])
  if len(chain_starts) > 1:
    for chain_break in chain_starts[1:]:
      plt.plot([chain_break,chain_break],[0,100],color="black",linewidth=1)
  if len(pae_plddt_per_model) < 6:
    plt.legend()
  plt.ylim(0, 100)
  plt.ylabel("Predicted LDDT")
  plt.xlabel("Positions")
  plt.savefig(f"{out_dir}/{name+('_' if name else '')}coverage_LDDT.png")
  ##################################################################

  ##################################################################
  if 'pae' in pae_plddt_per_model[list(pae_plddt_per_model.keys())[0]]:
    num_models = len(pae_plddt_per_model)
    horizontal = np.ceil(np.sqrt(num_models)).astype(int)
    vertical = np.ceil(num_models/horizontal).astype(int)
    fig = plt.figure(figsize=(3 * horizontal, 2.5 * vertical), dpi=300)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
      plt.subplot(vertical, horizontal, n + 1)
      plt.title(value["short_name"])
      plt.imshow(value["pae"], label=value["short_name"], cmap="bwr", vmin=0, vmax=30)
      if len(chain_starts) > 1:
        for chain_break in chain_starts[1:]:
          plt.plot([chain_break,chain_break],[0,len(indexes)],color="black",linewidth=1)
          plt.plot([0,len(indexes)],[chain_break,chain_break],color="black",linewidth=1) 
      plt.xlim(0, len(indexes))
      plt.ylim(len(indexes), 0)
      clb = plt.colorbar()
      clb.ax.set_title('Å')
    fig.tight_layout()
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}PAE.png")
  else:
    print('Unable to plot PAE. Try using "--model_preset=monomer_ptm" for your next Alphafold monomer job.')
  ##################################################################

parser = argparse.ArgumentParser(description='This script will generate plots containing the MSA, pLDDT distribution and \
				 Predicted Alignment Error (PAE) of a given AlphaFold output using the Pickle files (.pkl).')
required = parser.add_argument_group('required arguments')
required.add_argument('-i','--input_dir',dest='input_dir',metavar='<input_dir>',required=True,help='relative or absolute path to the input directory (AlphaFold output)')
parser.add_argument('-o','--output_dir',dest='output_dir',default='',metavar='<output_dir>',help='destination folder where files are generated')
parser.add_argument('-n','--name',dest='name',default='',metavar='<prefix>',help='add custom name as prefix to your images')
parser.add_argument('-m','--models',dest='models',default=0,type=int,metavar='<n>',help='limit the inspected pickles to n models')
parser.add_argument('--jsondump', action='store_true',help='skip the plotting and dump all relevant PAE and pLDDT information as human readable JSON file')
parser.add_argument('--jsonload',dest='json',default='',metavar='<file>',help='JSON file in the input directory to be read instead of pkl files')
parser.add_argument('-v', '--version', action='version', version='%(prog)s ('+__version__+') '+' by '+__author__)
groups_order = {
    'positional arguments': 0,
    'required arguments': 1,
    'optional arguments': 2
}
parser._action_groups.sort(key=lambda g: groups_order[g.title])
args = parser.parse_args()

feature_path = os.path.join(args.input_dir, 'features.pkl')
if not os.path.exists(feature_path):
  sys.exit(f'The file "{feature_path}" is mandatory. Please provide it at this specific location.')
else:
  feature_dict = pickle.load(open(f'{feature_path}','rb'))

if args.json:
  json_path = os.path.join(args.input_dir, args.json)
  if not os.path.exists(json_path):
    sys.exit(f'The file "{json_path}" is mandatory. Please provide it at this specific location.')
  else:
    pae_plddt_per_model = get_pae_plddt_from_json(json_path)
else:
  model_pkls = find_pkl_models(args.input_dir, args.models)
  pae_plddt_per_model = get_pae_plddt_from_pkl(model_pkls, args.input_dir)

if args.jsondump:
  generate_json_dump(pae_plddt_per_model, args.output_dir if args.output_dir else args.input_dir)
else:
  generate_output_images(feature_dict, args.output_dir if args.output_dir else args.input_dir, args.name, pae_plddt_per_model)
