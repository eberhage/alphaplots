__version_info__ = (1,0,1)
__version__ = '.'.join(map(str, __version_info__))
__author__ = 'Jan Eberhage, Institute for Biophysical Chemistry, Hannover Medical School (eberhage.jan@mh-hannover.de)'

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import argparse
import pickle

def get_pae_plddt(model_names, input_dir):
  out = {}
  for i,name in enumerate(model_names):
    newname = name.replace(input_dir+'result_','').replace('multimer_v2_pred_','').replace('pred_','').replace('.pkl','')
    print(f'Loading {name} as {newname}')
    d = pickle.load(open(name,'rb'))
    out[newname] = {'plddt': d['plddt']}
    if 'predicted_aligned_error' in d:
      out[newname]['pae'] = d['predicted_aligned_error']
  return out

def generate_output_images(feature_dict, out_dir, name, pae_plddt_per_model):
  print('Writing files to '+out_dir)
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
    plt.plot(value["plddt"], label=model_name)
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
      plt.title(model_name)
      plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
      if len(chain_starts) > 1:
        for chain_break in chain_starts[1:]:
          plt.plot([chain_break,chain_break],[0,len(indexes)],color="black",linewidth=1)
          plt.plot([0,len(indexes)],[chain_break,chain_break],color="black",linewidth=1) 
      plt.xlim(0, len(indexes))
      plt.ylim(len(indexes), 0)
      plt.colorbar()
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
parser.add_argument('-v', '--version', action='version', version='%(prog)s ('+__version__+') '+' by '+__author__)
groups_order = {
    'positional arguments': 0,
    'required arguments': 1,
    'optional arguments': 2
}
parser._action_groups.sort(key=lambda g: groups_order[g.title])
args = parser.parse_args()

if not os.path.exists(args.input_dir+'/features.pkl'):
  sys.exit(f'The file "{args.input_dir}features.pkl" is mandatory. Please provide it at this specific location.')
else:
  feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl','rb'))
model_names = []
for path, currentDirectory, files in os.walk(args.input_dir):
  for file in files:
    if file.startswith('result') and file.endswith('.pkl'):
      model_names.append(path+file)
model_names.sort()
print(f'Found {str(len(model_names))} models')
if args.models:
  print(f'Unpickling the first {args.models} models. Skipping the rest. Please wait.')
  model_names = model_names[:args.models]
else:
  print('Unpickling all models. Please wait.')

pae_plddt_per_model = get_pae_plddt(model_names, args.input_dir)
generate_output_images(feature_dict, args.output_dir if args.output_dir else args.input_dir, args.name, pae_plddt_per_model)
