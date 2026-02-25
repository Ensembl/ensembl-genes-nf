'''
This script read the fastqc output fastqc_data.txt and extract the top adapter sequence
for the sample and write it to a file called adapters.fa

Usage: python get_adapters.py -i <fastqc_data.txt> -a <adapters.tsv> -o <adapters.fa>

'''

import argparse
import pandas as pd
from io import StringIO


def get_adapter_dict(adapter_path: str) -> dict:
    '''
    parse the adapter file and return a dictionary of adapter sequences and names

    Input: 
        adapter file path

    Output:
        dictionary of adapter sequences and names
    '''
    df = pd.read_csv(adapter_path, sep='\t', header=None, comment='#')
    print(df)
    return df.set_index(0)[1].to_dict()


def get_adapter_module(fastqc_data_path: str) -> list:
    '''
    parse the fastqc data file and return the adapter module as a string 

    Input: 
        fastqc_data.txt path

    Output:
        adapter module as a list of strings
    '''
    adapter_module = []
    with open(fastqc_data_path, 'r') as f:
        for line in f:
            if line.startswith('>>Adapter'):
                adapter_module.append(f"%{line}")

            if adapter_module != []:
                if line.startswith('>>END_MODULE'):
                    break
                elif line.startswith('>>'):
                    continue
                else:
                    adapter_module.append(line)

    return adapter_module


def get_top_adapter(adapter_module: list, adapter_dict: dict) -> str:
    '''
    parse the adapter module and return the top adapter sequence

    Input:
        adapter module as a list of strings
        adapter_dict as a dictionary of adapter sequences and names

    Output:
        top adapter sequence as a string
    '''
    if len(adapter_module) > 1:

        module_string = ''.join(adapter_module)
        str = StringIO(module_string)
        adapter_df = pd.read_csv(str, sep='\t', comment='%').drop(columns=['#Position'])
        for col in adapter_df.columns:
            adapter_df[col] = adapter_df[col].astype(float)

        sums = adapter_df.sum(axis=0)
        top_adapter = sums.idxmax()
        with open(args.output, 'w') as f:
            f.write(f'>{top_adapter}\n')
            f.write(f'{adapter_dict[top_adapter].strip()}\n')
        return adapter_dict[top_adapter]
    else:
        with open(args.output, 'w') as f:
            f.write(f'\n')


def main(args):
    adapter_module = get_adapter_module(args.input) 
    adapter_dict = get_adapter_dict(args.adapters)
    top_adapter = get_top_adapter(adapter_module, adapter_dict)


    return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get top adapter sequence from fastqc_data.txt')
    parser.add_argument('-i', '--input', help='fastqc_data.txt file', required=True)
    parser.add_argument('-a', '--adapters', help='adapter sequence and names tab delimed file', required=True)
    parser.add_argument('-o', '--output', help='output file name', required=True)
    args = parser.parse_args()
    main(args)