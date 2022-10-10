import glob
import random

import json

with open('pdbs.jsonl', 'r') as json_file:
    json_list = list(json_file)

my_dict = {}
for json_str in json_list:
    result = json.loads(json_str)
    all_chain_list = [item[-1:] for item in list(result) if item[:9]=='seq_chain'] #['A','B', 'C',...]
    if '2mer' in result['name']:
        masked_chain_list = ['A','B']
    elif '3mer' in result['name']:
        masked_chain_list = ['A','B','C'] 
    elif '4mer' in result['name']:
        masked_chain_list = ['A','B','C','D']
    elif '5mer' in result['name']:
        masked_chain_list = ['A','B','C','D','E']
    elif '6mer' in result['name']:
        masked_chain_list = ['A','B','C','D','E','F']
    elif '7mer' in result['name']:
        masked_chain_list = ['A','B','C','D','E','F','G']
    elif '8mer' in result['name']:
        masked_chain_list = ['A','B','C','D','E','F','G','H']
    visible_chain_list = [] 
    my_dict[result['name']]= (masked_chain_list, visible_chain_list)

with open('pdbs_masked.jsonl', 'w') as f:
    f.write(json.dumps(my_dict) + '\n')


print('Finished')
# Output looks like this:
# {"5TTA": [["A"], ["B"]], "3LIS": [["A"], ["B"]]}
