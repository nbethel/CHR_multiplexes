import glob
import random
import numpy as np
import json
import itertools

#MODIFY this path
with open('pdbs.jsonl', 'r') as json_file:
    json_list = list(json_file)

my_dict = {}
for json_str in json_list:
    result = json.loads(json_str)
    all_chain_list = sorted([item[-1:] for item in list(result) if item[:9]=='seq_chain']) #A, B, C, ...
    tied_positions_list = []
    if True:#result['name'] == '3LIS':
        fin=open('pos/'+result['name']+'.pos')
        lline=fin.readline()
        fin.close()
        spl=lline.split(',')[:-1]
        spl2=result['name'].split('_')
        repL=float(spl2[0][1:])+float(spl2[1][1:])+float(spl2[2][1:])+float(spl2[3][1:])
        Tspl=spl
        for i in spl:
            T2spl=[]
            if i in Tspl:
                dupes=[int(i[:-1])]
                for j in Tspl:
                       
                    if int(j[:-1])%repL==int(i[:-1]) and j!=i:
                        dupes.append(int(j[:-1]))
                    elif j!=i:
                        T2spl.append(j)
                if True:
                        temp_dict = {}
                        temp_dict[all_chain_list[0]] = dupes
                        temp_dict[all_chain_list[1]] = dupes
                        if '3mer' in result['name']:
                            temp_dict[all_chain_list[2]] = dupes
                        elif '4mer' in result['name']:
                            temp_dict[all_chain_list[2]] = dupes
                            temp_dict[all_chain_list[3]] = dupes
                        elif '5mer' in result['name']:
                            temp_dict[all_chain_list[2]] = dupes
                            temp_dict[all_chain_list[3]] = dupes
                            temp_dict[all_chain_list[4]] = dupes
                        elif '6mer' in result['name']:
                            temp_dict[all_chain_list[2]] = dupes
                            temp_dict[all_chain_list[3]] = dupes
                            temp_dict[all_chain_list[4]] = dupes
                            temp_dict[all_chain_list[5]] = dupes
                        elif '7mer' in result['name']:
                            temp_dict[all_chain_list[2]] = dupes
                            temp_dict[all_chain_list[3]] = dupes
                            temp_dict[all_chain_list[4]] = dupes
                            temp_dict[all_chain_list[5]] = dupes
                            temp_dict[all_chain_list[6]] = dupes
                        elif '8mer' in result['name']:
                            temp_dict[all_chain_list[2]] = dupes
                            temp_dict[all_chain_list[3]] = dupes
                            temp_dict[all_chain_list[4]] = dupes
                            temp_dict[all_chain_list[5]] = dupes
                            temp_dict[all_chain_list[6]] = dupes
                            temp_dict[all_chain_list[7]] = dupes
                        tied_positions_list.append(temp_dict)
            Tspl=T2spl
        chain_length = len(result["seq_chain_A"])
        splNum=[]
        for i in spl:
            splNum.append(int(i[:-1]))
        for i in range(1,chain_length+1):
          if i not in splNum:
            temp_dict = {}
            temp_dict[all_chain_list[0]] = [i]  #all_chain_list[0] == "A" in this case
            temp_dict[all_chain_list[1]] = [i]  #all_chain_list[0] == "B"
            if '3mer' in result['name']:
                temp_dict[all_chain_list[2]] = [i]  #all_chain_list[0] == "B"
            elif '4mer' in result['name']:
                temp_dict[all_chain_list[2]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[3]] = [i]  #all_chain_list[0] == "B"
            elif '5mer' in result['name']:
                temp_dict[all_chain_list[2]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[3]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[4]] = [i]  #all_chain_list[0] == "B"
            elif '6mer' in result['name']:
                temp_dict[all_chain_list[2]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[3]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[4]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[5]] = [i]  #all_chain_list[0] == "B"
            elif '7mer' in result['name']:
                temp_dict[all_chain_list[2]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[3]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[4]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[5]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[6]] = [i]  #all_chain_list[0] == "B"
            elif '8mer' in result['name']:
                temp_dict[all_chain_list[2]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[3]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[4]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[5]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[6]] = [i]  #all_chain_list[0] == "B"
                temp_dict[all_chain_list[7]] = [i]
            tied_positions_list.append(temp_dict)
    else:
        tied_positions_list = []
    my_dict[result['name']] = tied_positions_list

#Write output to:    
with open('pdbs_tied.jsonl', 'w') as f:
    f.write(json.dumps(my_dict) + '\n')


print('Finished')
