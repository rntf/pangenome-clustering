# This program reads UTF-encoded files containing a list of bases 
# from a genome sequence and generates a sequence of feature values 
# using the RY rule. 
# Input values that do not correspond to a valid base are ignored. 
# Each line in the output files contains 
# the base, a whitespace, and the corresponding feature value. 

import os
path='input_folder'
files = os.listdir(path)
files_file = [f for f in files if os.path.isfile(os.path.join(path, f))]
# the output file are saved in folder ch
for i in files_file:
    with open(path+'/'+i) as fin, open(f'ch/step_{i}','w') as fout:
        tot=0
        while c:=fin.read(1):
            if c.upper() in 'AG':
                tot+=1
                fout.write(f'{c.upper()} {tot}\n')
            elif c.upper() in 'CT':
                tot-=1
                fout.write(f'{c.upper()} {tot}\n')
