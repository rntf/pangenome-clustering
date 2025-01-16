# This program reads a sequence of feature values 
# associated with corresponding bases and 
# performs clustering based on the mean 
# calculated using a sliding window.

import matplotlib.pyplot as plt
import numpy as np
import gc
import matplotlib
matplotlib.use('Agg')


# plot and text result
# output folders are ris and graph
def salva(i,nome):
    for key in list_of_n:
        with open(f'ris/ris_{nome}_{i}_{key}.txt', 'w') as fout:
            for k in sorted(clusters_per_l[key]):
                somma = sum(list(composition_per_l[key][k].values()))
                fout.write(
                    f'{k}\t{clusters_per_l[key][k]}\t{composition_per_l[key][k]["A"] / somma * 100:.2f}\t{composition_per_l[key][k]["C"] / somma * 100:.2f}\t{composition_per_l[key][k]["G"] / somma * 100:.2f}\t{composition_per_l[key][k]["T"] / somma * 100:.2f}\n')
        fig, ax = plt.subplots(figsize=(9, 6))

        ax.scatter(np.array(list(clusters_per_l[key])), np.array(list(clusters_per_l[key].values())))
        plt.gca().update(dict(xlabel='length', ylabel='P'))
        ax.set_yscale("log")
        ax.set_xscale("log")
        plt.savefig(f"graph/probabilita_{nome}_{i}_{key}.png")
        plt.close()
        fig, ax = plt.subplots(figsize=(9, 6))

        ax.plot(np.array(sorted(composition_per_l[key])[:MAX_FREQUENZE]),
                np.array([composition_per_l[key][l]["A"] / sum([composition_per_l[key][l][k] for k in composition_per_l[key][l]]) for l in
                          sorted(composition_per_l[key])[:MAX_FREQUENZE]]),
                label='A')
        ax.plot(np.array(sorted(composition_per_l[key])[:MAX_FREQUENZE]),
                np.array([composition_per_l[key][l]["C"] / sum([composition_per_l[key][l][k] for k in composition_per_l[key][l]]) for l in
                          sorted(composition_per_l[key])[:MAX_FREQUENZE]]),
                label='C')
        ax.plot(np.array(sorted(composition_per_l[key])[:MAX_FREQUENZE]),
                np.array([composition_per_l[key][l]["G"] / sum([composition_per_l[key][l][k] for k in composition_per_l[key][l]]) for l in
                          sorted(composition_per_l[key])[:MAX_FREQUENZE]]),
                label='G')
        ax.plot(np.array(sorted(composition_per_l[key])[:MAX_FREQUENZE]),
                np.array([composition_per_l[key][l]["T"] / sum([composition_per_l[key][l][k] for k in composition_per_l[key][l]]) for l in
                          sorted(composition_per_l[key])[:MAX_FREQUENZE]]),
                label='T')
        ax.legend()
        plt.gca().update(dict(xlabel='length', ylabel='P'))
        plt.savefig(f"graph/fequenze_{nome}_{i}_{key}.png")
        plt.cla()
        plt.clf()
        plt.close('all')
        # del ax
        # del fig
        # gc.collect()


# output file name generator
def nome(file):
    return file.split('.')[0].split('_',maxsplit=1)[1]


#dati tutti i file step in una cartella eseguo il clustering
INIZIO=1 #lower value of n
FINE=10 
PASSO=1 #step between n values
MAX_FREQUENZE=40
STOP=10**6
import os
path='input_folder'
files = os.listdir(path)
files_file = [f for f in files if os.path.isfile(os.path.join(path, f))]
for ff in files_file:
    # n: length of the moving window for calculating the average
    list_of_n=list(range(INIZIO+1,FINE,PASSO))+list(range(INIZIO*10,FINE*10,PASSO*10))+list(range(INIZIO*100,FINE*100+1,PASSO*100))
    # moving window for calculating the average 
    window=[0]*(list_of_n[-1]+1)
    # l: length of the clusters
    clusters_per_l= {}
    composition_per_l={}
    current_composition={}
    for key in list_of_n:
        clusters_per_l[key]={}
        composition_per_l[key]={}
        current_composition[key]={'A':0,'C':0,'G':0,'T':0}
    with open(path+'/'+ff) as fin:
        # reading of the first values in orter to fill in the buffers
        for i in range(list_of_n[-1]+1):
            row=fin.readline()
            row=row.split()
            window[i]=int(row[1])
        # the intersections correspond to a change of sign in the difference
        new_difference_per_n={}
        old_difference_per_n={}
        # starting difference between the average and the feature values
        for key in list_of_n:
            # if n is even the first and the last values of the windows are weighed half
            if key%2==0:
                new_difference_per_n[key]=(sum(window[-key:-1])+window[-key-1]/2+window[-1]/2)/key-window[-key//2-1]
            else:
                new_difference_per_n[key]=sum(window[-key:])/key-window[-key//2]
            old_difference_per_n[key]=new_difference_per_n[key]
        # main loop on the input data 
        for row in fin:
            #intermediate results saves
            stop=10000
            while i>stop:
                stop*=10
            if stop==i:
                salva(i,nome(ff))
            i+=1
            # each row in the input file has the base simbol, a white space, the value of the feature value
            # row[0]: base; row[1]: feature function value
            row=row.split()
            row[1]=int(row[1])
            window.pop(0)
            window.append(row[1])
            current_l_per_n={}
            # computation of the new difference for all the n values
            for key in list_of_n:
                if new_difference_per_n[key]:
                    old_difference_per_n[key]=new_difference_per_n[key]
                if key % 2 == 0:
                    new_difference_per_n[key]=(sum(window[-key:-1])+window[-key-1]/2+window[-1]/2)/key-window[-key//2-1]
                else:
                    new_difference_per_n[key]=sum(window[-key:])/key-window[-key//2]
                # l is computed from the composition
                current_l_per_n[key] = sum(list(current_composition[key].values()))
                # intersection
                if old_difference_per_n[key]*new_difference_per_n[key]<0 and current_l_per_n[key]:
                    if current_l_per_n[key] not in clusters_per_l[key]:
                        clusters_per_l[key][current_l_per_n[key]]=1
                        composition_per_l[key][current_l_per_n[key]]=current_composition[key]
                    else:
                        clusters_per_l[key][current_l_per_n[key]]+=1
                        for key2 in current_composition[key]:
                            composition_per_l[key][current_l_per_n[key]][key2]+=current_composition[key][key2]
                    current_composition[key] = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
                # update of the composition
                current_composition[key][row[0]]+=1
    salva(i,nome(ff))
