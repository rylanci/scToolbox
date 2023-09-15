import pandas as pd 
import os
import re
import torch
from sklearn.model_selection import train_test_split
from rfm import LaplaceRFM
import matplotlib.pyplot as plt
import csv

### This python script will be pointed towards a dir of matrices and metadata, choose files based off celltypes and comparisons,
### then run the rfm model on each pair. Lastly output will be written to file.

idir = "/projects/ps-epigen/users/rlan/Liver/RFM_Liver/cellclass_out/matrices/" 
odir = "/projects/ps-epigen/users/rlan/Liver/RFM_Liver/cellclass_out/model/" 

files = os.listdir(idir)
#print(files[1])

# extract first term from all files and take the unique terms. this will be our celltypes
celltypes = [f.split('_')[0] for f in files]
#print(celltypes)

# Extract Conditions 
compTable = pd.read_table("/projects/ps-epigen/users/rlan/Liver/manuscript_RNA/DESeq/comparisons.txt")
compList = list()
for index, row in compTable.iterrows():
    comp = row["Target"] + "_" + row["Reference"] 
    comp = comp.replace(" ", ".")
    compList.append(comp)

#print(compList)


# Extract matrix/metadata pairs from path an run model  
for c in celltypes:
    print(c)
    # Iterate over celltypes
    for j in compList:
        print(j)
        # Iterate over condition pairs 
        
        # compile regex terms and use to filter our list of files
        vmat = re.compile(f"{c}_{j}_rfm_mat.csv")
        vmeta = re.compile(f"{c}_{j}_rfm_meta.csv")
        matrix_name = list(filter(vmat.match, files))
        meta_name = list(filter(vmeta.match, files))
        matrix_name = ' '.join(matrix_name)
        meta_name = ' '.join(meta_name)

        #print(matrix_name)
        #print(meta_name)
        mat = pd.read_csv(f'{idir}{matrix_name}',index_col=0)
        meta = pd.read_csv(f'{idir}{meta_name}',index_col=0)


        X = torch.from_numpy(mat.T.to_numpy()).float() # .cuda()
        Y = torch.from_numpy(meta.conditions.map({j.split("_")[0]: 1, j.split("_")[1]: -1}).to_numpy()).float()[:,None]
        X_train, X_val, Y_train, Y_val = train_test_split(X, Y, test_size=0.1)

        print("Training Model")
        model = LaplaceRFM(bandwidth=10.)
        model.fit(
            (X_train, Y_train),  
            (X_val, Y_val), 
            loader=False, 
            iters=5,
            classif=True
        )

        plt.plot(model.M.diag()/model.M.diag().max())
        plt.xlabel('gene id')
        plt.title('gene importance')
        plt.savefig(f"{odir}{c}_{j}_rfm_plt.pdf")
        plt.close()

        topn = list(mat.index[model.M.diag().topk(100).indices])
        with open(f"{odir}{c}_{j}_rfm_top100.tsv", "w") as f:
            writer = csv.writer(f)
            writer.writerow(topn)



print("done")

