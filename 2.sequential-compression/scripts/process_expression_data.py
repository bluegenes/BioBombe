"""
extended/modified by N.T. Pierce from biobombe scripts https://github.com/greenelab/BioBombe
and https://github.com/greenelab/nf1_inactivation/blob/master/scripts/process_rnaseq.py
"""

import os
import sys
import random
import argparse
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import preprocessing

random.seed(42)

#wget https://osf.io/ek9nu/download -O haptophyta_orthogroup.quant.tsv

def read_counts(countfile):
    if '.tsv' in countfile or '.csv' in countfile:
        separator = '\t'
        if '.csv' in countfile:
            separator = ','
        try:
            counts = pd.read_csv(countfile, sep=separator, index_col=0, compression='infer')
        except Exception as e:
            sys.stderr.write(f"\n\tError: {countfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in countfile:
        try:
            counts = pd.read_excel(countfile, index_col=0, compression='infer')
        except Exception as e:
            sys.stderr.write(f"\n\tError: {countfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    return counts

def preprocess_data(countfile, out_file=None, scale = False, scale_method = "min_max", percTest = 0.1, mad= False, num_mad_genes = 8000, out_folder = ""):
    """
    Zero-one scale the expression data

    Output:
    Writes normalized matrix (if output=True) and mad genes to file
    (if mad=True); returns the normalized matrix if output=False
    """

    #################
    # Preprocessing #
    #################

    # load data to pd dataframe
    #expr_data = pd.read_csv(countfile, dtype='str', sep='\t', index_col=0)
    expr_data = read_counts(countfile)

    # Drop all row names with unidentified gene name
    #expr_data = data[-data.index.str.contains('?', regex=False)]

    expr_data = expr_data.sort_index() # sort by gene name

    if scale:
        if scale_method == "min_max":
           # Zero-one normalize
            min_max_scaler = preprocessing.MinMaxScaler()
            expr_scaled = min_max_scaler.fit_transform(expr_data.T)
            file_suffix = '.processed.zeroone.tsv.gz'
        elif scale_method == 'zscore':
            expr_scaled = preprocessing.scale(expr_data.T, axis=0)
            file_suffix = '.processed.zscore.tsv.gz'

        expr_norm = pd.DataFrame(expr_scaled, index=expr_data.columns, columns=expr_data.index) #.T  # transform back
    else:
        expr_norm = expr_data.T # transform
        file_suffix = '.processed.tsv.gz'

    # check that output folder exists, make it if not
    if out_folder:
        if not os.path.exists(out_folder):
            os.mkdirs(out_folder)

    # write scaled output

    if not out_file:
        out_file = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + file_suffix)

    expr_norm.to_csv(out_file, sep='\t', header=True, index=True, compression='gzip', float_format='%.3g')
    #print(min_max_scaler.data_max_)

    ############################
    # Split Test, Training sets
    ############################

    trainDF, testDF = train_test_split(expr_norm, test_size=percTest, random_state =42) #, shuffle=False

    # write training set to file
    train_file = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.train' + file_suffix)
    trainDF.to_csv(train_file, sep='\t', compression='gzip', float_format='%.3g')

    # write test set to file
    test_file = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.test' + file_suffix )
    testDF.to_csv(test_file, sep='\t',  header=True, index=True, compression='gzip', float_format='%.3g')

    ###########################################
    # Sort on Median Absolute Deviation "mad" #
    ###########################################

    if mad:

        mad_file = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.mad' + file_suffix) #+ str(num_mad_genes) + '.tsv.gz')
        mad_train = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.mad.train90' + file_suffix) #+ str(num_mad_genes) + '.tsv.gz')
        mad_test = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.mad.test10' + file_suffix) #+ str(num_mad_genes) + '.tsv.gz')

        mad_genes = expr_norm.mad(axis=0).sort_values(ascending=False)
        top_mad_genes = mad_genes.iloc[0:int(num_mad_genes), ].index
        mad_genes_df =expr_norm.loc[:, top_mad_genes] ##rnaseq_df.loc[:, top_mad_genes]
        #mad_genes_df.columns = ['gene_id', 'median_absolute_deviation']
        # write
        #mad_genes_df.to_csv(mad_file, sep='\t', index=False, compression='gzip', float_format='%.3g')

        # write full df
        mad_genes_df.to_csv(mad_file, sep='\t', header=True, index=True, compression='gzip', float_format='%.3g')

        #split mad --> test, training
        madtrainDF, madtestDF = train_test_split(mad_genes_df, test_size=percTest, random_state =42) #, shuffle=False
        madtrainDF.to_csv(mad_train, sep='\t', header=True, index=True, compression='gzip', float_format='%.3g')
        madtestDF.to_csv(mad_test, sep='\t',  header=True, index=True,compression='gzip', float_format='%.3g')


if __name__ == '__main__':
    #countfile = "haptophyta_orthogroup.quant.tsv"

    p = argparse.ArgumentParser()
    p.add_argument("countfile", help = "csv,tsv,or xls rnaseq gene expression table")
    p.add_argument("--scale", help="true/false: scale the expression data?", action="store_true", default=False)
    p.add_argument("--scale_method",help="method to scale expression matrix", default='min_max')
    p.add_argument("--percTest", help="percent of data to set aside as test set. Training set will be 1-percTest", default = 0.1, type=float)
    p.add_argument("--mad", help="subset mad genes", action="store_true", default=False)
    p.add_argument("--output_folder", default="")
    p.add_argument("--output_filename", default=None)
    p.add_argument("--num_mad_genes", help="number of highest median absolute deviation genes to output", type=int, default=8000)
    p.add_argument("--transform", help="transform dataset prior to processing", action="store_true")
    args = p.parse_args()
    sys.exit(preprocess_data(args.countfile, args.output_filename, args.scale, args.scale_method, args.percTest, args.mad, args.num_mad_genes, args.output_folder))
