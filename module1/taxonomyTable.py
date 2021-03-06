#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script creates the taxonomy table needed for loading the data into the 
phyloseq package.

The taxonomy table links every taxID to its taxonomy

@author: Adriel Latorre-Pérez
@company: Darwin Bioprospecting Excellence S.L.
@date: 19/10/2020
"""

import sys, os
from argparse import ArgumentParser # Para gestionar argumentos

def Arguments():
    """For input folder.
    """
    parser = ArgumentParser (description ="This script creates the taxonomy \
                            table needed for loading the data into the \
                            phyloseq package.\n \
                            The taxonomy table links every taxID to its \
                            taxonomy")
    parser.add_argument ('-i', '--input', dest='otu_table',
                           action ='store', required =True ,
                           help='Path to (pseudo)OTU table')
    
    parser.add_argument ('-t', '--taxonomy', dest='taxonomy_file',
                           action ='store', required =True ,
                           help="Path to database's taxonomy table")
    # Process the arguments
    try:
        args = parser.parse_args ()
        return args
    except:
        print('\nPlease, include the required arguments.')
        sys.exit()
    # end try

    
def taxnomyDIC(taxonomy_file):
    """
    It takes the database's taxonomy file and it creates a dictionary

    Parameters
    ----------
    taxonomy_file : PATH TO FILE
        PATH to databse's taxonomy file.

    Returns
    -------
    Taxonomy dictionary.

    """
    file = open (taxonomy_file)
    taxDIC = {}
    
    for line in file:
        taxID = line.strip().split("\t")[0]
        taxonomy = line.strip().split("\t")[1]
        taxonomy = ",".join(taxonomy.split("; "))
        taxDIC[taxID] = taxonomy
    
    file.close ()
    
    return taxDIC


def taxID2taxonomy(otu_table, taxDIC):
    """
    Reads the taxIDs from the (pseudo)OTU table and it writes the taxonomy
    of each ID.
    
    Parameters
    ----------
    otu_table : PATH TO FILE
        PATH to the (pseudo)OTU table generated by mergePAF.py.
    taxDIC : DICTIONARY
        Output of taxonomyDIC().

    Returns
    -------
    None. Writes output to stdout.

    """
    file = open (otu_table)
    
    print ("#OTU ID,Domain,Phylum,Class,Order,Family,Genus,Species")
    for line in file:
        if line[0] == "#":
            continue
        taxID = line.strip().split(",")[0]
        print (taxID + "," + taxDIC[taxID])
    
    file.close ()

if __name__ == "__main__":
   args = Arguments()
   taxDIC = taxnomyDIC(args.taxonomy_file)
   taxID2taxonomy(args.otu_table, taxDIC)
