#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Python script to generate Knowledge Graph csv file.
This script uses the GNU AGPLv3 license.
"""

# Import libraries
import pandas as pd
###############################################################################

# Author information
__author__ = "Kevin McCoy"
__copyright__ = "Copyright 2022, Laboratory for Pathology Dynamics"
__credits__ = ["Kevin McCoy", "David Kartchner"]
__license__ = "GNU AGPLv3"
__version__ = "0.1.0"
__maintainer__ = "Cassie Mitchell"
__email__ = "cassie.mitchell@bme.gatech.edu"
__status__ = "production"
__date__ = "2022-03-02"

###############################################################################


def load_input(input_fn):
    """Input raw csv data into pandas DataFrame."""

    return pd.read_csv(input_fn, header=None, usecols=range(12), encoding="ISO-8859-1")


def filter_input(df):
    """Filter data in Pandas DataFrame."""

    # Rename columns
    df.columns = ["PREDICATION_ID", "SENTENCE_ID", "PMID", "PREDICATE", 
                    "SUBJECT_CUI", "SUBJECT_NAME", "SUBJECT_SEMTYPE", 
                    "SUBJECT_NOVELTY", "OBJECT_CUI", "OBJECT_NAME", 
                    "OBJECT_SEMTYPE", "OBJECT_NOVELTY"]

    # Keep only relevant columns
    df = df[["PMID", "SUBJECT_CUI", "SUBJECT_SEMTYPE", "OBJECT_CUI", "OBJECT_SEMTYPE", 
                "PREDICATE"]]

    # Rename columns
    df.columns = ["pmid", "start_node", "start_type", "end_node", "end_type", 
                    "relation"]

    # Ignore duplicates in same paper
    df = df.drop_duplicates()

    # Keep only relevant columns
    df = df[["start_node", "start_type", "end_node", "end_type", "relation"]]

    return df


def make_generics_list(generic_concepts_fn):
    """Import list of generic concepts."""

    df = pd.read_csv(generic_concepts_fn, usecols=[1], header=None)

    return df[1].values.tolist()


def clean_input(df, generics_list):
    """Remove generic nodes and add weight column."""

    # Remove generic nodes
    df = df[(~df.start_node.isin(generics_list)) & (~df.end_node.isin(generics_list))]

    # Add weight of each predication
    df = df.groupby(list(df.columns)).size().reset_index().rename(columns={0:'weight'})

    df = df.reset_index()

    return df


def output(df, output_fn):
    """Output cleaned DataFrame to new csv."""

    df.to_csv(output_fn, index=False)


###############################################################################


if __name__ == "__main__":

    ## Change these variables
    # Filename of predications csv
    input_predications_fn = ""
    # Filename of generic concepts csv
    generic_concepts_fn = ""
    # Desired output csv filename
    output_fn = ""

    # Run script
    df = load_input(input_predications_fn)
    df = filter_input(df)
    generics_list = make_generics_list(generic_concepts_fn)
    df = clean_input(df, generics_list)
    output(df, output_fn)