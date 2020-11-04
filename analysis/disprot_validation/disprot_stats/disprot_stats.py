"""Plot various statistics associated with the entries in the DisProt database."""

import json
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import collections, functools, operator

# Use JSON module to load file
with open('../../../data/DisProt/JSON/DisProt.json') as file:
    data = json.load(file)['data']

def char_counter(string):
    dict = {}
    for char in string:
        if char in dict.keys():
            dict[char] = dict.get(char) + 1
        else:
            dict[char] = 1
    return dict

def difference_dict(tot_char, dis_char):
    or_char = {}
    for k in tot_char.keys():
        if k in dis_char.keys():
            or_char[k] = tot_char.get(k) - dis_char.get(k)
        else:
            or_char[k] = tot_char.get(k)
    return or_char

# Extract relevant fields from DisProt
rows = []
for record in data:
    row = {}
    for field in ['disprot_id', 'name', 'organism', 'length', 'released', 'disorder_content', 'sequence']:
        row[field] = record[field]
    row['taxonomy'] = record['taxonomy'][0]
    num_regions = 0
    num_residues = 0
    disordered_residues = []
    counts = {}
    seq = record['sequence']
    for region in record['disprot_consensus']['structural_state']:
        if region['type'] == 'D':
            placeholder_residues = []
            num_regions += 1
            num_residues += region['end'] - region['start'] + 1  # Add 1 since endpoints included
            for residue_index in range(region['start']-1,region['end']):
                placeholder_residues.append(record['sequence'][residue_index])
            disordered_residues.append(placeholder_residues)
    for region in disordered_residues:
        for residue in region:
            if residue in counts:
                counts[residue]+= 1
            else:
                counts[residue] = 1
    row['num_regions'] = num_regions
    row['num_residues'] = num_residues
    row['disordered_regions'] = disordered_residues
    row['disordered_counts'] = counts
    row['ordered_counts'] = difference_dict(char_counter(seq), counts)
    row['total_counts'] = char_counter(seq)
    rows.append(row)
df = pd.DataFrame(rows)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Distribution of entries across taxonomic groups
counts = df['taxonomy'].value_counts()
plt.bar(counts.index, counts.values)
plt.xlabel('Domain')
plt.ylabel('Number of DisProt entries')
plt.savefig('out/bar_taxa.png')
plt.close()

#Distribution of number of residues and protein length
with_outliers_residues = plt.hist(df.num_residues, ec='white')
plt.xlabel("Number of Residues")
plt.ylabel("Number of Disprot Entries")
plt.title("Number of Residues for Disprot Entries")
plt.savefig('out/len_res.png')
plt.close()

#Distribution of Number of Residues without Outliers
outliers = df[df['num_residues']>1500]
plt.hist(df.num_residues,range = [0,1500], ec='white')
plt.xlabel("Number of Residues")
plt.ylabel("Number of Disprot Entries")
plt.title("Number of Residues for Disprot Entries")
plt.savefig('out/len_out.png')
plt.close()

#Distribution of Disprot Entry Lengths
plt.hist(df.length, range = [0,4500], ec='white')
outliers = df[df['length']>4500]
plt.xlabel("Length")
plt.ylabel("Number of Disprot Entries")
plt.title("Length of Disprot Entries")
plt.savefig('out/length.png')
plt.close()

#Distribution of Disorder Content
plt.hist(df.disorder_content, ec='white')
plt.ylabel('Number of DisProt entries')
plt.xlabel("Disorder Content")
plt.title('Distribution of Disorder Content')
plt.savefig('out/dis_dist.png')
plt.close()

#Distribution of Length of Disordered Regions
#Disordered Length = Disordered Content * Total Length
df["length_of_disordered_regions"] = df.disorder_content * df.length
a = plt.hist(df.length_of_disordered_regions, ec='white')
plt.ylabel('Number of DisProt entries')
plt.xlabel("Length of Disordered Regions")
plt.title('Distribution of Length of Disordered Regions')
plt.savefig('out/dis_hist.png')
plt.close()

#Disordered Region Length Vs. Overall Length
#Is there an association between protein length and disordered content?
outliers_removed = df[df["length"] < 5000]
outliers_removed
colors = {'Eukaryota':'red', 'Bacteria':'green', 'Viruses':'blue', 'Archea':'yellow'}
outliers_removed.plot.scatter(x='length', y = 'length_of_disordered_regions')
# try 1) white lines 2) slightly transparent?
plt.xlabel("Length")
plt.ylabel("Total Length of Disordered Regions")
plt.title('Disordered Region Length vs Total Length')
plt.savefig('out/len_comp.png')
plt.close()

#Distribution of the number of Disordered Regions
plt.hist(df.num_regions, ec='white')
plt.xlabel("Number of Regions")
plt.ylabel('Number of DisProt entries')
plt.title('Distribution of the Number of Regions')
plt.savefig('out/num_reg.png')
plt.close()

# Violinplot of length by taxonomic origin
df_less_outliers = df[df['num_residues'] < 1000]
sns.violinplot(x='taxonomy', y='num_residues', data=df_less_outliers, palette='Set2')
plt.xlabel('Domain')
plt.ylabel('Number of Residues')
plt.title('Distribution of the Number of Residues by Taxonomic Origin')
plt.savefig('out/violin_taxa.png')
plt.close()

# get counts of amino acids totals
amino_total_dicts = df['total_counts'].to_list()
amino_acids_total = dict(functools.reduce(operator.add,
         map(collections.Counter, amino_total_dicts)))
amino_acids_total.pop('X')
amino_acids_total.pop('U')

# get counts of amino acids ordered
amino_ord_dicts = df['ordered_counts'].to_list()
amino_acids_ord = dict(functools.reduce(operator.add,
         map(collections.Counter, amino_ord_dicts)))
amino_acids_ord.pop('X')
amino_acids_ord.pop('U')

# get counts of amino acids disordered
amino_dis_dicts = df['disordered_counts'].to_list()
amino_acids_dis = dict(functools.reduce(operator.add,
         map(collections.Counter, amino_dis_dicts)))

# create data frame of amino acid data
amino_acid_lst = [amino_acids_total, amino_acids_ord, amino_acids_dis]
amino_acid_df = pd.DataFrame(amino_acid_lst)
amino_acid_df.rename(index={0:'total', 1:'ordered', 2:'disordered'}, inplace=True)

# Amino acid distribution over total
plt.bar(list(amino_acid_df.columns), list(amino_acid_df.loc['total',:]))
plt.xlabel('Amino Acids')
plt.ylabel('Number of Residues')
plt.title('Number of Residues per Amino Acid in All Regions')
plt.savefig('out/amino_total.png')
plt.close()

# Amino acid distribution over disordered
plt.bar(list(amino_acid_df.columns), list(amino_acid_df.loc['disordered',:]))
plt.xlabel('Amino Acids')
plt.ylabel('Number of Residues')
plt.title('Number of Residues per Amino Acid in Disordered Regions')
plt.savefig('out/amino_dis.png')
plt.close()

# Amino acid distribution over ordered
plt.bar(list(amino_acid_df.columns), list(amino_acid_df.loc['ordered',:]))
plt.xlabel('Amino Acids')
plt.ylabel('Number of Residues')
plt.title('Number of Residues per Amino Acid in Ordered Regions')
plt.savefig('out/amino_ord.png')
plt.close()

# Amino acid distribution by disorder
p1 = plt.bar(list(amino_acid_df.columns), list(amino_acid_df.loc['disordered',:]))
p2 = plt.bar(list(amino_acid_df.columns), list(amino_acid_df.loc['ordered',:]),
             bottom=list(amino_acid_df.loc['disordered',:]))
plt.xlabel('Amino Acids')
plt.ylabel('Number of Residues')
plt.title('Number of Residues per Amino Acid by Disorder')
plt.legend(['disordered', 'ordered'])
plt.savefig('out/amino_stack.png')
plt.close()

# Amino acid enrichment
total_aa_idr = amino_acid_df.apply(sum, axis=1)['disordered']
total_aa = amino_acid_df.apply(sum, axis=1)['total']
amino_enrich = (((amino_acid_df.loc['disordered', :] / total_aa_idr) * 10000) / ((amino_acid_df.loc['total', :] / total_aa) * 100)) - 100
plt.bar(list(amino_acid_df.columns), list(amino_enrich))
plt.xlabel('Amino Acids')
plt.ylabel('Enrichment')
plt.title('Amino Acid Enrichment in Disordered Regions')
plt.savefig('out/amino_enrich.png')
plt.close()

"""
DEPENDENCIES
../../../data/DisProt/JSON/DisProt.json
"""