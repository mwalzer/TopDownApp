import numpy as np
import pandas as pd
from datetime import datetime

tmp_file_loc = {
    'pro_single': 'TopPIC_tutorial/st_2_ms2_toppic_proteoform_single.tsv',
    'pro': 'TopPIC_tutorial/st_2_ms2_toppic_proteoform.tsv',
    'prsm_single': 'TopPIC_tutorial/st_2_ms2_toppic_prsm_single.tsv',
    'prsm': 'TopPIC_tutorial/st_2_ms2_toppic_prsm.tsv'
}
tmp_filename = "st_2.mzML"

def peek_toppic_res(file: str) -> int:
    with open(file, 'r') as f:
        lines = f.readlines()
    sl = 0
    for i,l in enumerate(lines):
        if '*** Parameters ***' in l:
            sl = i+1  # adjust for 0-start
    return sl+1  # add one more for a blank line

dfs_skiprows = {k: peek_toppic_res(v) for k,v in tmp_file_loc.items() }
l = min(dfs_skiprows.values())
dfs = {k: pd.read_csv(v, sep='\t', skiprows=l )  for k,v in tmp_file_loc.items() }

header = """
COM This is the mzTab to a minimal "Summary Top-Down Proteoform Identification Report"  
MTD mzTab-version   1.0.0 
MTD mzTab-mode  Summary 
MTD mzTab-type  Identification 
MTD title   {filename}
MTD description TopDown Proteomics experiment, note that PSMs are PrSMs 
MTD assay[1]-sample_ref sample[1] 
MTD assay[1]-ms_run_ref ms_run[1] 
MTD ms_run[1]-location  {filename} 
MTD ms_run[1]-format    [MS, MS:1000584, mzML file, ] 
MTD ms_run[1]-id_format [MS, MS:1000776 scan number only nativeID format, ] 
MTD software[1] [MS, MS:1003145, ThermoRawFileParser, ] 
MTD software[2] [MS, 1002714, TOPP FLASHDeconv, ]
MTD software[3] [MS, MS:1002901, TopPIC, ]
MTD psm_search_engine_score[1]  [MS, MS:1002928, TopPIC:spectral E-value, ]
MTD psm_search_engine_score[2]  [MS, MS:1002932, TopPIC:MIScore, ]
MTD protein_search_engine_score[1]  [MS:1002906, search engine specific score for proteoforms,]
"""

# TODO
# MTD software[1-n]-setting[1-n] ??? 
# MTD fixed_mod[1] [UNIMOD, UNIMOD:4, Carbamidomethyl, ] ..?
# MTD variable_mod[1] [UNIMOD, UNIMOD:35, Oxidation, ] ..?
# do we need full modification-choice flexibility in our app? default var M16 fix - 
# protein_search_engine_score is E-value too, but no CV?!
# MTD sample_processing[1-n] ??? 

#link psm-protein sub-tables how???
import re
def simpleparse_ncbi_tax(header: str) -> str:
    m = re.match(".*OX=(\d*) .*",header)
    if m:
        return m[1]
    else:
        return 'null'

def simpleparse_species(header: str) -> str:
    m = re.match(".*OS=([\w\s]*) .*",header)
    if m:
        return m[1]
    else:
        return 'null'

tmp_db_path = 'TopPIC_tutorial/TopPIC_tutorial_uniprot-st.fasta'
from Bio.SeqIO.FastaIO import SimpleFastaParser
with open(tmp_db_path, 'r') as h:
    pl = list(SimpleFastaParser(h))
    db_len = len(pl)
    db_tax = ','.join(list({simpleparse_ncbi_tax(e[0]) for e in pl}))
    db_species = ','.join(list({simpleparse_species(e[0]) for e in pl}))

rigged_cols_pr = { 
    'ambiguity_members': 'null',	
    'taxid': db_tax,
    'species': db_species,
    'database': 'UniProtKB',
    'database_version': datetime.today().strftime('%Y-%m') + ' (' + str(db_len) + ')',
    'search_engine': '[MS, MS:1002901, TopPIC, ]',
    'PRH': 'PRT'
}

def strip_seq(s:str) -> str: 
    rs = s.split('.')
    if len(rs) > 0:
        rs = '.'.join(rs[1:-1])  
    else:
        rs = rs[0]
    return rs

col_map_pr = {
    'Retention time':'opt_global_RT', 
    '#peaks':'opt_proteoform_peak_number', 
    'Charge':'opt_proteoform_charge', 
    'Protein accession': 'accession', 
    'Protein description': 'description', 
    '#unexpected modifications': 'opt_proteoform_unexpected_modifications', 
    '#variable PTMs': 'modifications', 
    '#matched peaks': 'opt_proteoform_peaks_matched', 
    '#matched fragment ions': 'opt_proteoform_fragments_matched', 
    'E-value': 'best_search_engine_score[1]',
    'Proteoform': 'opt_proteoform_sequence',
    'Spectrum ID': 'opt_proteoform_evidence_spectra',
    'Prsm ID': 'opt_PrSM_ID'
}

# TODO spectrum identifier for proteoforms more than one???
# TODO PrSM to proteoform rows via spectra ???
# TODO modifications

dfs['pro_single'].drop(columns=dfs['pro_single'].columns.difference(col_map_pr.keys()), inplace=True)
dfs['pro_single'].rename(columns=col_map_pr, errors="raise", inplace=True)
dfs['pro_single'] = dfs['pro_single'].join(pd.DataFrame({k:[v]*len(dfs['pro_single']) for k,v in rigged_cols_pr.items()}))
dfs['pro_single']["opt_proteoform_sequence"] = dfs['pro_single']["opt_proteoform_sequence"].apply(strip_seq)  

col_order = ["PRH","opt_proteoform_sequence"]
dfs['pro_single'] = dfs['pro_single'].reindex(columns=col_order+list((set(dfs['pro_single'].columns.to_list())- set(col_order))))

rigged_cols_sm = { 
    'database': 'UniProtKB',
    'database_version': datetime.today().strftime('%Y-%m') + ' (' + str(db_len) + ')',
    'search_engine': '[MS, MS:1002901, TopPIC, ]',
    'PSH': 'PSM',
    'unique':  False,
}

# TODO modifications
# TODO document PSM_ID == PrSM_ID in MTD
# TODO scan or spectrum id? also value = ms_run[1-n]:{SPECTRA_REF}

col_map_sm = {
    'Retention time':'retention_time', 
    'Proteoform': 'sequence',
    'Prsm ID': 'PSM_ID',
    'E-value': 'search_engine_score[1]',
    'MIScore': 'search_engine_score[2]',
    'Protein accession': 'accession', 
    'First residue': 'start',
    'Last residue': 'end',
    '#variable PTMs': 'modifications', 
    'Scan(s)': 'spectra_ref',
    '#peaks':'opt_prsm_peak_number', 
    'Charge':'charge', 
    'Precursor mass': 'opt_prsm_precursormass', 
    'Adjusted precursor mass': 'opt_prsm_adj_precursormass', 
    '#unexpected modifications': 'opt_prsm_unexpected_modifications', 
    '#matched peaks': 'opt_prsm_peaks_matched', 
    '#matched fragment ions': 'opt_prsm_fragments_matched', 
}

massaged_cols = {
    'pre': dfs['prsm_single']['Proteoform'].str.extract(r'^(\w)\.').replace(np.nan,'-'),
    'post': dfs['prsm_single']['Proteoform'].str.extract(r'.*\.(\w*)$').replace('','-')
}

# TODO what to do with:
# exp_mass_to_charge
# calc_mass_to_charge

dfs['prsm_single'].drop(columns=dfs['prsm_single'].columns.difference(col_map_sm.keys()), inplace=True)
dfs['prsm_single'].rename(columns=col_map_sm, errors="raise", inplace=True)
dfs['prsm_single'] = dfs['prsm_single'].join(
    pd.DataFrame({k:[v]*len(dfs['prsm_single']) for k,v in rigged_cols_sm.items()})
)
for k,v in massaged_cols.items():
    dfs['prsm_single'][k] = v
# dfs['prsm_single']["sequence"] = dfs['prsm_single']["sequence"].apply(strip_seq)  
col_order = ["PSH","sequence"]
dfs['prsm_single'] = dfs['prsm_single'].reindex(columns=col_order+list((set(dfs['prsm_single'].columns.to_list())- set(col_order))))


with open('/tmp/test.mzTab','w') as f:
    f.writelines( header.format(filename='file://tmp/st_2.mzML') )
    f.write('\n')
    f.write(dfs['pro_single'].to_csv(sep='\t', index=False))
    f.write('\n')
    f.write(dfs['prsm_single'].to_csv(sep='\t', index=False))
