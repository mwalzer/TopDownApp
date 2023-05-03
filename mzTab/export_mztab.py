import numpy as np
import pandas as pd
from datetime import datetime
import click
from click import command
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser
from os import path

# TODO
# nice to have: MTD software[1-n]-setting[1-n] ??? 
# do we need full modification-choice flexibility in our app?
# protein_search_engine_score is E-value too, but no CV?!
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
MTD sample_processing[1] [ERO, ERO:0000763, concentration calculation]|[NCIT, NCIT:C124326, Aliquotting, ]
MTD sample_processing[2] [CHMO, CHMO:0000524, liquid chromatography-mass spectrometry]
MTD fixed_mod[1] [UNIMOD, UNIMOD:4, Carbamidomethyl, ]
MTD variable_mod[1] [UNIMOD, UNIMOD:35, Oxidation, ]
COM PSM section entries are in fact PrSM
"""

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
INFO = '''
The selected toppic files has {n} different proteoforms registered, from {m} different PrSMs. 
'''
def print_help():
    """
    Print the help of the tool
    :return:
    """
    ctx = click.get_current_context()
    click.echo(ctx.get_help())
    ctx.exit()

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

# deprecated with click.File handle use
def peek_toppic_res_path(file: str) -> int:
    with open(file, 'r') as f:
        lines = f.readlines()
    sl = 0
    for i,l in enumerate(lines):
        if '*** Parameters ***' in l:
            sl = i+1  # adjust for 0-start
    return sl+1  # add one more for a blank line

def strip_seq(s:str) -> str: 
    rs = s.split('.')
    if len(rs) > 0:
        return '.'.join(rs[1:-1])  
    else:
        return rs[0]

def extract_peakfile_name(s:str, l:int) -> str:
    ps = pd.read_csv(s, sep='\t', skiprows=1, nrows=l-2, header=None, on_bad_lines='warn')
    ps[0] = ps[0].str.strip()
    ps.set_index(0, inplace=True)
    pl = ps.to_dict('index')
    tpfn = pl['Spectrum file:'][1]
    return path.splitext(path.basename(tpfn))[0] + '.mzML'

@click.command(short_help='toppic2mztab will export a summary style mzTab from the selected TopPic files and optionally linked h5 dataframes for further use.')
@click.option('-p', '--prsms_single', 'prsms_single', type=click.Path(exists=True,readable=True), 
    required=True, help="The PrSMs single run file from a single run TopPIC analysis")
@click.option('-f', '--proteoforms_single', 'proteoforms_single', type=click.Path(exists=True,readable=True),
    required=True, help="The Proteoforms single run file from a single run TopPIC analysis")
@click.option('-s', '--fasta', 'fasta', type=click.Path(exists=True,readable=True),
    required=True, help="The fasta file used from the same single run TopPIC analysis")
@click.argument('output_filepath', type=click.Path(writable=True) )  # help="The output destination path for the produced mzTab file")
@click.option('-5', '--hdf5', 'h5', type=click.Path(writable=True),
    help="The hdf5 destination path for the produced hdf5 dataframe (keys: proteoforms and prsms)")
def toppic2mztab(prsms_single, proteoforms_single, output_filepath, fasta, h5):
    """
    toppic2mztab will export a summary style mzTab from the selected TopPic files and 
    optionally linked h5 dataframes for further use. Warning: HDF5 may cause issues 
    between different versions.
    """
    if not any([prsms_single,proteoforms_single,fasta,output_filepath]):
        print_help()
    try:
        pn = peek_toppic_res_path(prsms_single)
        fn = peek_toppic_res_path(proteoforms_single)
        if pn!=fn:
            raise(BaseException("Input files appear to be from different runs, aborting!"))
        peakfile_name = extract_peakfile_name(prsms_single, pn)
        if peakfile_name!=extract_peakfile_name(proteoforms_single, fn):
            raise(BaseException("Input files appear to be from different runs, aborting!"))
        dfs = {'prsm_single': pd.read_csv(prsms_single, sep='\t', skiprows=pn ), 
                'pro_single':  pd.read_csv(proteoforms_single, sep='\t', skiprows=fn )}
        try:
            with open(fasta, 'r') as h:
                pl = list(SimpleFastaParser(h))
                db_len = len(pl)
                db_tax = ','.join(list({simpleparse_ncbi_tax(e[0]) for e in pl}))
                db_species = ','.join(list({simpleparse_species(e[0]) for e in pl}))
        except:
            raise(BaseException("Can't read the fasta file, aborting!"))
    except Exception as e:
        click.echo(e)
        print_help()

    ## Proteoform table section
    # map input table to mzTab-TDP where possible
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
    # rig columns that need special treatment but static
    rigged_cols_pr = { 
        'ambiguity_members': 'null',	
        'taxid': db_tax,
        'species': db_species,
        'database': 'UniProtKB',
        'database_version': datetime.today().strftime('%Y-%m') + ' (' + str(db_len) + ')',
        'search_engine': '[MS, MS:1002901, TopPIC, ]',
        'PRH': 'PRT'
    }
    # assemble proteoform mzTab table
    dfs['pro_single'].drop(columns=dfs['pro_single'].columns.difference(col_map_pr.keys()), inplace=True)
    dfs['pro_single'].rename(columns=col_map_pr, errors="raise", inplace=True)
    dfs['pro_single'] = dfs['pro_single'].join(pd.DataFrame({k:[v]*len(dfs['pro_single']) for k,v in rigged_cols_pr.items()}))
    dfs['pro_single']["opt_proteoform_sequence"] = dfs['pro_single']["opt_proteoform_sequence"].apply(strip_seq)  
    #sort the order of appearance for the cols
    col_order = ["PRH","opt_proteoform_sequence"]
    dfs['pro_single'] = dfs['pro_single'].reindex(columns=col_order+list((set(dfs['pro_single'].columns.to_list())- set(col_order))))

    ## PrSM table section
    # rig columns that need special treatment but static
    rigged_cols_sm = { 
        'database': 'UniProtKB',
        'database_version': datetime.today().strftime('%Y-%m') + ' (' + str(db_len) + ')',
        'search_engine': '[MS, MS:1002901, TopPIC, ]',
        'PSH': 'PSM',
        'unique':  False,
    }
    # map input table to mzTab-TDP where possible
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
    # adjust cols to mzTab spec
    massaged_cols = {
        'pre': dfs['prsm_single']['Proteoform'].str.extract(r'^(\w)\.').replace(np.nan,'-'),
        'post': dfs['prsm_single']['Proteoform'].str.extract(r'.*\.(\w*)$').replace('','-')
    }
    # clean the sequence column before assembly
    dfs['prsm_single']["Proteoform"] = dfs['prsm_single']["Proteoform"].apply(strip_seq)  

    dfs['prsm_single'].drop(columns=dfs['prsm_single'].columns.difference(col_map_sm.keys()), inplace=True)
    dfs['prsm_single'].rename(columns=col_map_sm, errors="raise", inplace=True)
    dfs['prsm_single'] = dfs['prsm_single'].join(
        pd.DataFrame({k:[v]*len(dfs['prsm_single']) for k,v in rigged_cols_sm.items()})
    )
    #sort the order of appearance for the cols
    for k,v in massaged_cols.items():
        dfs['prsm_single'][k] = v
    col_order = ["PSH","sequence"]
    dfs['prsm_single'] = dfs['prsm_single'].reindex(columns=col_order+list((set(dfs['prsm_single'].columns.to_list())- set(col_order))))

    print(INFO.format(m=len(dfs['prsm_single']), n=len(dfs['pro_single'])))

    # see 'header' for custom static mzTab file content documentation
    with open(output_filepath,'w') as f:
        f.writelines( header.format(filename=peakfile_name) )
        f.write('\n')
        f.write(dfs['pro_single'].to_csv(sep='\t', index=False))
        f.write('\n')
        f.write(dfs['prsm_single'].to_csv(sep='\t', index=False))

    if h5:
        dfs['prsm_single'].to_hdf(h5, key='prsms', mode='w')
        dfs['pro_single'].to_hdf(h5, key='proteoforms', mode='a')
    # read_hdf(h5_prsm, key='df')
    # make note that read will work guaranteed only on the same system, different environments may experience issues due protocol
    # see: https://github.com/DeepLabCut/DLCutils/issues/19 and https://stackoverflow.com/questions/63329657/python-3-7-error-unsupported-pickle-protocol-5


if __name__ == '__main__':
    toppic2mztab()
