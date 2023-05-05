import os
import re
import itertools
import numpy as np
import pandas as pd
from pyteomics import mzml
from typing import List, Dict, Tuple, Any
from dataclasses import dataclass, field

# TopDown isotope mass difference 55k u see OpenMS::Constants in kyowons branch
TIMD_CONST = 1.002371  
# TODO move into config file?!

@dataclass
class AppInputs:
    raw_path: str = ""
    fasta_path: str = ""

@dataclass
class WorkflowResults:
    run_name: str = "N/A"
    deconv_spectra: dict = field(default_factory=dict)
    annot_spectra: dict = field(default_factory=dict)
    vis_dict: dict = field(default_factory=dict)
    deconv_spectra_df:pd.DataFrame = pd.DataFrame()
    id_dfs: dict = field(default_factory=dict)
    id_mztabpath: str = ""


@dataclass
class SpecRef:
    spectrum_id: str
    tolerance: float
    massoffset: float
    chargemass: float
    isotoperangelimits: List[Tuple[int]]
    chargerangelimits: np.ndarray

@dataclass
class TargetRef:
    spectrum_id: str
    peak_index: int 
    mass: float
    charge: int
    mass_matches: np.ndarray
    intensity_matches: np.ndarray
    isotope_matches: np.ndarray


def parse_deconv_spectra_meta(spectrum: Dict[str,Any]) -> SpecRef:
  # value="tol=10;massoffset=0.000000;chargemass=1.007276;precursorscan=0;precursormass=0;peaks=1:1,0:1;1:1,0:1;1:1,0:1;cos=0.99441,0.995398,0.952226,;snr=7.9401,43.5187,3.78004,;qscore=0.817107,0.758096,0.651915,;qvalue=0.817107,0.758096,0.651915,"
  if "precursorscan" in spectrum['DeconvMassInfo'] and\
     "precursormass" in spectrum['DeconvMassInfo']:
    tolerance, massoffset, chargemass, precursorscan, precursormass, *peaknotes = spectrum['DeconvMassInfo'].split(';')
    precursorscan = int(precursorscan.split('=')[1])
    precursormass = float(precursormass.split('=')[1])
  else:
    tolerance, massoffset, chargemass, *peaknotes = spectrum['DeconvMassInfo'].split(';')
  tolerance = float(tolerance.split('=')[1])
  massoffset = float(massoffset.split('=')[1])
  chargemass = float(chargemass.split('=')[1])

  # TODO this is only necessary because of the ';' separator reuse
  # TODO best would be to avoid paired pairs altogether peaks=(zip(chargeranges, isotoperanges)), one chargerange=[min:max]
  peaknotes = peaknotes[::-1]
  token = peaknotes.pop()
  peaks = list()
  while token.startswith("peaks=") or not re.match("^[a-zA-Z]+=.*", token):
    peaks.append(token)
    token = peaknotes.pop()
  peaknotes.append(token)  # append last after condition test fail
  peaknotes.append('|'.join(peaks))
  peaknotes = peaknotes[::-1]

  notes = dict()
  for token in re.finditer(r"(([a-zA-Z]+=)([\d\:\.,|]+))", ';'.join(peaknotes)):
    # print(token)
    k=token.group(2).rstrip('=')
    v=token.group(3).rstrip(';,|')
    if k == "peaks":
      notes.update({k: v.split('|')})
    else:
      notes.update({k: [float(i) for i in v.split(',')]})   
  
  if "peaks" in notes:
    chargeranges, isotoperanges = [list(x) for x in zip(*[ tup.split(',') for tup in notes["peaks"]])]  #this is not correct
    chargerangelimits = [tuple(map(int, i.split(':'))) for i in chargeranges]
    isotoperangelimits = [tuple(map(int, j.split(':'))) for j in isotoperanges]
  # TODO above will fail if userParam is malformed (not exactly two \d\:\d per ';' split token)
  # TODO check len(spectrum) == len(chargerangelimits) == len(isotoperangelimits)
    specref = SpecRef(spectrum['id'],
                    tolerance,massoffset,chargemass,
                    chargerangelimits,isotoperangelimits)
  else:
    specref = SpecRef(None,None,None,None,None,None)
  return specref

def parse_source_spectra_meta(source_spectrum):
  par = [re.split('[,:]',x) for x in source_spectrum['DeconvMassPeakIndices'].split(';') if x]
  sublists = {m: list(map(int, i)) for m,*i in par}.values()
  indices = list(itertools.chain.from_iterable(sublists))
  indices.sort()
  return indices

def delta_ppm(m1,m2):
  return (m1 - m2) / m2 * 10**6

def get_match_window(mass,tolerance):
  """
  tolerance in ppm 
  mass in m/z
  returned is window border left and window border right
  """
  tol = tolerance/2
  md = ((tol / 10**6) * mass) 
  return mass-md, mass+md  

def calc_mz(givenmass, z, iso, massoffset, chargemass):
  if z>0:
    return (givenmass - massoffset + iso * TIMD_CONST)/z + chargemass
  else:
    return 0

def calc_range(givenmass, z, iso_min, iso_max, massoffset, chargemass):
  if z>0:
    return ((givenmass - massoffset + (iso_min-2) * TIMD_CONST)/z + chargemass, 
            (givenmass - massoffset + (iso_max+2) * TIMD_CONST)/z + chargemass )
  else:
    return (0,0)

def discharge_mz(givenmass, z, chargemass):
  return (givenmass - chargemass) * z

v_discharge_mz = np.vectorize(discharge_mz)
# TODO some peaks have 0 matches! => cannot call `vectorize` on size 0 inputs unless `otypes` is set

def get_source_peaks(range_l, range_r, source_spectrum, match_expect=None):
  target_idx = np.where(np.logical_and(source_spectrum["m/z array"] > range_l, source_spectrum["m/z array"] < range_r))[0]
  if not match_expect:
    match_expect = parse_source_spectra_meta(source_spectrum)
  isomatch_idx = np.intersect1d(target_idx, match_expect)
  return source_spectrum["m/z array"][target_idx], \
          source_spectrum["intensity array"][target_idx], \
          [1 if t in isomatch_idx else -1 for t in target_idx]

def acquire_targets_per_spectrum(deconv_spectrum: Dict[str,Any], source_spectrum: Dict[str,Any]):
  target_matches = list()
  specref = parse_deconv_spectra_meta(deconv_spectrum)
  isoref = parse_source_spectra_meta(source_spectrum)
  # TODO replace for-cascade with product 
  # see (http://stephantul.github.io/python/2019/07/20/product/
  # or https://note.nkmk.me/en/python-itertools-product/)
  for i in range(0,len(deconv_spectrum["m/z array"])):  # "mass array" rather
    for z in range(specref.chargerangelimits[i][0], specref.chargerangelimits[i][1]+1):
      range_l, range_r = calc_range(deconv_spectrum["m/z array"][i], z, 
                                    specref.isotoperangelimits[i][0], 
                                    specref.isotoperangelimits[i][1], 
                                    specref.massoffset, specref.chargemass)
      target_range_mz, target_range_int, target_range_iso = get_source_peaks(range_l, range_r, source_spectrum, isoref)
      if target_range_mz.size > 0:
        target_range_mass = v_discharge_mz(target_range_mz, z, specref.chargemass)
      else:
        target_range_mass = target_range_mz
      target_matches.append(
          TargetRef(deconv_spectrum["id"],
            i, deconv_spectrum["m/z array"][i], z,
            target_range_mass, target_range_int, target_range_iso)
      )
  return target_matches

def load_mzml(base_dir, common_prefix=None, suffix_pair=("deconv.mzML","annot.mzML")):
  if not common_prefix:
    spec_paths = [f for f in os.listdir(base_dir) if f.endswith('.mzML')]
    common_prefix = os.path.commonprefix(spec_paths)
  # print(common_prefix)
  with mzml.read(os.path.join(base_dir, common_prefix + suffix_pair[0])) as reader:
    deconv_spectra = {spectrum['id']: spectrum for spectrum in reader}

  with mzml.read(os.path.join(base_dir + common_prefix + suffix_pair[1])) as reader:
    annot_spectra = {spectrum['id']:spectrum for spectrum in reader}
  
  run_name = os.path.basename(common_prefix).rstrip('_')
	
  vis_dict = dict()  # is a dict of spectrum id ('controllerType=0 .. scan=101') to list of corresponding TargetRefs
  for spectrum_id in set(deconv_spectra.keys()).intersection(set(annot_spectra.keys())):
    vis_dict[spectrum_id] = acquire_targets_per_spectrum(deconv_spectra[spectrum_id], annot_spectra[spectrum_id])
  
  return run_name, deconv_spectra, annot_spectra, vis_dict

def load_ids(base_dir):
    h5_paths = [f for f in os.listdir(base_dir) if f.endswith('.h5')]
    if len(h5_paths) > 1:
        print("Warning about the presence of more than one h5 file in target. Picking the first one (alphabetically).")
    try:
        h5_path = os.path.join(base_dir, h5_paths[0])
        proteoforms = pd.read_hdf(h5_path, key='proteoforms')
        prsms = pd.read_hdf(h5_path, key='prsms')
        ids = {'proteoforms': proteoforms, 'prsms':prsms}
    except:
        ids = None
      
    mztab_paths = [f for f in os.listdir(base_dir) if f.endswith('.mzTab')]
    if len(mztab_paths) > 1:
        print("Warning about the presence of more than one mzTab file in target. Picking the first one (alphabetically).")
    try:
        mztab_path = os.path.join(base_dir, mztab_paths[0])
    except:
        mztab_path = None
    return ids, mztab_path

