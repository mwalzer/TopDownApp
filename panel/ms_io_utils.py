import os
import re
import itertools
import numpy as np
import pandas as pd
from pyteomics import mzml
from typing import List, Dict, Tuple, Any
from dataclasses import dataclass

# TopDown isotope mass difference 55k u see OpenMS::Constants in kyowons branch
TIMD_CONST = 1.002371  
# TODO move into config file?!

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
  tolerance, massoffset, chargemass, *peaknotes = spectrum['DeconvMassInfo'].split(';')
  tolerance = float(tolerance.split('=')[1])
  massoffset = float(massoffset.split('=')[1])
  chargemass = float(chargemass.split('=')[1])
  peaknotes[0] = peaknotes[0].split('=')[1]  # remove 'header' before and up to '='
  peaknotes = [i.split(',') for i in peaknotes if i]  # 'if i' necessary because of trailing ','
  chargeranges, isotoperanges = list(map(list, zip(*peaknotes)))
  chargerangelimits = [tuple(map(int, i.split(':'))) for i in chargeranges]
  isotoperangelimits = [tuple(map(int, j.split(':'))) for j in isotoperanges]
  # TODO check len(spectrum) == len(chargerangelimits) == len(isotoperangelimits)
  specref = SpecRef(spectrum['id'],
                    tolerance,massoffset,chargemass,
                    chargerangelimits,isotoperangelimits)
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
    return 0

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

def load_mzml(base_dir):
  spec_paths = [f for f in os.listdir(base_dir) if f.endswith('.mzML')]
  comprefina = os.path.commonprefix(spec_paths)
  print(comprefina)
  with mzml.read(os.path.join(base_dir, comprefina + "deconv.mzML")) as reader:
    deconv_spectra = {spectrum['id']: spectrum for spectrum in reader}

  with mzml.read(os.path.join(base_dir + comprefina + "annot.mzML")) as reader:
    annot_spectra = {spectrum['id']:spectrum for spectrum in reader}
  
  run_name = os.path.basename(comprefina).rstrip('_')
	
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

