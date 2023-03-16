import os
from typing import List, Set, Dict, Tuple, Union, Any
from dataclasses import dataclass
import panel as pn  # Panel is a simple, flexible and enterprise-ready data app framework

from pyteomics import mzml, auxiliary
import numba
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import numpy  as np
from scipy.optimize import curve_fit
from collections import defaultdict

from alphapept.constants import averagine_aa, averagine_avg
from alphapept.constants import isotopes as averagine_iso
from alphapept.chem import get_average_formula, mass_to_dist

# import hvplot.pandas  # Adds .hvplot and .interactive methods to Pandas dataframes
import holoviews as hv
from holoviews import opts
from holoviews.operation.datashader import datashade, shade, dynspread, spread, rasterize
import holoviews.plotting.mpl
import holoviews.plotting.bokeh
import matplotlib.pyplot as plt
import datashader as ds
from mpl_toolkits.mplot3d.axis3d import Axis
from matplotlib.cm import get_cmap

import colorcet as cc
from bokeh.palettes import Magma, Inferno, Plasma, Viridis, Cividis
from holoviews.plotting.util import bokeh_palette_to_palette
from bokeh.models import HoverTool
# import random, pandas as pd, numpy as np, holoviews as hv, datashader as ds, colorcet as cc
# from hv.datashader import datashade, shade, dynspread, spread, rasterize

import nextflow # pip name nextflowpy  (inofficial btw)

pn.extension(loading_spinner='dots', loading_color='#4d8060', sizing_mode="stretch_width", template="fast")
pd.set_option("display.precision", 0)
hv.extension('bokeh', 'matplotlib')  # https://panel.holoviz.org/reference/panes/HoloViews.html#switching-backends

PALETTE = [
    "#ff6f69",
    "#ffcc5c",
    "#88d8b0",
]
ACCENT_BASE_COLOR = PALETTE[2]
TIMD_CONST = 1.002371  # TopDown isotope mass difference 55k u see OpenMS::Constants in kyowons branch
MAX_SIZE_MB = 150
BASE_DIR = "/tmp/results/"

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
  return (givenmass - massoffset + iso * TIMD_CONST)/z + chargemass

def calc_range(givenmass, z, iso_min, iso_max, massoffset, chargemass):
  return ((givenmass - massoffset + (iso_min-2) * TIMD_CONST)/z + chargemass, 
          (givenmass - massoffset + (iso_max+2) * TIMD_CONST)/z + chargemass )

def discharge_mz(givenmass, z, chargemass):
  return (givenmass - chargemass) * z

v_discharge_mz = np.vectorize(discharge_mz)
# TODO some peaks have 0 matches! => cannot call `vectorize` on size 0 inputs unless `otypes` is set

def get_source_peaks(range_l, range_r, source_spectrum):
  target_idx = np.where(np.logical_and(source_spectrum["m/z array"] > range_l, source_spectrum["m/z array"] < range_r))
  return source_spectrum["m/z array"][target_idx], \
          source_spectrum["intensity array"][target_idx], \
          np.ones(len(target_idx[0]))*-1

def acquire_targets_per_spectrum(deconv_spectrum: Dict[str,Any], source_spectrum: Dict[str,Any]):
  target_matches = list()
  specref = parse_deconv_spectra_meta(deconv_spectrum)
  # TODO replace for-cascade with product 
  # see (http://stephantul.github.io/python/2019/07/20/product/
  # or https://note.nkmk.me/en/python-itertools-product/)
  for i in range(0,len(deconv_spectrum["m/z array"])):
    for z in range(specref.chargerangelimits[i][0], specref.chargerangelimits[i][1]+1):
      range_l, range_r = calc_range(deconv_spectrum["m/z array"][i], z, 
                                    specref.isotoperangelimits[i][0], 
                                    specref.isotoperangelimits[i][1], 
                                    specref.massoffset, specref.chargemass)
      target_range_mz, target_range_int, target_range_iso = get_source_peaks(range_l, range_r, source_spectrum)
      if target_range_mz.size > 0:
        target_range_mass = v_discharge_mz(target_range_mz, z, specref.chargemass)
      else:
        target_range_mass = target_range_mz
      for e in range(specref.isotoperangelimits[i][0], specref.isotoperangelimits[i][1]+1):
          target_mz_iso = calc_mz(source_spectrum["m/z array"][i],z,e, specref.massoffset, specref.chargemass)
          mz_window_l, mz_window_r = get_match_window(target_mz_iso, specref.tolerance)
          np.put(target_range_iso, np.where(np.logical_and(target_range_mz > mz_window_l, target_range_mz < mz_window_r)), e)
      target_matches.append(
          TargetRef(deconv_spectrum["id"],
                    i, deconv_spectrum["m/z array"][i], z,
                    target_range_mass, target_range_int, target_range_iso)
      )
  return target_matches

def load_mzml():
	spec_paths = [f for f in os.listdir(BASE_DIR) if f.endswith('.mzML')]
	comprefina = os.path.commonprefix(spec_paths)
	with mzml.read(BASE_DIR + comprefina + "deconv.mzML") as reader:
		deconv_spectra = [spectrum for spectrum in reader]	

	with mzml.read(BASE_DIR + comprefina + "annot.mzML") as reader:
		annot_spectra = [spectrum for spectrum in reader]

	vis_dict = dict()
	for s_a, s_o in zip(deconv_spectra, annot_spectra):
		vis_dict[s_a['id']] = acquire_targets_per_spectrum(s_a,s_o)
	return deconv_spectra, annot_spectra, vis_dict

def load_ids():
    h5_paths = [f for f in os.listdir(BASE_DIR) if f.endswith('.h5')]
    if len(h5_paths) > 1:
        print("Warning of the presence of more than one h5 file in target. Picking the first one (alphabetically).")
    try:
        h5_path = os.path.join(BASE_DIR, h5_paths[0])
        proteoforms = pd.read_hdf(h5_path, key='proteoforms')
        prsms = pd.read_hdf(h5_path, key='prsms')
    except:
        return None, None
    return proteoforms, prsms

def plot_2d_spectra(sidx, annot_spectra, deconv_spectra):
  if not sidx:
    fig_d = hv.Spikes(pd.DataFrame({},columns=['mass','intensity'])).opts(color='green', title="Deconvolved Spectrum {}".format('spectrum["id"]'))
    fig_o = hv.Spikes(pd.DataFrame({},columns=['m/z','intensity'])).opts(color='blue', title="Original Spectrum {}".format('spectrum["id"]'))
    fig_s = fig_o.opts(aspect=3, padding=0.1) + fig_d.opts(aspect=3, padding=0.1)
    fig_s.opts(aspect_weight=True, tight=False, fig_inches=300, fig_size=3).cols(1)
    return fig_s
  spectrum = deconv_spectra[sidx-2]
  ori_spectrum = annot_spectra[sidx-2]
  peak_coord = pd.DataFrame(np.concatenate([spectrum["m/z array"][np.newaxis].T, spectrum["intensity array"][np.newaxis].T], axis=1), columns = ['mass','intensity'])
  fig_d = hv.Spikes(peak_coord).opts(color='green', title="Deconvolved Spectrum {}".format(spectrum["id"]))
  peak_coord = pd.DataFrame(np.concatenate([ori_spectrum["m/z array"][np.newaxis].T, ori_spectrum["intensity array"][np.newaxis].T], axis=1), columns = ['m/z','intensity'])
  fig_o = hv.Spikes(peak_coord).opts(color='blue', title="Original Spectrum {}".format(spectrum["id"]))
  fig_s = fig_o.opts(aspect=3, padding=0.1) + fig_d.opts(aspect=3, padding=0.1)
  # fig_s = fig_o.opts(aspect=2, padding=0.10) + fig_d.opts(aspect=2, padding=0.1)
  fig_s.opts(aspect_weight=True, tight=False, fig_inches=300, fig_size=3).cols(1)
  # fig_s.opts(aspect_weight=True, tight=False).cols(1)
  return fig_s


def plot_peak_map(spectra):
    # palette_inv = palette.reversed()
    # p_df = pd.DataFrame([ [s['m/z array'],s['intensity array'], s['id'], s['ms level']] for s in spectra if s['ms level']==1 ], columns = ['m/z array','intensity array', 'id', 'ms level'])
    # df = p_df.explode(['m/z array', 'intensity array'])

    N = int(10e6)
    x_r = (0,100)
    y_r = (100,2000)
    z_r = (0,10e8)
    x = np.random.randint(x_r[0]*1000,x_r[1]*1000,size=(N, 1))
    y = np.random.randint(y_r[0]*1000,y_r[1]*1000,size=(N, 1))
    z = np.random.randint(z_r[0]*1000,z_r[1]*1000,size=(N, 1))
    z2 = np.ones((N,1)).astype(int)
    df = pd.DataFrame(np.column_stack([x,y,z,z2]), columns=['x','y','z','z2'])
    df[['x','y','z']] = df[['x','y','z']].div(1000, axis=0)

    p=hv.Points(df,['x','y'], ['z','z2'])
    palette = get_cmap('viridis')
    fig_map=rasterize(p, aggregator=ds.sum("z2"),x_range=(0,100))
    fig_map.opts(xlim=(0,100), ylim=(100,2000), cmap=palette)
    return fig_map.hist(dimension='y',weight_dimension='x_y z2',num_bins = 2000,normed=True)


# interactive_2d = pn.bind(plot_2d_spectra, sidx=spec_selec, deconv_spectra=deconv_spectra, annot_spectra=annot_spectra)
# interactive_text = pn.bind(message, sidx=spec_selec)
# interactive_3d = pn.bind(plot_3d_spectrum, sidx=spec_selec, pidx=peak_selec, deconv_spectra=deconv_spectra, vis_dict=vis_dict)

def handle_incoming_peakfile(event):
  analysis_pipeline = nextflow.Pipeline(
    "/home/walzer/ms-tools/TopDown/topdown_local_v2.nf", 
    config="/opt/app/nf.conf"
  )
  rfp = "/home/walzer/ms-tools/TopDown/st_1.raw"
  # ?! this is a blocking process - oh why?!
  apr = analysis_pipeline.run(params={
      "raw_file": rfp, 
      "mods": "/home/walzer/ms-tools/TopDown/common_mods.txt", 
      "fasta": "/home/walzer/ms-tools/TopDown/TopPIC_tutorial/TopPIC_tutorial_uniprot-st.fasta"
  })
      # "raw_file": "/home/walzer/ms-tools/TopDown/TopPIC_tutorial/st_1.raw", 

  deconv_spectra, annot_spectra, vis_dict = load_mzml()


def tdv_app():

  pn.state.template.param.update(
      site="TopDownViz",
      title="Turn topdown MS runs into deconvolved spectra visualisations with identifications.",
      accent_base_color=ACCENT_BASE_COLOR,
      header_background=ACCENT_BASE_COLOR,
  )

  # Sidebar setup 
  # ===
  drop_field = pn.widgets.FileInput(name="Input raw file", accept='.raw,.Raw,.RAW').servable(area="sidebar")
  drop_field.param.watch(print, 'value')
  progress = pn.widgets.Progress(active=False).servable(area="sidebar")
  # drop_field.param.watch(reset_load_indicator, "value")
  # drop_field.param.watch(handle_incoming_peakfile, "value")
  pn.pane.Markdown("""Some TEXT.""").servable(area="sidebar")

  deconv_spectra_panel = pn.widgets.Tabulator(deconv_spectra_df).servable(area="sidebar", height=300)  # heigth arg orig Tabulator() param

  # Main window setup
  # ===
  pn.pane.Markdown("""
  The raw file you chose is being processed. 
  """).servable()


# TODO for much later: a help/report page
def markdown_app():
    return '# This is (going to be) the TopDownVis app help.'

def update_deconv_spectra_panel():
  global deconv_spectra
  global deconv_spectra_df
  if deconv_spectra:
    deconv_spectra_df = pd.DataFrame({'spectrum level':[s['ms level'] for s in deconv_spectra]})
  else: 
    deconv_spectra_df = pd.DataFrame({},columns=['index','ms level', 'RT'])
  
if __name__ == '__main__':
  # scope global vars - these must be replaced with a session-id'ed element in some fast (redis?)db when scaling up
  global deconv_spectra
  global deconv_spectra_df
  global vis_dict
   
  deconv_spectra = list()
  annot_spectra = list()
  vis_dict = dict()
  deconv_spectra_df = pd.DataFrame()
  update_deconv_spectra_panel()

  # https://panel.holoviz.org/user_guide/Server_Configuration.html
  # https://panel.holoviz.org/reference/widgets/FileInput.html#server-context
  pn.serve(
      {'help': markdown_app,
      'TopDownVis': tdv_app},
      # Increase the maximum websocket message size allowed by Bokeh
      websocket_max_message_size=MAX_SIZE_MB*1024*1014,
      # Increase the maximum buffer size allowed by Tornado
      http_server_kwargs={'max_buffer_size': MAX_SIZE_MB*1024*1014}
  )
