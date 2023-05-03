import holoviews as hv
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d.axis3d import Axis
from matplotlib.cm import get_cmap
from holoviews.operation.datashader import rasterize
import datashader as ds
from alphapept.constants import averagine_aa, averagine_avg
from alphapept.constants import isotopes as averagine_iso
from alphapept.chem import get_average_formula, mass_to_dist
from collections import defaultdict

def plot_2d_spectra(sidx, annot_spectra, deconv_spectra):
  if sidx not in range(0,len(deconv_spectra)):
    fig_d = hv.Spikes(pd.DataFrame({},columns=['mass','intensity'])).opts(color='green', title="Deconvolved Spectrum {}".format('spectrum["id"]'))
    fig_o = hv.Spikes(pd.DataFrame({},columns=['m/z','intensity'])).opts(color='blue', title="Original Spectrum {}".format('spectrum["id"]'))
    fig_s = fig_o.opts(aspect=3, padding=0.1) + fig_d.opts(aspect=3, padding=0.1)
    fig_s.opts(aspect_weight=True, tight=False, fig_inches=300, fig_size=3).cols(1)
    return fig_s
  spectrum = list(deconv_spectra.values())[sidx]
  ori_spectrum = annot_spectra[spectrum['id']]
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

# Remove margins in matplotlib 3D plots
# https://stackoverflow.com/a/42648316/3319796
def _get_coord_info_new(self, renderer):
    mins, maxs, cs, deltas, tc, highs = self._get_coord_info_old(renderer)
    correction = deltas * [1.0/4 + 6.0/11,
                           1.0/4 + 6.0/11,
                           1.0/4]
    mins += correction
    maxs -= correction
    return mins, maxs, cs, deltas, tc, highs

def dummy_3d_fig():
    fig = hv.Path3D(list()).opts(color="blue", azimuth=40, elevation=20)
    fig.opts(ylabel="Mass",
              xlabel="Charge",
              zlabel="Intensity",
    ).opts(title='Peaks of [ID]', invert_xaxis=True, fig_inches=6,)  #ax.yaxis.get_major_locator().set_params(integer=True)
    return fig

def plot_3d_spectrum(pidx, sidx, deconv_spectra, vis_dict):
  hv.extension('matplotlib')

  if not hasattr(Axis, "_get_coord_info_old"):
      Axis._get_coord_info_old = Axis._get_coord_info  
  Axis._get_coord_info = _get_coord_info_new

  if sidx not in range(0,len(deconv_spectra)): 
    return dummy_3d_fig()

  spectrum = list(deconv_spectra.values())[sidx]
  target_matches = vis_dict[spectrum['id']]
  # print(target_matches)

  #each peak has its matches, choice comes through pidx
  choice_peak_matches = list(filter(lambda t: t.peak_index==pidx, target_matches))

  if not choice_peak_matches: 
    return dummy_3d_fig()

  choice_peak_mass = choice_peak_matches[-1].mass

  plot_title = ' '.join([spectrum["id"].split(' ')[-1], 'precursor mass=', str(choice_peak_mass), 'selected peak index=', str(pidx)])

  #for each charge there is an element in target_matches with resp. peak index that now gets formed into vis_peaks
  hv_vis_peaks = list()
  vis_int_per_z = defaultdict(list)
  for t in choice_peak_matches:
    for idx, m in enumerate(t.mass_matches):
      hv_vis_peaks.append(
          {('x', 'y', 'z'): [[t.charge,m,0],[t.charge,m,t.intensity_matches[idx]]], 
          'type': 'noise' if t.isotope_matches[idx]<0 else 'isomatch'}
      )
      vis_int_per_z[t.charge].append(t.intensity_matches[idx])
  
  print(vis_int_per_z)
  maxz = max(vis_int_per_z.keys())
  minz = min(vis_int_per_z.keys())
  for c in range(minz,maxz+1):
    if c in vis_int_per_z:
      vis_int_per_z[c] = max(vis_int_per_z[c])
  maxi = max(vis_int_per_z.values())
  # mini = min(vis_int_per_z.values())  # only the min of max's
  print("maxz",maxz,"maxi",maxi,"minz",minz)

  axis_charge_min = max(0,minz-1)
  axis_charge_max = maxz+1

  avrgn_masses, avrgn_ints = mass_to_dist(choice_peak_mass, averagine_aa, averagine_iso)
  # print("avrgn", avrgn_masses, avrgn_ints)

  #for each charge there is also the averagine model peaks for the selected peak's precursor_mass 
  # which we collected in averagine_peaks with intensities 
  averagine_peaks = list() 
  for c,cmaxi in vis_int_per_z.items():
    # c size of avrgn_masses zip avrgn_ints
    averagine_peaks.append(
        {('x', 'y', 'z'): [[c, m, i*cmaxi] for m,i in zip(avrgn_masses, avrgn_ints)]}
    )

  # print("peaks", len(vis_peaks), len(averagine_peaks))
  if not hv_vis_peaks: 
    return dummy_3d_fig()

  explicit_cmapping = {'noise': 'lightcoral', 'isomatch':'mediumblue'}  # https://holoviews.org/user_guide/Styling_Plots.html
  fig = hv.Path3D(hv_vis_peaks, vdims='type').opts(azimuth=40, elevation=20)\
    .opts(color='type', cmap=explicit_cmapping)\
    .opts(xlim=(axis_charge_min,axis_charge_max),xticks=list(range(axis_charge_min,axis_charge_max+1)))
  # fog = hv.Scatter3D(iso_traces).opts(c='grey', s=.1, azimuth=40, elevation=20, alpha=0.1)
  fug = hv.Path3D(averagine_peaks).opts(color='black', linestyle='dashed', linewidth=.6, azimuth=40, elevation=20)  #, alpha=0.1
  fig_fin = (fug*fig).opts(
              ylabel="Mass",
              xlabel="Charge",
              zlabel="Intensity",
  ).opts(title=plot_title, invert_xaxis=True, fig_inches=6,)  #ax.yaxis.get_major_locator().set_params(integer=True)

  # server debug
  # f = '/tmp/results/'+'_'.join([str(s) for s in [pidx, sidx,len(deconv_spectra),len(vis_dict)]]) + '.png'
  # print(f)
  # hv.save(fig_fin, f, backend='matplotlib')

  return fig_fin