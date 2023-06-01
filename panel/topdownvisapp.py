import numpy as np
import matplotlib.pyplot as plt
import panel as pn
import re
import nextflow
import pandas as pd
from typing import List,Dict,Any
from ms_io_utils import load_mzml, load_ids, AppInputs, WorkflowResults
import panel_plot_utils as ppu
import holoviews as hv

# font_url="https://fonts.googleapis.com/css?family=Open+Sans"
# # template=pn.template.FastListTemplate(font_url="https://fonts.googleapis.com/css?family=Open+Sans")
# # TypeError: 'FastListTemplate' object is not callable
# print(pn.template.FastListTemplate.font_url)
# pn.template.FastListTemplate.font_url = font_url

pn.extension(loading_spinner='dots', 
    loading_color='#4d8060', 
    sizing_mode="stretch_width", 
    template="fast",
    theme='default',
    # theme_toggle=False,  # ?! AttributeError: 'theme_toggle' is not a valid config parameter.
)

# pn.template.FastListTemplate.font_url = font_url
# print(pn.template.FastListTemplate.font_url)

pd.set_option("display.precision", 0)
hv.extension('matplotlib')  # https://panel.holoviz.org/reference/panes/HoloViews.html#switching-backends

PALETTE = [
    "#ff6f69",
    "#ffcc5c",
    "#88d8b0",
]
ACCENT_BASE_COLOR = PALETTE[2]
TIMD_CONST = 1.002371  # TopDown isotope mass difference 55k u see OpenMS::Constants in kyowons branch
MAX_SIZE_MB = 150
BASE_DIR = "/tmp/results/"
MODS_IN = {'Meth': 'common_mods_Meth.txt',
    'Ox': 'common_mods_Ox.txt',
    'OxMeth': 'common_mods_OxMeth.txt',
    'OxMethAcet': 'common_mods_OxMethAcet.txt',
    'OxMethAcetPhos': 'common_mods_OxMethAcetPhos.txt',
    'OxMethPhos': 'common_mods_OxMethPhos.txt'}

input_path_store = AppInputs()
workflow_result_store = WorkflowResults()

def fileselect_callback(*events):
    event = events[-1]
    if not event.new:
        fileinputselections.object = default_fipm
    else:
        reraw = re.compile(".*raw", re.IGNORECASE)
        raws_selected = list(filter(reraw.match, event.new))
        refasta = re.compile(".*fasta", re.IGNORECASE)
        fastas_selected = list(filter(refasta.match, event.new))
        fileinputselections.object = untangle_inputs(
            raws_selected, fastas_selected, 
            list(set(event.new) - set(raws_selected) - set(fastas_selected)))

def untangle_inputs(raw, fasta, rest):
    global input_path_store
    selection_message = "You selected "
    if not raw:
        selection_message += "_no_ raw file (**required**)"
        cond1 = False
    else:
        selection_message += "`{}` as your raw file".format(raw[0])
        if len(raw)>1:
            selection_message += " (mind you though that you selected multiple raw files, I pick the _first_ regardeless)"
        cond1 = True
    if not fasta:
        selection_message += " and _no_ fasta file (**required**)"
        cond2 = False
    else:
        selection_message += " and `{}` as your fasta file".format(fasta[0])
        if len(fasta)>1:
            selection_message += " (mind you though that you selected multiple raw files, I pick the _first_ regardeless)"
        cond2 = True
    if rest:
        selection_message += ", but I dont know what to do with these ({}), hence will _ignore_.".format(','.join(rest))
    
    if cond1 and cond2:
        trigger_wf_btn.disabled = False
        input_path_store = AppInputs(raw_path=raw[0],fasta_path=fasta[0])
    else:
        trigger_wf_btn.disabled = True
        input_path_store = AppInputs()
    
    return selection_message

def update_resultoverview():
    overview = """Result Overview on {}"""
    if not workflow_result_store.run_name:
        resultoverview.object = default_ro
    else:
        resultoverview.object = overview.format(workflow_result_store.run_name)
    mztab_filepath = workflow_result_store.id_mztabpath
    mztab_button = pn.widgets.FileDownload(
        file=mztab_filepath, button_type='success', auto=False,
        embed=False, name="(Right-click using 'Save as' dialog)"
    )
    sidebar[1] = mztab_button
    id_tbl.value = workflow_result_store.id_dfs['prsms'][['sequence','search_engine_score[2]','search_engine','charge','accession', 'opt_prsm_precursormass','opt_prsm_fragments_matched','spectra_ref']]

def update_spectra_tbl():
    spec_df = pd.DataFrame({
        # 'Idx': [s['index'] for s in workflow_result_store.deconv_spectra],
        'MS level': [s['ms level'] for s in workflow_result_store.deconv_spectra.values()],
        'Scan': [s['id'] for s in workflow_result_store.deconv_spectra.values()],
        'RT': [s['scanList']['scan'][0]['scan start time'] for s in workflow_result_store.deconv_spectra.values()],
        '# peaks': [len(s['intensity array']) for s in workflow_result_store.deconv_spectra.values()],
    }) #columns=['MS level', 'Scan', 'RT', '# peaks'])
    spec_df.Scan = spec_df.Scan.str.extract(r'scan=(\d+)')
    spectra_tbl.value = spec_df

def trigger_wf_func(event):
    global workflow_result_store
    print('Clicked wf_btn {0} times'.format(trigger_wf_btn.clicks))
    print('I will use the files: {} and {}'.format(input_path_store.raw_path, input_path_store.fasta_path))
    print('I will include modifications from: {}'.format(mods_radio_group.value))
    input_panel.visible = False
    result_panel.visible = True

    print("Analysing data")
    # TODO move hardcode-paths to config of configs
    analysis_pipeline = nextflow.Pipeline(
        "/opt/app/wf/topdown_local_app.nf", 
        config="/opt/app/config/nf.config"
    )

    with pn.param.set_values(result_panel, loading=True):
        # ?! this is a blocking process - oh why?!
        apr = analysis_pipeline.run(params={
            "raw_file": input_path_store.raw_path, 
            "fasta_file": input_path_store.fasta_path,
            "mods": "/opt/app/modconf/" + mods_radio_group.value, 
        })
        # might be switched off for debugging
        print("Loading spectra")
        workflow_result_store = WorkflowResults(*load_mzml(BASE_DIR))
        print("Loading ids")
        workflow_result_store.id_dfs, workflow_result_store.id_mztabpath = load_ids(BASE_DIR)
        print("Workflow finished.")
        update_resultoverview()
        update_spectra_tbl()

def update_peak_tbl(sidx):
    with pn.param.set_values(peak_tbl, loading=True):
        ssidx = spectra_tbl.value.Scan.iloc[sidx]
        spectra_name.object = spec_detail_msg(sidx,ssidx)
        sidx_sst = list(workflow_result_store.deconv_spectra.values())[sidx]['scanList']['scan'][0]['scan start time']
        peak_tbl.value = pd.DataFrame({
            'RT': [sidx_sst]*len(list(workflow_result_store.deconv_spectra.values())[sidx]['intensity array']),
            'Mass': list(workflow_result_store.deconv_spectra.values())[sidx]['m/z array'],
            'Intensity': list(workflow_result_store.deconv_spectra.values())[sidx]['intensity array'],
        })
        # reset peak selected and plotted
        peak_name.object= "No peak selected."
        result_panel[0][2] = ppu.plot_3d_spectrum(None, None, [], {})

def spec_detail_msg(sidx,ssidx):
    return f"""Selected spectrum: index **{sidx}** (scan={ssidx})."""

def peak_detail_msg(pidx,sidx):
    return f"""Selected peak: index **{pidx}** (spectrum index {sidx})."""

def update_spec_fig(sidx):
    global result_panel
    # print('update_spec_fig', sidx) 

    with pn.param.set_values(result_panel[1][3], loading=True):
        spectra_fig = ppu.plot_2d_spectra(sidx, 
            workflow_result_store.annot_spectra, 
            workflow_result_store.deconv_spectra)
        result_panel[1][3] = spectra_fig
        #TODO replace result_panel[N] (2d) with named element replace with .object?

default_fipm = """You still need: 
1. a fasta file as identification sequences basis
2. a raw file as peak inputs.
"""
default_ro = """No peak file loaded, no identifications available."""

fileinputselections = pn.pane.Markdown(object=default_fipm)

fise = pn.widgets.FileSelector('/opt/')
fise.param.watch(fileselect_callback, 'value', onlychanged=False)
mods_radio_group = pn.widgets.RadioButtonGroup(
    name='Modifications Choice', options=MODS_IN, 
    orientation='vertical', button_type='warning')

trigger_wf_btn = pn.widgets.Button(name='\u25b6 Start Workflow', width=50, disabled=True)
trigger_wf_btn.on_click(trigger_wf_func)

resultoverview = pn.pane.Markdown(object=default_ro)
current_spec_idx = None
current_peak_idx = None

mztab_button = pn.widgets.Button(name='download mzTab', width=50, disabled=True)

spectra_tbl = pn.widgets.Tabulator(
        pd.DataFrame({},columns=['MS level', 'Scan', 'RT', '# peaks']), 
        # height=400,
        disabled=True,
        pagination='local',
        page_size=10,
        header_filters=True,
        widths={'index':'5%', 'MS level':'20%', 'Scan':'25%', 'RT':'25%', '# peaks':'25%'},
    )

def spectra_tbl_click(*events):
    event = events[-1]
    # print('spectra_tbl_click',f'Clicked cell in {event.column!r} column, row {event.row!r} with value {event.value!r}')
    update_peak_tbl(event.row)  # spectra_tbl.value.Index.iloc[event.row]
    update_spec_fig(event.row)
    global current_spec_idx
    current_spec_idx = event.row

spectra_tbl.on_click(spectra_tbl_click) 
spectra_name = pn.pane.Markdown(object="No spectrum selected.")

peak_tbl = pn.widgets.Tabulator(
        pd.DataFrame({},columns=['RT', 'Mass', 'Intensity']),
        height=300,
        disabled=True,
        sizing_mode='stretch_width',
        pagination='local',
)

peak_name = pn.pane.Markdown(object="No peak selected.")

def peak_tbl_click(*events):
    global result_panel
    global current_peak_idx
    if events:
        event = events[-1]
        current_peak_idx = event.row
    rpp = result_panel[0][2]
    peak_name.object = peak_detail_msg(current_peak_idx, current_spec_idx)
    with pn.param.set_values(rpp, loading=True):
        spectra_fig = ppu.plot_3d_spectrum(pidx=current_peak_idx, 
            sidx=current_spec_idx, 
            vis_dict=workflow_result_store.vis_dict,
            deconv_spectra=workflow_result_store.deconv_spectra,
            azimuth=azi_slider.value,
            elevation=ele_slider.value)
        rpp = spectra_fig
        result_panel[0][2] = rpp
        #TODO replace result_panel[M] (3d) with named element replace with .object?

peak_tbl.on_click(peak_tbl_click) 

id_tbl = pn.widgets.Tabulator(
        pd.DataFrame({},columns=['RT', 'mass', 'sequence']),
        # height=300,
        pagination='local',
        page_size=10,
        disabled=True,
        widths={'index': '5%', 'sequence': '15%'}, 
        sizing_mode='stretch_width',
)

def id_tbl_click(*events):
    event = events[-1]
    # print('id_tbl_click',f'Clicked cell in {event.column!r} column, row {event.row!r} with value {event.value!r}')
    # print('and in spectra_ref col that is', id_tbl.value.spectra_ref.iloc[event.row])  # I need col spectra_ref
    update_peak_tbl(id_tbl.value.spectra_ref.iloc[event.row])  # spectra_tbl.value.Index.iloc[event.row]
    update_spec_fig(id_tbl.value.spectra_ref.iloc[event.row])
    global current_spec_idx
    current_spec_idx = id_tbl.value.spectra_ref.iloc[event.row]

id_tbl.on_click(id_tbl_click) 

azi_slider = pn.widgets.IntSlider(name='3D peak plot Azimuth', start=0, end=180, step=5, value=40, value_throttled=(1, 45))
ele_slider = pn.widgets.IntSlider(name='3D peak plot Elevation', start=0, end=90, step=5, value=20, value_throttled=(1, 20))

def update_peak_plot_aspect(value):
    peak_tbl_click()

pn.bind(update_peak_plot_aspect, value=azi_slider.param.value_throttled, watch=True)
pn.bind(update_peak_plot_aspect, value=ele_slider.param.value_throttled, watch=True)

# ===
# LAYOUT
# ===

sidebar = pn.Column(
    pn.Row(spectra_tbl),
    pn.Row(mztab_button),
    pn.Row(azi_slider),
    pn.Row(ele_slider),
).servable(target='sidebar')

input_panel = pn.Column(
    pn.Row('# File stuffs'),
    pn.Row(
        pn.Column(fileinputselections, width=400),
        pn.Spacer(width=25),
        pn.Column(mods_radio_group, width=200),
        pn.Spacer(width=25),
        pn.Column(trigger_wf_btn, width=200)
    ),
    pn.Row(fise),
    sizing_mode='stretch_width',
).servable(target='main')

col_right = pn.Column(
    pn.Row('### Spectrum Details'),
    pn.Row(resultoverview),
    pn.Row(spectra_name),
    pn.Row(ppu.plot_2d_spectra(None, [], [])),
    pn.Row(peak_tbl),
    width=300
)
col_left = pn.Column(
    pn.Row('### Deconvolved Peaks'),
    pn.Row(peak_name),
    pn.Row(ppu.plot_3d_spectrum(None, None, [], {})),
    pn.Row('### Identified PrSM'),
    pn.Row(id_tbl),
    width=700
)
result_panel = pn.Row(
    col_left,
    col_right,
    visible = False,
).servable(target='main')


pn.state.template.param.update(
    site="TopDownViz",
    title="From TopDown raw data to deconvolved spectra visualisations with identifications",
    accent_base_color=ACCENT_BASE_COLOR,
    header_background=ACCENT_BASE_COLOR,
)
