##### Edit here

# Give the ion as the element letters followed by the valence number without + or -
ion = 'Er3'
q_start = 0
q_step = 0.02
q_end = 10
scale = 220
background = 140

##### End edits

import numpy as np
import mslice.cli as mc
from mslice.models.workspacemanager.workspace_provider import add_workspace
import mslice.app as mslice_app
from mantid.simpleapi import CreateSampleWorkspace, ScaleX, LoadInstrument, SofQW, MagFormFactorCorrection, \
    CreateSimulationWorkspace, CopyLogs, AddSampleLog, ConvertToMD, CreateMDHistoWorkspace, DeleteWorkspace, mtd
from mslice.workspace import wrap_workspace
from mslice.models.axis import Axis
ws = CreateSampleWorkspace(binWidth=0.1, XMin=0, XMax=1000, XUnit='DeltaE')
ws = ScaleX(ws, -15, "Add")
LoadInstrument(ws, InstrumentName='MARI', RewriteSpectraMap = True)
ws = SofQW(ws, [q_start-q_step/2., q_step, q_end+q_step/2.], 'Direct', 1000)
Q = ws.getAxis(1).extractValues()
for i in range(len(Q)-1):
    qv = ( (Q[i]+Q[i+1])*0.5 ) / 4 / np.pi
    y = ws.dataY(i)
    y *= np.exp(-16*qv*qv)
ws_corr = MagFormFactorCorrection(ws, IonName=ion, FormFactorWorkspace='FFCalc')
ws_out = (mtd['FFCalc'] * mtd['FFCalc']) * scale + background
DeleteWorkspace(mtd['FFCalc'])
DeleteWorkspace(ws)
name = ion + '_F2'
xdim = ws_out.getDimension(0)
extents = " ,".join(map(str, (xdim.getMinimum(), xdim.getMaximum())))
_tmpws = CreateSimulationWorkspace(Instrument='MAR', BinParams=[-1, 1, 1], UnitX='DeltaE', OutputWorkspace=name,
                                   EnableLogging=False)
CopyLogs(ws_out, _tmpws, EnableLogging=False)
AddSampleLog(_tmpws, LogName='Ei', LogText='3.', LogType='Number', EnableLogging=False)
_tmpws1 = ConvertToMD(_tmpws, EnableLogging=False, StoreInADS=False, PreprocDetectorsWS='-',
                     QDimensions='|Q|', dEAnalysisMode='Direct')
xx = ws_out.extractX()
ws_out = CreateMDHistoWorkspace(SignalInput=ws_out.extractY(), ErrorInput=ws_out.extractE(), Dimensionality=1,
                                Extents=extents, NumberOfBins=xdim.getNBins(), Names=name, Units='DeltaE',
                                EnableLogging=False, OutputWorkspace='__MSL'+name)
ws_out.copyExperimentInfos(_tmpws1)
DeleteWorkspace(_tmpws, EnableLogging=False)
DeleteWorkspace(_tmpws1, EnableLogging=False)
ws_out = wrap_workspace(ws_out, name)
cut_axis = Axis('|Q|', np.min(xx), np.max(xx), np.mean(np.diff(xx)), 'meV')
int_axis = Axis('DeltaE', -1, 1, 2, 'meV')
ws_out.axes = [cut_axis, int_axis]
add_workspace(ws_out, name)
mslice_app.MAIN_WINDOW.workspace_presenter.update_displayed_workspaces()
mslice_app.presenters.get_slice_plotter_presenter().update_displayed_workspaces()
mslice_app.presenters.get_cut_plotter_presenter().update_main_window()
