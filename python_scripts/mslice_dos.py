##### Edit here

workspace_name = 'aluminium_Ei60'
temperature = 300
workspace_name = 'silicon_Ei80_300K'; temperature = 300
workspace_name = 'zinc_Ei50_300K'; temperature = 300
workspace_name = 'zinc_Ei50_100K'; temperature = 100
q_integration_range = 'Qmax/2, Qmax'
e_bin = 'Emax/10, Emax/100, Emax*0.9'
mean_square_displacement = 0

##### End edits

import numpy as np
import mslice.cli as mc
from mslice.models.workspacemanager.workspace_provider import add_workspace
import mslice.app as mslice_app
from mantid.simpleapi import ComputeIncoherentDOS, CreateSimulationWorkspace, CopyLogs, AddSampleLog, \
    ConvertToMD, CreateMDHistoWorkspace, DeleteWorkspace
from mslice.workspace import wrap_workspace
from mslice.models.axis import Axis
ws = mc.get_workspace_handle(workspace_name)
wss = mc.Slice(ws)
ws_out = ComputeIncoherentDOS(wss.raw_ws, Temperature=temperature, MeanSquareDisplacement=mean_square_displacement,
    QSumRange=q_integration_range, EnergyBinning=e_bin, StoreInADS=False)
name = workspace_name + '_DOS'
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
cut_axis = Axis('DeltaE', np.min(xx), np.max(xx), np.mean(np.diff(xx)), 'meV')
q_ax = wss.axes[0]
if isinstance(q_integration_range, str):
    try:
        out = [eval(qstr, None, {'Qmax':q_ax.end, 'Qmin':q_ax.start}) for qstr in q_integration_range.split(',')]
    except NameError:
        raise ValueError("Malformed QSumRange. Only the variables 'Qmin', 'Qmax' are allowed.")
    except:
        raise SyntaxError('Syntax error in QSumRange')
else:
    out = q_integration_range
int_axis = Axis('MomentumTransfer', out[0], out[1], out[1]-out[0], 'meV')
ws_out.axes = [cut_axis, int_axis]
add_workspace(ws_out, name)
mslice_app.MAIN_WINDOW.workspace_presenter.update_displayed_workspaces()
mslice_app.presenters.get_slice_plotter_presenter().update_displayed_workspaces()
mslice_app.presenters.get_cut_plotter_presenter().update_main_window()