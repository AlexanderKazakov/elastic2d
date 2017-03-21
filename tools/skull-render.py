#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import glob
vtkFiles = glob.glob("snapshots/vtk/*.vtu")
vtkFiles.sort()
vtkFiles = vtkFiles[:10]

# create a new 'XML Unstructured Grid Reader'
mesh1core00snap0 = XMLUnstructuredGridReader(FileName=vtkFiles)
mesh1core00snap0.CellArrayStatus = []
mesh1core00snap0.PointArrayStatus = ['Velocity', 'pressure', 'material_index']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [887, 496]

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Clip'
clip1 = Clip(Input=mesh1core00snap0)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'pressure']
clip1.Value = 137721.0
clip1.InsideOut = 0
clip1.Crinkleclip = 0

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [-0.35504913330078125, 1.255530834197998, 140.08167266845703]
clip1.ClipType.Normal = [1.0, 0.0, 0.0]
clip1.ClipType.Offset = 0.0

# # show data in view
clip1Display = Show(clip1, renderView1)
# # trace defaults for the display properties.
clip1Display.Representation = 'Surface'
# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'pressure'))
# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

animationScene1.GoToLast()
# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)
animationScene1.GoToFirst()

# current camera placement for renderView1
renderView1.CameraPosition = [-57.1526875876323, 1.255530834197998, 140.08167266845703]
renderView1.CameraFocalPoint = [-0.35504913330078125, 1.255530834197998, 140.08167266845703]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 26.042496856192816

# save animation images/movie
import shutil
shutil.rmtree('snapshots/pictures', ignore_errors=True)
import os
os.mkdir('snapshots/pictures')
WriteAnimation('snapshots/pictures/p.png', Magnification=1, FrameRate=15.0, Compression=True)
