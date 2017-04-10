#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import glob

meshIds = [1, 2, 3, 4, 5]
gridReaders = []
for i in meshIds:
	fileNames = 'snapshots/vtk/mesh' + str(i) + 'core00snap*.vtu'
	vtkFiles = glob.glob(fileNames)
	vtkFiles.sort()
        step_through_files = 1
	vtkFiles = [vtkFiles[i] for i in range(1, len(vtkFiles), step_through_files)]
	print(vtkFiles)
	gridReader = XMLUnstructuredGridReader(FileName=vtkFiles)
	gridReader.PointArrayStatus = ['index_of_node', 'Velocity', 'pressure', 'material_index']
	gridReaders.append(gridReader)
	
# get animation scene
animationScene1 = GetAnimationScene()
# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# set a specific view size
renderView1.ViewSize = [1024, 1024]

# create a new 'Group Datasets'
groupDatasets1 = GroupDatasets(Input=gridReaders)

# create a new 'Clip'
clip1 = Clip(Input=groupDatasets1)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'pressure']
clip1.Value = 58656.0

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.0, 0.0, 146.0]
clip1.ClipType.Normal = [0.0, 0.0, -1.0]

# create a new 'Clip'
clip2 = Clip(Input=clip1)
clip2.ClipType = 'Plane'
clip2.Scalars = ['POINTS', 'pressure']
clip2.Value = 58655.5
clip2.ClipType.Origin = [0.28154797922160146, -0.19209043385484362, 136.3104248046875]
clip2.ClipType.Normal = [0.40254987371824447, -0.9153980550391319, 0.0]

# show data in view
clip2Display = Show(clip2, renderView1)

# set scalar coloring
ColorBy(clip2Display, ('POINTS', 'pressure'))

# show color bar/color legend
clip2Display.SetScalarBarVisibility(renderView1, True)


# reset view to fit data
renderView1.ResetCamera()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [30.530873875287273, 40.39820072237624, 162.1000319302619]
renderView1.CameraFocalPoint = [-0.35504913330078125, -3.669614315032959, 136.3104248046875]
renderView1.CameraViewUp = [-0.18298316009790655, -0.397873885547437, 0.899006971229906]
renderView1.CameraParallelScale = 15.444833020027161


# get color transfer function/color map for 'pressure'
pressureLUT = GetColorTransferFunction('pressure')

# Rescale transfer function
pressureLUT.RescaleTransferFunction(-0.0005, 0.0005)

# save animation images/movie
import shutil
shutil.rmtree('snapshots/pictures', ignore_errors=True)
import os
os.mkdir('snapshots/pictures')

# for animation
WriteAnimation('snapshots/pictures/skull.ogv', FrameRate=15.0, Compression=True)
# for batch of pictures
#WriteAnimation('snapshots/pictures/p.png')

shutil.copy('snapshots/pictures/skull.ogv', os.path.expanduser('~'))
