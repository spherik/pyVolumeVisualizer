# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 11:39:10 2016

@author: Antoni Gurgui
"""

import os
import numpy
import vtk
import scipy.io as sio
from vtk.util import numpy_support
import sys

class VTKVolumeVisualizer:
    def __init__(self):
        """ Class constructor. 
        Inputs:
            * mri_volume: RxCxDx1 volume with 1 component per voxel corresponding to MRI acquisition 
            * dti_volume: RxCxDx3 volume with 3 components per voxel corresponding to vector information
            * mask: RxCxDx1 volume with boolean values to mask the input volumes
        """
    
    
        # Create renderer, render window, interactor and assign TrackBallCamera as
        # interaction style
        self.renderer = vtk.vtkRenderer()
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.AddRenderer(self.renderer)
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self.iren.SetRenderWindow(self.renderWindow)
        self.iren.AddObserver("KeyPressEvent", self.Keypress)
            
        

        self.planes = vtk.vtkPlanes()
        
        self.lut = vtk.vtkLookupTable()
        
        self.has_mask = False


    def Keypress(self, obj, event):
        key = obj.GetKeySym()
        if key == "q":
            obj.InvokeEvent("DeleteAllObjects")
            sys.exit()
        elif key=="l":
            self.longAxis.SetEnabled(not(self.longAxis.GetEnabled()))
        elif key=="s":
            self.shortAxis.SetEnabled(not(self.shortAxis.GetEnabled()))
        self.iren.Render()
    
    # When interaction starts, the requested frame rate is increased.
    def StartInteraction(self,obj, event):
        self.renderWindow.SetDesiredUpdateRate(10)

    # When interaction ends, the requested frame rate is decreased to
    # normal levels. This causes a full resolution render to occur.
    def EndInteraction(self, obj, event):
        self.renderWindow.SetDesiredUpdateRate(0.001)

    # The implicit function vtkPlanes is used in conjunction with the
    # volume ray cast mapper to limit which portion of the volume is
    # volume rendered.
    def ClipVolumeRender(self, obj, event):
        obj.GetPlanes(self.planes)
        self.mapper.SetClippingPlanes(self.planes)

    def SetImageValueRange(self, min_value, max_value):
        self.lut.SetTableRange(min_value,max_value) 

    def SetMask(self, matImage):
        VTK_DATA = numpy_support.numpy_to_vtk(num_array=matImage.ravel(), deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
        self.mask_volume = vtk.vtkImageData()
        self.mask_volume.SetDimensions(matImage.shape[2], matImage.shape[1], matImage.shape[0]) #set dimensions as necessary
        self.mask_volume.SetOrigin(0,0,0) #set origin as necessary
        self.mask_volume.SetSpacing(1, 1, 1) #set spacing as necessary  
        self.mask_volume.GetPointData().SetScalars(VTK_DATA)

        self.has_mask = True
        self.mask = vtk.vtkImageMask()
        self.mask.SetMaskInputData(self.mask_volume)
        self.mask.SetInputData(self.image_volume)
        self.mask.SetMaskedOutputValue(255.0)

    def SetMRIImage(self, matImage):
        ###################################################
        ### Volume plane cuts
        ###################################################
    
        VTK_DATA = numpy_support.numpy_to_vtk(num_array=matImage.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        self.image_volume = vtk.vtkImageData()
        self.image_volume.SetDimensions(matImage.shape[2], matImage.shape[1], matImage.shape[0]) #set dimensions as necessary
        self.image_volume.SetOrigin(0,0,0) #set origin as necessary
        self.image_volume.SetSpacing(1, 1, 1) #set spacing as necessary  
        self.image_volume.GetPointData().SetScalars(VTK_DATA)
    
        
        
    def CreateVisualization(self):
        # Cell picker
        cellpicker = vtk.vtkCellPicker()
        cellpicker.SetTolerance(0.005)
        
        # Look-up table to color
        print(self.image_volume.GetScalarRange())
        self.lut.SetTableRange(-255,255)    # image intensity range
        self.lut.SetAlphaRange(1.0,1.0)
        self.lut.SetValueRange(0,1)  # from black to white
        self.lut.SetHueRange(0.66667, 0.0)
        self.lut.SetSaturationRange(1,1)
        self.lut.SetNumberOfColors(256)
        self.lut.Build()
        
        # Short Axis image plane
        self.shortAxis = vtk.vtkImagePlaneWidget()
        self.shortAxis.DisplayTextOn()
        if self.has_mask:
            self.shortAxis.SetInputConnection(self.mask.GetOutputPort())
        else:   
            self.shortAxis.SetInputData(self.image_volume)
        self.shortAxis.SetPlaneOrientationToXAxes()
        self.shortAxis.SetSliceIndex(self.image_volume.GetExtent()[2]/2) # put the plane in the middle
        self.shortAxis.SetPicker(cellpicker)
        self.shortAxis.SetKeyPressActivationValue('s')
        self.shortAxis.SetRightButtonAction(1)
        self.shortAxis.SetLookupTable(self.lut)
        
        # Short Axis image plane
        self.longAxis = vtk.vtkImagePlaneWidget()
        self.longAxis.DisplayTextOn()
        if self.has_mask:
            self.longAxis.SetInputConnection(self.mask.GetOutputPort())
        else:
            self.longAxis.SetInputData(self.image_volume)
        self.longAxis.SetPlaneOrientationToZAxes()
        self.longAxis.SetSliceIndex(self.image_volume.GetExtent()[1]/2) # put the plane in the middle
        self.longAxis.SetPicker(cellpicker)
        self.longAxis.SetKeyPressActivationValue('l')
        self.longAxis.SetRightButtonAction(1)
        self.longAxis.SetLookupTable(self.lut)
        
        # Add planes to the scene (no visible by default)
        self.shortAxis.SetInteractor(self.iren);
        self.longAxis.SetInteractor(self.iren);


        # Alpha channel function
        alphaChannelFunc = vtk.vtkPiecewiseFunction()
        alphaChannelFunc.AddPoint(-255.0, 0.0)
        #alphaChannelFunc.AddPoint(-90.01, 0.1)
        alphaChannelFunc.AddPoint(-90.0, 1.0)
        alphaChannelFunc.AddPoint(90.0, 1.0)
        #alphaChannelFunc.AddPoint(90.01, 0.1)
        alphaChannelFunc.AddPoint(255.0, 0.0)

        # Color transfer functions
        color1 = vtk.vtkColorTransferFunction()
        color1.SetColorSpaceToRGB()
        color1.AddRGBPoint(-255.0, 0.0, 0.0, 0.0);
        color1.AddRGBPoint(-90.0, 1.0,0.0,0.0)
        color1.AddRGBPoint(0.0, 0.0,1.0,0.0)
        color1.AddRGBPoint(90.0, 0.0,0.0,1.0)
        color1.AddRGBPoint(255.0, 0.0, 0.0, 0.0);

        funcOpacityGradient = vtk.vtkPiecewiseFunction()
        funcOpacityGradient.AddPoint(-255.0, 0.0)
        #alphaChannelFunc.AddPoint(-90.01, 0.1)
        funcOpacityGradient.AddPoint(-90.0, 1.0)
        funcOpacityGradient.AddPoint(90.0, 1.0)
        #alphaChannelFunc.AddPoint(90.01, 0.1)
        funcOpacityGradient.AddPoint(255.0, 0.0)

        # Volume rendering property
        propVolume = vtk.vtkVolumeProperty()
        propVolume.ShadeOff()
        propVolume.SetInterpolationTypeToLinear()
        propVolume.SetScalarOpacity(alphaChannelFunc)
        propVolume.SetGradientOpacity(funcOpacityGradient)
        propVolume.SetColor(color1)

        # lookupTable = vtk.vtkLookupTable();
        # lookupTable.SetNumberOfTableValues(256);
        # lookupTable.SetRange(-90.0, 90.0);
        # lookupTable.SetRampToLinear()
        # lookupTable.SetValueRange(0, 1)
        # lookupTable.SetAlphaRange(1.0, 1.0)
        # lookupTable.SetHueRange(0, 1)
        # lookupTable.SetSaturationRange(0, 0)
        # lookupTable.Build();

        # scalarValuesToColors = vtk.vtkImageMapToColors();
        # scalarValuesToColors.SetLookupTable(lookupTable);
        # scalarValuesToColors.PassAlphaToOutputOn();
        # scalarValuesToColors.SetInputData(self.image_volume);

        self.mapper = vtk.vtkOpenGLGPUVolumeRayCastMapper()
        if self.has_mask:
            self.mapper.SetInputConnection(self.mask.GetOutputPort())
        else:
            self.mapper.SetInputData(self.image_volume)
        #############################################

        newvol = vtk.vtkVolume()
        newvol.SetMapper(self.mapper)
        newvol.SetProperty(propVolume)
        
        # Add actor to renderer        
        self.renderer.AddActor(newvol)
        
        # Box widget
        self.boxWidget = vtk.vtkBoxWidget()
        self.boxWidget.SetInteractor(self.iren)
        self.boxWidget.SetPlaceFactor(1.0)
        self.boxWidget.RotationEnabledOff()
        # Place the interactor initially. The output of the reader is used to
        # place the box widget.
        self.boxWidget.SetInputData(self.image_volume)
        self.boxWidget.PlaceWidget()
        self.boxWidget.InsideOutOn()
        self.boxWidget.On()
        self.boxWidget.AddObserver("StartInteractionEvent", self.StartInteraction)
        self.boxWidget.AddObserver("InteractionEvent", self.ClipVolumeRender)
        self.boxWidget.AddObserver("EndInteractionEvent", self.EndInteraction)

        outlineProperty = self.boxWidget.GetOutlineProperty()
        outlineProperty.SetRepresentationToWireframe()
        outlineProperty.SetAmbient(1.0)
        outlineProperty.SetAmbientColor(1, 1, 1)
        outlineProperty.SetLineWidth(3)

        selectedOutlineProperty = self.boxWidget.GetSelectedOutlineProperty()
        selectedOutlineProperty.SetRepresentationToWireframe()
        selectedOutlineProperty.SetAmbient(1.0)
        selectedOutlineProperty.SetAmbientColor(1,0, 0)
        selectedOutlineProperty.SetLineWidth(3)
        
        # Add axes to the scene
        self.axes = vtk.vtkAxesActor()
        self.axes.SetTotalLength(100,100,100)
        self.renderer.AddActor(self.axes)

    def Start(self):
        self.renderer.SetBackground(0.0,0.0,0.0)
        self.renderWindow.SetSize(800,800)

        # Finally, execute pipeline and render
        self.iren.Initialize()
        
        self.renderer.ResetCamera()
        self.renderer.GetActiveCamera().Zoom(-1.5)
        self.renderWindow.Render()
        
        self.iren.Start() 

class VTKDTIVectorVisualizer:
    def __init__(self):
        """ Class constructor. 
        Inputs:
            * mri_volume: RxCxDx1 volume with 1 component per voxel corresponding to MRI acquisition 
            * dti_volume: RxCxDx3 volume with 3 components per voxel corresponding to vector information
            * mask: RxCxDx1 volume with boolean values to mask the input volumes
        """
    
    
        # Create renderer, render window, interactor and assign TrackBallCamera as
        # interaction style
        self.renderer = vtk.vtkRenderer()
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.AddRenderer(self.renderer)
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self.iren.SetRenderWindow(self.renderWindow)
        self.iren.AddObserver("KeyPressEvent", self.Keypress)

        # Create polydata to hold vector info
        self.polydata_v1 = vtk.vtkPolyData()
        self.polydata_v2 = vtk.vtkPolyData()
        self.polydata_v3 = vtk.vtkPolyData()    
        
        
    def Keypress(self, obj, event):
        key = obj.GetKeySym()
        if key == "q":
            obj.InvokeEvent("DeleteAllObjects")
            sys.exit()
        elif key == "1":
            self.ptMask.SetInputData(self.polydata_v1)
    
        elif key =="2":
            self.ptMask.SetInputData(self.polydata_v2)
            
        elif key =="3":
            self.ptMask.SetInputData(self.polydata_v3)
        elif key == "plus":
            self.ptMask.SetOnRatio(self.ptMask.GetOnRatio()+5)
            print(self.ptMask.GetOnRatio())
        elif key == "minus":
            if(self.ptMask.GetOnRatio()>0):
                self.ptMask.SetOnRatio(self.ptMask.GetOnRatio()-5)
            print(self.ptMask.GetOnRatio())
        elif key=="l":
            self.longAxis.SetEnabled(not(self.longAxis.GetEnabled()))
        elif key=="s":
            self.shortAxis.SetEnabled(not(self.shortAxis.GetEnabled()))
        self.iren.Render()

    def SetDTIPolydatas(self, dti_volume_v1, dti_volume_v2, dti_volume_v3, mask):
        
        points = vtk.vtkPoints()
        vertices = vtk.vtkCellArray()
        
        normals_v1 = vtk.vtkFloatArray()
        normals_v1.Allocate(1,1)
        normals_v1.SetNumberOfComponents(3)
        
        normals_v2 = vtk.vtkFloatArray()
        normals_v2.Allocate(1,1)
        normals_v2.SetNumberOfComponents(3)
        
        normals_v3 = vtk.vtkFloatArray()
        normals_v3.Allocate(1,1)
        normals_v3.SetNumberOfComponents(3)
        
        colors_v1 = vtk.vtkUnsignedCharArray();
        colors_v1.SetName('colors')
        colors_v1.SetNumberOfComponents(3)
        
        colors_v2 = vtk.vtkUnsignedCharArray();
        colors_v2.SetName('colors')
        colors_v2.SetNumberOfComponents(3)
        
        colors_v3 = vtk.vtkUnsignedCharArray();
        colors_v3.SetName('colors')
        colors_v3.SetNumberOfComponents(3)
        
        
        print('Copying data')
        for i in range(mask.shape[2]):
        	for j in range(mask.shape[1]):
        		for k in range(mask.shape[0]):
        			if(mask[k,j,i]):
        				#normal = [0.0,0.0,0.0]
        				
        				# Insert point
        				id = points.InsertNextPoint(i,j,k)
        
        				# Insert point cell to represent it
        				vertices.InsertNextCell(id)
        
            			# Radial local frame vectors
        				#normals_v1.InsertNextTuple(numpy.flipud(dti_volume_v1[k,j,i,:]))
        				normals_v1.InsertNextTuple(dti_volume_v1[k,j,i,:])
            
        				# Radial component color
        				colors_v1.InsertNextTuple(numpy.multiply(255,numpy.absolute(dti_volume_v1[k,j,i,:])))
        
        				# Circular component normal
        				normals_v2.InsertNextTuple(dti_volume_v2[k,j,i,:])
        
        				# circular component color
        				colors_v2.InsertNextTuple(numpy.multiply(255,numpy.absolute(dti_volume_v2[k,j,i,:])))
        				
        				# longitudinal component normal
        				normals_v3.InsertNextTuple(dti_volume_v3[k,j,i,:])
        
        				# longitudinal component color
        				colors_v3.InsertNextTuple(numpy.multiply(255,numpy.absolute(dti_volume_v3[k,j,i,:])))
        
        print('Done')

        # Fill radial component polydata
        self.polydata_v1.SetPoints(points);
        self.polydata_v1.SetVerts(vertices);
        self.polydata_v1.GetPointData().SetNormals(normals_v1);
        self.polydata_v1.GetPointData().SetScalars(colors_v1);
        
        # Fill circular component polydata
        self.polydata_v2.SetPoints(points);
        self.polydata_v2.SetVerts(vertices);
        self.polydata_v2.GetPointData().SetNormals(normals_v2);
        self.polydata_v2.GetPointData().SetScalars(colors_v2);
        
        # Fill longitudinal component polydata
        self.polydata_v3.SetPoints(points);
        self.polydata_v3.SetVerts(vertices);
        self.polydata_v3.GetPointData().SetNormals(normals_v3);
        self.polydata_v3.GetPointData().SetScalars(colors_v3);
        
        self.CreateVisualization()

    def SetMRIImage(self, matImage):
        ###################################################
        ### Volume plane cuts
        ###################################################
    
        VTK_DATA = numpy_support.numpy_to_vtk(num_array=matImage.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        image_volume = vtk.vtkImageData()
        image_volume.SetDimensions(matImage.shape[2], matImage.shape[1], matImage.shape[0]) #set dimensions as necessary
        image_volume.SetOrigin(0,0,0) #set origin as necessary
        image_volume.SetSpacing(1, 1, 1) #set spacing as necessary  
        image_volume.GetPointData().SetScalars(VTK_DATA)
    
        # Cell picker
        cellpicker = vtk.vtkCellPicker()
        cellpicker.SetTolerance(0.005)
        
        # Look-up table to color
        lut3D = vtk.vtkLookupTable()
        lut3D.SetTableRange(-6.283,6.283)	# image intensity range
        lut3D.SetAlphaRange(1.0,1.0)
        lut3D.SetValueRange(0,1)	# from black to white
        lut3D.SetSaturationRange(0,0)
        lut3D.SetRange(0,1)
        lut3D.Build()
        
        # Look-up table to color
        lut = vtk.vtkLookupTable()
        lut.SetTableRange(0,255)	# image intensity range
        lut.SetAlphaRange(1.0,1.0)
        lut.SetValueRange(0,1)	# from black to white
        lut.SetSaturationRange(0,0)
        lut.SetRange(0,1)
        lut.Build()
        
        # Short Axis image plane
        self.shortAxis = vtk.vtkImagePlaneWidget()
        self.shortAxis.DisplayTextOn()
        self.shortAxis.SetInputData(image_volume)
        self.shortAxis.SetPlaneOrientationToXAxes()
        self.shortAxis.SetSliceIndex(image_volume.GetExtent()[2]/2) # put the plane in the middle
        self.shortAxis.SetPicker(cellpicker)
        self.shortAxis.SetKeyPressActivationValue('s')
        self.shortAxis.SetRightButtonAction(1)
        self.shortAxis.SetLookupTable(lut3D)
        
        # Short Axis image plane
        self.longAxis = vtk.vtkImagePlaneWidget()
        self.longAxis.DisplayTextOn()
        self.longAxis.SetInputData(image_volume)
        self.longAxis.SetPlaneOrientationToZAxes()
        self.longAxis.SetSliceIndex(image_volume.GetExtent()[1]/2) # put the plane in the middle
        self.longAxis.SetPicker(cellpicker)
        self.longAxis.SetKeyPressActivationValue('l')
        self.longAxis.SetRightButtonAction(1)
        self.longAxis.SetLookupTable(lut)
        
        # Add planes to the scene (no visible by default)
        self.shortAxis.SetInteractor(self.iren);
        self.longAxis.SetInteractor(self.iren);
        
    def CreateVisualization(self):
        # Resample 
        # We don't really need one vector per voxel
        self.ptMask = vtk.vtkMaskPoints()
        self.ptMask.SetInputData(self.polydata_v1)
        self.ptMask.SetOnRatio(20)
        self.ptMask.RandomModeOff()
        self.ptMask.Update()
        
        # Glyph3D
        # We will represent vectors in 3D using arrows
        arrowSource = vtk.vtkArrowSource()
        
        
        glyph = vtk.vtkGlyph3D();
        glyph.SetInputConnection(self.ptMask.GetOutputPort())
        glyph.SetSourceConnection(arrowSource.GetOutputPort())
        glyph.ScalingOn()
        glyph.SetScaleFactor(3.0)	# Scale respect color value
        glyph.OrientOn()
        glyph.SetVectorModeToUseNormal()		# Orient the arrow using 'normals' 
        glyph.SetScaleModeToDataScalingOff()	# Normals (and colors) might not be normalized. Don't use them to scale.
        glyph.SetColorModeToColorByScalar()		# Color using the values in 'color'
        glyph.Update();

        # Mapper for glyphs
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(glyph.GetOutputPort())
        
        # Actor to represent glyphs
        actorVolume = vtk.vtkActor()
        actorVolume.SetMapper(mapper)
        
        # Add actor to renderer        
        self.renderer.AddActor(actorVolume)
        
        # Add axes to the scene
        self.axes = vtk.vtkAxesActor()
        self.axes.SetTotalLength(100,100,100)
        self.renderer.AddActor(self.axes)
        
    def CreateIntegrationVisualization(self):
        rk = vtk.vtkRungeKutta45()
        # Create source for streamtubes
        streamer = vtk.vtkStreamTracer()
        streamer.SetInputData(self.polydata_v1)
        streamer.SetSourceConnection(self.ptMask.GetOutputPort())
        streamer.SetMaximumPropagation(500)
        streamer.SetIntegrationStepUnit(2)
        streamer.SetMinimumIntegrationStep(0.972)
        streamer.SetMaximumIntegrationStep(20)
        streamer.SetInitialIntegrationStep(10)
        streamer.SetIntegrationDirection(0)
        streamer.SetIntegrator(rk)
        streamer.SetRotationScale(1.0)
        streamer.SetMaximumError(1.0)
        aa = vtk.vtkAssignAttribute()
        aa.SetInputConnection(streamer.GetOutputPort())
        aa.Assign("Normals","NORMALS","POINT_DATA")
        rf1 = vtk.vtkRibbonFilter()
        rf1.SetInputConnection(aa.GetOutputPort())
        rf1.SetWidth(0.1)
        rf1.VaryWidthOff()
        mapStream = vtk.vtkPolyDataMapper()
        mapStream.SetInputConnection(rf1.GetOutputPort())
        #mapStream.SetScalarRange(-1,1)
        streamActor = vtk.vtkActor()
        streamActor.SetMapper(mapStream)
        self.renderer.AddActor(streamActor)
        
    def Start(self):
        self.renderer.SetBackground(0.0,0.0,0.0)
        self.renderWindow.SetSize(800,800)

        # Finally, execute pipeline and render
        self.iren.Initialize()
        
        self.renderer.ResetCamera()
        self.renderer.GetActiveCamera().Zoom(-1.5)
        self.renderWindow.Render()
        
        self.iren.Start() 

class VTKDTIVolumeVisualizer:
    def __init__(self):
        """ Class constructor. 
        Inputs:
            * mri_volume: RxCxDx1 volume with 1 component per voxel corresponding to MRI acquisition 
            * dti_volume: RxCxDx3 volume with 3 components per voxel corresponding to vector information
            * mask: RxCxDx1 volume with boolean values to mask the input volumes
        """
    
    
        # Create renderer, render window, interactor and assign TrackBallCamera as
        # interaction style
        self.renderer = vtk.vtkRenderer()
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.AddRenderer(self.renderer)
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self.iren.SetRenderWindow(self.renderWindow)
        self.iren.AddObserver("KeyPressEvent", self.Keypress)
           
        # Add axes to the scene
        self.axes = vtk.vtkAxesActor()
        self.axes.SetTotalLength(100,100,100)
        self.renderer.AddActor(self.axes)

        # Create polydata to hold vector info
        self.image_v1 = vtk.vtkImageData()
        self.image_v2 = vtk.vtkImageData()
        self.image_v3 = vtk.vtkImageData()    
        
        self.planes = vtk.vtkPlanes()

    def Keypress(self, obj, event):
        key = obj.GetKeySym()
        print(key)
        if key == "q":
            obj.InvokeEvent("DeleteAllObjects")
            sys.exit()
        elif key == "1":
            self.mapper.SetInputData(self.image_v1)
    
        elif key =="2":
            self.mapper.SetInputData(self.image_v2)
            
        elif key =="3":
            self.mapper.SetInputData(self.image_v3)

        elif key == "i":
            self.boxWidget.SetEnabled(not(self.boxWidget.GetEnabled()))

        self.iren.Render()

    # When interaction starts, the requested frame rate is increased.
    def StartInteraction(self,obj, event):
        self.renderWindow.SetDesiredUpdateRate(10)

    # When interaction ends, the requested frame rate is decreased to
    # normal levels. This causes a full resolution render to occur.
    def EndInteraction(self, obj, event):
        self.renderWindow.SetDesiredUpdateRate(0.001)

    # The implicit function vtkPlanes is used in conjunction with the
    # volume ray cast mapper to limit which portion of the volume is
    # volume rendered.
    
    def ClipVolumeRender(self, obj, event):
        obj.GetPlanes(self.planes)
        self.mapper.SetClippingPlanes(self.planes)

    def SetDTIImages(self, dti_image_v1, dti_image_v2, dti_image_v3, mask):
        
        self.image_v1.SetDimensions(mask.shape[2], mask.shape[1], mask.shape[0]) #set dimensions as necessary
        self.image_v1.SetOrigin(0,0,0) #set origin as necessary
        self.image_v1.SetSpacing(1, 1, 1) #set spacing as necessary
        self.image_v1.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 4);

        self.image_v2 = vtk.vtkImageData()
        self.image_v2.SetDimensions(mask.shape[2], mask.shape[1], mask.shape[0]) #set dimensions as necessary
        self.image_v2.SetOrigin(0,0,0) #set origin as necessary
        self.image_v2.SetSpacing(1, 1, 1) #set spacing as necessary
        self.image_v2.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 4);

        self.image_v3 = vtk.vtkImageData()
        self.image_v3.SetDimensions(mask.shape[2], mask.shape[1], mask.shape[0]) #set dimensions as necessary
        self.image_v3.SetOrigin(0,0,0) #set origin as necessary
        self.image_v3.SetSpacing(1, 1, 1) #set spacing as necessary
        self.image_v3.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 4);

        print('Copying data')
        for i in range(mask.shape[2]):
            for j in range(mask.shape[1]):
                for k in range(mask.shape[0]):
                    self.image_v1.SetScalarComponentFromFloat(i,j,k,0, dti_image_v1[k,j,i,0])
                    self.image_v1.SetScalarComponentFromFloat(i,j,k,1, dti_image_v1[k,j,i,1])
                    self.image_v1.SetScalarComponentFromFloat(i,j,k,2, dti_image_v1[k,j,i,2])
                    self.image_v1.SetScalarComponentFromFloat(i,j,k,3, mask[k,j,i])

                    self.image_v2.SetScalarComponentFromFloat(i,j,k,0, dti_image_v2[k,j,i,0])
                    self.image_v2.SetScalarComponentFromFloat(i,j,k,1, dti_image_v2[k,j,i,1])
                    self.image_v2.SetScalarComponentFromFloat(i,j,k,2, dti_image_v2[k,j,i,2])
                    self.image_v2.SetScalarComponentFromFloat(i,j,k,3, mask[k,j,i])

                    self.image_v3.SetScalarComponentFromFloat(i,j,k,0, dti_image_v3[k,j,i,0])
                    self.image_v3.SetScalarComponentFromFloat(i,j,k,1, dti_image_v3[k,j,i,1])
                    self.image_v3.SetScalarComponentFromFloat(i,j,k,2, dti_image_v3[k,j,i,2])
                    self.image_v3.SetScalarComponentFromFloat(i,j,k,3, mask[k,j,i])
        print('Done')   
        
        
        self.CreateVisualization()

    def SetMRIImage(self, matImage):
        ###################################################
        ### Volume plane cuts
        ###################################################
    
        VTK_DATA = numpy_support.numpy_to_vtk(num_array=matImage.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
        image_volume = vtk.vtkImageData()
        image_volume.SetDimensions(matImage.shape[2], matImage.shape[1], matImage.shape[0]) #set dimensions as necessary
        image_volume.SetOrigin(0,0,0) #set origin as necessary
        image_volume.SetSpacing(1, 1, 1) #set spacing as necessary  
        image_volume.GetPointData().SetScalars(VTK_DATA)
    
        # Cell picker
        cellpicker = vtk.vtkCellPicker()
        cellpicker.SetTolerance(0.005)
        
        # Look-up table to color
        lut3D = vtk.vtkLookupTable()
        lut3D.SetTableRange(-6.283,6.283)   # image intensity range
        lut3D.SetAlphaRange(1.0,1.0)
        lut3D.SetValueRange(0,1)    # from black to white
        lut3D.SetSaturationRange(0,0)
        lut3D.SetRange(0,1)
        lut3D.Build()
        
        # Look-up table to color
        lut = vtk.vtkLookupTable()
        lut.SetTableRange(0,255)    # image intensity range
        lut.SetAlphaRange(1.0,1.0)
        lut.SetValueRange(0,1)  # from black to white
        lut.SetSaturationRange(0,0)
        lut.SetRange(0,1)
        lut.Build()
        
        # Short Axis image plane
        self.shortAxis = vtk.vtkImagePlaneWidget()
        self.shortAxis.DisplayTextOn()
        self.shortAxis.SetInputData(image_volume)
        self.shortAxis.SetPlaneOrientationToXAxes()
        self.shortAxis.SetSliceIndex(image_volume.GetExtent()[2]/2) # put the plane in the middle
        self.shortAxis.SetPicker(cellpicker)
        self.shortAxis.SetKeyPressActivationValue('s')
        self.shortAxis.SetRightButtonAction(1)
        self.shortAxis.SetLookupTable(lut3D)
        
        # Short Axis image plane
        self.longAxis = vtk.vtkImagePlaneWidget()
        self.longAxis.DisplayTextOn()
        self.longAxis.SetInputData(image_volume)
        self.longAxis.SetPlaneOrientationToZAxes()
        self.longAxis.SetSliceIndex(image_volume.GetExtent()[1]/2) # put the plane in the middle
        self.longAxis.SetPicker(cellpicker)
        self.longAxis.SetKeyPressActivationValue('l')
        self.longAxis.SetRightButtonAction(1)
        self.longAxis.SetLookupTable(lut)
        
        # Add planes to the scene (no visible by default)
        self.shortAxis.SetInteractor(self.iren);
        self.longAxis.SetInteractor(self.iren);
        
    def CreateVisualization(self):
        
        # Alpha channel function
        alphaChannelFunc = vtk.vtkPiecewiseFunction()
        alphaChannelFunc.AddPoint(0, 0.0)
        alphaChannelFunc.AddPoint(1.0, 0.8)
        alphaChannelFunc.AddPoint(255.0, 1.0)

        funcOpacityGradient = vtk.vtkPiecewiseFunction()
        funcOpacityGradient.AddPoint(0.0, 0.0)
        #alphaChannelFunc.AddPoint(90.01, 0.1)
        funcOpacityGradient.AddPoint(255.0, 0.0)


        # Color transfer functions
        color1 = vtk.vtkColorTransferFunction()
        color1.AddRGBPoint(0.0, 0.0, 0.0, 0.0);
        color1.AddRGBPoint(255.0, 1.0, 0.0, 0.0);

        color2 = vtk.vtkColorTransferFunction()
        color2.AddRGBPoint(0.0, 0.0, 0.0, 0.0);
        color2.AddRGBPoint(255.0, 0.0, 1.0, 0.0);

        color3 = vtk.vtkColorTransferFunction()
        color3.AddRGBPoint(0.0, 0.0, 0.0, 0.0);
        color3.AddRGBPoint(255.0, 0.0, 0.0, 1.0);

        # Volume rendering property
        propVolume = vtk.vtkVolumeProperty()
        propVolume.ShadeOff()
        propVolume.SetInterpolationTypeToLinear()
        propVolume.IndependentComponentsOff()
        #propVolume.SetGradientOpacity(funcOpacityGradient)
        propVolume.SetScalarOpacity(0,alphaChannelFunc)
        propVolume.SetScalarOpacity(1,alphaChannelFunc)
        propVolume.SetScalarOpacity(2,alphaChannelFunc)
        propVolume.SetScalarOpacity(3,alphaChannelFunc)
        propVolume.SetColor(0, color1)
        propVolume.SetColor(1, color2)
        propVolume.SetColor(2, color3)
        propVolume.SetColor(3, color3)

        self.mapper = vtk.vtkOpenGLGPUVolumeRayCastMapper()
        self.mapper.SetInputData(self.image_v1)
        #############################################

        newvol = vtk.vtkVolume()
        newvol.SetMapper(self.mapper)
        newvol.SetProperty(propVolume)
        
        # Add actor to renderer        
        self.renderer.AddActor(newvol)
        
        # Box widget
        self.boxWidget = vtk.vtkBoxWidget()
        self.boxWidget.SetInteractor(self.iren)
        self.boxWidget.SetPlaceFactor(1.0)
        self.boxWidget.RotationEnabledOff()
        self.boxWidget.ScalingEnabledOff()
        # Place the interactor initially. The output of the reader is used to
        # place the box widget.
        self.boxWidget.SetInputData(self.image_v1)
        self.boxWidget.PlaceWidget()
        self.boxWidget.InsideOutOn()
        self.boxWidget.On()
        self.boxWidget.AddObserver("StartInteractionEvent", self.StartInteraction)
        self.boxWidget.AddObserver("InteractionEvent", self.ClipVolumeRender)
        self.boxWidget.AddObserver("EndInteractionEvent", self.EndInteraction)

        outlineProperty = self.boxWidget.GetOutlineProperty()
        outlineProperty.SetRepresentationToWireframe()
        outlineProperty.SetAmbient(1.0)
        outlineProperty.SetAmbientColor(1, 1, 1)
        outlineProperty.SetLineWidth(3)

        selectedOutlineProperty = self.boxWidget.GetSelectedOutlineProperty()
        selectedOutlineProperty.SetRepresentationToWireframe()
        selectedOutlineProperty.SetAmbient(1.0)
        selectedOutlineProperty.SetAmbientColor(1,0, 0)
        selectedOutlineProperty.SetLineWidth(3)
        
    def CreateIntegrationVisualization(self):
        rk = vtk.vtkRungeKutta45()
        # Create source for streamtubes
        streamer = vtk.vtkStreamTracer()
        streamer.SetInputData(self.polydata_v1)
        streamer.SetSourceConnection(self.ptMask.GetOutputPort())
        streamer.SetMaximumPropagation(500)
        streamer.SetIntegrationStepUnit(2)
        streamer.SetMinimumIntegrationStep(0.972)
        streamer.SetMaximumIntegrationStep(20)
        streamer.SetInitialIntegrationStep(10)
        streamer.SetIntegrationDirection(0)
        streamer.SetIntegrator(rk)
        streamer.SetRotationScale(1.0)
        streamer.SetMaximumError(1.0)
        aa = vtk.vtkAssignAttribute()
        aa.SetInputConnection(streamer.GetOutputPort())
        aa.Assign("Normals","NORMALS","POINT_DATA")
        rf1 = vtk.vtkRibbonFilter()
        rf1.SetInputConnection(aa.GetOutputPort())
        rf1.SetWidth(0.1)
        rf1.VaryWidthOff()
        mapStream = vtk.vtkPolyDataMapper()
        mapStream.SetInputConnection(rf1.GetOutputPort())
        #mapStream.SetScalarRange(-1,1)
        streamActor = vtk.vtkActor()
        streamActor.SetMapper(mapStream)
        self.renderer.AddActor(streamActor)
        
    def Start(self):
        self.renderer.SetBackground(0.0,0.0,0.0)
        self.renderWindow.SetSize(800,800)

        # Finally, execute pipeline and render
        self.iren.Initialize()
        
        self.renderer.ResetCamera()
        self.renderer.GetActiveCamera().Zoom(-1.5)
        self.renderWindow.Render()
        
        self.iren.Start() 