import os
import numpy
import vtk
import scipy.io as sio
from vtk.util import numpy_support
import sys
import VTKDTIVisualizer
visualizer = VTKDTIVisualizer.VTKDTIVectorVisualizer()

##############################################        
# MRI volume
##############################################        
matFilename = '/Users/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/I.mat'
#matFilename = '/home/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/I.mat'
#matFilename = 'I.mat'
mat_contents = sio.loadmat(matFilename)
matImage = mat_contents['I'].copy('C')
visualizer.SetMRIImage(matImage)
##############################################        
# DTI volume
##############################################        
matFilename = '/Users/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/LocalReferenceFrame.mat'
#matFilename = '/home/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/LocalReferenceFrame.mat'
#matFilename = 'LocalReferenceFrame.mat'
mat_contents = sio.loadmat(matFilename)
radial_vectors = mat_contents['radial_vectors'].copy(order='C')
circ_vectors = mat_contents['circ_vectors'].copy(order='C') 
long_vectors = mat_contents['long_vectors'].copy(order='C')

############################################## 
# Mask volume
##############################################        
matFilename = '/Users/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/mask.mat'
#matFilename = '/home/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/mask.mat'
#matFilename = 'mask.mat'
mat_contents = sio.loadmat(matFilename)
mask = mat_contents['mask2'].copy(order='C')

visualizer.SetDTIPolydatas(radial_vectors, circ_vectors, long_vectors, mask)
visualizer.Start()