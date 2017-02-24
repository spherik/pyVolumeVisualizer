#!/usr/bin/env python

import os
import numpy
import vtk
import scipy.io as sio
from vtk.util import numpy_support
import sys
import VTKDTIVisualizer

matFilename = '/Users/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/DTI.mat'
#matFilename = '/home/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/DTI.mat'
#matFilename = 'LocalReferenceFrame.mat'
mat_contents = sio.loadmat(matFilename)
vectorField = mat_contents['VectorField']
#vectorField[:,:,:,:,0] = -vectorField[:,:,:,:,0]
#vectorField[:,:,:,:,2] = -vectorField[:,:,:,:,2]
radial_vectors = numpy.multiply(255,numpy.absolute(vectorField[:,:,:,0,:])).copy(order='C')
circ_vectors = numpy.multiply(255,numpy.absolute(vectorField[:,:,:,1,:])).copy(order='C')
long_vectors = numpy.multiply(255,numpy.absolute(vectorField[:,:,:,2,:])).copy(order='C')

print('Fliping vectors')
# Swap x and z components
for i in range(radial_vectors.shape[0]):
    for j in range(radial_vectors.shape[1]):
        for k in range(radial_vectors.shape[2]):
            radial_vectors[i,j,k,:] = numpy.flipud(radial_vectors[i,j,k,:])
            circ_vectors[i,j,k,:] = numpy.flipud(circ_vectors[i,j,k,:])
            long_vectors[i,j,k,:] = numpy.flipud(long_vectors[i,j,k,:])


# Load mask
matFilename = '/Users/spherik/ownCloudZoidberg/Hearts/Datasets/Oxford/Toni/mask.mat'
mat_contents = sio.loadmat(matFilename)
mask = numpy.multiply(255,mat_contents['mask2']).copy(order='C')

mat_contents_shape = radial_vectors.shape

visualizer = VTKDTIVisualizer.VTKDTIVolumeVisualizer()
visualizer.SetDTIImages(radial_vectors,circ_vectors,long_vectors,mask)
visualizer.Start()