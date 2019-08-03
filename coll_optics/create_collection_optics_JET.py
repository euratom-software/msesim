import xml.etree.ElementTree
import numpy as np
from xml.dom import minidom

import matplotlib.pyplot as plt

def create_file(x_length, y_length, n_channels):

    grid_x, grid_y = np.meshgrid(x_length, y_length)

    names = []

    for i in range(n_channels-1):
        names.append('ch' + str([i+1][0]))

    xml.etree.ElementTree.SubElement(tree.find('coll'), 'bundleID',
                                         attrib={'description': '', 'type': 'string', 'value':', '.join(names)})

    channels = np.arange(1,26,1)

    for i in range(np.shape(grid_x)[1]):
        xml.etree.ElementTree.SubElement(tree.find('coll'), 'ch'+str(channels[i])+'l', attrib={'description':'', 'type':'str', 'value':','.join(map(str, grid_x[:,i]))})
        xml.etree.ElementTree.SubElement(tree.find('coll'), 'ch'+str(channels[i])+'m', attrib={'description':'', 'type':'str', 'value':','.join(map(str, grid_y[:,i]))})

    xmlstr = minidom.parseString(xml.etree.ElementTree.tostring(root)).toprettyxml(indent="   ")
    with open(fname, "w") as f:
        f.write(xmlstr)

#width and height of whole fiber array in mm
fiber_height = 36
fiber_width =  7.8
n_channels = 25
n_fibers = 6 #number fibers per channel
fiber_radius_x = 1.5/2 #slightly wider distance between the fibers to account for the struts between them
fiber_radius_y = 1.3/2

x_length = np.linspace((-1*fiber_height/2)+fiber_radius_x,(fiber_height/2)-fiber_radius_x, n_channels)
y_length = np.linspace( (-1*fiber_width/2)+fiber_radius_y, (fiber_width/2)-fiber_radius_y, n_fibers)

# Open the collection optics xml file template:
tree = xml.etree.ElementTree.parse('xml_template.xml')
root = tree.getroot()

fname = "JET_MSE_collection_optics.xml"


