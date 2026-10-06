import xml.etree.ElementTree
import numpy as np
from xml.dom import minidom

# Open the collection optics xml file template:
tree = xml.etree.ElementTree.parse('xml_template.xml')
root = tree.getroot()

x_length = np.linspace(-8.3,8.3,32)
y_length = np.linspace(-7,7,32)+0.0675 #np.zeros((len(x_length)))+0.0675

fname = "DIIID_imse_2d_lenspos.xml"

grid_x, grid_y = np.meshgrid(x_length, y_length)

pixel_names = []
attrib_x = dict()
attrib_y = dict()
pixel_number = np.arange(0, len(x_length) ** 2 + 1, 1)
coords = []

for i in range(len(grid_x)):
    for j in range(len(grid_y)):
        coords.append(([grid_x[i,j], grid_y[i,j]]))

names = []

for i in range(len(coords)):
    names.append('p' + str([i][0]))

xml.etree.ElementTree.SubElement(tree.find('coll'), 'bundleID', attrib={'description': '', 'type': 'string', 'value': ', '.join(names)})

for i in range(len(coords)):
    xml.etree.ElementTree.SubElement(tree.find('coll'), 'p'+str([i][0])+'l', attrib={'description':'', 'type':'str', 'value':str(coords[i][0])})
    xml.etree.ElementTree.SubElement(tree.find('coll'), 'p'+str([i][0])+'m', attrib={'description':'', 'type':'str', 'value':str(coords[i][1])})


xmlstr = minidom.parseString(xml.etree.ElementTree.tostring(root)).toprettyxml(indent="   ")
with open(fname, "w") as f:
    f.write(xmlstr)
