import xml.etree.ElementTree as ET
import os

"""
Splits chemsep1 and chemsep2 into
individual compound files.
"""

def build(filename):
    # Parse the XML file
    tree = ET.parse(filename)
    root = tree.getroot()
    for child in root:
        child_elem = child.find("CompoundID")
        child_name = child_elem.get("value").lower()
        child_elem.set('value', child_name)
        if child_name is not None:
            tree = ET.ElementTree(child)
            tree.write(f"data/data_files/{child_name}.xml")

def run():
    # Remove all files in the directory
    data_directory = "data/data_files"
    for fname in os.listdir(data_directory):
        file_path = os.path.join(data_directory, fname)
        os.remove(file_path)
    
    # Build from XML files
    build("data/chemsep/chemsep1.xml")
    build("data/chemsep/chemsep2.xml")

run()
