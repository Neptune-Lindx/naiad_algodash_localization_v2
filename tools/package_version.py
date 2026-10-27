#!/usr/bin/env python3

import sys
import xml.etree.ElementTree as XTree

pkg_file = '../package.xml'

if len(sys.argv) > 1:
    pkg_file = sys.argv[1]

try:
    tree = XTree.parse(pkg_file)
    root = tree.getroot()
    item = root.find('version')    
    if item != None:
        print(item.text, end='')
        sys.exit(0)
    else:
        print("***Tag 'version' not found", file=sys.stderr)
except Exception as e :
    print("***XML parsed failed:", e, file=sys.stderr)

sys.exit(1)
