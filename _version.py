#__version__ = "0.0.3"
# Update version from VERSION file into module
import os
with open(os.path.dirname(__file__)+'/VERSION', 'r') as fversion:
    __version__ = fversion.readline().rstrip()
