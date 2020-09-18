#__version__ = "0.0.3"
# Update version from VERSION file into module
with open('VERSION', 'r') as fversion:
    __version__ = fversion.readline().rstrip()
