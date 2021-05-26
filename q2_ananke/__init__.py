from ._version import get_versions
from .dereplicate import dereplicate
#from .compute import compute

__version__ = get_versions()['version']
del get_versions

__all__ = ['dereplicate']
