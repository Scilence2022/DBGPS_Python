from setuptools import setup
from Cython.Build import cythonize



setup(
    name='DBGPS',
    ext_modules=cythonize(["utils.py", "crc16pure.py", "deBruijnGraph.py", "glass.py", "DNAdroplet.py", "droplet.py", "fountain.py", "DNAfountain.py"]),

)