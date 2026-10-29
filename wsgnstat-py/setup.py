"""Build the shared wsgnstat C back-end at INSTALL time, as an Extension.

Declaring the sources as an Extension hands the compile to pip, inside build isolation, and lets
real platform wheels be built, after which users need no compiler at all. Same pattern as
arcstat-py and exactcond-py.

init.c is excluded deliberately: it is R's symbol-registration glue and references R headers.
"""
import glob
import os

from setuptools import Extension, setup

# Paths must be RELATIVE to setup.py: setuptools rejects absolute source paths.
SRC = os.path.join("src", "wsgnstat")

sources = sorted(
    s for s in glob.glob(os.path.join(SRC, "*.c"))
    if os.path.basename(s) != "init.c"
)

# Leading underscore: it carries no PyInit and is loaded with ctypes.CDLL by _core.py.
ext = Extension(
    name="wsgnstat._libwsgnstat",
    sources=sources,
    include_dirs=[SRC],
    libraries=[] if os.name == "nt" else ["m"],
    optional=False,
)

setup(ext_modules=[ext])
