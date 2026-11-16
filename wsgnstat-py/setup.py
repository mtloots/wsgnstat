"""Build the shared wsgnstat C back-end at INSTALL time, as an Extension.

Declaring the sources as an Extension hands the compile to pip, inside build isolation, and lets
real platform wheels be built, after which users need no compiler at all. Same pattern as
arcstat-py and exactcond-py.

init.c is excluded deliberately: it is R's symbol-registration glue and references R headers.

OpenMP. The back-end runs the independent trials inside a fit on several threads when it is
compiled with OpenMP (the ncores argument of fit_wsgn and profile_wsgn). The build asks for it
where the compiler supports it and falls back to a serial library, silently apart from one line of
output, when it does not, so the build never fails for lack of OpenMP. WSGN_OPENMP controls this:
unset or "auto" tries OpenMP on Linux and other Unix systems and stays serial on macOS and Windows;
"1" tries it everywhere it can (on macOS this needs Homebrew's libomp, found under
WSGN_LIBOMP_PREFIX, /opt/homebrew/opt/libomp or /usr/local/opt/libomp); "0" forces a serial
library. macOS is serial by default because a wheel that bundles libomp collides with the other
copy of libomp that scikit-learn, PyTorch and others load into the same process, and delocate
refuses a libomp built for a newer macOS than the wheel targets.
"""
import copy
import glob
import os
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.errors import CompileError, ExecError, LinkError, PlatformError

# Paths must be RELATIVE to setup.py: setuptools rejects absolute source paths.
SRC = os.path.join("src", "wsgnstat")

sources = sorted(
    s for s in glob.glob(os.path.join(SRC, "*.c"))
    if os.path.basename(s) != "init.c"
)


def _openmp_trials():
    """Candidate (compile args, link args, include dirs, library dirs) for OpenMP, best first."""
    mode = os.environ.get("WSGN_OPENMP", "auto").strip().lower()
    if mode in ("0", "no", "false", "off", "none"):
        return []
    forced = mode in ("1", "yes", "true", "on")
    if sys.platform == "win32":
        return []
    if sys.platform == "darwin":
        if not forced:
            return []
        prefixes = [os.environ.get("WSGN_LIBOMP_PREFIX", ""), "/opt/homebrew/opt/libomp",
                    "/usr/local/opt/libomp", os.environ.get("CONDA_PREFIX", "")]
        trials = []
        for p in prefixes:
            if p and os.path.isdir(os.path.join(p, "include")):
                trials.append((["-Xpreprocessor", "-fopenmp"], ["-lomp"],
                               [os.path.join(p, "include")], [os.path.join(p, "lib")]))
        return trials
    return [(["-fopenmp"], ["-fopenmp"], [], [])]


class BuildExtOpenMP(build_ext):
    """Try the extension with OpenMP first and fall back to a serial build."""

    def build_extension(self, ext):
        for cargs, largs, incs, libs in _openmp_trials():
            trial = copy.copy(ext)
            trial.include_dirs = list(ext.include_dirs) + incs
            trial.library_dirs = list(ext.library_dirs) + libs
            trial.extra_compile_args = list(ext.extra_compile_args) + cargs
            trial.extra_link_args = list(ext.extra_link_args) + largs
            try:
                build_ext.build_extension(self, trial)
            except (CompileError, LinkError, ExecError, PlatformError, OSError) as exc:
                print("wsgnstat: OpenMP build failed (%s); trying without." % (exc,))
                continue
            print("wsgnstat: built with OpenMP.")
            return
        print("wsgnstat: built without OpenMP (serial); ncores has no effect.")
        build_ext.build_extension(self, ext)


# Leading underscore: it carries no PyInit and is loaded with ctypes.CDLL by _core.py.
ext = Extension(
    name="wsgnstat._libwsgnstat",
    sources=sources,
    include_dirs=[SRC],
    libraries=[] if os.name == "nt" else ["m"],
    optional=False,
)

setup(ext_modules=[ext], cmdclass={"build_ext": BuildExtOpenMP})
