#### LOG FOR LAMMPS BRANCH OF aiGCMC PROJECT
## 10/ 15/ 20
Started up my Xenial Puppy Linux machine and downloaded a tarball for LAMMPS. I am using cmake to configure with the following flags:
-C ../cmake/presets/minimal.cmake -D BUILD_SHARED_LIBS=on -D LAMMPS_EXCEPTIONS=on \
-D PKG_PYTHON=on -D PKG_MC=on -D PKG_REPLICA=on \
-DPYTHON_LIBRARY=/root/miniconda3/lib/libpython3.7m.so \
-DPYTHON_INCLUDE_DIR=/root/miniconda3/include/python3.7m/ \
-DPYTHON_EXECUTABLE=/root/miniconda3/bin/python

cmake --install . doesn't work for me, make install does

/////////////////////////////////////////////////////////

None of that worked. It wouldn't create a shared library liblammps.so. Of course, I could create a static library liblammps.a in a 
static cmake build. Inspired by Axel, I decide to just do a traditional make build. To make sure the python install works, I replace
the Makefile.lammps with the Makefile.lammps.python3 to have it search for python3 since the Puppy also has python 2.7. 
make yes-mc
make mode=shlib serial
make install-python

I ran the python interpreter
from lammps import lammps
l = lammps()
l.command("info config")

It all worked! Next I will try to run the mc.py code and in.mc to test it all together in a real script. Unfortunately putting it on mode-shlibs
means that it does not create an executable, so I ran it again without that mode. To disable track pad scrolling in bash
xinput set-prop 6 "Device Enabled" 0

## 10/ 19/ 20
So I tried the mc.py python code with the preset in.mc, which perturbs a hexagonal lattice. This works well, so I took my old in.sic script
and got rid of the extra stuff to make it work as the input for mc.py. I need to have MANYBODY installed to use the tersoff potential.
Also this time I added
make yes-manybody
make yes-python
It seems that lammps and python work even if the python package isn't installed, so I'll see if installing it affects anything even.
Ok the installation didn't work with 'make yes-python'. Reverted it to no, and it works
I have decided that having 'mc.uvt_new_structure()' and 'mc.uvt_new_structure()' is confusing and frustrating to switch all calls of the 
method from one to another. I think it is simpler to have a child class MC_NP which overrides the method, so the calls are the same.