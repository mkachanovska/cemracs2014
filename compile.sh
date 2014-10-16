gfortran -c tdp.f90
f2py --overwrite-signature -m sub -h sub.pyf tdp.f90
f2py -c --fcompiler=gnu95 sub.pyf  tdp.f90
