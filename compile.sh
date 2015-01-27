gfortran -c tdpn.f90
f2py --overwrite-signature -m sub -h sub.pyf tdpn.f90
f2py -c --fcompiler=gnu95 sub.pyf  tdpn.f90
