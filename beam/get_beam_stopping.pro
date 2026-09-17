; Get the excited n=2 population                                                                                                                                                                            

files    = ['/home/adas/adas/adf22/bmp97#h/bmp97#h_2_h1.dat', $
            '/home/adas/adas/adf22/bmp97#h/bmp97#h_2_c6.dat' ]

energy = 63000 ;eV
te = 800 ;eV
dens = 2e19 ;m^-3 gets converted to cm^3 in the method
fraction = 1

read_adf21, files=files, energy=energy, te=te,  dens=dens/1e6,   $
            fraction=fraction,  data=n2pop
