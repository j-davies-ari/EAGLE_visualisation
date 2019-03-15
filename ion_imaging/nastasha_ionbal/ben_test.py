import ionbal_benstripped as ionbal
import numpy as np

lT = np.linspace(3.0,6.0,31)

lnH = np.asarray(lT)*0-2.0

lT = np.asarray(lT)

lnH = np.linspace(-6.0,-2.0,41)

lT = np.asarray(lnH)*0+4.0

bal = ionbal.find_ionbal(3.0,"c4",lnH, lT)

print "bal= ", bal
print "logTK= ", lT
print "lognHcm3= ", lnH

