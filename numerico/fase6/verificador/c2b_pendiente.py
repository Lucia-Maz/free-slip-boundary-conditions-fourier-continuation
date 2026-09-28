"""C2b: N=3 spurious pair with a frozen slope s0: Re(lam) vs the estimate nu*(2 s0^2/eta0^2 - k^2)."""
import math, sys
import numpy as np
sys.dont_write_bytecode = True
import c2_continuo as C
print(" N  k     eta0    s0     Re(lam) root   estimate nu(2s0^2/eta0^2-k^2)   Im(lam)")
for (k, eta0, s0) in [(1.414, -0.01, 0.02), (1.414, -0.01, 0.05), (1.414, -0.01, -0.05), (1.414, -0.01, 0.1),
                      (1.414, -0.02, 0.05), (1.414, -0.04, 0.05), (3.0, -0.01, 0.05), (5.0, -0.01, 0.05),
                      (8.0, -0.01, 0.05), (1.414, -0.005, 0.01), (1.414, -0.04, 0.02)]:
    for sgn in (1, -1):
        mm = (1 - sgn * 1j) / abs(eta0)
        lam0 = C.NU * (mm * mm - k * k)
        lam, q = C.newton(lam0, k, eta0, s0, 3)
        est = C.NU * (2 * s0**2 / eta0**2 - k * k)
        print(" 3 %5.2f %+.3f %+.3f   %+10.5f        %+10.5f                  %+9.3f" % (k, eta0, s0, lam.real, est, lam.imag))
# N=2 with slope: does the unstable root survive for all these?
for (k, eta0, s0) in [(1.414, -0.01, 0.1), (10.0, -0.01, 0.1), (1.414, -0.04, 0.1)]:
    lam, q = C.newton(C.NU * (1 / eta0**2 - k * k), k, eta0, s0, 2)
    print(" 2 %5.2f %+.3f %+.3f   Re lam %+.4f  (nu/eta0^2 - nu k^2 = %.4f)" % (k, eta0, s0, lam.real, C.NU * (1 / eta0**2 - k * k)))
