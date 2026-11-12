"""PVC check script from the 2026-11-11 design session: compat (multiscatter 1.2.11) and native modes. Reference for milestone M4; not project code."""
import sys, subprocess
import numpy as np
from scipy.special import erfcx

U1, U2, V1SQ, V2SQ = 0.3, 0.5, 16.0, 0.16
YW, YG0 = 0.89, 0.038

def drange(i, n, r):
    if i <= 0: return abs(r[1]-r[0])
    if i >= n-1: return abs(r[n-1]-r[n-2])
    return abs(r[i+1]-r[i-1])*0.5

def coangle_var(s2, z2, c, R, s2max, native):
    g2 = z2 + s2/R**2 - 2*c/R
    if s2max > 0 and np.real(s2) > np.real(s2max):
        cg = c - s2/R
        rho2 = cg*cg/(s2*(g2 if native else z2))
        g2 = g2*(rho2*s2max/s2 + 1 - rho2)
    return g2

def F_drop(g2, Th2): return (1-U1-U2) + U1/(1+g2*V1SQ/Th2) + U2/(1+g2*V2SQ/Th2)
def F_ice(g2):
    g = np.sqrt(g2)*0.5/YG0
    return (1-YW) + YW*(1 - g*np.sqrt(np.pi)*erfcx(g))

def pvc(r, ext, rad, S, ext_air, fd, fi, lam, alt, rt, rr, native=False):
    n = len(r); rt2 = rt*rt
    k_air = 3/(8*np.pi) if lam <= 1e-6 else 0.0
    kappa = 1/(1-np.exp(-rr*rr/rt2))
    R0 = abs(alt - 1.5*r[0] + 0.5*r[1])
    s2 = R0*R0*rt2; z2 = rt2; c = R0*rt2; T = 1.0
    s2r, z2r, cr, Tr = s2, z2, c, T
    Tu = 1.0; Tu_prev = 1.0; M_prev = 0.0; Fp_prev = 1.0
    out = np.zeros((4, n), dtype=complex)   # particulate MS, air MS, particulate SS, air SS
    for i in range(n):
        dr = drange(i, n, r); Th2 = (lam/(np.pi*rad[i]))**2
        b = ext[i]/S[i]; ba = ext_air[i]*k_air
        od = 2*(ext[i]+ext_air[i])*dr
        R = abs(alt-r[i]) + 0.5*dr; s2max = (R*rr)**2
        s2u = R*R*rt2; cu = R*rt2
        Tu = Tu*np.exp(-od)
        s2 = s2 + z2*dr*dr + 2*c*dr + ext[i]*Th2*dr**3/3
        c = c + z2*dr + ext[i]*Th2*dr*dr/2
        z2 = z2 + ext[i]*Th2*dr
        T = T*np.exp(-(ext[i]+2*ext_air[i])*dr)
        s2r = s2r + z2r*dr*dr + 2*cr*dr + ext[i]*Th2*dr**3/3
        cr = cr + z2r*dr + ext[i]*Th2*dr*dr/2
        z2r = z2r + ext[i]*Th2*dr
        Tr = Tr*np.exp(-(ext[i]+2*ext_air[i])*dr)
        Ta = Tr - Tu
        if np.real(Ta) < 0: Ta = 0.0*Ta
        M = 0.0; Fp = 1.0; Fp_next = 1.0
        ss = Tu_prev*(1-np.exp(-od))/od if np.real(od) > 0 else Tu_prev
        if np.real(Ta) > np.real(Tr)*1e-12:
            s2a = (Tr*s2r - Tu*s2u)/Ta; z2a = (Tr*z2r - Tu*rt2)/Ta; ca = (Tr*cr - Tu*cu)/Ta
            Tb = T - Tr
            if np.real(s2a) > s2max:
                rho2 = ca*ca/(s2a*z2a); f = s2max/s2a
                s2a = s2max + 0*s2a; Ta = Ta*f; ca = ca*f; z2a = z2a*(rho2*f + 1 - rho2)
                Tr = Tu + Ta; s2r = (Ta*s2a + Tu*s2u)/Tr
                z2r = (Ta*z2a + Tu*rt2)/Tr; cr = (Ta*ca + Tu*cu)/Tr
            aniso = fd[i] > 0 or fi[i] > 0
            aniso_n = i < n-1 and (fd[i+1] > 0 or fi[i+1] > 0)
            if aniso or aniso_n:
                Tg = T - Tu
                s2g = (T*s2 - Tu*s2u)/Tg; z2g = (T*z2 - Tu*rt2)/Tg; cg = (T*c - Tu*cu)/Tg
                g2 = coangle_var(s2g, z2g, cg, R, s2max, native)
                def Fpart(j):
                    Thj = (lam/(np.pi*rad[j]))**2
                    return fd[j]*F_drop(g2, Thj) + fi[j]*F_ice(g2) + 1 - fd[j] - fi[j]
                if aniso: Fp = Fpart(i)
                if aniso_n: Fp_next = Fpart(i+1)
            if np.real(Tb) > 0:
                s2b = (T*s2 - Tr*s2r)/Tb
                if np.real(s2b) <= np.real(s2a): s2b = s2a
                M = ((1-np.exp(-s2max/s2a))*Ta + (1-np.exp(-s2max/s2b))*Tb)*kappa/Tu
            else:
                M = kappa*(1-np.exp(-s2max/s2a))*Ta/Tu
            def gate(Fm, Fp_):   # gate mean with (1+F M) exponential between half-levels
                Xm = 1 + Fm*M_prev; Xp = 1 + Fp_*M
                return Tu_prev*(Xm - Xp*np.exp(-od))/(od - np.log(Xp/Xm))
            if native:           # class-specific factors: particulate F_p, air F=1
                out[0, i] = gate(Fp_prev, Fp)*b; out[1, i] = gate(1.0, 1.0)*ba
            else:                # 1.2.11: backscatter-weighted mean factor applied to both
                def avg(F, j):
                    bj = ext[j]/S[j]; baj = ext_air[j]*k_air
                    return (bj*F + baj)/(bj + baj)
                fac = gate(avg(Fp_prev, i), avg(Fp, i))
                out[0, i] = fac*b; out[1, i] = fac*ba
        elif np.real(od) > 0:
            out[0, i] = ss*b; out[1, i] = ss*ba
        out[2, i] = ss*b; out[3, i] = ss*ba
        Tu_prev = Tu; M_prev = M; Fp_prev = Fp_next
    return out

def halo_case(alpha=0.03, a=8e-6, nz=120, dz=15.0, ztop=1500.0, thick=300.0, alt=10800.0):
    z = 1800.0 - dz*np.arange(nz)
    cloud = (z <= ztop) & (z > ztop-thick)
    ext = np.where(cloud, alpha, 0.0); rad = np.full(nz, a); S = np.full(nz, 18.8)
    ext_air = 1.3e-5*np.exp(-z/8000.0)
    fd = np.where(cloud, 1.0, 0.0); fi = np.zeros(nz)
    return z, ext, rad, S, ext_air, fd, fi, cloud

if __name__ == "__main__":
    exe = sys.argv[1]
    lam, alt, rt, rr = 532e-9, 10800.0, 0.4e-3/np.sqrt(2), 0.5e-3
    z, ext, rad, S, ext_air, fd, fi, cloud = halo_case()
    lines = [f"{len(z)} {lam} {alt} {float(rt):.17g} {rr}"]
    for k in range(len(z)):
        lines.append(f"{z[k]:.17g} {ext[k]:.17g} {rad[k]:.17g} {S[k]:.17g} {ext_air[k]:.17g} 1 0.85 1 {fd[k]:.17g} {fi[k]:.17g}")
    inp = "\n".join(lines) + "\n"
    def run(flags):
        o = subprocess.run([exe, "-quiet", "-algorithms", "fast", "none"] + flags, input=inp,
                           capture_output=True, text=True).stdout
        return np.array([[float(x) for x in l.split()] for l in o.strip().splitlines() if not l.startswith('#')])
    ref_h = run(["-hsrl"]); ref_t = run([])
    c = pvc(z, ext, rad, S, ext_air, fd, fi, lam, alt, rt, rr, native=False).real
    nat = pvc(z, ext, rad, S, ext_air, fd, fi, lam, alt, rt, rr, native=True).real
    rel = lambda a, b: np.max(np.abs(a-b)/np.maximum(np.abs(b), 1e-300))
    print("C output columns:", ref_h.shape, ref_t.shape)
    print("compat vs C, particulate (-hsrl): max rel diff = %.2e" % rel(c[0], ref_h[:, 4]))
    print("compat vs C, air        (-hsrl): max rel diff = %.2e" % rel(c[1], ref_h[:, 5]))
    print("compat vs C, total              : max rel diff = %.2e" % rel(c[0]+c[1], ref_t[:, 4]))
    ic = np.where(cloud)[0]
    print("native vs compat in cloud: particulate max rel = %.3e, air max rel = %.3e"
          % (rel(nat[0, ic], c[0, ic]), rel(nat[1, ic], c[1, ic])))
    # multiple-scattering enhancement and Hu-type depolarization (liquid class), native mode
    dr = np.array([drange(i, len(z), z) for i in range(len(z))])
    IT = np.cumsum(nat[0]*dr); IS = np.cumsum(nat[2]*dr)
    Iperp = 0.5*(IT - np.sqrt(IT*IS))
    bperp = np.diff(np.concatenate([[0.0], Iperp]))/dr
    bpar = nat[0] - bperp
    print("\n depth_m  tau   beta_MS/beta_SS  A_s   delta_acc  delta_bin(particulate)")
    for k in ic[:20:2]:
        depth = 1500.0 - z[k] + 7.5; tau = 0.03*depth
        As = IS[k]/IT[k]; dacc = (1-np.sqrt(As))/(1+np.sqrt(As))
        print(" %6.1f  %4.2f   %8.3f      %5.3f   %6.3f     %6.3f" % (depth, tau, nat[0, k]/nat[2, k], As, dacc, bperp[k]/bpar[k]))
    # complex-step directional derivative vs central finite difference (d/d ln ext in cloud)
    t = np.where(cloud, ext, 0.0)                    # tangent: uniform relative change of cloud extinction
    h = 1e-30
    cs = pvc(z, ext + 1j*h*t, rad, S, ext_air, fd, fi, lam, alt, rt, rr, native=True)[0].imag/h
    e = 1e-6
    fdiff = (pvc(z, ext*(1+e*cloud), rad, S, ext_air, fd, fi, lam, alt, rt, rr, native=True)[0].real
             - pvc(z, ext*(1-e*cloud), rad, S, ext_air, fd, fi, lam, alt, rt, rr, native=True)[0].real)/(2*e)
    print("\ncomplex-step vs central FD (d beta_hat / d ln alpha_cloud), cloud bins: max rel diff = %.2e"
          % rel(cs[ic], fdiff[ic]))
