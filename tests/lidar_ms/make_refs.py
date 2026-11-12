"""Generate the frozen multiscatter 1.2.11 references used by the lidar_ms regression tests.

Usage: python make_refs.py <multiscatter executable> <multiscatter examples dir> <output dir>
Build note: link the executable with -lm after -lmultiscatter (the shipped Makefile has them reversed).
"""
import sys, subprocess, shutil, os
import numpy as np

exe, src_examples, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
TAGS = {
    "single":   ["single", "none"],
    "pvc":      ["fast", "none"],
    "pvc_hsrl": ["fast", "none", "-hsrl"],
    "h06":      ["original", "none"],
    "explicit2": ["explicit", "none", "-explicit", "2"],
    "pvc_tdts": ["fast", "lidar"],
    "tdts_only": ["none", "lidar"],
    "single_tdts_nolobe": ["single", "radar"],
    "tdts_only_nolobe": ["none", "radar"],
}
# HALO geometry: divergence 0.8 mrad (full angle, 1/e^2), field of view 1 mrad (full angle), instrument at 10.8 km
lam, alt, rt, rr = 532e-9, 10800.0, 0.4e-3/np.sqrt(2), 0.5e-3

def write_in(name, z, ext, rad, S, air, fd, fi, comment):
    with open(os.path.join(out, name + ".in"), "w") as f:
        f.write(f"# lidar_ms regression input: {comment}; HALO geometry\n")
        f.write(f"{len(z)} {lam} {alt} {float(rt):.17g} {rr}\n")
        for k in range(len(z)):
            f.write(f"{z[k]:.17g} {ext[k]:.17g} {rad[k]:.17g} {S[k]:.17g} {air[k]:.17g} 1 0.85 1 {fd[k]:.17g} {fi[k]:.17g}\n")

z = 1800.0 - 15.0*np.arange(120); air = 1.3e-5*np.exp(-z/8000.0); tiny = np.full(120, 1e-9)
cl = (z <= 1500) & (z > 1200); one = np.ones(120); zero = np.zeros(120)
write_in("halo_liquid_30", z, 0.03*cl, 8e-6*one, 18.8*one, air, 1.0*cl, zero, "liquid 30 km-1, a = 8 um, droplet near-180 function")
write_in("halo_liquid_10", z, 0.01*cl, 5e-6*one, 18.8*one, air, 1.0*cl, zero, "liquid 10 km-1, a = 5 um, droplet near-180 function")
write_in("halo_liquid_30_iso", z, 0.03*cl, 8e-6*one, 18.8*one, tiny, zero, zero, "liquid 30 km-1, a = 8 um, isotropic near 180, negligible air")
write_in("halo_liquid_10_iso", z, 0.01*cl, 5e-6*one, 18.8*one, tiny, zero, zero, "liquid 10 km-1, a = 5 um, isotropic near 180, negligible air")
write_in("halo_thin_slab", z, 1e-4*cl, 8e-6*one, 18.8*one, np.full(120, 1e-12), zero, zero, "one-way optical depth 0.03, isotropic near 180")
z2 = 5000.0 - 15.0*np.arange(300); air2 = 1.3e-5*np.exp(-z2/8000.0)
ice = (z2 <= 4500) & (z2 > 2500); liq = (z2 <= 2500) & (z2 > 2200)
write_in("halo_ice_over_liquid", z2, 1e-3*ice + 0.03*liq, np.where(ice, 50e-6, 8e-6), np.where(ice, 30.0, 18.8), air2,
         1.0*liq, np.zeros(300), "ice of optical depth 2 (a = 50 um, isotropic near 180) above liquid 30 km-1")

SHIPPED = ["i3rc_lidar_exp1", "i3rc_lidar_exp2", "i3rc_lidar_exp3", "i3rc_lidar_exp5", "i3rc_lidar_exp6", "hogan2008fig8",
           "calipso_cirrus_aerosol_wa_example", "ground_based_cirrus_example", "benchmark_profile50_aniso"]
for f in SHIPPED:
    shutil.copy(os.path.join(src_examples, f + ".in"), os.path.join(out, f + ".in"))

PLAN = {
    "halo_liquid_30": ["single", "pvc", "pvc_hsrl", "h06"],
    "halo_liquid_10": ["single", "pvc", "pvc_hsrl", "h06"],
    "halo_liquid_30_iso": ["single", "pvc", "h06"],
    "halo_liquid_10_iso": ["single", "pvc", "h06"],
    "halo_thin_slab": ["single", "pvc", "explicit2"],
    "halo_ice_over_liquid": ["single", "pvc", "pvc_hsrl"],
    "i3rc_lidar_exp1": ["single", "single_tdts_nolobe", "tdts_only_nolobe"],
    "i3rc_lidar_exp2": ["single", "single_tdts_nolobe", "tdts_only_nolobe"],
    "i3rc_lidar_exp3": ["single", "single_tdts_nolobe", "tdts_only_nolobe"],
    "i3rc_lidar_exp5": ["single", "pvc", "pvc_tdts", "tdts_only"],
    "i3rc_lidar_exp6": ["single", "pvc", "pvc_tdts", "tdts_only"],
    "hogan2008fig8": ["single", "pvc", "pvc_tdts", "tdts_only"],
    "calipso_cirrus_aerosol_wa_example": ["pvc", "pvc_tdts", "tdts_only"],
    "ground_based_cirrus_example": ["single", "pvc"],
    "benchmark_profile50_aniso": ["pvc", "pvc_hsrl", "pvc_tdts", "tdts_only"],
}
n = 0
for name, tags in PLAN.items():
    for tag in tags:
        opts = ["-algorithms", *TAGS[tag]]
        p = subprocess.run([exe, "-quiet", *opts, name + ".in"], capture_output=True, text=True, timeout=60, cwd=out)
        with open(os.path.join(out, f"{name}__{tag}.ref"), "w") as f:
            f.write(f"# multiscatter 1.2.11 (Apache-2.0), options: {' '.join(opts)}\n")
            f.write("# columns: gate, apparent range above ground (m), extinction (m-1), radius (um), "
                    "apparent backscatter (m-1 sr-1)[, air backscatter (m-1 sr-1) with -hsrl]\n")
            f.write(p.stdout)
        n += 1
print(n, "reference files written to", out)
