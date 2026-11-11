#!/usr/bin/env python3
"""Mean annual surface temperature on Mercury, as a lat/lon grid.

Mercury sits in a 3:2 spin-orbit resonance on an eccentric orbit (e=0.2056), so insolation
depends on longitude as well as latitude: the hot poles at 0 and 180 deg E face the Sun at
perihelion, the warm poles at 90 and 270 deg at aphelion. Obliquity (~2 arcmin) is ignored.

A 1-D heat equation is integrated over the 176-day solar day for every grid site at once.
After each day the whole column is shifted by the residual mean surface flux imbalance, which
drives the profile to radiative steady state in tens of days instead of the ~70 solar days
free relaxation would need. At steady state the deep temperature equals the mean surface
temperature, and that is the boundary condition the vapor-diffusion front sees at depth.
"""
import numpy as np

AU, S0, SIG = 1.495978707e11, 1361.0, 5.670374419e-8
A_BOND, EMIS = 0.09, 0.95
P_ORB, P_ROT = 87.9691, 58.6461
P_SOLAR = 1.0 / (1.0 / P_ROT - 1.0 / P_ORB)
ECC, SMA = 0.205630, 0.3870993 * AU
RHO, CP, K = 1600.0, 700.0, 0.002        # thermal inertia sqrt(k rho cp) = 47, lunar-like
KAPPA = K / (RHO * CP)
DSKIN = np.sqrt(KAPPA * P_SOLAR * 86400.0 / np.pi)

def _flux(lat, lon, nt):
    """Absorbed flux over one solar day. lat, lon in degrees, shape (nsite,)."""
    t = np.linspace(0, P_SOLAR, nt, endpoint=False)
    M = 2 * np.pi * t / P_ORB
    E = M.copy()
    for _ in range(80):
        E = E - (E - ECC * np.sin(E) - M) / (1 - ECC * np.cos(E))
    nu = 2 * np.arctan2(np.sqrt(1 + ECC) * np.sin(E / 2), np.sqrt(1 - ECC) * np.cos(E / 2))
    r = SMA * (1 - ECC * np.cos(E))
    lam = (2 * np.pi * t / P_ROT - nu) % (2 * np.pi)
    hour = np.deg2rad(lon)[:, None] - lam[None, :]
    cosi = np.cos(np.deg2rad(lat))[:, None] * np.cos(hour)
    F = (1 - A_BOND) * S0 * (AU / r)[None, :] ** 2 * np.clip(cosi, 0, None)
    return t * 86400.0, F

def mean_annual_T(lat, lon, nz=70, nt=900, ndays=45, verbose=False):
    lat, lon = np.atleast_1d(lat).astype(float), np.atleast_1d(lon).astype(float)
    t, F = _flux(lat, lon, nt)
    dt, dz = t[1] - t[0], DSKIN / 5.0
    assert dt < 0.45 * dz * dz / KAPPA, "unstable timestep"
    ns = lat.size
    T = np.full((ns, nz), 300.0)
    c = KAPPA * dt / dz**2
    Tmax = np.zeros(ns)
    for day in range(ndays):
        net = np.zeros(ns); Tmax = np.zeros(ns)
        for n in range(nt):
            T[:, 1:-1] += c * (T[:, 2:] - 2 * T[:, 1:-1] + T[:, :-2])
            T[:, -1] = T[:, -2]
            Ts, T1 = T[:, 0].copy(), T[:, 1]
            for _ in range(8):                      # implicit radiative surface balance
                f = EMIS * SIG * Ts**4 - F[:, n] + K * (Ts - T1) / dz
                Ts -= f / (4 * EMIS * SIG * Ts**3 + K / dz)
            T[:, 0] = Ts
            Tmax = np.maximum(Tmax, Ts)
            net += F[:, n] - EMIS * SIG * Ts**4
        net /= nt
        shift = net / (4 * EMIS * SIG * np.maximum(T[:, -1], 50.0) ** 3)
        T += np.clip(shift, -40, 40)[:, None]
        if verbose and day % 10 == 9:
            print(f"    day {day+1:3d}  max |imbalance| {np.abs(net).max():8.3f} W/m2")
    return T[:, -1], Tmax, np.abs(net).max()

def build_grid(dlat=5.0, dlon=5.0, **kw):
    """T on |lat| 0-90, lon 0-90. T is symmetric in latitude (no obliquity) and has
    period 180 deg in longitude, mirrored about the hot pole, so this quadrant is complete."""
    la = np.arange(0, 90 + dlat / 2, dlat)
    lo = np.arange(0, 90 + dlon / 2, dlon)
    LA, LO = np.meshgrid(la, lo, indexing="ij")
    Td, Tx, res = mean_annual_T(LA.ravel(), LO.ravel(), **kw)
    return la, lo, Td.reshape(LA.shape), Tx.reshape(LA.shape), res

def lookup(la, lo, G, lat, lon):
    """Bilinear interpolation, folding an arbitrary lat/lon into the computed quadrant."""
    aa = np.abs(np.asarray(lat, float))
    x = np.abs(((np.asarray(lon, float) % 180.0) + 180.0) % 180.0)
    x = np.where(x > 90.0, 180.0 - x, x)
    i = np.clip(np.searchsorted(la, aa) - 1, 0, len(la) - 2)
    j = np.clip(np.searchsorted(lo, x) - 1, 0, len(lo) - 2)
    u = (aa - la[i]) / (la[i + 1] - la[i])
    v = (x - lo[j]) / (lo[j + 1] - lo[j])
    return ((1 - u) * (1 - v) * G[i, j] + u * (1 - v) * G[i + 1, j]
            + (1 - u) * v * G[i, j + 1] + u * v * G[i + 1, j + 1])

if __name__ == "__main__":
    print(f"solar day {P_SOLAR:.2f} d | skin depth {DSKIN:.3f} m | "
          f"thermal inertia {np.sqrt(K*RHO*CP):.0f} J m-2 K-1 s-0.5")
    la, lo, G, X, res = build_grid(verbose=True)
    np.savez("thermal_grid.npz", la=la, lo=lo, G=G, X=X)
    print(f"converged: max residual surface flux imbalance {res:.4f} W/m2\n")
    print("  mean annual surface temperature (K)")
    print("   lat |   0E(hot)    30E     60E    90E(warm)")
    for a in (0, 15, 30, 45, 60, 75, 85):
        r = [lookup(la, lo, G, a, x) for x in (0, 30, 60, 90)]
        print(f"   {a:3d}  | " + "  ".join(f"{float(v):8.1f}" for v in r))
