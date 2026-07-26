import numpy as np
import matplotlib.pyplot as plt
from .green_functions import DiscrGF, AvDiscrGF

def construir_placa(A, N, z):
    """
    Construye la malla de centros para una placa de lado A con N divisiones a una altura z.
    """
    Delta = A / N
    coords = np.linspace(-A/2 + Delta/2, A/2 - Delta/2, N)
    xmat, ymat = np.meshgrid(coords, coords)
    xvec = xmat.flatten(order='F')
    yvec = ymat.flatten(order='F')
    zvec = np.full_like(xvec, z)
    return xvec, yvec, zvec, Delta

def solve_dual_plates(a_top=1.0, a_bot=0.5, N_bot=20, h_sep=0.3,
                      tresh=10.0, flag_galerkin=False, verbose=True):
    """
    Resuelve el sistema de dos placas paralelas de tamaños arbitrarios (a_top y a_bot) a una distancia h_sep.

    Retorna:
    - Dict con matriz de capacitancia, cargas, coordenadas y métricas de reciprocidad.
    """
    N_top = round(a_top / a_bot * N_bot)
    assert abs(N_top - a_top / a_bot * N_bot) < 1e-9, "a_top/a_bot * N_bot debe ser entero."

    x_top, y_top, z_top, d_top = construir_placa(a_top, N_top, +h_sep/2)
    x_bot, y_bot, z_bot, d_bot = construir_placa(a_bot, N_bot, -h_sep/2)
    Delta = d_top
    assert abs(d_top - d_bot) < 1e-12, "El tamaño de celda debe ser idéntico en ambas placas."

    x = np.concatenate([x_top, x_bot])
    y = np.concatenate([y_top, y_bot])
    z = np.concatenate([z_top, z_bot])

    n_top_cells = N_top**2
    n_bot_cells = N_bot**2
    ndim = n_top_cells + n_bot_cells

    is_top = np.zeros(ndim, dtype=bool)
    is_top[:n_top_cells] = True

    XC = x[:, np.newaxis] - x[np.newaxis, :]
    YC = y[:, np.newaxis] - y[np.newaxis, :]
    ZC = z[:, np.newaxis] - z[np.newaxis, :]

    rh = np.sqrt(XC**2 + YC**2)
    R = np.sqrt(rh**2 + ZC**2)

    gf_cerca = AvDiscrGF(XC, YC, ZC, Delta, Delta) if flag_galerkin else DiscrGF(XC, YC, ZC, Delta, Delta)

    with np.errstate(divide='ignore'):
        gf_lejos = 1.0 / R
    mom = np.where(rh > tresh * Delta, gf_lejos, gf_cerca)

    def resolver_escenario(v_top, v_bot):
        exc = np.where(is_top, v_top, v_bot).reshape(-1, 1)
        charge = np.linalg.solve(mom, exc)
        return charge, float(np.sum(charge[is_top])), float(np.sum(charge[~is_top]))

    charge_A, C11, C21 = resolver_escenario(1.0, 0.0)
    charge_B, C12, C22 = resolver_escenario(0.0, 1.0)

    err_recip = 100 * abs(C12 - C21) / (0.5 * abs(C12 + C21) + 1e-15)

    if verbose:
        motor = "Galerkin" if flag_galerkin else "PM"
        print(f"[{motor}] Dos Placas: N_bot={N_bot}, N_top={N_top}, Delta={Delta:.4f}, dim={ndim}")
        print(f"  Capacitancias: C11={C11:.5f}, C22={C22:.5f}, C12={C12:.5f}, C21={C21:.5f}")
        print(f"  Reciprocidad (C12 vs C21): Dif = {err_recip:.4f}%")

    return {
        "N_top": N_top, "N_bot": N_bot, "Delta": Delta, "ndim": ndim,
        "C11": C11, "C22": C22, "C12": C12, "C21": C21, "err_recip": err_recip,
        "x": x, "y": y, "z": z, "is_top": is_top,
        "charge_A": charge_A, "charge_B": charge_B,
        "a_top": a_top, "a_bot": a_bot, "h_sep": h_sep,
        "flag_galerkin": flag_galerkin
    }

def plot_dual_plates_3d(res):
    """
    Grafica la densidad de carga 3D en el sistema de dos placas.
    """
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    x, y, z = res['x'], res['y'], res['z']
    is_top = res['is_top']
    charge_A = res['charge_A']
    Delta = res['Delta']

    rho = (charge_A / (Delta**2)).flatten()

    sc = ax.scatter(x, y, z, c=rho, cmap='plasma', s=20, alpha=0.8)
    fig.colorbar(sc, ax=ax, label=r'Densidad de Carga $\rho_s$ (C/m$^2$)')

    ax.set_title(f"Sistema de Dos Placas (A_top={res['a_top']}, A_bot={res['a_bot']}, h={res['h_sep']})",
                 fontsize=12, fontweight='bold')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')

    return fig
