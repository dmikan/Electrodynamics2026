import numpy as np
import matplotlib.pyplot as plt
from .green_functions import DiscrGF, AvDiscrGF

def solve_single_plate(aside=1.0, bside=1.0, M=20, N=20, flag_galerkin=False,
                       scenario='constant', tresh=10.0):
    """
    Resuelve el problema electrostático para una placa aislada plana mediante el Método de los Momentos.

    Parámetros:
    - aside, bside: dimensiones de la placa en X e Y (m)
    - M, N: número de divisiones por lado (M x N celdas)
    - flag_galerkin: False para Point-Matching, True para Galerkin
    - scenario: 'constant' (potencial unitario U=1) o 'point_charge' (carga puntual sobre el centro)
    - tresh: umbral de distancia para aproximación en espacio libre (1/R)

    Retorna:
    - result: dict con campos calculados ('charge', 'cap_norm', 'rho_mat', 'xmat', 'ymat', etc.)
    """
    a, b = aside / M, bside / N
    x = np.linspace(a/2, aside - a/2, M)
    y = np.linspace(b/2, bside - b/2, N)
    xmat, ymat = np.meshgrid(x, y)

    xvec = xmat.flatten(order='F')
    yvec = ymat.flatten(order='F')
    ndim = M * N

    XC = np.abs(xvec[:, np.newaxis] - xvec[np.newaxis, :])
    YC = np.abs(yvec[:, np.newaxis] - yvec[np.newaxis, :])
    ZC = 0.0

    R = np.sqrt(XC**2 + YC**2 + ZC**2)
    rh = np.sqrt(XC**2 + YC**2)

    threshold_val = tresh * max(a, b)

    if flag_galerkin:
        gf_exact = AvDiscrGF(XC, YC, ZC, a, b)
    else:
        gf_exact = DiscrGF(XC, YC, ZC, a, b)

    with np.errstate(divide='ignore'):
        gf_far = 1.0 / R
    mom = np.where(rh > threshold_val, gf_far, gf_exact)

    # Vector de excitación
    if scenario == 'constant':
        unorm = 1.0
        excv = unorm * np.ones((ndim, 1))
    elif scenario == 'point_charge':
        unorm = 1.0
        xq, yq, zq = aside/2.0, bside/2.0, 0.5
        excv = 1.0 / np.sqrt((xq - xvec)**2 + (yq - yvec)**2 + zq**2)
        excv = excv.reshape(-1, 1)
    else:
        raise ValueError(f"Escenario no válido: {scenario}")

    charge = np.linalg.solve(mom, excv)
    charge_tot = float(np.sum(charge))
    cap_norm = charge_tot / unorm if scenario == 'constant' else charge_tot

    rho_s = charge / (a * b)
    rho_mat = rho_s.reshape((N, M), order='F')

    return {
        'aside': aside, 'bside': bside, 'M': M, 'N': N,
        'flag_galerkin': flag_galerkin, 'scenario': scenario,
        'charge': charge, 'charge_tot': charge_tot, 'cap_norm': cap_norm,
        'rho_mat': rho_mat, 'xmat': xmat, 'ymat': ymat
    }

def plot_single_plate_density(res):
    """
    Genera figura 3D de la densidad de carga de una placa aislada.
    """
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(res['xmat'], res['ymat'], res['rho_mat'],
                           cmap='viridis', edgecolor='black', lw=0.1)

    metodo = "Galerkin" if res['flag_galerkin'] else "Point-Matching"
    ax.set_title(f"Distribución de Densidad de Carga en Placa ({metodo})", fontsize=12, fontweight='bold')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('Densidad de Carga Normalizada')
    fig.colorbar(surf, shrink=0.5, aspect=10, label=r'Densidad $\rho_s$')

    return fig
