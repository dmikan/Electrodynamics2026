import numpy as np
import matplotlib.pyplot as plt
from .green_functions import DiscrGF, AvDiscrGF

def solve_cylinder(R=0.5, H=1.0, M=32, N=24, flag_galerkin=False,
                   pos_q1=(0.0, 0.0, 0.5), pos_q2=(0.8, -0.5, 0.5)):
    """
    Resuelve la solución MoM para un cilindro conductor cerrado ante fuentes puntuales.

    Parámetros:
    - R: Radio del cilindro
    - H: Altura del cilindro
    - M: Divisiones angulares
    - N: Divisiones verticales
    - flag_galerkin: 0 para PM, 1 para Galerkin
    - pos_q1: Posición de carga puntual interna
    - pos_q2: Posición de carga puntual externa

    Retorna:
    - Dict con mallas, vectores de carga inducida y parámetros geométricos.
    """
    a_cell = (2 * np.pi * R) / M
    b_cell = H / N

    theta = np.linspace(0, 2*np.pi, M, endpoint=False)
    z = np.linspace(b_cell/2, H - b_cell/2, N)
    theta_mat, z_mat = np.meshgrid(theta, z)

    xvec = (R * np.cos(theta_mat)).flatten()
    yvec = (R * np.sin(theta_mat)).flatten()
    zvec = z_mat.flatten()

    dx = xvec[:, np.newaxis] - xvec[np.newaxis, :]
    dy = yvec[:, np.newaxis] - yvec[np.newaxis, :]
    dz = zvec[:, np.newaxis] - zvec[np.newaxis, :]
    R_dist = np.sqrt(dx**2 + dy**2 + dz**2)

    mom = 1.0 / (R_dist + np.finfo(float).eps)

    if flag_galerkin:
        self_val = AvDiscrGF(0, 0, 0, a_cell, b_cell)
    else:
        self_val = DiscrGF(0, 0, 0, a_cell, b_cell)

    np.fill_diagonal(mom, self_val)

    # Excitaciones
    xq1, yq1, zq1 = pos_q1
    dist_q1 = np.sqrt((xvec - xq1)**2 + (yvec - yq1)**2 + (zvec - zq1)**2)
    excv_int = 1.0 / dist_q1

    xq2, yq2, zq2 = pos_q2
    dist_q2 = np.sqrt((xvec - xq2)**2 + (yvec - yq2)**2 + (zvec - zq2)**2)
    excv_ext = 1.0 / dist_q2

    # Solución de sistema
    charge_int = np.linalg.solve(mom, -excv_int)
    charge_ext = np.linalg.solve(mom, -excv_ext)

    return {
        'R': R, 'H': H, 'M': M, 'N': N, 'a_cell': a_cell, 'b_cell': b_cell,
        'theta_mat': theta_mat, 'z_mat': z_mat,
        'xvec': xvec, 'yvec': yvec, 'zvec': zvec,
        'charge_int': charge_int, 'charge_ext': charge_ext,
        'pos_q1': pos_q1, 'pos_q2': pos_q2,
        'flag_galerkin': flag_galerkin
    }

def plot_cylinder_density(res):
    """
    Genera gráficas 3D de la densidad de carga inducida en la superficie del cilindro.
    """
    R, H, M, N = res['R'], res['H'], res['M'], res['N']
    theta_mat, z_mat = res['theta_mat'], res['z_mat']
    a_cell, b_cell = res['a_cell'], res['b_cell']

    X_plot = R * np.cos(theta_mat)
    Y_plot = R * np.sin(theta_mat)
    Z_plot = z_mat

    X_plot = np.hstack((X_plot, X_plot[:, :1]))
    Y_plot = np.hstack((Y_plot, Y_plot[:, :1]))
    Z_plot = np.hstack((Z_plot, Z_plot[:, :1]))

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), subplot_kw={'projection': '3d'})
    metodo = "Galerkin" if res['flag_galerkin'] else "Point-Matching"

    cases = [
        (res['charge_int'], res['pos_q1'], "Carga Interna", 'viridis_r'),
        (res['charge_ext'], res['pos_q2'], "Carga Externa", 'viridis_r')
    ]

    for i, (charge, pos_q, title, cmap) in enumerate(cases):
        rho_mat = (charge / (a_cell * b_cell)).reshape((N, M))
        rho_v = np.hstack((rho_mat, rho_mat[:, :1]))

        ax = axes[i]
        norm = plt.Normalize(rho_v.min(), rho_v.max())
        surf = ax.plot_surface(X_plot, Y_plot, Z_plot,
                               facecolors=plt.colormaps[cmap](norm(rho_v)),
                               shade=False, edgecolor='black', lw=0.1, alpha=0.8)
        ax.scatter(*pos_q, color='black', s=100, label='Fuente (q=1)')
        ax.set_title(f"{title}\n({metodo})", fontsize=12, fontweight='bold')
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                     shrink=0.5, label=r'Densidad $\rho_s$ (C/m$^2$)')

    return fig

def potential_cylinder(res):
    """
    Calcula y gráfica los contornos del potencial total V(x,y) en el plano z = 0.5.
    """
    xvec, yvec, zvec = res['xvec'], res['yvec'], res['zvec']
    charge_int, charge_ext = res['charge_int'], res['charge_ext']
    pos_q1, pos_q2 = res['pos_q1'], res['pos_q2']
    R = res['R']

    xo = np.linspace(-1.5, 1.5, 100)
    yo = np.linspace(-1.5, 1.5, 100)
    Xo, Yo = np.meshgrid(xo, yo)
    Zo = 0.5

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    casos = [
        (charge_int, pos_q1, "Carga Interna (0, 0, 0.5)"),
        (charge_ext, pos_q2, "Carga Externa (0.8, 0, 0.5)")
    ]

    for i, (charge, pos_q, titulo) in enumerate(casos):
        xq, yq, zq = pos_q
        dist_q = np.sqrt((Xo - xq)**2 + (Yo - yq)**2 + (Zo - zq)**2)
        V_q = 1.0 / (dist_q + 1e-6)

        dist_cells = np.sqrt((Xo[:, :, np.newaxis] - xvec)**2 +
                             (Yo[:, :, np.newaxis] - yvec)**2 + (Zo - zvec)**2)
        V_ind = np.sum(charge / (dist_cells + 1e-6), axis=2)

        Potencial = V_q + V_ind

        ax = axes[i]
        niveles = np.linspace(-0.5, 0.5, 21)
        cp = ax.contourf(Xo, Yo, Potencial, levels=niveles, cmap='RdBu_r', extend='both')
        ax.contour(Xo, Yo, Potencial, levels=niveles, colors='black', linewidths=0.5, alpha=0.5)

        circle = plt.Circle((0, 0), R, fill=False, color='black', linestyle='--', linewidth=2, label='Cilindro')
        ax.add_artist(circle)
        ax.scatter(xq, yq, color='yellow', edgecolor='black', s=100, label='Fuente (q=1)', zorder=5)

        ax.set_title(titulo, fontsize=12, fontweight='bold')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
        ax.legend()
        fig.colorbar(cp, ax=ax, label=r'Potencial Normalizado ($4\pi\varepsilon_0V$)')

    return fig
