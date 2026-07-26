import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from src.green_functions import DiscrGF, AvDiscrGF
from src.single_plate import solve_single_plate, plot_single_plate_density
from src.dual_plates import solve_dual_plates, plot_dual_plates_3d
from src.cylinder import solve_cylinder, plot_cylinder_density, potential_cylinder
from src.utils import ensure_output_dir, save_figure

def plot_green_functions_coplanar(num_points=200):
    """
    Genera la comparación de Funciones de Green para celdas coplanares (h=0, yc=0).
    """
    a, b = 1.0, 1.0
    xc_vals = np.linspace(0.0, 2.0, num_points)
    psi_a = np.array([AvDiscrGF(x, 0.0, 0.0, a, b) for x in xc_vals])
    gamma_a = np.array([DiscrGF(x, 0.0, 0.0, a, b) for x in xc_vals])

    xc_g = np.linspace(0.27, 2.0, num_points)
    g_a = 1.0 / xc_g

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xc_vals, psi_a, label=r'$\Psi$ (Galerkin)', color='#3366cc', lw=2)
    ax.plot(xc_vals, gamma_a, label=r'$\Gamma$ (PM)', color='#cc3333', lw=2)
    ax.plot(xc_g, g_a, label=r'$G$ (Espacio libre)', color='black', linestyle='--', lw=1.8)

    ax.set_title(r'Celdas Coplanares ($h = 0$, $y_c = 0$)', fontsize=12, fontweight='bold')
    ax.set_xlabel(r'$x_c$ = Distancia horizontal entre centros', fontsize=11)
    ax.set_ylabel('Valor GF Normalizado', fontsize=11)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 4)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', frameon=True)
    return fig

def plot_green_functions_parallel(num_points=200):
    """
    Genera la comparación de Funciones de Green para celdas paralelas (xc=0, yc=0).
    """
    a, b = 1.0, 1.0
    zc_vals = np.linspace(0.0, 2.0, num_points)
    psi_b = np.array([AvDiscrGF(0.0, 0.0, z, a, b) for z in zc_vals])
    gamma_b = np.array([DiscrGF(0.0, 0.0, z, a, b) for z in zc_vals])

    zc_g = np.linspace(0.27, 2.0, num_points)
    g_b = 1.0 / zc_g

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(zc_vals, psi_b, label=r'$\Psi$ (Galerkin)', color='#3366cc', lw=2)
    ax.plot(zc_vals, gamma_b, label=r'$\Gamma$ (PM)', color='#cc3333', lw=2)
    ax.plot(zc_g, g_b, label=r'$G$ (Espacio libre)', color='black', linestyle='--', lw=1.8)

    ax.set_title(r'Celdas Paralelas ($x_c = 0$, $y_c = 0$)', fontsize=12, fontweight='bold')
    ax.set_xlabel(r'$z_c$ = Distancia vertical ($h$) entre centros', fontsize=11)
    ax.set_ylabel('Valor GF Normalizado', fontsize=11)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 4)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', frameon=True)
    return fig

def plot_green_functions_surface_3d(h_val=0.3, a=1.0, b=1.0, grid_size=30):
    """
    Genera superficie 3D de la Función de Green para un h especificado.
    """
    x = np.linspace(0.0, 2.0, grid_size)
    y = np.linspace(0.0, 2.0, grid_size)
    X, Y = np.meshgrid(x, y)

    Z_psi = np.zeros_like(X)
    Z_gamma = np.zeros_like(X)

    for i in range(grid_size):
        for j in range(grid_size):
            Z_psi[i, j] = AvDiscrGF(X[i, j], Y[i, j], h_val, a, b)
            Z_gamma[i, j] = DiscrGF(X[i, j], Y[i, j], h_val, a, b)

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')

    surf_psi = ax.plot_surface(X, Y, Z_psi, cmap='viridis', alpha=0.85, edgecolor='none')
    ax.plot_wireframe(X, Y, Z_gamma, color='red', linewidth=0.5, rstride=2, cstride=2)

    ax.set_title(f'Función de Green Galerkin vs PM ($h = {h_val}$)', fontsize=12, fontweight='bold')
    ax.set_xlabel(r'$x_c$ (Distancia en X)', fontsize=10)
    ax.set_ylabel(r'$y_c$ (Distancia en Y)', fontsize=10)
    ax.set_zlabel('Valor GF Normalizado', fontsize=10)

    fig.colorbar(surf_psi, ax=ax, shrink=0.5, aspect=10, pad=0.1, label=r'Valor $\Psi$')
    return fig

def main():
    print("==========================================================")
    print("      ELECTRODYNÁMICA 2026 - MÉTODO DE LOS MOMENTOS (MoM)  ")
    print("==========================================================")
    print("Ejecutando simulaciones y generando gráficas principales...\n")

    output_dir = ensure_output_dir()

    # 1. Comparación de Funciones de Green
    print("[1/5] Generando comparación de Funciones de Green...")
    fig1 = plot_green_functions_coplanar()
    save_figure(fig1, "01_green_functions_coplanar.png", output_dir)
    plt.close(fig1)

    fig2 = plot_green_functions_parallel()
    save_figure(fig2, "02_green_functions_parallel.png", output_dir)
    plt.close(fig2)

    fig3 = plot_green_functions_surface_3d(h_val=0.3)
    save_figure(fig3, "03_green_functions_surface_h03.png", output_dir)
    plt.close(fig3)

    # 2. Solución Placa Única Conductora
    print("\n[2/5] Calculando solución para Placa Conductora Aislada...")
    res_plate = solve_single_plate(aside=1.0, bside=1.0, M=20, N=20, flag_galerkin=True)
    print(f"   -> Capacitancia normalizada calculada (N=20): {res_plate['cap_norm']:.4f}")
    fig4 = plot_single_plate_density(res_plate)
    save_figure(fig4, "04_single_plate_charge_density.png", output_dir)
    plt.close(fig4)

    # 3. Solución Sistema de Dos Placas
    print("\n[3/5] Calculando solución para Sistema de Dos Placas (A_top=1.0, A_bot=0.5, h=0.3)...")
    res_dual = solve_dual_plates(a_top=1.0, a_bot=0.5, N_bot=20, h_sep=0.3, flag_galerkin=True, verbose=True)
    fig5 = plot_dual_plates_3d(res_dual)
    save_figure(fig5, "05_dual_plates_charge_density.png", output_dir)
    plt.close(fig5)

    # 4. Solución Cilindro Conductor
    print("\n[4/5] Calculando solución para Cilindro Conductor Cerrado...")
    res_cyl = solve_cylinder(R=0.5, H=1.0, M=32, N=24, flag_galerkin=False)
    print(f"   -> Carga total inducida por carga interna (0,0,0.5): {np.sum(res_cyl['charge_int']):.4f}")
    print(f"   -> Carga total inducida por carga externa (0.8,-0.5,0.5): {np.sum(res_cyl['charge_ext']):.4f}")

    fig6 = plot_cylinder_density(res_cyl)
    save_figure(fig6, "06_cylinder_charge_density.png", output_dir)
    plt.close(fig6)

    # 5. Mapas de Contorno de Potencial en Cilindro
    print("\n[5/5] Generando mapas de contorno de potencial V(x,y) para Cilindro Conductor...")
    fig7 = potential_cylinder(res_cyl)
    save_figure(fig7, "07_cylinder_potential_contours.png", output_dir)
    plt.close(fig7)

    print("\n==========================================================")
    print(f" ¡PROCESO COMPLETADO EXITOSAMENTE!")
    print(f" Todas las gráficas se guardaron en la carpeta: {os.path.abspath(output_dir)}")
    print("==========================================================")

if __name__ == "__main__":
    main()
