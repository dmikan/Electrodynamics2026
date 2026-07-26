import os
import matplotlib.pyplot as plt

def ensure_output_dir(output_dir="output"):
    """
    Garantiza que la carpeta de salidas exista.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    return output_dir

def save_figure(fig, filename, output_dir="output", dpi=300):
    """
    Guarda una figura de matplotlib en la carpeta de salidas especificada.
    """
    ensure_output_dir(output_dir)
    filepath = os.path.join(output_dir, filename)
    fig.tight_layout()
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
    print(f"[OK] Gráfica guardada exitosamente en: {filepath}")
    return filepath
