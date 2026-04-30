import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
import os

def load_matrix(matrix_path):
    return np.loadtxt(matrix_path, dtype=int)

def load_path(path_path):
    path = set()
    with open(path_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 2:
                path.add((int(parts[0]), int(parts[1])))
    return path

def plot_heatmap(matrix_path, path_path, output_img=None):
    print(f"Reading matrix: {matrix_path}")
    print(f"Reading path: {path_path}")

    H = load_matrix(matrix_path)
    path = load_path(path_path)

    print(f"Matrix size: {H.shape}")
    print(f"Path points: {len(path)}")
    print(f"Max score: {H.max()}")

    if output_img is None:
        base_name = os.path.splitext(os.path.basename(matrix_path))[0]
        output_img = os.path.join("Heatmaps", base_name + ".png")

    fig, ax = plt.subplots(figsize=(14, 12))

    im = ax.imshow(H, cmap="YlOrRd", aspect="auto", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Score", fraction=0.03, pad=0.04)

    if path:
        path_coords = [(j, i) for (i, j) in path if i < H.shape[0] and j < H.shape[1]]
        if path_coords:
            xs, ys = zip(*path_coords)
            ax.scatter(xs, ys, s=10, c='steelblue', alpha=0.8, linewidths=0)

    ax.set_xlabel("Sequence 2 (positions)", fontsize=11)
    ax.set_ylabel("Sequence 1 (positions)", fontsize=11)
    ax.set_title(
        f"Smith-Waterman Score Matrix\n"
        f"Size={H.shape} | Max score={H.max()}",
        fontsize=12, pad=15
    )

    patch_path = mpatches.Patch(color='steelblue', alpha=0.8, label='Optimal traceback path')
    patch_high = mpatches.Patch(color='darkred', label='High score')
    patch_low = mpatches.Patch(color='lightyellow', label='Low score / zero')
    ax.legend(handles=[patch_low, patch_high, patch_path], loc='lower right', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_img, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_img}")
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python ShowMatrix.py <_matrix.txt> <_path.txt> [output.png]")
        sys.exit(1)

    matrix_file = sys.argv[1]
    path_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) >= 4 else None

    plot_heatmap(matrix_file, path_file, output_file)