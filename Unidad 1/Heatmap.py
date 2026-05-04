import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def load_matrix(matrix_path):
    return np.loadtxt(matrix_path, delimiter=",", dtype=int)

def plot_heatmap(matrix_path, header1, header2):
    print(f"Reading matrix: {matrix_path}")

    H = load_matrix(matrix_path)
    print(f"Matrix size: {H.shape}")

    output_img = os.path.join("Heatmaps", header1 + "_vs_" + header2 + "_alignment_matrix.png")

    fig, ax = plt.subplots(figsize=(14, 12))

    im = ax.imshow(H, cmap="YlOrRd", aspect="auto", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Score", fraction=0.03, pad=0.04)

    ax.set_xlabel("Sequence 2 (positions)", fontsize=11)
    ax.set_ylabel("Sequence 1 (positions)", fontsize=11)
    ax.set_title(
        f"Needleman-Wunsch Score Matrix\nSize={H.shape} | Max={H.max()} | Min={H.min()}",
        fontsize=12, pad=15
    )

    plt.tight_layout()
    plt.savefig(output_img, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_img}")
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python ShowMatrix.py <matrix.csv> <header1> <header2>")
        sys.exit(1)

    plot_heatmap(sys.argv[1], sys.argv[2], sys.argv[3])