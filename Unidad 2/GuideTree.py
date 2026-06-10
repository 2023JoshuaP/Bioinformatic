import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def parse_results(filename):
    sequences = {}
    merges = []

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    in_tree_section = False

    while i < len(lines):
        line = lines[i].strip()
        line_upper = line.upper()

        if "GUIDE TREE DATA" in line_upper:
            in_tree_section = True
            i += 1
            continue

        if in_tree_section:
            if line_upper == "SEQUENCES":
                i += 1
                while i < len(lines):
                    current = lines[i].strip().upper()
                    if current == "MERGES":
                        i -= 1
                        break
                    parts = lines[i].strip().split()
                    if len(parts) >= 2:
                        try:
                            sequences[int(parts[0])] = parts[1]
                        except ValueError:
                            pass
                    i += 1

            elif line_upper == "MERGES":
                i += 1
                while i < len(lines):
                    current = lines[i].strip().upper()
                    if current == "END":
                        break
                    parts = lines[i].strip().split()
                    if len(parts) == 3:
                        try:
                            merges.append((int(parts[0]), int(parts[1]), float(parts[2])))
                        except ValueError:
                            pass
                    i += 1

            elif line_upper == "END":
                break

        i += 1

    return sequences, merges


def get_leaf_order(sequences, merges):
    clusters = {idx: [idx] for idx in sequences.keys()}

    for idx_a, idx_b, _ in merges:
        merged = clusters[idx_a] + clusters[idx_b]
        clusters[idx_a] = merged
        clusters[idx_b] = []

    for idx in sequences.keys():
        if clusters[idx]:
            return clusters[idx]
    return list(sequences.keys())


def draw_tree(sequences, merges):
    fig, ax = plt.subplots(figsize=(14, 9))
    ax.set_title("UPGMA Guide Tree", fontsize=16, fontweight='bold', pad=20)
    ax.axis('off')

    n = len(sequences)
    colors = ['#E63946', '#2A9D8F', '#E9C46A', '#F4A261', '#6A4C93', '#457B9D']

    # Ordenar hojas segun merges
    leaf_order = get_leaf_order(sequences, merges)
    y_positions = {idx: (n - i) * 3.0 for i, idx in enumerate(leaf_order)}

    # Normalizar distancias al ancho util [0.5, 0.95]
    all_dists = [d for _, _, d in merges]
    min_d, max_d = min(all_dists), max(all_dists)

    def norm(d):
        if max_d == min_d:
            return 0.7
        return 0.5 + (d - min_d) / (max_d - min_d) * 0.45

    LEAF_X = 0.32

    # Etiquetas hojas
    for idx, name in sequences.items():
        ax.text(LEAF_X - 0.01, y_positions[idx], name,
                ha='right', va='center', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#AED6F1',
                          edgecolor='#2E86C1', linewidth=1.5))

    node_y = dict(y_positions)
    node_x = {idx: LEAF_X for idx in sequences}

    for step, (idx_a, idx_b, dist) in enumerate(merges):
        color = colors[step % len(colors)]
        x_merge = norm(dist)
        y_a = node_y[idx_a]
        y_b = node_y[idx_b]
        x_a = node_x[idx_a]
        x_b = node_x[idx_b]
        y_new = (y_a + y_b) / 2.0

        # Rama A horizontal
        ax.plot([x_a, x_merge], [y_a, y_a], color=color, lw=2.5, solid_capstyle='round')
        # Rama B horizontal
        ax.plot([x_b, x_merge], [y_b, y_b], color=color, lw=2.5, solid_capstyle='round')
        # Vertical
        ax.plot([x_merge, x_merge], [y_b, y_a], color=color, lw=2.5, solid_capstyle='round')
        # Nodo
        ax.plot(x_merge, y_new, 'o', color=color, markersize=9, zorder=5,
                markeredgecolor='white', markeredgewidth=1.5)
        # Etiqueta
        ax.text(x_merge + 0.012, y_new + 0.25,
                f"S{step+1}  d={dist:.2f}",
                fontsize=8, color='white', ha='left', va='bottom', zorder=6,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=color,
                          edgecolor='none', alpha=0.9))

        node_y[idx_a] = y_new
        node_x[idx_a] = x_merge

    # Leyenda
    patches = [
        mpatches.Patch(color=colors[i % len(colors)],
                       label=f"S{i+1}: {sequences[merges[i][0]]} + {sequences[merges[i][1]]}")
        for i in range(len(merges))
    ]
    ax.legend(handles=patches, loc='lower right', fontsize=9,
              framealpha=0.9, title="Merge steps")

    ax.set_xlim(0.0, 1.1)
    ax.set_ylim(0.5, (n + 1) * 3.0)

    plt.tight_layout()
    plt.savefig("guide_tree.png", dpi=150, bbox_inches='tight')
    plt.show()
    print("Saved: guide_tree.png")


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python GuideTree.py <results.txt>")
        sys.exit(1)

    sequences, merges = parse_results(sys.argv[1])
    draw_tree(sequences, merges)