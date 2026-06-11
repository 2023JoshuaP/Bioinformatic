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
    # Sort merges by distance so leaf ordering respects actual tree topology
    for _, (idx_a, idx_b, _) in sorted(enumerate(merges), key=lambda t: t[1][2]):
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

    # Y positions
    leaf_order = get_leaf_order(sequences, merges)
    y_positions = {idx: (n - i) * 3.0 for i, idx in enumerate(leaf_order)}

    # X scale: leaves at 0 (left), root at 1.0 (right)
    all_dists = [d for _, _, d in merges]
    max_d = max(all_dists)

    LABEL_X = 0.0
    ROOT_X  = 1.0

    def dist_to_x(d):
        if max_d == 0:
            return ROOT_X * 0.5
        return (d / max_d) * ROOT_X

    # Draw leaf labels
    for idx, name in sequences.items():
        ax.text(LABEL_X - 0.01, y_positions[idx], name,
                ha='right', va='center', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#AED6F1',
                          edgecolor='#2E86C1', linewidth=1.5))

    node_y = dict(y_positions)
    node_x = {idx: LABEL_X for idx in sequences}

    # Sort merges by distance (ascending) so x positions are always non-decreasing.
    # Keep orig_step for label / color so S1/S2/… stay correct.
    sorted_merges = sorted(enumerate(merges), key=lambda t: t[1][2])

    for orig_step, (idx_a, idx_b, dist) in sorted_merges:
        color   = colors[orig_step % len(colors)]
        x_merge = dist_to_x(dist)
        y_a     = node_y[idx_a]
        y_b     = node_y[idx_b]
        x_a     = node_x[idx_a]
        x_b     = node_x[idx_b]
        y_new   = (y_a + y_b) / 2.0

        # Horizontal arms from each child to the merge bar
        ax.plot([x_a, x_merge], [y_a, y_a], color=color, lw=2.5, solid_capstyle='round')
        ax.plot([x_b, x_merge], [y_b, y_b], color=color, lw=2.5, solid_capstyle='round')
        # Vertical bar
        ax.plot([x_merge, x_merge], [y_b, y_a], color=color, lw=2.5, solid_capstyle='round')
        # Dot at midpoint
        ax.plot(x_merge, y_new, 'o', color=color, markersize=9, zorder=5,
                markeredgecolor='white', markeredgewidth=1.5)
        # Label
        ax.text(x_merge + 0.015, y_new + 0.2,
                f"S{orig_step+1}  d={dist:.2f}",
                fontsize=8, color='white', ha='left', va='bottom', zorder=6,
                bbox=dict(boxstyle='round,pad=0.3', facecolor=color,
                          edgecolor='none', alpha=0.9))

        # Update merged cluster representative
        node_y[idx_a] = y_new
        node_x[idx_a] = x_merge

    # Legend (in original step order)
    patches = [
        mpatches.Patch(color=colors[i % len(colors)],
                       label=f"S{i+1}: {sequences[merges[i][0]]} + {sequences[merges[i][1]]}")
        for i in range(len(merges))
    ]
    ax.legend(handles=patches, loc='lower right', fontsize=9,
              framealpha=0.9, title="Merge steps")

    ax.set_xlim(-0.45, ROOT_X + 0.15)
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