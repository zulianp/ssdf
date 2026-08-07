#!/usr/bin/env python3

import numpy as np

def generate_sparse_rhombic_lattice(radius, height, num_layers, points_per_circle):
    # Z-levels for vertical layers
    z_vals = np.linspace(0, height, num_layers)
    theta = np.linspace(0, 2 * np.pi, points_per_circle, endpoint=False)

    # Create circular layers
    layers = []
    for z in z_vals:
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        z_layer = np.full_like(x, z)
        layers.append(np.stack([x, y, z_layer], axis=1))  # shape: (points_per_circle, 3)

    line0 = []
    line1 = []

    for i in range(num_layers):
        current = layers[i]

        # Horizontal connections (around circumference)
        for j in range(points_per_circle):
            p1 = current[j]
            p2 = current[(j + 1) % points_per_circle]
            # lines.append((p1, p2))
            line0.append(p1)
            line1.append(p2)

        # Diagonal connections to next layer
        if i < num_layers - 1:
            next_layer = layers[i + 1]
            for j in range(points_per_circle):
                if i % 2 == 0:
                    # Diagonal to next point in next layer
                    p1 = current[j]
                    p2 = next_layer[(j + 1) % points_per_circle]
                else:
                    # Diagonal to previous point in next layer
                    p1 = current[j]
                    p2 = next_layer[(j - 1) % points_per_circle]
                # lines.append((p1, p2))
                line0.append(p1)
                line1.append(p2)

    return line0, line1

# Parameters
radius = 0.5
height = 4.8
num_layers = 12            # Fewer layers → lower vertical frequency
points_per_circle = 20      # Fewer points per ring → lower circumferential frequency

# Generate lattice
l0, l1 = generate_sparse_rhombic_lattice(radius, height, num_layers, points_per_circle)

n = len(l0)

l0.extend(l1)
p = np.array(l0)

i0 = np.array(range(0, n)).astype(np.int32)
i1 = np.array(range(n, 2*n)).astype(np.int32)

x = p[:,0].astype(np.float32)
y = p[:,1].astype(np.float32)
z = p[:,2].astype(np.float32) + 0.1

i0.tofile("i0.int32")
i1.tofile("i1.int32")

x.tofile("x.float32")
y.tofile("y.float32")
z.tofile("z.float32")
