"Author: Mayukhmala Jana"
import matplotlib.pyplot as plt
import numpy as np
import gc


def clear_all():
    # Clear local variables
    local_vars = locals()
    for var in local_vars:
        if var != "clear_all":
            del local_vars[var]

    # Run garbage collection to free up memory
    gc.collect()


def drawNucContoursByClass_SA2(M, I, centroids, classes, colors, tickness=0.5):
    # clear_all()
    numCent = len(centroids)
    boundaries = find_boundaries(M)

    # plt.figure()
    plt.imshow(I, cmap="gray", origin="upper")
    plt.axis("off")

    for i in range(len(boundaries)):
        b = boundaries[i]
        # print(len(boundaries))
        w, x, y, z = min(b[:, 1]), max(b[:, 1]), min(b[:, 0]), max(b[:, 0])
        # print(w, x, y, z)
        for j in range(numCent):
            if (
                centroids[j, 0] > y
                and centroids[j, 0] < z
                and centroids[j, 1] > w
                and centroids[j, 1] < x
            ):
                index = int(classes[j])
                if 1 <= index <= len(colors):
                    if index == 1:
                        plt.plot(
                            b[:, 1],
                            b[:, 0],
                            color=colors[index - 1],
                            linewidth=tickness,
                        )  # color=(1, 1, 0), linewidth=tickness)   # Yellow for class 0
                    elif index == 2:
                        plt.plot(
                            b[:, 1],
                            b[:, 0],
                            color=colors[index - 1],
                            linewidth=tickness,
                        )  # color=(0, 0, 1), linewidth=tickness)  # Blue for class 1
                    else:
                        print(f"Unexpected class index: {index}")
                else:
                    print(f"Index {index} is out of range for the 'colors' list.")

    plt.show()


def find_boundaries(M):
    boundaries = []
    rows, cols = M.shape

    for i in range(rows):
        for j in range(cols):
            if M[i, j] == 1:
                boundary = trace_boundary(M, i, j)
                if boundary.size > 0:  # Check if the array is not empty
                    boundaries.append(boundary)

    return boundaries


def trace_boundary(M, i, j):
    rows, cols = M.shape
    boundary = []

    stack = [(i, j)]
    while stack:
        current_i, current_j = stack.pop()
        if (
            0 <= current_i < rows
            and 0 <= current_j < cols
            and M[current_i, current_j] == 1
        ):
            M[current_i, current_j] = 0  # Mark as visited
            boundary.append((current_i, current_j))
            stack.extend(
                [
                    (current_i, current_j + 1),
                    (current_i, current_j - 1),
                    (current_i + 1, current_j),
                    (current_i - 1, current_j),
                ]
            )

    return np.array(boundary)
