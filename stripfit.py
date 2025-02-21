import math
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-((x - mean) / sigma)**2 / 2)

def fit_gaussian(hist, bins):
    bin_centers = (bins[1:] + bins[:-1]) / 2
    popt, pcov = curve_fit(gaussian, bin_centers, hist)
    return popt

class Strip:
    def __init__(self, a, b, c, w, L, mod):
        self.a = a
        self.b = b
        self.c = c
        self.w = w
        self.L = L
        self.mod = mod

class Point3D:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError("Index out of range")

def extract_numbers(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        blocks = []
        block = []
        numbers = []
        i = 0
        while i < len(lines):
            if lines[i].strip() == "start":
                i += 1
                block = []
                numbers = []
                block.append([float(x) for x in lines[i].strip().split(',')])
                i += 1
            while lines[i].strip() != "stop":
                numbers.append([float(x) for x in lines[i].strip().split(',')])
                i += 1
            block.append(numbers)
            blocks.append(block)
            i += 1
        return blocks

def linetorectangle(params, rectangle):
    x0, y0, z0, ux, uy, uz = params
    a, b, c, W, L, modu = rectangle.a, rectangle.b, rectangle.c, rectangle.w, rectangle.L, rectangle.mod

    # Ensure the direction vector is normalized
    norm = math.sqrt(ux * ux + uy * uy + uz * uz)
    ux /= norm
    uy /= norm
    uz /= norm

    # Parametric line: p(t) = (x0 + t*ux, y0 + t*uy, z0 + t*uz)
    # Project the rectangle center onto the line
    t = (a - x0) * ux + (b - y0) * uy + (c - z0) * uz

    # Closest point on the line
    px = x0 + t * ux
    py = y0 + t * uy
    pz = z0 + t * uz

    # Determine the dimensions of the rectangle based on modu
    if modu == 0.0 or modu == 2.0:
        # Rectangle is 79.2 long in x and 1 long in y
        half_length_x = 79.2 / 2
        half_length_y = 1
        half_length_z = 0.001  # Assuming the rectangle is 1 unit thick in z
        dx = px - a
        dy = py - b
        dz = pz - c
    elif modu == 1.0 or modu == 3.0:
        # Rectangle is 79.2 long in y and 1 long in x
        half_length_x = 1
        half_length_y = 79.2 / 2
        half_length_z = 0.001  # Assuming the rectangle is 1 unit thick in z
        dx = px - a
        dy = py - b
        dz = pz - c

    # Clamp the closest point to the rectangle's bounds
    clamped_x = max(-half_length_x, min(dx, half_length_x))
    clamped_y = max(-half_length_y, min(dy, half_length_y))
    clamped_z = max(-half_length_z, min(dz, half_length_z))

    # Calculate the vector from the closest point on the line to the clamped point on the rectangle
    closest_point_on_rectangle = (a + clamped_x, b + clamped_y, c + clamped_z)
    closest_point_on_line = (px, py, pz)

    # Calculate the distance between these two points
    distance = math.sqrt((closest_point_on_rectangle[0] - closest_point_on_line[0]) ** 2 +
                         (closest_point_on_rectangle[1] - closest_point_on_line[1]) ** 2 +
                         (closest_point_on_rectangle[2] - closest_point_on_line[2]) ** 2)

    return distance

def FitFunction(params, strips):
    x0, y0, z0, ux, uy, uz = params

    # Add penalties for violating the constraints
    penalty = 0.0
    if x0 < -2:
        penalty += 1e9 * (abs(x0) - 2)
    elif x0 > 2:
        penalty += 1e9 * (x0 - 2)
    if y0 < -2:
        penalty += 1e9 * (abs(y0) - 2)
    elif y0 > 2:
        penalty += 1e9 * (y0 - 2)
    if z0 < -49.26:
        penalty += 1e9 * (abs(z0) + 49.26)
    elif z0 > 94.74:
        penalty += 1e9 * (z0 - 94.74)

    total_sum = 0.0

    for strip in strips:
        total_sum += linetorectangle(params, strip)
    return total_sum + penalty

def FitRect(strips, first):
    def fit_function(params):
        return FitFunction(params, strips)

    initial_params = np.array([first[0], first[1], first[2], 0, 0, 1])
    bounds = [(-2, 2), (-2, 2), (-49.26, 94.74), (-1, 1), (-1, 1), (-1, 1)]
    result = minimize(fit_function, initial_params, method='Powell', bounds=bounds)

    centroid = result.x[:3]
    direction = result.x[3:]

    return direction, centroid

def findClosestPointOnLine(centroid, direction, point):
    # Unpack components
    xc, yc, zc = centroid[0], centroid[1], centroid[2]
    dx, dy, dz = direction[0], direction[1], direction[2]
    a, b, c = point[0], point[1], point[2]
    dist = 99999
    tmin = 0

    for t in [x * 0.01 for x in range(-2000, 2000)]:
        x = xc + t * dx
        y = yc + t * dy
        z = zc + t * dz
        diff = math.sqrt((x - a) ** 2 + (y - b) ** 2 + (z - c) ** 2)

        if diff < dist:
            tmin = t
            dist = diff

    # Compute the closest point
    closestPoint = Point3D(
        xc + tmin * dx,  # x-coordinate
        yc + tmin * dy,  # y-coordinate
        zc + tmin * dz  # z-coordinate
    )

    return closestPoint

filename = 'box4.csv'  # replace with your file name
result = extract_numbers(filename)
reco_z_err = []
count = 0

for block in result:
    strips = []
    origin = Point3D(block[0][2], block[0][3], block[0][4])
    count = block[0][0]

    if count > 100:
        break

    for line in block[1]:
        print(line)
        if line[0] == 0.0:
            continue
        strip = Strip(line[0], line[1], line[2], line[3], line[4], line[5])
        strips.append(strip)

    first = [0, 0, 1]
    min_diff = 1e6

    for i in range(len(strips)):
        fitvec, fitcent = FitRect([strips[i]], first)
        reco_v = findClosestPointOnLine(fitcent, fitvec, origin)
        reco_z_err.append(reco_v[2] - origin[2])

hist, bins = np.histogram(reco_z_err, bins=200, range=(-100, 100))
bin_centers = (bins[1:] + bins[:-1]) / 2
popt, pcov = curve_fit(gaussian, bin_centers, hist)

min_idx = np.digitize(popt[1] - abs(3 * popt[2]), bins)
max_idx = np.digitize(popt[1] + abs(3 * popt[2]), bins)

integral = 0

for i in range(max_idx - min_idx):
    integral += hist[min_idx + i]

print("mean = ", popt[1], " +/-", np.sqrt(pcov[1, 1]), " [mm]")
print("sigma = ", popt[2], " +/-", np.sqrt(pcov[2, 2]), " [mm]")
print("integral =", integral, ", 3sig =", integral / count)

x = np.linspace(bin_centers[0], bin_centers[-1], 100)
y = gaussian(x, *popt)

plt.hist(reco_z_err, bins=200, range=(-100, 100), alpha=0.5, label='Histogram')
plt.plot(x, y, label='Fitted Gaussian')
plt.legend()
plt.show()

#Understanding geometry of reconstruction incorrectly its lines not surfaces
