import math
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.decomposition import PCA

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-((x - mean) / sigma)**2 / 2)

def fit_gaussian(hist, bins):
    bin_centers = (bins[1:] + bins[:-1]) / 2
    popt, pcov = curve_fit(gaussian, bin_centers, hist)
    return popt

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

def linetopoint(params, point):
    x0, y0, z0, ux, uy, uz = params

    P = np.array([point.x, point.y, point.z])
    O = np.array([x0, y0, z0])
    U = np.array([ux, uy, uz])

    OP = P - O

    cross_product = np.cross(OP, U)
    distance = np.linalg.norm(cross_product) / np.linalg.norm(U)

    return distance

def FitFunction(params, strips):
    x0, y0, z0, ux, uy, uz = params

    # Add penalties for violating the constraints
    penalty = 0.0
    if x0 < -3:
        penalty += 1e6 * (abs(x0) - 3)
    elif x0 > 3:
        penalty += 1e6 * (x0 - 3)
    if y0 < -3:
        penalty += 1e6 * (abs(y0) - 3)
    elif y0 > 3:
        penalty += 1e6 * (y0 - 3)
    if z0 < -49.26:
        penalty += 1e6 * (abs(z0) + 49.26)
    elif z0 > 94.74:
        penalty += 1e6 * (z0 - 94.74)

    total_sum = 0.0

    for strip in strips:
        total_sum += linetopoint(params, strip)
    return total_sum + penalty

def FitPoints(points):
    # Use PCA to estimate the initial direction vector
    points_array = np.array([[p.x, p.y, p.z] for p in points])
    pca = PCA(n_components=1)
    pca.fit(points_array)
    direction = pca.components_[0]
    centroid = pca.mean_

    # Ensure the initial parameters are within the bounds
    x0, y0, z0 = centroid
    ux, uy, uz = direction

    # Clip the initial parameters to the bounds
    x0 = np.clip(x0, -3, 3)
    y0 = np.clip(y0, -3, 3)
    z0 = np.clip(z0, -49.26, 94.74)
    ux = np.clip(ux, -1, 1)
    uy = np.clip(uy, -1, 1)
    uz = np.clip(uz, -1, 1)

    initial_params = np.array([x0, y0, z0, ux, uy, uz])
    bounds = [(-3, 3), (-3, 3), (-49.26, 94.74), (-1, 1), (-1, 1), (-1, 1)]
    result = minimize(FitFunction, initial_params, args=(points,), method='Powell', bounds=bounds)

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

filename = 'build/10mm_space.csv'  # replace with your file name
result = extract_numbers(filename)
reco_z_err = []
count = 0

for block in result:
    points = []
    origin = Point3D(block[0][2], block[0][3], block[0][4])
    count = block[0][0]
    print(count)

    if count > 1000:
        break

    for line in block[1]:
        if line[0] == 0.0:
            continue
        point = Point3D(line[0], line[1], line[2])
        points.append(point)

    fitvec, fitcent = FitPoints(points)
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
print("integral =", integral, ", 5sig =", integral / count)

x = np.linspace(bin_centers[0], bin_centers[-1], 200)
y = gaussian(x, *popt)

plt.hist(reco_z_err, bins=200, range=(-100, 100), alpha=0.5, label='Histogram')
plt.plot(x, y, label='Fitted Gaussian')
plt.legend()
plt.show()
