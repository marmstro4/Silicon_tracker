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

class Cylinder:
    def __init__(self, a, b, c, r, L, mod):
        self.a = a
        self.b = b
        self.c = c
        self.r = r
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

class Cylinder:
    def __init__(self, a, b, c, r, L, mod):
        self.a = a
        self.b = b
        self.c = c
        self.r = r
        self.L = L
        self.mod = mod

def linetocylinder(params, cylinder):
    x0, y0, z0, ux, uy, uz = params

    # Ensure the direction vector is normalized
    norm = math.sqrt(ux * ux + uy * uy + uz * uz)
    ux /= norm
    uy /= norm
    uz /= norm

    # Parametric line: p(t) = (x0 + t*ux, y0 + t*uy, z0 + t*uz)
    # Project the cylinder center onto the line
    a, b, c, r, L, modu = cylinder.a, cylinder.b, cylinder.c, cylinder.r, cylinder.L, cylinder.mod

    t = (a - x0) * ux + (b - y0) * uy + (c - z0) * uz

    # Closest point on the line
    px = x0 + t * ux
    py = y0 + t * uy
    pz = z0 + t * uz

    if modu == 0 or modu == 2:
        # Distance along the cylinder axis (x-axis)
        dx = px - a
        if ((dx<a-(L-vac_r-outer_end)) or (dx>a+(vac_r-inner_end))):
            # If out of cylinder bounds, penalize
            return 1e9

        # Perpendicular distance to the cylinder's surface
        dy = py - b
        dz = pz - c
        radial_distance = math.sqrt(dy * dy + dz * dz)

    elif modu == 1 or modu == 3:
        # Distance along the cylinder axis (x-axis)
        dy = py - b
        if ((dy<a-(L-vac_r-outer_end)) or (dy>b+(vac_r-inner_end))):
            # If out of cylinder bounds, penalize
            return 1e9

        # Perpendicular distance to the cylinder's surface
        dx = px - a
        dz = pz - c
        radial_distance = math.sqrt(dx * dx + dz * dz)

    # Return the squared difference from the radius
    return (radial_distance - r) * (radial_distance - r)

def FitFunction(params, cylinders):
    x0, y0 = params[:2]

    # Add penalties for violating the constraints
    penalty = 0.0
    if x0 < -0.3:
        penalty += 1e9 * (abs(x0) - 0.3)
    elif x0 > 0.3:
        penalty += 1e9 * (x0 - 0.3)
    if y0 < -0.3:
        penalty += 1e9 * (abs(y0) - 0.3)
    elif y0 > 0.3:
        penalty += 1e9 * (y0 - 0.3)

    total_sum = 0.0
    for cylinder in cylinders:
        total_sum += linetocylinder(params, cylinder)
    return total_sum + penalty

def FitRadii(cylinders, first):
    def fit_function(params):
        return FitFunction(params, cylinders)

    initial_params = np.array([first[0], first[1], first[2], 0, 0, 1])
    bounds = [(-0.2, 0.2), (-0.2, 0.2), (-7.5, 7.5), (-1, 1), (-1, 1), (-1, 1)]
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

filename = 'shift_4mm_4layers_gap.csv'  # replace with your file name
result = extract_numbers(filename)
reco_z_err = []
vac_r = 6.35
inner_end = 1
outer_end = 3

count = 0

for block in result:
    cylinders = []
    origin = Point3D(block[0][1],block[0][2],block[0][3])
    girth = block[0][7]
    girth = 17.55
    print("block ",block[0][0])
    count = block[0][0]

    #if count>1000:
        #break

    for line in block[1]:
        cylinder = Cylinder(line[0], line[1], line[2], line[3], girth, line[4])
        cylinders.append(cylinder)

    first = [0,0,1]
    min_diff = 1e6;

    for i in range(len(cylinders)):
        diff = math.sqrt(cylinders[i].a**2 + cylinders[i].b**2 + cylinders[i].c**2)
        if diff < min_diff:
            min_diff = diff
            norm = math.sqrt(cylinders[i].a**2+cylinders[i].a**2+cylinders[i].c**2)
            first = [cylinders[i].a, cylinders[i].b, cylinders[i].c]

    fitvec, fitcent = FitRadii([cylinders[i]],first)
    reco_v = findClosestPointOnLine(fitvec,fitcent,origin)
    reco_z_err.append(reco_v[2]-origin[2])


hist, bins = np.histogram(reco_z_err, bins=200, range=(-1, 1))
bin_centers = (bins[1:] + bins[:-1]) / 2
popt, pcov = curve_fit(gaussian, bin_centers, hist)

min = np.digitize(popt[1]-abs(3*popt[2]), bins)
max = np.digitize(popt[1]+abs(3*popt[2]), bins)

integral = 0

for i in range(max-min):
    integral = integral + hist[min + i]

print("mean = ",popt[1]*10," +/-",np.sqrt(pcov[1,1])*10," [mm]")
print("sigma = ",popt[2]*10," +/-",np.sqrt(pcov[2,2])*10," [mm]")
print("integral =",integral,", 3sig =",integral/count)


x = np.linspace(bin_centers[0], bin_centers[-1], 100)
y = gaussian(x, *popt)

plt.hist(reco_z_err, bins=200, range=(-1, 1), alpha=0.5, label='Histogram')
plt.plot(x, y, label='Fitted Gaussian')
plt.legend()
plt.show()





