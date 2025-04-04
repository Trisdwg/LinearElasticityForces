#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
 
class Mesh:
    def __init__(self, fname):
        self.fname = fname
        with open(fname, "r") as f:
            self.nnodes = int(f.readline().split(" ")[3])
            self.nodes = np.zeros((self.nnodes, 2))
            for i in range(self.nnodes):
                line = f.readline()
                parts = [s.strip() for s in line.split(':')]
                self.nodes[i] = [float(val.strip()) for val in parts[1].split()]
            line = f.readline()
            while("triangles" not in line and "quads" not in line):
                line = f.readline()
            self.nlocal = 3 if "triangles" in line else 4
            self.nelem = int(line.split(" ")[3])
            self.elem = np.zeros((self.nelem, self.nlocal), dtype=np.int32)
            for i in range(self.nelem):
                line = f.readline()
                parts = [s.strip() for s in line.split(':')]
                self.elem[i] = [int(val.strip()) for val in parts[1].split()]
 
    def plotfield(self, field, displace=None, *args, **kwargs):
        coord = self.nodes if displace is None else np.array(displace) + self.nodes
        field = np.array(field).ravel()
        if field.size == self.nnodes:
            a = plt.tripcolor(coord[:, 0], coord[:, 1], self.elem, field, shading="gouraud", *args, **kwargs)
        elif field.size == self.nelem:
            a = plt.tripcolor(coord[:, 0], coord[:, 1], self.elem, facecolors=field, shading="flat", *args, **kwargs)
        else:
            raise ValueError("Le champ doit être un tableau 1D dont la taille est égale à nnodes ou nelem.")
        return a
 
    def plot(self, displace=None, **kwargs):
        coord = self.nodes if displace is None else np.array(displace) + self.nodes
        plt.triplot(coord[:, 0], coord[:, 1], self.elem, **kwargs)
 
def load_stress_field(csv_fname):
    data = np.loadtxt(csv_fname, delimiter=",", skiprows=1)
    sigma_x = data[:, 1]
    sigma_y = data[:, 2]
    tau_xy  = data[:, 3]
    stress_norm = np.sqrt(sigma_x**2 + sigma_y**2 + 2 * tau_xy**2)
    return stress_norm
 
if __name__ == "__main__":
    mesh = Mesh("data/elasticity.txt")
    stress_norm = load_stress_field("data/stress.csv")
    plt.figure(figsize=(8,6))
    tpc = mesh.plotfield(stress_norm, cmap="jet", vmin=0, vmax=3e5)
    plt.colorbar(tpc, label="Norme des contraintes")
    mesh.plot(color="k", lw=0.5)
    plt.title("Champ de contraintes sur la structure")
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.gca().set_aspect("equal")
    plt.grid(True, alpha=0.3)
    plt.show()