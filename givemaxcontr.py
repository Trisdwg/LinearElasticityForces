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
 
def load_von_mises_field(csv_fname):
    """
    Charge les contraintes et calcule la contrainte équivalente de von Mises 2D.
    """
    data = np.loadtxt(csv_fname, delimiter=",", skiprows=1)
    sigma_x = data[:, 1]
    sigma_y = data[:, 2]
    tau_xy = data[:, 3]
 
    # Von Mises (2D déformations planes)
    sigma_vm = np.sqrt(sigma_x**2 - sigma_x * sigma_y + sigma_y**2 + 3 * tau_xy**2)
    return sigma_vm
 
if __name__ == "__main__":
    # Constante matériau
    limite_acier = 240e6  # 240 MPa
 
    mesh = Mesh("data/mesh2,7k.txt")
    sigma_vm = load_von_mises_field("build/stress.csv")
 
    vmax = np.max(sigma_vm) * 1.1  # on élargit un peu pour bien voir
    print(f"Contrainte von Mises max : {np.max(sigma_vm):.2e} Pa")
 
    