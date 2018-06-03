# Filename: create_disk.py
# Author: Bailey Miller
# Description: Creates voxel data to model a protoplanetary disk.

import math
from subprocess import call

class Point:
  def __init__(self, x = 0, y = 0, z = 0):
    self.x = x
    self.y = y
    self.z = z

  def r(self):
    return math.sqrt(pow(self.x,2) + pow(self.y,2))

  def h(self):
    return self.z

  def __subtract__(self,p):
    return Point(self.x - p.x, self.y - p.y, self.z - p.z)

  def __str__(self):
    return "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"

class ProtoplanetaryDisk:
  def __init__(self, params = {}):
    default_params = {
      'mass' : 0.00057,
      'height' : 16.2,
      'R_in' : 0.1,
      'R_out' : 125,
      'alpha' : -1.77,
      'beta' : 1.19,
      'inclination' : 0,
      'pbrt_file' : '.disk.pbrt',
      'pbrt_path' : '../pbrt-v3/build/pbrt',
      'samples' : 128
    }
   
    for key, value in default_params.iteritems():
      if key not in params:
        params[key] = value

    self.pbrt_file = params['pbrt_file']
    self.pbrt_path = params['pbrt_path']
    self.samples = params['samples']

    self.mass = params['mass']
    self.R_in = params['R_in']
    self.R_out = params['R_out']
    self.alpha = params['alpha']
    self.beta = params['beta']
    self.inclination = params['inclination']

    self.height_o = params['height']
    self.sigma_o = (self.alpha + 2) 
    self.sigma_o *= pow(self.R_out,self.alpha) * self.mass
    self.sigma_o /= 2 * math.pi
    self.sigma_o /= (pow(self.R_out, self.alpha + 2)\
                     - pow(self.R_in, self.alpha + 2)) 

  def surface_density(self, p):
    return self.sigma_o * pow(p.r() / self.R_out, self.alpha)

  def height(self, p):
    return self.height_o * pow(p.r() / self.R_out, self.beta)

  def density(self, p):
    if p.r() > self.R_out or p.r() < self.R_in:
      return 0 
    elif abs(p.h()) > self.height(p):
      return 0
    else:
      return self.surface_density(p) / self.height(p) 

  def create_voxel_grid(self, nx = 10, ny = 10, nz = 10, 
                        write_to_file = False, create_pbrt_file = True):
    # voxel grid will be in a unity cube.
    # scale for x and y is self.R_out * 2
    # scale for z is max of height at R_out and R_in
    # (0,0,0) --> (0.5, 0.5, 0.5)
    max_height = max(self.height(Point(self.R_out, 0, 0)), 
                     self.height(Point(self.R_in, 0, 0)))
    voxel = [0] * (nx + 1) * (ny + 1) * (nz + 1)
    
    min_nonzero_density = float("inf")
    for k in range(nz + 1):
      for i in range(nx + 1):
        for j in range(ny + 1):
          x = float(i) / nx
          y = float(j) / ny
          z = float(k) / nz
          p = Point(
                   x * 2 * self.R_out - self.R_out, 
                   y * 2 * self.R_out - self.R_out, 
                   z * 2 * max_height - max_height)
          density = self.density(p)
          voxel[self.three_to_one(i, j, k, 
                                  nx + 1, ny + 1, nz + 1)] = density
          if not density == 0:
            min_nonzero_density = min(min_nonzero_density, density)
   
    for i in range(len(voxel)):
      voxel[i] /= min_nonzero_density

    if write_to_file:
      f = open("disk", 'w')
      f.write("[")
      for i in range(len(voxel)):
        f.write("{0}".format(voxel[i]))
      f.write("]")
      f.close()
    
    disk_str = "["
    for i in range(len(voxel)):
      disk_str +=  "  " + str(voxel[i]) + " "
    disk_str += "]"

    name = "disk_" + str(self.inclination) + "_" + str(self.height_o) + "_" + str(self.alpha) + "_" + str(self.beta) + ".exr"
    self.write_pbrt_file(disk_str, nx + 1, ny + 1, nz + 1, max_height / self.R_out, self.inclination, name, self.samples)

  def three_to_one(self, x, y, z, nx, ny, nz):
    return int(x + nx * (y + ny * z))

  def write_pbrt_file(self, voxel, nx, ny, nz, height = 0.5, 
                      inclination = 0, name = "out.exr", samples = 128,
                      sigma_s = 1, sigma_a = 1, 
                      xres=100, yres=100):
    pbrt = "LookAt\n\
      4 0 0 # Eye position\n\
      0 0 0       # look at point\n\
      0 0 1       # up vector\n\
    Camera \"perspective\" \"float fov\" 15\n\
    Sampler \"random\" \"integer pixelsamples\" {6}\n\
    Integrator \"volpath\"\n\
    Film \"image\"\
      \"string filename\" \"{10}\"\n\
      \"integer xresolution\" [{7}]\n\
      \"integer yresolution\" [{8}]\n\
    WorldBegin\n\
      # Star position\n\
      AttributeBegin\n\
        LightSource \"point\"\n\
          \"rgb I\" [250 250 250]\n\
      AttributeEnd\n\
      # Replace with dust eventually.\n\
      Rotate {11} 0 1 0\n\
      AttributeBegin\n\
        Material \"\"\n\
        MakeNamedMedium \"mygrid\"\n\
          \"string type\"  [\"heterogeneous\"]\n\
          \"point p0\" [-0.5 -0.5 -{9}]\n\
          \"point p1\" [0.5 0.5 {9}]\n\
          \"rgb sigma_s\" [{1} {1} {1}]\n\
          \"rgb sigma_a\" [{2} {2} {2}]\n\
          \"integer nx\" [{3}]\n\
          \"integer ny\" [{4}]\n\
          \"integer nz\" [{5}]\n\
          \"float density\" {0}\n\
          \"float g\" 0\n\
        MediumInterface  \"mygrid\" \"\"\n\
        Shape \"sphere\" \"float radius\" 1\n\
      AttributeEnd\n\
    WorldEnd\n\
    ".format(voxel, sigma_s, sigma_a, nx, ny, nz, samples, xres, yres, height, name, inclination)
    
    f = open(self.pbrt_file, 'w')
    f.write(pbrt)
    f.close()

  def render(self):
    call([self.pbrt_path, self.pbrt_file])

def create_edge_on_disk():
  disk = ProtoplanetaryDisk({'samples' : 1024})
  disk.create_voxel_grid()
  disk.render()

def create_multiple_inclinations():
  num_disks = 20
  for i in range(num_disks + 1):
    inclination = float(i) / num_disks * 90
    disk = ProtoplanetaryDisk({
      "inclination": inclination
    })
    disk.create_voxel_grid()
    disk.render()

def create_multiple_inclinations_top_to_side():
  num_disks = 20
  for i in range(num_disks + 1):
    inclination = 90 + float(i) / num_disks * 90
    disk = ProtoplanetaryDisk({
      "inclination": inclination
    })
    disk.create_voxel_grid()
    disk.render()

def create_multiple_height():
  num_disks = 6
  start_height = 10
  for i in range(num_disks):
    height = start_height + 5  * i
    disk = ProtoplanetaryDisk({
      "height": height
    })
    disk.create_voxel_grid()
    disk.render()

def create_multiple_alpha():
  alpha_base = -2.01
  num_disks = 20
  for i in range(num_disks):
    alpha = alpha_base + (i * 0.1)
    disk = ProtoplanetaryDisk({
      "alpha": alpha
    })
    disk.create_voxel_grid()
    disk.render()

def create_multiple_alpha_top_view():
  alpha_base = -2.01
  num_disks = 20
  for i in range(num_disks):
    alpha = alpha_base + (i * 0.1)
    disk = ProtoplanetaryDisk({
      "alpha": alpha,
      "inclination" : 90
    })
    disk.create_voxel_grid()
    disk.render()

def create_multiple_beta():
  beta_base = 1
  num_disks = 10
  for i in range(num_disks):
    beta = beta_base + (i * 0.1)
    disk = ProtoplanetaryDisk({
      "beta": beta
    })
    disk.create_voxel_grid()
    disk.render()


if __name__ == "__main__":
  #create_multiple_inclinations()
  create_multiple_inclinations_top_to_side()
  #create_multiple_height()
  #create_multiple_alpha_top_view()
  #create_multiple_beta()
  #create_edge_on_disk()
