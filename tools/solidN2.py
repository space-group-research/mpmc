#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#  @2013, Christian Cioce
#  Space Research Group
#  Department of Chemistry
#  University of South Florida
#  Email: ccioce@mail.usf.edu
#
#  This script has been used in modeling the crystalline phases of N2
#  and CO2, which has resulted in the following publications:
#
#  + C.R. Cioce, K. McLaughlin, J.L. Belof, B. Space, "A Polarizable and 
#    Transferable PHAST N2 Potential for use in Materials Simulation", 
#    J. Chem. Theory Comput., 2013, 9 (12), 5550–5557
#
#  - AND -
#
#  + A.L. Mullen, T. Pham, K.A. Forrest, C.R. Cioce, K. McLaughlin, B. Space, 
#    "A Polarizable and Transferable PHAST CO2 Potential for Materials Simulation", 
#    J. Chem. Theory Comput., 2013, 9 (12), 5421–5429
#
#  license: GNU LGPL
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  NOTE: The author kindly requests acknowledgement in any further use of
#        this script. Though optional, it is the right thing to do!

"""This script REQUIRES the input molecule to be diatomic. Homo/Heteronuclear independent."""
"""This version REQUIRES that the initial unit cell be cubic. For generalization, add conditionals to check min and max coordinates for each dimension independently."""

"""
For alphaN2 (Pa3(bar)):

PRIMITIVE VECTORS:
A1 = aX		a = 5.644 (recip: 0.17717931)
A2 = aY
A3 = aZ

BASIS VECTORS:
B1 =       uA1 + uA2 + uA3        =       uaX + uaY + uaZ        	(N) (8c)
B2 = (1/2-u)A1 - uA2 + (1/2+u)A3  = (1/2-u)aX - uaY + (1/2+u)aZ  	(N) (8c)
B3 = -uA1 + (1/2+u)A2 + (1/2-u)A3 = -uaX + (1/2+u)aY + (1/2-u)aZ 	(N) (8c)
B4 = (1/2+u)A1 + (1/2-u)A2 - uA3  = (1/2+u)aX + (1/2-u)aY - uaZ  	(N) (8c)
B5 =      -uA1 - uA2 - uA3        =      -uaX - uaY - uaZ        	(N) (8c)
B6 = (1/2+u)A1 + uA2 + (1/2-u)A3  = (1/2+u)aX + uaY + (1/2-u)aZ  	(N) (8c)
B7 = uA1 + (1/2-u)A2 + (1/2+u)A3  = uaX + (1/2-u)aY + (1/2+u)aZ  	(N) (8c)
B8 = (1/2-u)A1 + (1/2+u)A2 + uA3  = (1/2-u)aX + (1/2+u)aY + uaZ   	(N) (8c)
"""
"""
For gammaN2 (P4_2/mnm):

PRIMITIVE VECTORS:
A1 = aX		a = 3.957 (recip: 0.25271670)
A2 = aY
A3 = cZ		c = 5.109 (recip: 0.19573302)

BASIS VECTORS:
B1 =           uA1 + uA2             =          uaX + uaY		(N) (4f)
B2 =          -uA1 - uA2             =         -uaX - uaY		(N) (4f)
B3 = (1/2+u)A1 + (1/2-u)A2 + (1/2)A3 = (1/2+u)aX + (1/2-u)aY + (1/2)cZ	(N) (4f)
B4 = (1/2-u)A1 + (1/2+u)A2 + (1/2)A3 = (1/2-u)aX + (1/2+u)aY + (1/2)cZ	(N) (4f)
"""

import sys
import os
import time
from numpy import *
from scitools import *
from math import *
from subprocess import call
from scipy import *
from pylab import *
from shutil import move
import scitools.filetable as ft
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

########################################################################
# search a dictionary for key or value using named functions or a class
# tested with Python25   by Ene Uran    01/19/2008
# http://www.daniweb.com/software-development/python/code/217019
def find_key(dic, val):
    """return the key of dictionary dic given the value"""
    return [k for k, v in symbol_dic.iteritems() if v == val][0]

def find_value(dic, key):
    """return the value of dictionary dic given the key"""
    return dic[key]

class Lookup(dict):
    """
    a dictionary which can lookup value by key, or keys by value
    """
    def __init__(self, items=[]):
        """items can be a list of pair_lists or a dictionary"""
        dict.__init__(self, items)

    def get_key(self, value):
        """find the key(s) as a list given a value"""
        return [item[0] for item in self.items() if item[1] == value]

    def get_value(self, key):
        """find the value given a key"""
        return self[key]
########################################################################

def rmsd(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] numpy arrays"""

  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)

  n_vec = shape(crds1)[0]
  correlation_matrix = dot(transpose(crds1), crds2)
  v, s, w = linalg.svd(correlation_matrix)
  is_reflection = (linalg.det(v) * linalg.det(w)) < 0.0 
  if is_reflection:
    s[-1] = - s[-1]
  E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  return sqrt(rmsd_sq)

# Alpha Loop
def alpha_loop(mixrule, pol_on, n_cells, abc_min, abc_max, abc_del, x, breal, parentdir, trajfile, input_file, job_name, pqr_input, pqr_output, pqr_restart, energy_output, traj_output, dipole_output, field_output, sitelist, masslist, qlist, pollist, epslist, siglist, site, totalelems, totalrows):
  """Generates alpha-N2 crystals"""

  for a_current in arange(abc_min, abc_max, abc_del):
      trajfile.write("REMARK Trajectory File created in Python by Chris Cioce\n")
      trajfile.write("REMARK a = %f A\n" % (a_current))
      print "****************************************************"
      print "Current a=b=c = %f Angstroms" % a_current
      start_local = time.time()
      # Keep N cartesian coordinate static (specific code here...generalize)	<-- Determined by generating N2E sites and looking at xyz of N closest to basis origin
      N_fixed = 0.316965	# Reciprocal space coordinate
      u = N_fixed / float(a_current)

      # Generate Unit Cell from Basis Vectors (bv)
      d = zeros(x.size)
      e = zeros(x.size)
      f = zeros(x.size)

      # N(1) : B_4
      bvx = (0.5+u) * a_current
      bvy = (0.5-u) * a_current
      bvz =     -u  * a_current
      d.itemset(0,bvx)
      e.itemset(0,bvy)
      f.itemset(0,bvz)

      # N(2) : B_8
      bvx = (0.5-u) * a_current
      bvy = (0.5+u) * a_current
      bvz =      u  * a_current
      d.itemset(1,bvx)
      e.itemset(1,bvy)
      f.itemset(1,bvz)

      # N(3) : B_3
      bvx =     -u  * a_current
      bvy = (0.5+u) * a_current
      bvz = (0.5-u) * a_current
      d.itemset(2,bvx)
      e.itemset(2,bvy)
      f.itemset(2,bvz)

      # N(4) : B_7
      bvx =      u  * a_current
      bvy = (0.5-u) * a_current
      bvz = (0.5+u) * a_current
      d.itemset(3,bvx)
      e.itemset(3,bvy)
      f.itemset(3,bvz)

      # N(5) : B_1
      bvx =      u  * a_current
      bvy =      u  * a_current
      bvz =      u  * a_current
      d.itemset(4,bvx)
      e.itemset(4,bvy)
      f.itemset(4,bvz)

      # N(6) : B_5
      bvx =     -u  * a_current
      bvy =     -u  * a_current
      bvz =     -u  * a_current
      d.itemset(5,bvx)
      e.itemset(5,bvy)
      f.itemset(5,bvz)

            # N(7) : B_2
      bvx = (0.5-u) * a_current
      bvy =     -u  * a_current
      bvz = (0.5+u) * a_current
      d.itemset(6,bvx)
      e.itemset(6,bvy)
      f.itemset(6,bvz)

      # N(8) : B_6
      bvx = (0.5+u) * a_current
      bvy =      u  * a_current
      bvz = (0.5-u) * a_current
      d.itemset(7,bvx)
      e.itemset(7,bvy)
      f.itemset(7,bvz)

      x = d.copy()
      y = e.copy()
      z = f.copy()

      # ^^^^^^^^^^^^^^^^^^^^ STEP 4: Generate N2G sites (by midpoint formula) ^^^^^^^^^^^^^^^^^^^^ #
      a = zeros(x.size*0.5)
      b = zeros(x.size*0.5)
      c = zeros(x.size*0.5)
      j = 0

      for i in range(0, x.size, 2):
          xcom = ((x.item(i) + x.item(i+1)) * 0.5)
          ycom = ((y.item(i) + y.item(i+1)) * 0.5)
          zcom = ((z.item(i) + z.item(i+1)) * 0.5)
          a.itemset(j,round(xcom,6))
          b.itemset(j,round(ycom,6))
          c.itemset(j,round(zcom,6))
          j += 1

      # ^^^^^^^^^^^^^^^^^^^^ STEP 6: Generate N2N sites ^^^^^^^^^^^^^^^^^^^^ #
      blGN = 0.738					# LB nonpolar (default)
      if pol_on:
          blGN = 0.745					# Bond length from the G to the N site (for LB polar model)
      if mixrule == "WH":
          blGN = 0.789
      if mixrule == "WH" and pol_on:
          blGN = 0.791
      if S == "+S1" and mixrule == "LB":
          blGN = 0.790983
      if S == "+S1" and mixrule == "LB" and pol_on:
          blGN = 0.790811
      if S == "+S2" and mixrule == "LB":
          blGN = 0.788201
      if S == "+S2" and mixrule == "LB" and pol_on:
          blGN = 0.789864
      if S == "+S3" and mixrule == "LB":
          blGN = 0.788072
      if S == "+S3" and mixrule == "LB" and pol_on:
          blGN = 0.790914
      if S == "+S4" and mixrule == "LB":
          blGN = 0.788475
      if S == "+S4" and mixrule == "LB" and pol_on:
          blGN = 0.785869
      if S == "+S5" and mixrule == "LB":
          blGN = 0.788258
      if S == "+S5" and mixrule == "LB" and pol_on:
          blGN = 0.792170
      if S == "+S6" and mixrule == "LB":
          blGN = 0.787394
      if S == "+S6" and mixrule == "LB" and pol_on:
          blGN = 0.789624
      if S == "+S7" and mixrule == "LB":
          blGN = 0.786035
      if S == "+S7" and mixrule == "LB" and pol_on:
          blGN = 0.789760
      if S == "+S8" and mixrule == "LB":
          blGN = 0.788488
      if S == "+S8" and mixrule == "LB" and pol_on:
          blGN = 0.792476
      if S == "+S9" and mixrule == "LB":
          blGN = 0.788909
      if S == "+S9" and mixrule == "LB" and pol_on:
          blGN = 0.791495
      if S == "+S10" and mixrule == "LB":
          blGN = 0.788053
      if S == "+S10" and mixrule == "LB" and pol_on:
          blGN = 0.794897
      if S == "+S1" and mixrule == "WH":
          blGN = 0.784579
      if S == "+S1" and mixrule == "WH" and pol_on:
          blGN = 0.793307
      if S == "+S2" and mixrule == "WH":
          blGN = 0.787110
      if S == "+S2" and mixrule == "WH" and pol_on:
          blGN = 0.790668
      if S == "+S3" and mixrule == "WH":
          blGN = 0.783514
      if S == "+S3" and mixrule == "WH" and pol_on:
          blGN = 0.785596
      if S == "+S4" and mixrule == "WH":
          blGN = 0.788218
      if S == "+S4" and mixrule == "WH" and pol_on:
          blGN = 0.792883
      if S == "+S5" and mixrule == "WH":
          blGN = 0.785007
      if S == "+S5" and mixrule == "WH" and pol_on:
          blGN = 0.790113
      if S == "+S6" and mixrule == "WH":
          blGN = 0.785943
      if S == "+S6" and mixrule == "WH" and pol_on:
          blGN = 0.785050
      if S == "+S7" and mixrule == "WH":
          blGN = 0.784082
      if S == "+S7" and mixrule == "WH" and pol_on:
          blGN = 0.787391
      if S == "+S8" and mixrule == "WH":
          blGN = 0.784705
      if S == "+S8" and mixrule == "WH" and pol_on:
          blGN = 0.788553
      if S == "+S9" and mixrule == "WH":
          blGN = 0.788614
      if S == "+S9" and mixrule == "WH" and pol_on:
          blGN = 0.785022
      if S == "+S10" and mixrule == "WH":
          blGN = 0.786695
      if S == "+S10" and mixrule == "WH" and pol_on:
          blGN = 0.787558

      t = (blGN / (breal*0.5))
      N2N = zeros(x.size*3)
      N2N.shape = (x.size,3)
      j = 0
      k = 0
      for i in range(a.size):
          dx = a.item(i) - x.item(j)
          dy = b.item(i) - y.item(j)
          dz = c.item(i) - z.item(j)
          xnew = a.item(i) - t*dx
          ynew = b.item(i) - t*dy
          znew = c.item(i) - t*dz

          N2N.itemset((k,0),round(xnew,6))
          N2N.itemset((k,1),round(ynew,6))
          N2N.itemset((k,2),round(znew,6))
          k += 1

          xnew = a.item(i) + t*dx
          ynew = b.item(i) + t*dy
          znew = c.item(i) + t*dz

          N2N.itemset((k,0),round(xnew,6))
          N2N.itemset((k,1),round(ynew,6))
          N2N.itemset((k,2),round(znew,6))

          k += 1
          j += 2

      # Refill three_site with current x,y,z and a,b,c positions
      three_site = zeros(totalelems)
      three_site.shape = (totalrows,3)
      k = 0     # 0 --> 12 Total counter
      l = 0     # 0 -->  8 x,y,z counter
      for i in range(a.size):
          three_site.itemset((k,0),a.item(i))
          three_site.itemset((k,1),b.item(i))
          three_site.itemset((k,2),c.item(i))
          k += 1
          for j in range(2):
              three_site.itemset((k,0),x.item(l))
              three_site.itemset((k,1),y.item(l))
              three_site.itemset((k,2),z.item(l))
              l += 1
              k += 1
      totalrows = x.size + a.size + x.size              # 2nd x.size is b/c there are as many "N" sites are there are "E" sites
      totalelems = totalrows * 3
      five_site = zeros(totalelems)
      three_site = three_site.ravel()

      k = 0     # Total counter
      l = 0     # Counter for three_site
      m = 0     # Counter for N2N
      for i in range(a.size):
          for j in range(9):
              five_site.itemset(k,three_site.item(l))
              l += 1
              k += 1
          for j in range(6):
              five_site.itemset(k,N2N.item(m))
              m += 1
              k += 1

      five_site.shape = (totalrows,3)
      #print "\nDEBUG: 5-SITE N2 MODEL ~> EXPERIMENTAL DISTANCE OF %f ANGSTROMS" % breal
      #print five_site
      #print

      # Write unit cell to file
      for i in range(len(a)):
          for j in range(len(sitelist)):
              site.append(sitelist[j])

      abc_new = str(a_current)
      ac = 'unit_cell-' + abc_new + '.xyz'
      ucell = open(ac, 'w')
      ucell.write("%d\n\n" % len(five_site))
      for i in range(len(five_site)):
          ucell.write("%s\t%f\t%f\t%f\n" % (site[i],five_site.item(i,0), five_site.item(i,1), five_site.item(i,2)))
      ucell.close()
      site = []

      # ^^^^^^^^^^^^^^^^^^^^ STEP 7: Minor Translate ^^^^^^^^^^^^^^^^^^^^ #
      coor_min_x = float(five_site[:, 0].min())
      coor_min_y = float(five_site[:, 1].min())
      coor_min_z = float(five_site[:, 2].min())

      if (coor_min_x != 0 or coor_min_y != 0 or coor_min_z != 0):

          if coor_min_x < 0:        # Minor Translate: Ensure that a vertex is at (0,0,0)
              five_site[:, 0] += abs(coor_min_x)
          else:
              five_site[:, 0] -= coor_min_x

          if coor_min_y < 0:
              five_site[:, 1] += abs(coor_min_y)
          else:
              five_site[:, 1] -= coor_min_y

          if coor_min_z < 0:
              five_site[:, 2] += abs(coor_min_z)
          else:
              five_site[:, 2] -= coor_min_z

      # ^^^^^^^^^^^^^^^^^^^^ STEP 7: Generate SuperCell ^^^^^^^^^^^^^^^^^^^^ #
      # Determine initial a (ASSUMES A CUBIC LATTICE, SO a=b=c) XXX
      x_min = float(a.min())
      x_max = float(a.max())

      #print "DEBUG: AFTER MINOR TRANSLATE"
      #print five_site
      #print

      abc_orig = abs(x_min) + abs(x_max)
      delta = abc_orig * 2.0
      abc_orig = delta

      #if(abc_orig != float(a_current)):
      #    print "Error: Calculated lattice parameter a does not match the request (%f != %f). Exiting..." % (abc_orig, a_current)
      #    sys.exit(1)

      scelltmp_1 = open('scell_tmp.xyz', 'w')
      # Write initial cell to tmp file
      for i in range(len(five_site)):
          scelltmp_1.write("%f\t%f\t%f\n" % (five_site.item(i,0), five_site.item(i,1), five_site.item(i,2)))

      croot = (n_cells**(1/3.0))
      tstcroot = str(croot)
      tstcroot = tstcroot.split('.')
      if(len(tstcroot)==2 and tstcroot[1]=='0'):
          croot = int(round(croot,0))

          # Replicate unit cell in x-dimension
          for i in range(croot-1):
              five_site[:, 0] += delta
              for j in range(len(five_site)):
                  scelltmp_1.write("%f\t%f\t%f\n" % (five_site.item(j,0), five_site.item(j,1), five_site.item(j,2)))

          scelltmp_1.close()

          # Replicate supercell in y-dimension
          scelltmp_2 = open('scell_tmp_2.xyz', 'w')
          scelltmp_1 = open('scell_tmp.xyz', 'r')
          scelltmp = ft.read(scelltmp_1)
          scelltmp_1.close()
          for i in range(len(scelltmp)):
              scelltmp_2.write("%f\t%f\t%f\n" % (scelltmp.item(i,0), scelltmp.item(i,1), scelltmp.item(i,2)))

          for i in range(croot-1):
              scelltmp[:, 1] += delta
              for j in range(len(scelltmp)):
                  scelltmp_2.write("%f\t%f\t%f\n" % (scelltmp.item(j,0), scelltmp.item(j,1), scelltmp.item(j,2)))

          scelltmp_2.close()

          # Replicate supercell in z-dimension
          scelltmp_3 = open('scell_tmp_3.xyz', 'w')
          scelltmp_2 = open('scell_tmp_2.xyz', 'r')
          scelltmp = ft.read(scelltmp_2)
          scelltmp_2.close()
          for i in range(len(scelltmp)):
              scelltmp_3.write("%f\t%f\t%f\n" % (scelltmp.item(i,0), scelltmp.item(i,1), scelltmp.item(i,2)))

          for i in range(croot-1):
              scelltmp[:, 2] += delta
              for j in range(len(scelltmp)):
                  scelltmp_3.write("%f\t%f\t%f\n" % (scelltmp.item(j,0), scelltmp.item(j,1), scelltmp.item(j,2)))

          scelltmp_3.close()

          # Read fully replicated supercell into array
          scelltmp_3 = open('scell_tmp_3.xyz', 'r')
          scelltmp = ft.read(scelltmp_3)
          scelltmp_3.close()

          os.remove('scell_tmp.xyz')
          os.remove('scell_tmp_2.xyz')
          os.remove('scell_tmp_3.xyz')

      else:
          print "Not a cube root"
          sys.exit(1)

      # ^^^^^^^^^^^^^^^^^^^^ STEP 8: Center SuperCell at origin ^^^^^^^^^^^^^^^^^^^^ #
      coor_min_x = float(scelltmp[:, 0].min())
      coor_min_y = float(scelltmp[:, 1].min())
      coor_min_z = float(scelltmp[:, 2].min())
      coor_max_x = float(scelltmp[:, 0].max())
      coor_max_y = float(scelltmp[:, 1].max())
      coor_max_z = float(scelltmp[:, 2].max())

      check_x = (coor_min_x + coor_max_x) == 0.0
      check_y = (coor_min_y + coor_max_y) == 0.0
      check_z = (coor_min_z + coor_max_z) == 0.0

      if not (check_x or check_y or check_z):
          #print "Centering supercell at origin.\n"
          half_max_x = coor_max_x * 0.5
          half_max_y = coor_max_y * 0.5
          half_max_z = coor_max_z * 0.5

          if coor_max_x > 0:
              scelltmp[:, 0] -= half_max_x
          else:
              scelltmp[:, 0] += abs(half_max_x)

          if coor_max_y > 0:
              scelltmp[:, 1] -= half_max_y
          else:
              scelltmp[:, 1] += abs(half_max_y)

          if coor_max_z > 0:
              scelltmp[:, 2] -= half_max_z
          else:
              scelltmp[:, 2] += abs(half_max_z)

      # DEBUG 
      #print
      #for i in range(len(scelltmp)):
      #    print "%f\t%f\t%f" % (scelltmp.item(i,0), scelltmp.item(i,1), scelltmp.item(i,2))

      # ^^^^^^^^^^^^^^^^^^^^ STEP 8: Scale lattice parameters ^^^^^^^^^^^^^^^^^^^^ #
      nmols = int(len(scelltmp)/5.0)            # Based on a 5-site model
      for i in range(nmols):
          for j in range(len(sitelist)):
              site.append(sitelist[j])
              molnum.append(i+1)
              mass.append(masslist[j])
              charge.append(qlist[j])
              polar.append(pollist[j])
              eps.append(epslist[j])
              sig.append(siglist[j])

      scell = scelltmp.copy()

      abclist.append(a_current)
      basis = a_current * croot         # XXX This assumes that the number of requested supercells has a cubic root
      print "Basis = %f" % basis

      # ^^^^^^^^^^^^^^^^^^^^ STEP 8: Create .pqr file for MPMC ^^^^^^^^^^^^^^^^^^^^ #
      if not os.path.exists(abc_new):
          os.mkdir(abc_new)

      os.chdir(abc_new)
      dircurr = os.getcwd()
      # Move unit cell file to pwd
      move(parentdir + '/' + ac, dircurr + '/' + ac)

      coord_file = open(pqr_input, 'w')
      raw_coords = open('coords.dat', 'w')
      # NOTE: MAX NUMBER OF DECIMALS IN COORDS ALLOWED IS 3 (%8.3f, ideally %11.6f, or %21.16f).
      for i in range(len(scell)):
          raw_coords.write("%11.6f\t%11.6f\t%11.6f\n" % (scell.item(i,0),scell.item(i,1),scell.item(i,2)))
          coord_file.write("ATOM  %5d %-4.45s %-3.3s %-1.1s %4d   %11.6f%11.6f%11.6f%9.5f%9.5f%9.5f%9.5f%9.5f%9.5f\n" % (i+1,site[i],'N2','M',molnum[i],scell.item(i,0),scell.item(i,1),scell.item(i,2),mass[i],charge[i],polar[i],eps[i],sig[i],0.0))

      coord_file.write("END")
      raw_coords.close()
      coord_file.close()

      # RMSD
      crds_list.append(scell)           # Store [nx3] numpy array in list
      if len(crds_list) == 2:           # We have 2 structures to compute RMSD
          # Do RMSD
          crds1 = crds_list[0]
          crds2 = crds_list[1]
          assert(crds1.shape[1] == 3)
          assert(crds1.shape == crds2.shape)
          n_vec = shape(crds1)[0]
          correlation_matrix = dot(transpose(crds1), crds2)
          v, s, w = linalg.svd(correlation_matrix)
          is_reflection = (linalg.det(v) * linalg.det(w)) < 0.0
          if is_reflection:
            s[-1] = - s[-1]
          E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
          rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
          rmsd_sq = max([rmsd_sq, 0.0])
          rmsd = sqrt(rmsd_sq)
          rmsd_list.append(rmsd)
          print "RMSD: %f" % rmsd
          crds_list.pop(0)              # Remove old coordinates, sets current coords to list position 0

      # Create input file
      seed = random_integers(100000000,999999999)
      infile = open(input_file, 'w')
      infile.write('\n')
      infile.write('job_name\t\t%s\n' % job_name)
      infile.write('\n')
      infile.write('ensemble\t\ttotal_energy\n')
      if mixrule == "WH":
          infile.write('waldmanhagler\t\ton\n')
      infile.write('\n')
      infile.write('numsteps\t\t1\n')
      infile.write('corrtime\t\t1\n')
      infile.write('\n')
      infile.write('rd_lrc\t\t\ton\n')
      infile.write('rd_crystal\t\ton\n')
      infile.write('rd_crystal_order\t10\n')
      infile.write('\n')
      infile.write('basis1\t\t\t%f\t0.0\t\t0.0\n' % basis)
      infile.write('basis2\t\t\t0.0\t\t%f\t0.0\n' % basis)
      infile.write('basis3\t\t\t0.0\t\t0.0\t\t%f\n' % basis)
      infile.write('\n')
      if pol_on:
          infile.write('polarization\t\ton\n')
          infile.write('polar_damp_type\t\texponential\n')
          infile.write('polar_damp\t\t2.1304\n')
          infile.write('\n')
          infile.write('polar_ewald\t\ton\n')
          infile.write('ewald_kmax\t\t10\n')
          #infile.write('polar_wolf_alpha\t0.13\n')
          infile.write('polar_palmo\t\ton\n')
          infile.write('polar_gs_ranked\t\ton\n')
          infile.write('polar_iterative\t\ton\n')
          infile.write('polar_max_iter\t\t10\n')
          infile.write('polar_gamma\t\t1.03\n')
          infile.write('\n')
      infile.write('wrapall\t\t\ton\n')
      #infile.write('pqr_input\t\t%s\n' % pqr_input)
      #infile.write('pqr_output\t\t%s\n' % pqr_output)
      infile.write('pqr_restart\t\t%s\n' % pqr_restart)
      #infile.write('energy_output\t\t%s\n' % energy_output)
      infile.write('traj_output\t\t%s\n' % traj_output)
      if pol_on:
          #infile.write('dipole_output\t\t%s\n' % dipole_output)
          #infile.write('field_output\t\t%s\n' % field_output)
          infile.write('\n')
      infile.write('\n')
      infile.close()

      # Execute MPMC
      #call("mpmc " + input_file + " >& runlog.log", shell=True)

      os.chdir(parentdir)

      # Print timing statistics


# Gamma Loop
def gamma_loop(mixrule, pol_on, n_cells, abc_min, abc_max, abc_del, x, breal, parentdir, trajfile, input_file, job_name, pqr_input, pqr_output, pqr_restart, energy_output, traj_output, dipole_output, field_output, sitelist, masslist, qlist, pollist, epslist, siglist, site, totalelems, totalrows):
  """Generates gamma-N2 crystals"""

  # First, statically set c lattice
  a2c_ratio = 3.957/5.109	# volume ratio between the a and c lattice constants, and is always our starting point

  # Must decide the initial parameters, depending on +/-
  if inc_dec == "-":
      a_current = abc_max
      abc_max = 3.957
      c_current = a_current / a2c_ratio
  else:
      a_current = 3.957     # via experiment
      c_current = 5.109     # via experiment

  for a_current in arange(a_current, abc_max, abc_del):
      trajfile.write("REMARK Trajectory File created in Python by Chris Cioce\n")
      trajfile.write("REMARK a=b = %f A, c = %f A\n" % (a_current, c_current))
      print "****************************************************"
      print "Current a=b = %f Ang / c = %f Ang" % (a_current, c_current)
      start_local = time.time()
      # Keep N cartesian coordinate static (specific code here...generalize)	<-- Determined by generating N2E sites and looking at xyz of N closest to basis origin
      N_fixed = 0.388202	# Reciprocal space coordinate
      u = N_fixed / float(a_current)

      # Generate Unit Cell from Basis Vectors (bv)
      d = zeros(x.size)
      e = zeros(x.size)
      f = zeros(x.size)

      # N(1) : B_1
      bvx = u * a_current
      bvy = bvx
      bvz = 0.0
      d.itemset(0,bvx)
      e.itemset(0,bvy)
      f.itemset(0,bvz)

      # N(2) : B_2
      bvx = -bvx
      bvy =  bvx
      bvz =  0.0
      d.itemset(1,bvx)
      e.itemset(1,bvy)
      f.itemset(1,bvz)

      # N(3) : B_3
      bvx = (0.5+u) * a_current
      bvy = (0.5-u) * a_current
      bvz =  0.5    * c_current
      d.itemset(2,bvx)
      e.itemset(2,bvy)
      f.itemset(2,bvz)

      # N(4) : B_4
      bvx = (0.5-u) * a_current
      bvy = (0.5+u) * a_current
      bvz =  0.5    * c_current
      d.itemset(3,bvx)
      e.itemset(3,bvy)
      f.itemset(3,bvz)

      x = d.copy()
      y = e.copy()
      z = f.copy()

      # ^^^^^^^^^^^^^^^^^^^^ STEP 4: Generate N2G sites (by midpoint formula) ^^^^^^^^^^^^^^^^^^^^ #
      a = zeros(x.size*0.5)
      b = zeros(x.size*0.5)
      c = zeros(z.size*0.5)
      j = 0

      for i in range(0, x.size, 2):
          xcom = ((x.item(i) + x.item(i+1)) * 0.5)
          ycom = ((y.item(i) + y.item(i+1)) * 0.5)
          zcom = ((z.item(i) + z.item(i+1)) * 0.5)
          a.itemset(j,round(xcom,6))
          b.itemset(j,round(ycom,6))
          c.itemset(j,round(zcom,6))
          j += 1

      # ^^^^^^^^^^^^^^^^^^^^ STEP 6: Generate N2N sites ^^^^^^^^^^^^^^^^^^^^ #
      blGN = 0.738					# LB nonpolar (default)
      if pol_on:
          blGN = 0.745					# Bond length from the G to the N site (for LB polar model)
      if mixrule == "WH":
          blGN = 0.789
      if mixrule == "WH" and pol_on:
          blGN = 0.791
      if S == "+S1" and mixrule == "LB":
          blGN = 0.790983
      if S == "+S1" and mixrule == "LB" and pol_on:
          blGN = 0.790811
      if S == "+S2" and mixrule == "LB":
          blGN = 0.788201
      if S == "+S2" and mixrule == "LB" and pol_on:
          blGN = 0.789864
      if S == "+S3" and mixrule == "LB":
          blGN = 0.788072
      if S == "+S3" and mixrule == "LB" and pol_on:
          blGN = 0.790914
      if S == "+S4" and mixrule == "LB":
          blGN = 0.788475
      if S == "+S4" and mixrule == "LB" and pol_on:
          blGN = 0.785869
      if S == "+S5" and mixrule == "LB":
          blGN = 0.788258
      if S == "+S5" and mixrule == "LB" and pol_on:
          blGN = 0.792170
      if S == "+S6" and mixrule == "LB":
          blGN = 0.787394
      if S == "+S6" and mixrule == "LB" and pol_on:
          blGN = 0.789624
      if S == "+S7" and mixrule == "LB":
          blGN = 0.786035
      if S == "+S7" and mixrule == "LB" and pol_on:
          blGN = 0.789760
      if S == "+S8" and mixrule == "LB":
          blGN = 0.788488
      if S == "+S8" and mixrule == "LB" and pol_on:
          blGN = 0.792476
      if S == "+S9" and mixrule == "LB":
          blGN = 0.788909
      if S == "+S9" and mixrule == "LB" and pol_on:
          blGN = 0.791495
      if S == "+S10" and mixrule == "LB":
          blGN = 0.788053
      if S == "+S10" and mixrule == "LB" and pol_on:
          blGN = 0.794897
      if S == "+S1" and mixrule == "WH":
          blGN = 0.784579
      if S == "+S1" and mixrule == "WH" and pol_on:
          blGN = 0.793307
      if S == "+S2" and mixrule == "WH":
          blGN = 0.787110
      if S == "+S2" and mixrule == "WH" and pol_on:
          blGN = 0.790668
      if S == "+S3" and mixrule == "WH":
          blGN = 0.783514
      if S == "+S3" and mixrule == "WH" and pol_on:
          blGN = 0.785596
      if S == "+S4" and mixrule == "WH":
          blGN = 0.788218
      if S == "+S4" and mixrule == "WH" and pol_on:
          blGN = 0.792883
      if S == "+S5" and mixrule == "WH":
          blGN = 0.785007
      if S == "+S5" and mixrule == "WH" and pol_on:
          blGN = 0.790113
      if S == "+S6" and mixrule == "WH":
          blGN = 0.785943
      if S == "+S6" and mixrule == "WH" and pol_on:
          blGN = 0.785050
      if S == "+S7" and mixrule == "WH":
          blGN = 0.784082
      if S == "+S7" and mixrule == "WH" and pol_on:
          blGN = 0.787391
      if S == "+S8" and mixrule == "WH":
          blGN = 0.784705
      if S == "+S8" and mixrule == "WH" and pol_on:
          blGN = 0.788553
      if S == "+S9" and mixrule == "WH":
          blGN = 0.788614
      if S == "+S9" and mixrule == "WH" and pol_on:
          blGN = 0.785022
      if S == "+S10" and mixrule == "WH":
          blGN = 0.786695
      if S == "+S10" and mixrule == "WH" and pol_on:
          blGN = 0.787558

      t = (blGN / (breal*0.5))
      N2N = zeros(x.size*3)
      N2N.shape = (x.size,3)
      j = 0 
      k = 0 
      for i in range(a.size):
          dx = a.item(i) - x.item(j)
          dy = b.item(i) - y.item(j)
          dz = c.item(i) - z.item(j)
          xnew = a.item(i) - t*dx
          ynew = b.item(i) - t*dy
          znew = c.item(i) - t*dz

          N2N.itemset((k,0),round(xnew,6))
	  N2N.itemset((k,1),round(ynew,6))
	  N2N.itemset((k,2),round(znew,6))
	  k += 1

          xnew = a.item(i) + t*dx
          ynew = b.item(i) + t*dy
          znew = c.item(i) + t*dz

	  N2N.itemset((k,0),round(xnew,6))
	  N2N.itemset((k,1),round(ynew,6))
	  N2N.itemset((k,2),round(znew,6))

          k += 1
          j += 2

      # Refill three_site with current x,y,z and a,b,c positions
      three_site = zeros(totalelems)
      three_site.shape = (totalrows,3)
      k = 0	# 0 --> 12 Total counter
      l = 0	# 0 -->  8 x,y,z counter
      for i in range(a.size):
	  three_site.itemset((k,0),a.item(i))
	  three_site.itemset((k,1),b.item(i))
	  three_site.itemset((k,2),c.item(i))
	  k += 1
	  for j in range(2):
              three_site.itemset((k,0),x.item(l))
	      three_site.itemset((k,1),y.item(l))
	      three_site.itemset((k,2),z.item(l))
	      l += 1
	      k += 1

      totalrows = x.size + a.size + x.size		# 2nd x.size is b/c there are as many "N" sites are there are "E" sites
      totalelems = totalrows * 3
      five_site = zeros(totalelems)
      three_site = three_site.ravel()
      
      k = 0	# Total counter
      l = 0	# Counter for three_site
      m = 0	# Counter for N2N
      for i in range(a.size):
          for j in range(9):
              five_site.itemset(k,three_site.item(l))
	      l += 1
              k += 1
          for j in range(6):
              five_site.itemset(k,N2N.item(m))
	      m += 1
              k += 1

      five_site.shape = (totalrows,3)
      #print "\nDEBUG: 5-SITE N2 MODEL ~> EXPERIMENTAL DISTANCE OF %f ANGSTROMS" % breal
      #print five_site
      #print

      # Write unit cell to file
      for i in range(len(a)):
          for j in range(len(sitelist)):
              site.append(sitelist[j])

      abc_new = str(a_current)
      c_updated = str(c_current)	#@CRC
      ac = 'unit_cell-' + abc_new + '-' + c_updated + '.xyz'
      ucell = open(ac, 'w')
      ucell.write("%d\n\n" % len(five_site))
      for i in range(len(five_site)):
          ucell.write("%s\t%f\t%f\t%f\n" % (site[i],five_site.item(i,0), five_site.item(i,1), five_site.item(i,2)))
      ucell.close()
      site = []

      # ^^^^^^^^^^^^^^^^^^^^ STEP 7: Minor Translate ^^^^^^^^^^^^^^^^^^^^ #
      coor_min_x = float(five_site[:, 0].min())
      coor_min_y = float(five_site[:, 1].min())
      coor_min_z = float(five_site[:, 2].min())

      if (coor_min_x != 0 or coor_min_y != 0 or coor_min_z != 0): 

          if coor_min_x < 0:        # Minor Translate: Ensure that a vertex is at (0,0,0)
              five_site[:, 0] += abs(coor_min_x)
          else:
              five_site[:, 0] -= coor_min_x

          if coor_min_y < 0:
              five_site[:, 1] += abs(coor_min_y)
          else:
              five_site[:, 1] -= coor_min_y

          if coor_min_z < 0:
              five_site[:, 2] += abs(coor_min_z)
          else:
              five_site[:, 2] -= coor_min_z

      # ^^^^^^^^^^^^^^^^^^^^ STEP 7: Generate SuperCell ^^^^^^^^^^^^^^^^^^^^ #
      # Determine initial a, b and c
      x_min = float(a.min())
      x_max = float(a.max())
      y_min = float(b.min())
      y_max = float(b.max())
      z_min = float(c.min())
      z_max = float(c.max())

      #print "DEBUG: AFTER MINOR TRANSLATE"
      #print five_site
      #print

      a_orig = abs(x_min) + abs(x_max)
      b_orig = abs(y_min) + abs(y_max)
      c_orig = abs(z_min) + abs(z_max)
      delta_a = a_orig * 2.0
      delta_b = b_orig * 2.0
      delta_c = c_orig * 2.0
      a_orig = delta_a
      b_orig = delta_b
      c_orig = delta_c

      scelltmp_1 = open('scell_tmp.xyz', 'w')
      # Write initial cell to tmp file
      for i in range(len(five_site)):
          scelltmp_1.write("%f\t%f\t%f\n" % (five_site.item(i,0), five_site.item(i,1), five_site.item(i,2)))

      croot = (n_cells**(1/3.0))
      tstcroot = str(croot)
      tstcroot = tstcroot.split('.')
      if(len(tstcroot)==2 and tstcroot[1]=='0'):
          croot = int(round(croot,0))

          # Replicate unit cell in x-dimension
          for i in range(croot-1):
	      five_site[:, 0] += delta_a
	      for j in range(len(five_site)):
	          scelltmp_1.write("%f\t%f\t%f\n" % (five_site.item(j,0), five_site.item(j,1), five_site.item(j,2)))

          scelltmp_1.close()

          # Replicate supercell in y-dimension
          scelltmp_2 = open('scell_tmp_2.xyz', 'w')
          scelltmp_1 = open('scell_tmp.xyz', 'r')
          scelltmp = ft.read(scelltmp_1)
          scelltmp_1.close()
          for i in range(len(scelltmp)):
	      scelltmp_2.write("%f\t%f\t%f\n" % (scelltmp.item(i,0), scelltmp.item(i,1), scelltmp.item(i,2)))

          for i in range(croot-1):
	      scelltmp[:, 1] += delta_b
	      for j in range(len(scelltmp)):
	          scelltmp_2.write("%f\t%f\t%f\n" % (scelltmp.item(j,0), scelltmp.item(j,1), scelltmp.item(j,2)))

          scelltmp_2.close()

          # Replicate supercell in z-dimension
          scelltmp_3 = open('scell_tmp_3.xyz', 'w')
          scelltmp_2 = open('scell_tmp_2.xyz', 'r')
          scelltmp = ft.read(scelltmp_2)
          scelltmp_2.close()
          for i in range(len(scelltmp)):
	      scelltmp_3.write("%f\t%f\t%f\n" % (scelltmp.item(i,0), scelltmp.item(i,1), scelltmp.item(i,2)))

          for i in range(croot-1):
	      scelltmp[:, 2] += delta_c
	      for j in range(len(scelltmp)):
	          scelltmp_3.write("%f\t%f\t%f\n" % (scelltmp.item(j,0), scelltmp.item(j,1), scelltmp.item(j,2)))

          scelltmp_3.close()

          # Read fully replicated supercell into array
          scelltmp_3 = open('scell_tmp_3.xyz', 'r')
          scelltmp = ft.read(scelltmp_3)
          scelltmp_3.close()

          os.remove('scell_tmp.xyz')
          os.remove('scell_tmp_2.xyz')
          os.remove('scell_tmp_3.xyz')

      else:
          print "Not a cube root"
          sys.exit(1)

      # ^^^^^^^^^^^^^^^^^^^^ STEP 8: Center SuperCell at origin ^^^^^^^^^^^^^^^^^^^^ #
      coor_min_x = float(scelltmp[:, 0].min())
      coor_min_y = float(scelltmp[:, 1].min())
      coor_min_z = float(scelltmp[:, 2].min())
      coor_max_x = float(scelltmp[:, 0].max())
      coor_max_y = float(scelltmp[:, 1].max())
      coor_max_z = float(scelltmp[:, 2].max())

      check_x = (coor_min_x + coor_max_x) == 0.0 
      check_y = (coor_min_y + coor_max_y) == 0.0 
      check_z = (coor_min_z + coor_max_z) == 0.0 

      if not (check_x or check_y or check_z):
          #print "Centering supercell at origin.\n"
          half_max_x = coor_max_x * 0.5 
          half_max_y = coor_max_y * 0.5 
          half_max_z = coor_max_z * 0.5 

          if coor_max_x > 0:
              scelltmp[:, 0] -= half_max_x
          else:
              scelltmp[:, 0] += abs(half_max_x)

          if coor_max_y > 0:
              scelltmp[:, 1] -= half_max_y
          else:
              scelltmp[:, 1] += abs(half_max_y)

          if coor_max_z > 0:
              scelltmp[:, 2] -= half_max_z
          else:
              scelltmp[:, 2] += abs(half_max_z)

      # DEBUG 
      #print
      #for i in range(len(scelltmp)):
      #    print "%f\t%f\t%f" % (scelltmp.item(i,0), scelltmp.item(i,1), scelltmp.item(i,2))

      # ^^^^^^^^^^^^^^^^^^^^ STEP 8: Scale lattice parameters ^^^^^^^^^^^^^^^^^^^^ #
      nmols = int(len(scelltmp)/5.0)		# Based on a 5-site model
      for i in range(nmols):
          for j in range(len(sitelist)):
              site.append(sitelist[j])
              molnum.append(i+1)
              mass.append(masslist[j])
              charge.append(qlist[j])
              polar.append(pollist[j])
              eps.append(epslist[j])
              sig.append(siglist[j])

      scell = scelltmp.copy()
    
      abclist.append(a_current)
      basis_c  = c_current * croot
      basis_ab = a_current * croot		# XXX This assumes that the number of requested supercells has a cubic root
      print "Basis_ab = %f / c = %f" % (basis_ab, basis_c)
    
      # ^^^^^^^^^^^^^^^^^^^^ STEP 8: Create .pqr file for MPMC ^^^^^^^^^^^^^^^^^^^^ #
      if not os.path.exists(abc_new):
	  os.mkdir(abc_new)

      os.chdir(abc_new)
      dircurr = os.getcwd()
      # Move unit cell file to pwd
      move(parentdir + '/' + ac, dircurr + '/' + ac)

      coord_file = open(pqr_input, 'w')
      raw_coords = open('coords.dat', 'w')
      # NOTE: MAX NUMBER OF DECIMALS IN COORDS ALLOWED IS 3 (%8.3f, ideally %11.6f, or %21.16f).
      for i in range(len(scell)):
	  raw_coords.write("%11.6f\t%11.6f\t%11.6f\n" % (scell.item(i,0),scell.item(i,1),scell.item(i,2)))
          coord_file.write("ATOM  %5d %-4.45s %-3.3s %-1.1s %4d   %11.6f%11.6f%11.6f%9.5f%9.5f%9.5f%9.5f%9.5f%9.5f\n" % (i+1,site[i],'N2','M',molnum[i],scell.item(i,0),scell.item(i,1),scell.item(i,2),mass[i],charge[i],polar[i],eps[i],sig[i],0.0))
    
      coord_file.write("END")
      raw_coords.close()
      coord_file.close()

      # RMSD
      crds_list.append(scell)		# Store [nx3] numpy array in list
      if len(crds_list) == 2:		# We have 2 structures to compute RMSD
          # Do RMSD
	  crds1 = crds_list[0]
	  crds2 = crds_list[1]
	  assert(crds1.shape[1] == 3)
	  assert(crds1.shape == crds2.shape)
	  n_vec = shape(crds1)[0]
	  correlation_matrix = dot(transpose(crds1), crds2)
	  v, s, w = linalg.svd(correlation_matrix)
	  is_reflection = (linalg.det(v) * linalg.det(w)) < 0.0
	  if is_reflection:
	    s[-1] = - s[-1]
	  E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
	  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
	  rmsd_sq = max([rmsd_sq, 0.0])
	  rmsd = sqrt(rmsd_sq)
	  rmsd_list.append(rmsd)
	  print "RMSD: %f" % rmsd
          crds_list.pop(0)		# Remove old coordinates, sets current coords to list position 0
      
      # Create input file
      seed = random_integers(100000000,999999999)
      infile = open(input_file, 'w')
      infile.write('\n')
      infile.write('job_name\t\t%s\n' % job_name)
      infile.write('\n')
      infile.write('ensemble\t\ttotal_energy\n')
      if mixrule == "WH":
          infile.write('waldmanhagler\t\ton\n')
      infile.write('\n')
      infile.write('numsteps\t\t1\n')
      infile.write('corrtime\t\t1\n')
      infile.write('\n')
      infile.write('rd_lrc\t\t\ton\n')
      infile.write('rd_crystal\t\ton\n')
      infile.write('rd_crystal_order\t10\n')
      infile.write('\n')
      infile.write('basis1\t\t\t%f\t0.0\t\t0.0\n' % basis_ab)
      infile.write('basis2\t\t\t0.0\t\t%f\t0.0\n' % basis_ab)
      infile.write('basis3\t\t\t0.0\t\t0.0\t\t%f\n' % basis_c)
      infile.write('\n')
      if pol_on:
	  infile.write('polarization\t\ton\n')
	  infile.write('polar_damp_type\t\texponential\n')
	  infile.write('polar_damp\t\t2.1304\n')
	  infile.write('\n')
	  infile.write('polar_ewald\t\ton\n')
	  infile.write('ewald_kmax\t\t10\n')
	  #infile.write('polar_wolf_alpha\t0.13\n')
	  infile.write('polar_palmo\t\ton\n')
          infile.write('polar_gs_ranked\t\ton\n')
	  infile.write('polar_iterative\t\ton\n')
	  infile.write('polar_max_iter\t\t10\n')
          infile.write('polar_gamma\t\t1.03\n')
	  infile.write('\n')
      infile.write('wrapall\t\t\ton\n')
      #infile.write('pqr_input\t\t%s\n' % pqr_input)
      #infile.write('pqr_output\t\t%s\n' % pqr_output)
      infile.write('pqr_restart\t\t%s\n' % pqr_restart)
      #infile.write('energy_output\t\t%s\n' % energy_output)
      infile.write('traj_output\t\t%s\n' % traj_output)
      if pol_on:
	  #infile.write('dipole_output\t\t%s\n' % dipole_output)
	  #infile.write('field_output\t\t%s\n' % field_output)
	  infile.write('\n')
      infile.write('\n')
      infile.close()

      # Execute MPMC
      #call("mpmc " + input_file + " >& runlog.log", shell=True)

      os.chdir(parentdir)

      # Print timing statistics

      c_current = (a_current + abc_del) / a2c_ratio


if __name__ == '__main__':
  print
  # ^^^^^^^^^^^^^^^^^^^^ STEP 1: Load in original coordinates ^^^^^^^^^^^^^^^^^^^^ #
  if (len(sys.argv) > 1):
      if (len(sys.argv) < 9):
  	  print "ERROR: Insufficient number of arguments. Usage: $ python solidN2.py [alpha | gamma] [np | pol] [LB | WH] [# UCs] [min_dist] [max_dist] [increment] [+ | -] [+SX | -SX]"
	  sys.exit(1)

      phase   = sys.argv[1]
      pol_on  = sys.argv[2]
      mixrule = sys.argv[3]
      n_cells = int(sys.argv[4])
      abc_min = float(sys.argv[5])
      abc_max = float(sys.argv[6])
      abc_del = float(sys.argv[7])
      inc_dec = sys.argv[8]
      S       = sys.argv[9]

      if phase == "alpha":
          print "Alpha-N2 crystal will be generated."
      elif phase == "gamma":
          print "Gamma-N2 crystal will be generated."
      else:
          print "Specified phase is not compatible. Sorry, try again."
          sys.exit(1)

      if mixrule == "LB":
          print "Lorentz-Berthelot mixing rules will be used."
      elif mixrule == "WH":
          print "Waldman-Hagler mixing rules will be used."
      else:
          print "Specified mixing rule is not compatible. Sorry, try again."
          sys.exit(1)

      if pol_on == "pol":
	  pol_on = 1    # Turn on polarization
	  print "Polar job requested."
      else:
	  pol_on = None # Turn off polarization
	  print "Nonpolar job requested."
  else:
      print "ERROR: Insufficient number of arguments. Usage: $ python solidN2.py [alpha | gamma] [np | pol] [LB | WH] [# UCs] [min_dist] [max_dist] [increment] [+ | -] [+SX | -SX]"
      sys.exit(1)

  ncycles = (round(abc_max,0) - abc_min) / float(abc_del)               # Number of cycles: not perfect, but won't mess up important things, just timing stats

  if phase == "alpha":
      ifile = open('/home/ccioce/N2/Solid/script/alphaN2.UC', 'r')
  elif phase == "gamma":
      ifile = open('/home/ccioce/N2/Solid/script/gammaN2.UC', 'r')

  x, y, z = ft.read_columns(ifile)
  ifile.close()

  # DEBUG: Print initial unit cell to stdout
  #print "\nINITIAL UNIT CELL:"
  #for i in xrange(len(x)):
  #    print x.item(i), y.item(i), z.item(i)
  #print

  # ^^^^^^^^^^^^^^^^^^^^ STEP 2: Form distance list from combinatorics (nCr), as to correctly pair neighbors 
  # (obtained XYZ of unit cell may not have neighbors in correct order) ^^^^^^^^^^^^^^^^^^^^ #
  distlist = []
  numdists = int((factorial(x.size) / (2*factorial(x.size-2))))
  distarr = zeros(numdists)
  k = 1
  for i in range(x.size-1):
      for j in range(i+1, x.size):
	  dx = x.item(i) - x.item(j)
	  dy = y.item(i) - y.item(j)
	  dz = z.item(i) - z.item(j)
	  r2 = dx*dx + dy*dy + dz*dz
	  d = sqrt(r2)
	  d = round(d,6)				# Limit distance to 6 decimal places
	  distarr.itemset(k-1,d)
	  distlist.append([((i+1),(j+1)),d])
	  #print "%d : %d-->%d = %f" % (k, i+1, j+1, d)
	  k += 1

  distarr.sort()
  bondlen = round(distarr.item(0),6)		# Limit bondlength match to 6 decimals (not as strict as full 16 decimals matching)
  look2 = Lookup(distlist)
  pairslist = look2.get_key(bondlen)
  numatoms = len(pairslist)*2

  if numatoms != x.size:				# Error check
      print "\nERROR in determining the correct pairs from bond lengths. Exiting..."
      sys.exit(1)

  pairs = zeros(numatoms, int)
  k = 0
  for i in range(len(pairslist)):
      for j in range(2):
	  pairs.itemset(k,pairslist[i][j])
	  k += 1

  pairs.shape = (len(pairslist),2)

  # Armed with the list of atom pairs, re-arrange input file so as to list pairs on neighboring lines
  d = zeros(x.size)
  e = zeros(x.size)
  f = zeros(x.size)
  k = 0
  for i in range(len(pairs)):
      for j in range(2):
	  atnum = pairs[i][j]
	  d.itemset(k,x.item(atnum-1))
	  e.itemset(k,y.item(atnum-1))
	  f.itemset(k,z.item(atnum-1))
	  k += 1

  x = d.copy()
  y = e.copy()
  z = f.copy()

  # DEBUG
  #print "Rearranged Unit Cell:"
  #for i in range(len(x)):
  #    print "%f\t%f\t%f" % (x.item(i), y.item(i), z.item(i))
  #print

  # ^^^^^^^^^^^^^^^^^^^^ STEP 4: Generate N2G sites (by midpoint formula) ^^^^^^^^^^^^^^^^^^^^ #
  a = zeros(x.size*0.5)
  b = zeros(x.size*0.5)
  if phase == "alpha":
      c = zeros(x.size*0.5)
  elif phase == "gamma":
      c = zeros(z.size*0.5)
  j = 0

  for i in range(0, x.size, 2):
      xcom = ((x.item(i) + x.item(i+1)) * 0.5)
      ycom = ((y.item(i) + y.item(i+1)) * 0.5)
      zcom = ((z.item(i) + z.item(i+1)) * 0.5)
      a.itemset(j,round(xcom,6))
      b.itemset(j,round(ycom,6))
      c.itemset(j,round(zcom,6))
      j += 1

  # Print COM coordinates to stdout
  #for i in xrange(a.shape[0]):
  #    print a.item(i), b.item(i), z.item(i)
  
  # Print N2G + N2E coordinates (xyz) to stdout
  #print "\nDEBUG: REARRANGED UNIT CELL + COM SITES:"
  #k = 0
  #for i in range(a.size):
  #    print a.item(i), b.item(i), c.item(i)
  #    for j in range(2):
  #	print x.item(k), y.item(k), z.item(k)
  #	k += 1
  #print

  # ^^^^^^^^^^^^^^^^^^^^ STEP 5: Generate N2E sites (crystallographic data yields bond lengths less than 0.549 A) ^^^^^^^^^^^^^^^^^^^^ #
  # XXX Should have a check here that if the provided unit cell has the correct bond length, don't need to execute this step.
  breal = 1.098
  bhalf = breal * 0.5
  t = (bhalf / (bondlen*0.5))
  print "Bond length (via input)  = %f Angstroms" % bondlen
  print "Bond length (experiment) = %f Angstroms" % breal
  print "Optimal t for G-->E: %f" % t
  print

  N2E = zeros(x.size*3)
  j = 0
  k = 0
  for i in range(a.size):
      dx = a.item(i) - x.item(j)
      dy = b.item(i) - y.item(j)
      dz = c.item(i) - z.item(j)
      xnew = a.item(i) - t*dx
      ynew = b.item(i) - t*dy
      znew = c.item(i) - t*dz

      N2E.itemset(k,round(xnew,6))
      k += 1
      N2E.itemset(k,round(ynew,6))
      k += 1
      N2E.itemset(k,round(znew,6))
      k += 1

      xnew = a.item(i) + t*dx
      ynew = b.item(i) + t*dy
      znew = c.item(i) + t*dz

      N2E.itemset(k,round(xnew,6))
      k += 1
      N2E.itemset(k,round(ynew,6))
      k += 1
      N2E.itemset(k,round(znew,6))
      k += 1
      j += 2

  N2E.shape = (x.size,3)

  totalrows = x.size + a.size
  totalelems = totalrows * 3
  three_site = zeros(totalelems)

  l = 0
  m = 0
  for i in range(a.size):
      three_site.itemset(m,a.item(i))
      m += 1
      three_site.itemset(m,b.item(i))
      m += 1
      three_site.itemset(m,c.item(i))
      m += 1
      for j in range(2):
          for k in range(3):
	      three_site.itemset(m,N2E.item(l,k))
	      m += 1
	  l += 1

  three_site.shape = (totalrows,3)
  # DEBUG
  #print
  #print "3-Site (Scaled Nitrogens for appropriate bond length):"
  #for i in range(len(three_site)):
  #    print "%d: %f\t%f\t%f" % (i+1, three_site.item(i,0), three_site.item(i,1), three_site.item(i,2))
  #print
  
  # Job Preparation (don't need to loop this)
  parentdir = os.getcwd()
  trajfile = open(parentdir + '/traj.pqr', 'w')

  input_file    = 'n2.input'
  job_name      = 'n2'
  pqr_input     = 'n2.initial.pqr'      # This should ALWAYS be set. There are dependencies on this keyword within this script!
  pqr_output    = 'n2.final.pqr'
  pqr_restart   = 'n2.restart.pqr'
  energy_output = 'n2.energy.dat'
  traj_output   = 'n2.traj.pqr'
  dipole_output = 'dipole.dat'
  field_output  = 'field.dat'

  sitelist      = [    'N2G',    'N2E',    'N2E',    'N2N',    'N2N' ]
  masslist      = [  0.00000, 14.00670, 14.00670,  0.00000,  0.00000 ]
  qlist         = [  1.04742, -0.52371, -0.52371,  0.00000,  0.00000 ]
  pollist       = [      0.0,      0.0,      0.0,      0.0,      0.0 ]

  # Set LB LJ parameters (default)
  epslist       = [ 17.60293,  0.00000,  0.00000, 18.12772, 18.12772 ]
  siglist       = [  3.44522,  0.00000,  0.00000,  3.15125,  3.15125 ]

  if pol_on: # Replace LB alpha and LJ parameters with polar ones
      pollist   = [  1.49645,  0.45510,  0.45510,      0.0,      0.0 ]
      epslist   = [ 20.63650,  0.00000,  0.00000, 16.14200, 16.14200 ]
      siglist   = [  3.42344,  0.00000,  0.00000,  3.16141,  3.16141 ]

  # ...but if WH mixing rules are requested, then use these LJ parameters
  if mixrule == "WH":
      epslist   = [ 32.01121,  0.00000,  0.00000, 13.09958, 13.09958 ]
      siglist   = [  3.39190,  0.00000,  0.00000,  3.09821,  3.09821 ]

  # ...if both WH and polarization are requested, use these parameters
  if mixrule == "WH" and pol_on:
      print "Both WH and pol_on are set..."
      pollist   = [  1.49645,  0.45510,  0.45510,      0.0,      0.0 ]
      epslist   = [ 32.10001,  0.00000,  0.00000, 12.72315, 12.72315 ]
      siglist   = [  3.39450,  0.00000,  0.00000,  3.10175,  3.10175 ]

  # LB Parameters for the S configuration
  if S == "+S1" and mixrule == "LB":
      print "Using S1 / NP / LB parameters..."
      epslist   = [ 25.03351,  0.00000,  0.00000, 15.93472, 15.93472 ]
      siglist   = [  3.45163,  0.00000,  0.00000,  3.06412,  3.06412 ]
  if S == "+S1" and mixrule == "LB" and pol_on:
      print "Using S1 / Pol / LB parameters..."
      epslist   = [ 27.02224,  0.00000,  0.00000, 14.55166, 14.55166 ]
      siglist   = [  3.43565,  0.00000,  0.00000,  3.08409,  3.08409 ]

  if S == "+S2" and mixrule == "LB":
      print "Using S2 / NP / LB parameters..."
      epslist       = [ 26.01446,  0.00000,  0.00000, 15.44499, 15.44499 ]
      siglist       = [  3.44119,  0.00000,  0.00000,  3.07386,  3.07386 ]
  if S == "+S2" and mixrule == "LB" and pol_on:
      print "Using S2 / Pol / LB parameters..."
      epslist       = [ 28.19592,  0.00000,  0.00000, 14.07209, 14.07209 ]
      siglist       = [  3.42602,  0.00000,  0.00000,  3.09286,  3.09286 ]

  if S == "+S3" and mixrule == "LB":
      print "Using S3 / NP / LB parameters..."
      epslist       = [ 26.67119,  0.00000,  0.00000, 15.16878, 15.16878 ]
      siglist       = [  3.43807,  0.00000,  0.00000,  3.07685,  3.07685 ]
  if S == "+S3" and mixrule == "LB" and pol_on:
      print "Using S3 / Pol / LB parameters..."
      epslist       = [ 27.26533,  0.00000,  0.00000, 14.47356, 14.47356 ]
      siglist       = [  3.43441,  0.00000,  0.00000,  3.08518,  3.08518 ]

  if S == "+S4" and mixrule == "LB":
      print "Using S4 / NP / LB parameters..."
      epslist       = [ 27.00444,  0.00000,  0.00000, 14.94928, 14.94928 ]
      siglist       = [  3.43335,  0.00000,  0.00000,  3.08084,  3.08084 ]
  if S == "+S4" and mixrule == "LB" and pol_on:
      print "Using S4 / Pol / LB parameters..."
      epslist       = [ 25.84594,  0.00000,  0.00000, 15.07049, 15.07049 ]
      siglist       = [  3.44659,  0.00000,  0.00000,  3.08442,  3.08442 ]

  if S == "+S5" and mixrule == "LB":
      print "Using S5 / NP / LB parameters..."
      epslist       = [ 25.64425,  0.00000,  0.00000, 15.53200, 15.53200 ]
      siglist       = [  3.44416,  0.00000,  0.00000,  3.07293,  3.07293 ]
  if S == "+S5" and mixrule == "LB" and pol_on:
      print "Using S5 / Pol / LB parameters..."
      epslist       = [ 24.09181,  0.00000,  0.00000, 16.06115, 16.06115 ]
      siglist       = [  3.46325,  0.00000,  0.00000,  3.06249,  3.06249 ]

  if S == "+S6" and mixrule == "LB":
      print "Using S6 / NP / LB parameters..."
      epslist       = [ 24.42768,  0.00000,  0.00000, 16.19522, 16.19522 ]
      siglist       = [  3.45331,  0.00000,  0.00000,  3.06721,  3.06721 ]
  if S == "+S6" and mixrule == "LB" and pol_on:
      print "Using S6 / Pol / LB parameters..."
      epslist       = [ 28.10306,  0.00000,  0.00000, 14.07191, 14.07191 ]
      siglist       = [  3.42862,  0.00000,  0.00000,  3.09229,  3.09229 ]

  if S == "+S7" and mixrule == "LB":
      print "Using S7 / NP / LB parameters..."
      epslist       = [ 25.92003,  0.00000,  0.00000, 15.45715, 15.45715 ]
      siglist       = [  3.44168,  0.00000,  0.00000,  3.07802,  3.07802 ]
  if S == "+S7" and mixrule == "LB" and pol_on:
      print "Using S7 / Pol / LB parameters..."
      epslist       = [ 27.34094,  0.00000,  0.00000, 14.44778, 14.44778 ]
      siglist       = [  3.43299,  0.00000,  0.00000,  3.08472,  3.08472 ]

  if S == "+S8" and mixrule == "LB":
      print "Using S8 / NP / LB parameters..."
      epslist       = [ 26.17654,  0.00000,  0.00000, 15.31917, 15.31917 ]
      siglist       = [  3.44049,  0.00000,  0.00000,  3.07569,  3.07569 ]
  if S == "+S8" and mixrule == "LB" and pol_on:
      print "Using S8 / Pol / LB parameters..."
      epslist       = [ 27.26834,  0.00000,  0.00000, 14.49144, 14.49144 ]
      siglist       = [  3.43809,  0.00000,  0.00000,  3.08059,  3.08059 ]

  if S == "+S9" and mixrule == "LB":
      print "Using S9 / NP / LB parameters..."
      epslist       = [ 24.75185,  0.00000,  0.00000, 15.98923, 15.98923 ]
      siglist       = [  3.45095,  0.00000,  0.00000,  3.06702,  3.06702 ]
  if S == "+S9" and mixrule == "LB" and pol_on:
      print "Using S9 / Pol / LB parameters..."
      epslist       = [ 26.96095,  0.00000,  0.00000, 14.74373, 14.74373 ]
      siglist       = [  3.43808,  0.00000,  0.00000,  3.08184,  3.08184 ]

  if S == "+S10" and mixrule == "LB":
      print "Using S10 / NP / LB parameters..."
      epslist       = [ 26.65155,  0.00000,  0.00000, 15.07068, 15.07068 ]
      siglist       = [  3.43539,  0.00000,  0.00000,  3.07904,  3.07904 ]
  if S == "+S10" and mixrule == "LB" and pol_on:
      print "Using S10 / Pol / LB parameters..."
      epslist       = [ 27.05609,  0.00000,  0.00000, 14.68874, 14.68874 ]
      siglist       = [  3.44079,  0.00000,  0.00000,  3.07444,  3.07444 ]

  # WH Parameters for the S configuration
  if S == "+S1" and mixrule == "WH":
      print "Using S1 / NP / WH parameters..."
      epslist       = [ 27.61722,  0.00000,  0.00000, 15.19242, 15.19242 ]
      siglist       = [  3.41827,  0.00000,  0.00000,  3.07809,  3.07809 ]
  if S == "+S1" and mixrule == "WH" and pol_on:
      print "Using S1 / Pol / WH parameters..."
      epslist       = [ 31.03292,  0.00000,  0.00000, 13.47803, 13.47803 ]
      siglist       = [  3.40469,  0.00000,  0.00000,  3.08985,  3.08985 ]

  if S == "+S2" and mixrule == "WH":
      print "Using S2 / NP / WH parameters..."
      epslist       = [ 28.10847,  0.00000,  0.00000, 15.10849, 15.10849 ]
      siglist       = [  3.41665,  0.00000,  0.00000,  3.07546,  3.07546 ]
  if S == "+S2" and mixrule == "WH" and pol_on:
      print "Using S2 / Pol / WH parameters..."
      epslist       = [ 29.21930,  0.00000,  0.00000, 14.16880, 14.16880 ]
      siglist       = [  3.41148,  0.00000,  0.00000,  3.08606,  3.08606 ]

  if S == "+S3" and mixrule == "WH":
      print "Using S3 / NP / WH parameters..."
      epslist       = [ 27.26005,  0.00000,  0.00000, 15.30942, 15.30942 ]
      siglist       = [  3.42032,  0.00000,  0.00000,  3.07774,  3.07774 ]
  if S == "+S3" and mixrule == "WH" and pol_on:
      print "Using S3 / Pol / WH parameters..."
      epslist       = [ 28.19577,  0.00000,  0.00000, 14.51448, 14.51448 ]
      siglist       = [  3.41610,  0.00000,  0.00000,  3.08761,  3.08761 ]

  if S == "+S4" and mixrule == "WH":
      print "Using S4 / NP / WH parameters..."
      epslist       = [ 29.00184,  0.00000,  0.00000, 14.66607, 14.66607 ]
      siglist       = [  3.41348,  0.00000,  0.00000,  3.07822,  3.07822 ]
  if S == "+S4" and mixrule == "WH" and pol_on:
      print "Using S4 / Pol / WH parameters..."
      epslist       = [ 29.53328,  0.00000,  0.00000, 14.14841, 14.14841 ]
      siglist       = [  3.41320,  0.00000,  0.00000,  3.08045,  3.08045 ]

  if S == "+S5" and mixrule == "WH":
      print "Using S5 / NP / WH parameters..."
      epslist       = [ 27.28212,  0.00000,  0.00000, 15.35816, 15.35816 ]
      siglist       = [  3.42104,  0.00000,  0.00000,  3.07580,  3.07580 ]
  if S == "+S5" and mixrule == "WH" and pol_on:
      print "Using S5 / Pol / WH parameters..."
      epslist       = [ 29.74861,  0.00000,  0.00000, 13.93183, 13.93183 ]
      siglist       = [  3.41115,  0.00000,  0.00000,  3.08681,  3.08681 ]

  if S == "+S6" and mixrule == "WH":
      print "Using S6 / NP / WH parameters..."
      epslist       = [ 27.66681,  0.00000,  0.00000, 15.24137, 15.24137 ]
      siglist       = [  3.41800,  0.00000,  0.00000,  3.07565,  3.07565 ]
  if S == "+S6" and mixrule == "WH" and pol_on:
      print "Using S6 / Pol / WH parameters..."
      epslist       = [ 29.83458,  0.00000,  0.00000, 13.73245, 13.73245 ]
      siglist       = [  3.40462,  0.00000,  0.00000,  3.09956,  3.09956 ]

  if S == "+S7" and mixrule == "WH":
      print "Using S7 / NP / WH parameters..."
      epslist       = [ 27.59946,  0.00000,  0.00000, 15.24748, 15.24748 ]
      siglist       = [  3.41814,  0.00000,  0.00000,  3.07852,  3.07852 ]
  if S == "+S7" and mixrule == "WH" and pol_on:
      print "Using S7 / Pol / WH parameters..."
      epslist       = [ 29.66824,  0.00000,  0.00000, 13.95069, 13.95069 ]
      siglist       = [  3.40969,  0.00000,  0.00000,  3.09263,  3.09263 ]

  if S == "+S8" and mixrule == "WH":
      print "Using S8 / NP / WH parameters..."
      epslist       = [ 27.63584,  0.00000,  0.00000, 15.19941, 15.19941 ]
      siglist       = [  3.41651,  0.00000,  0.00000,  3.07763,  3.07763 ]
  if S == "+S8" and mixrule == "WH" and pol_on:
      print "Using S8 / Pol / WH parameters..."
      epslist       = [ 30.22857,  0.00000,  0.00000, 13.60798, 13.60798 ]
      siglist       = [  3.40488,  0.00000,  0.00000,  3.09790,  3.09790 ]

  if S == "+S9" and mixrule == "WH":
      print "Using S9 / NP / WH parameters..."
      epslist       = [ 29.28207,  0.00000,  0.00000, 14.54166, 14.54166 ]
      siglist       = [  3.41080,  0.00000,  0.00000,  3.08026,  3.08026 ]
  if S == "+S9" and mixrule == "WH" and pol_on:
      print "Using S9 / Pol / WH parameters..."
      epslist       = [ 28.76589,  0.00000,  0.00000, 14.14972, 14.14972 ]
      siglist       = [  3.41006,  0.00000,  0.00000,  3.09524,  3.09524 ]

  if S == "+S10" and mixrule == "WH":
      print "Using S10 / NP / WH parameters..."
      epslist       = [ 29.14182,  0.00000,  0.00000, 14.51428, 14.51428 ]
      siglist       = [  3.41029,  0.00000,  0.00000,  3.08428,  3.08428 ]
  if S == "+S10" and mixrule == "WH" and pol_on:
      print "Using S10 / Pol / WH parameters..."
      epslist       = [ 29.13527,  0.00000,  0.00000, 14.15119, 14.15119 ]
      siglist       = [  3.41287,  0.00000,  0.00000,  3.08987,  3.08987 ]

  site          = []
  molnum        = []
  mass          = []
  charge        = []
  polar         = []
  eps           = []
  sig           = []
  total_e       = []
  total_e_norm  = []
  coul_e        = []
  coul_e_norm   = []
  rd_e          = []
  rd_e_norm     = []
  pol_e         = []
  pol_e_norm    = []
  abclist       = []
  rmsd_list     = []
  crds_list     = []
  cell_volume   = []
  
  # Main Loop over all a's
  cc = 0		# Set coordinates counter (for RMSD)

  # Time MPMC
  times = zeros(0)      # Declare numpy array which will store each local time value for averaging
  start_global = time.time()

  # Call the appropriate loop
  if phase == "alpha":
      al = alpha_loop(mixrule, pol_on, n_cells, abc_min, abc_max, abc_del, x, breal, parentdir, trajfile, input_file, job_name, pqr_input, pqr_output, pqr_restart, energy_output, traj_output, dipole_output, field_output, sitelist, masslist, qlist, pollist, epslist, siglist, site, totalelems, totalrows)
      print al
  elif phase == "gamma":
      gl = gamma_loop(mixrule, pol_on, n_cells, abc_min, abc_max, abc_del, x, breal, parentdir, trajfile, input_file, job_name, pqr_input, pqr_output, pqr_restart, energy_output, traj_output, dipole_output, field_output, sitelist, masslist, qlist, pollist, epslist, siglist, site, totalelems, totalrows)
      print gl
  else:
      print "Unsure as to how you managed to get to this point without prior erroring, but whatev...I have no loop for you, as your phase is not 'alpha' or 'gamma', so try again!"
      sys.exit(1)


      # Execute MPMC
      call("mpmc " + input_file + " >& runlog.log", shell=True)

      # Append final coordinates to trajectory file
      infinal = open(pqr_output, 'r')
      finalpqr = infinal.read()
      infinal.close()
      stop = finalpqr.find('CONECT')
      trajpqr = finalpqr[:stop]
      trajfile.write('%s' % trajpqr)
      trajfile.write('ENDMDL\n')

      # Process Energies
      norm_fac = n_cells * 2
      enfile = open(energy_output, 'r')
      enfile.readline()
      step, tote, coul, rd, pol, vdw, kin, tempe, partnum, sr, vol = ft.read_columns(enfile)
      enfile.close()
      total_e.append(tote.item(0))
      total_e_norm.append(tote.item(0)/norm_fac)
      coul_e.append(coul.item(0))
      coul_e_norm.append(coul.item(0)/norm_fac)
      rd_e.append(rd.item(0))
      rd_e_norm.append(rd.item(0)/norm_fac)
      pol_e.append(pol.item(0))
      pol_e_norm.append(pol.item(0)/norm_fac)
      cell_volume.append(vol.item(0))
     
      os.chdir(parentdir)

      # Print timing statistics
      end = time.time()
      elapsed_global = end - start_global
      elapsed_local  = end - start_local
      times = append(times,elapsed_local)               # Running list of local times
      avg = float(average(times))                       # Average local time (time/MPMC job)
      pred_time = avg * ncycles                         # Predicted total time
      eta = pred_time - elapsed_global                  # Predicted ETA to job completion
      print "Global time = %08.2f s / %07.2f m / %5.2f h" % (elapsed_global, elapsed_global/60.0, elapsed_global/3600.0)
      print "Local  time = %08.2f s / %07.2f m / %5.2f h" % (elapsed_local, elapsed_local/60.0, elapsed_local/3600.0)
      print "    ETA     = %08.2f s / %07.2f m / %5.2f h" % (eta, eta/60.0, eta/3600.0)

     #################################################	#@CRC
     # UPDATE c_current!!!				#@CRC
      c_current = (a_current + abc_del) / a2c_ratio	#@CRC
     #################################################	#@CRC

  # Print Energy & Stats
  print
  esummary      = open('energy_total_summary.dat', 'w')
  esummary_norm = open('energy_total_summary_NORM.dat', 'w')
  csummary      = open('energy_coul_summary.dat', 'w')
  csummary_norm = open('energy_coul_summary_NORM.dat', 'w')
  rsummary      = open('energy_rd_summary.dat', 'w')
  rsummary_norm = open('energy_rd_summary_NORM.dat', 'w')
  psummary      = open('energy_polar_summary.dat', 'w')
  psummary_norm = open('energy_polar_summary_NORM.dat', 'w')
  rmsdsum       = open('rmsd_summary.dat', 'w')
  volsum        = open('volumes.dat', 'w')

  for i in range(len(abclist)):
      print abclist[i], total_e[i]
      esummary.write("%f\t%f\n" % (abclist[i], total_e[i]))
      esummary_norm.write("%f\t%f\n" % (abclist[i], total_e_norm[i]))
      csummary.write("%f\t%f\n" % (abclist[i], coul_e[i]))
      csummary_norm.write("%f\t%f\n" % (abclist[i], coul_e_norm[i]))
      rsummary.write("%f\t%f\n" % (abclist[i], rd_e[i]))
      rsummary_norm.write("%f\t%f\n" % (abclist[i], rd_e_norm[i]))
      psummary.write("%f\t%f\n" % (abclist[i], pol_e[i]))
      psummary_norm.write("%f\t%f\n" % (abclist[i], pol_e_norm[i]))
      volsum.write(  "%f\t%f\n" % (abclist[i], cell_volume[i]))

  trajfile.close()
  esummary.close()
  csummary.close()
  rsummary.close()
  psummary.close()
  esummary_norm.close()
  csummary_norm.close()
  rsummary_norm.close()
  psummary_norm.close()

  for i in range(len(rmsd_list)):
      rmsdsum.write("%f\n" % rmsd_list[i])

  rmsdsum.close()
