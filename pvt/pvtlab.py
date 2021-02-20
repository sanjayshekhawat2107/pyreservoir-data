"""
Programs for Processing PVT Lab Reports

@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

def linear_interpolate(p, p1, prop):
  """
  Linear interpolation to assign properties from PVT data to respective
  reservoir pressures in production data

  Input:

  p = pressure in PVT data
  p1 = pressure in production data
  prop = property in PVT data to be interpolated to production data 
         (has same array length as "p")

  Output:

  prop_interpolated = the interpolated property values for the production data
  """
  import numpy as np

  prop_interpolated = []
  for i in range(len(p1)):
    for j,k in zip(range(1,len(p)), range(len(p)-1)):
      if p1[i] < p[j-1] and p1[i] > p[k+1]:
        # interpolating if value in prod data is between two values in the PVT data
        prop_plus, p_plus = prop[j-1], p[j-1] 
        prop_min, p_min = prop[k+1], p[k+1]
        prop_int = ((p1[i] - p_min) * (prop_plus - prop_min) / (p_plus - p_min)) + prop_min
        prop_interpolated.append(prop_int)
      if p1[i] == p[j-1]:
        # using the value in the PVT data if equals to value in the prod data
        prop_int = prop[j-1]
        prop_interpolated.append(prop_int)
      if p1[i] == p[k+1]:
        # using the value in the PVT data if equals to value in the prod data
        prop_int = prop[k+1]
        prop_interpolated.append(prop_int)
  prop_interpolated = np.array(prop_interpolated)
  return prop_interpolated

def cvd_condensate(z, z2, temp, p, Gp, Np, Vo):
    """
    Calculate volatile oil-gas ratio of condensate from Constant-Volume Depletion (CVD) Study
    Walsh and Towler (1995)

    Inputs
    z: measured gas-phase compressibility factor (array)
    z2: measured two-phase compressibility factor (array)
    p: measured pressure (array)
    Gp: gas produced in the PVT cell
    Np: condensate produced in the PVT cell
    Vo: condensate volume in the PVT cell
    """

    z_j = z; Gp_j = Gp; Np_j = Np; z2_j = z2; Vo_j = Vo

    # calculate gas FVF (Bg)
    Bg = (0.00503676 * z_j * (temp + 460)) / p  # in RB/scf

    # initial gas FVF
    Bgi = Bg[0]

    # initial Gfg
    Gfgi = Gp_j[-1]  # in scf

    # initial Nfo
    Nfoi = Np_j[-1]

    # calculate initial Vtg (Vtg1)
    Vtg1 = Gfgi * Bgi  # in res bbl

    # initial values for Eq 10.14
    ntj_nt1 = 1
    delta_ngj_to_nt1 = 0

    # empty arrays for appending
    Vtoj_arr = []
    ntj_nt1_arr = []
    Vtoj_Vtgj_arr = []
    Vtgj_arr = []
    delta_Vtgj_arr = []
    ngj_nt1_arr = []
    delta_ngj_to_ngj_arr = []
    delta_ngj_to_nt1_arr = []
    delta_Gpj_arr = []
    delta_Npj_arr = []
    Gfgj_arr = []
    Nfgj_arr = []
    Gj_arr = []
    Nj_arr = []
    Gfoj_arr = []
    Nfoj_arr = []
    Boj_arr = []
    Bgj_arr = []
    Rsj_arr = []
    Rvj_arr = []

    for i in range(len(j) - 1):

        # Eq 10.13
        Vtoj = Vo_j[i] * Vtg1
        Vtoj_arr.append(Vtoj)

        # Eq 10.14
        ntj_nt1 = ntj_nt1 - delta_ngj_to_nt1
        ntj_nt1_arr.append(ntj_nt1)

        # Eq 10.15
        Vtoj_Vtgj = ((Vtg1 * z2_j[i] * p[0]) / (z2_j[0] * p[i])) * (ntj_nt1)
        Vtoj_Vtgj_arr.append(Vtoj_Vtgj)

        # Eq 10.16
        Vtgj = Vtoj_Vtgj - Vtoj
        Vtgj_arr.append(Vtgj)

        # Eq 10.17
        delta_Vtgj = Vtoj_Vtgj - Vtg1
        delta_Vtgj_arr.append(delta_Vtgj)

        # Eq 10.18
        ngj_nt1 = (Vtgj * z_j[0] * p[i]) / (z_j[i] * Vtg1 * p[0])
        ngj_nt1_arr.append(ngj_nt1)

        # Eq 10.19
        delta_ngj_to_ngj = delta_Vtgj / Vtgj
        delta_ngj_to_ngj_arr.append(delta_ngj_to_ngj)

        # Eq 10.20
        delta_ngj_to_nt1 = delta_ngj_to_ngj * ngj_nt1
        delta_ngj_to_nt1_arr.append(delta_ngj_to_nt1)

        if i == 0:
            # Eq 10.21
            delta_Gpj = Gp_j[i] - 0

            # Eq 10.22
            delta_Npj = Np_j[i] - 0

            # Eq 10.23
            Gj = Gfgi - delta_Gpj

            # Eq 10.24
            Nj = Nfoi - delta_Gpj

            # Eq 10.25
            Gfgj = Gfgi

            # Eq 10.26
            Nfgj = Nfoi

        if i > 0:
            # Eq 10.21
            delta_Gpj = Gp_j[i] - Gp_j[i - 1]

            # Eq 10.22
            delta_Npj = Np_j[i] - Np_j[i - 1]

            # Eq 10.23
            Gj = Gj - delta_Gpj_arr[-1]

            # Eq 10.24
            Nj = Nj - delta_Npj_arr[-1]

            # Eq 10.25
            Gfgj = (Vtgj * delta_Gpj) / delta_Vtgj

            # Eq 10.26
            Nfgj = (Vtgj * delta_Npj) / delta_Vtgj

        delta_Gpj_arr.append(delta_Gpj)
        delta_Npj_arr.append(delta_Npj)
        Gj_arr.append(Gj)
        Nj_arr.append(Nj)
        Gfgj_arr.append(Gfgj)
        Nfgj_arr.append(Nfgj)

        # Eq 10.27
        Gfoj = Gj - Gfgj
        Gfoj_arr.append(Gfoj)

        # Eq 10.28
        Nfoj = Nj - Nfgj
        Nfoj_arr.append(Nfoj)

        # Eq 10.29
        Boj = Vtoj / Nfoj
        Boj_arr.append(Boj)

        # Eq 10.30
        Bgj = Vtgj / Gfgj
        Bgj_arr.append(Bgj)

        # Eq 10.31
        Rsj = Gfoj / Nfoj
        Rsj_arr.append(Rsj)

        # Eq 10.32
        Rvj = (Nfgj / Gfgj) * 1E+06  # result in STB/scf
        Rvj_arr.append(Rvj)

    Rv = Rvj_arr
    return(Rv)
