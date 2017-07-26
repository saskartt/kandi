import sys
import numpy as np
import netCDF4 as nc

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createTableRow(qty,label):
  qty = qty.astype(str)
  qty = np.insert(qty, 0, label)
  return qty

def openDataSet(filename):
  try:
    ds = nc.Dataset(filename)
  except RuntimeError:
    raise IOError("Input file {} not found!".format(filename))
  return ds

def averageProfilesWS(statdomain, t_inds, h_inds, ds):
  # Mean values
  u = np.mean(ds.variables['u_0'+statdomain][t_inds,:],axis=0)
  # Interpolate values
  u = np.interp(h_inds, ds.variables['zu_0'+statdomain][:], u)
  str_u = createTableRow(u,'u')

  v = np.mean(ds.variables['v_0'+statdomain][t_inds,:],axis=0)
  v = np.interp(h_inds, ds.variables['zv_0'+statdomain][:], v)
  str_v = createTableRow(v,'v')

  w = np.mean(ds.variables['w_0'+statdomain][t_inds,:],axis=0)
  w = np.interp(h_inds, ds.variables['zw_0'+statdomain][:], w)
  str_w = createTableRow(w,'w')

  U = np.linalg.norm([u,v,w],axis=0)
  str_U = createTableRow(U,'U')
  return [str_u,str_v,str_w,str_U]

def averageProfilesVariances(statdomain, t_inds, h_inds, ds):
  # Variances
  var_u = np.mean(ds.variables['u*2_0'+statdomain][t_inds,:],axis=0)
  var_u = np.interp(h_inds, ds.variables['zu*2_0'+statdomain][:], var_u)
  var_u = createTableRow(var_u,'u*2')

  var_v = np.mean(ds.variables['v*2_0'+statdomain][t_inds,:],axis=0)
  var_v = np.interp(h_inds, ds.variables['zv*2_0'+statdomain][:], var_v)
  var_v = createTableRow(var_v,'v*2')

  var_w = np.mean(ds.variables['w*2_0'+statdomain][t_inds,:],axis=0)
  var_w = np.interp(h_inds, ds.variables['zw*2_0'+statdomain][:], var_w)
  var_w = createTableRow(var_w,'w*2')

  return [var_u,var_v,var_w]

def averageProfilesMomentumFluxes(statdomain,t_inds, h_inds, ds):
  flx_u = np.mean(ds.variables['wu_0'+statdomain][t_inds,:],axis=0)
  flx_u = np.interp(h_inds, ds.variables['zwu_0'+statdomain][:], flx_u)
  flx_u_str = createTableRow(flx_u, "wu")

  flx_v = np.mean(ds.variables['wv_0'+statdomain][t_inds,:],axis=0)
  flx_v = np.interp(h_inds, ds.variables['zwv_0'+statdomain][:], flx_v)
  flx_v_str = createTableRow(flx_v, "wv")

  flx_uv = flx_u + flx_v
  flx_uv = createTableRow(flx_uv, "u'w'+ v'w'")

  return [flx_u_str, flx_v_str, flx_uv]

def averageProfilesTKE(statdomain, t_inds, h_inds, ds):
  # Resolved perturbation energy and subgrid-scale TKE
  tke_e = np.mean(ds.variables['e*_0'+statdomain][t_inds,:],axis=0)
  tke_e = np.interp(h_inds, ds.variables['ze*_0'+statdomain][:], tke_e)
  tke_e_sg = np.mean(ds.variables['e_0'+statdomain][t_inds,:],axis=0)
  tke_e_sg = np.interp(h_inds, ds.variables['ze_0'+statdomain][:], tke_e_sg)
  tke_e = createTableRow(tke_e+tke_e_sg, "TKE")

  return [tke_e]

def compareProfilesWS(statdomain, t_inds, ct_inds, h_inds, ds, cds):
  u = np.mean(ds.variables['u_0'+statdomain][t_inds,:],axis=0)
  u = np.interp(h_inds, ds.variables['zu_0'+statdomain][:], u)
  cu = np.mean(cds.variables['u_0'+statdomain][ct_inds,:],axis=0)
  cu = np.interp(h_inds, ds.variables['zu_0'+statdomain][:], cu)
  cmp_u = createTableRow(np.divide(u,cu)*100.,'u (%)')

  v = np.mean(ds.variables['v_0'+statdomain][t_inds,:],axis=0)
  v = np.interp(h_inds, ds.variables['zv_0'+statdomain][:], v)
  cv = np.mean(cds.variables['v_0'+statdomain][ct_inds,:],axis=0)
  cv = np.interp(h_inds, ds.variables['zv_0'+statdomain][:], cv)
  cmp_v = createTableRow(np.divide(v,cv)*100.,'v (%)')

  w = np.mean(ds.variables['w_0'+statdomain][t_inds,:],axis=0)
  w = np.interp(h_inds, ds.variables['zw_0'+statdomain][:], w)
  cw = np.mean(cds.variables['w_0'+statdomain][ct_inds,:],axis=0)
  cw = np.interp(h_inds, ds.variables['zv_0'+statdomain][:], cw)
  cmp_w = createTableRow(np.divide(w,cw)*100.,'w (%)')

  U = np.linalg.norm([u,v,w],axis=0)
  cU = np.linalg.norm([cu,cv,cw],axis=0)
  cmp_U = createTableRow(np.divide(U,cU)*100.,'U (%)')

  return [cmp_u, cmp_v, cmp_w, cmp_U]

def compareProfilesVariances(statdomain, t_inds, ct_inds, h_inds, ds, cds):
  var_u = np.mean(ds.variables['u*2_0'+statdomain][t_inds,:],axis=0)
  var_u = np.interp(h_inds, ds.variables['zu*2_0'+statdomain][:], var_u)
  var_cu = np.mean(cds.variables['u*2_0'+statdomain][ct_inds,:],axis=0)
  var_cu = np.interp(h_inds, ds.variables['zu*2_0'+statdomain][:], var_cu)
  cmp_var_u = createTableRow(np.divide(var_u,var_cu)*100.,'u*2 (%)')

  var_v = np.mean(ds.variables['v*2_0'+statdomain][t_inds,:],axis=0)
  var_v = np.interp(h_inds, ds.variables['zv*2_0'+statdomain][:], var_v)
  var_cv = np.mean(cds.variables['v*2_0'+statdomain][ct_inds,:],axis=0)
  var_cv = np.interp(h_inds, ds.variables['zv*2_0'+statdomain][:], var_cv)
  cmp_var_v = createTableRow(np.divide(var_v,var_cv)*100.,'v*2 (%)')

  var_w = np.mean(ds.variables['w*2_0'+statdomain][t_inds,:],axis=0)
  var_w = np.interp(h_inds, ds.variables['zw*2_0'+statdomain][:], var_w)
  var_cw = np.mean(cds.variables['w*2_0'+statdomain][ct_inds,:],axis=0)
  var_cw = np.interp(h_inds, ds.variables['zw*2_0'+statdomain][:], var_cw)
  cmp_var_w = createTableRow(np.divide(var_w,var_cw)*100.,'w*2 (%)')

  return [cmp_var_u, cmp_var_v, cmp_var_w]

def compareProfilesMomentumFluxes(statdomain, t_inds, ct_inds, h_inds, ds, cds):
  flx_u = np.mean(ds.variables['wu_0'+statdomain][t_inds,:],axis=0)
  flx_u = np.interp(h_inds, ds.variables['zwu_0'+statdomain][:], flx_u)
  flx_cu = np.mean(cds.variables['wu_0'+statdomain][ct_inds,:],axis=0)
  flx_cu = np.interp(h_inds, ds.variables['zwu_0'+statdomain][:], flx_cu)
  cmp_flx_u = createTableRow(np.divide(flx_u,flx_cu)*100.,'wu (%)')

  flx_v = np.mean(ds.variables['wv_0'+statdomain][t_inds,:],axis=0)
  flx_v = np.interp(h_inds, ds.variables['zwv_0'+statdomain][:], flx_v)
  flx_cv = np.mean(cds.variables['wv_0'+statdomain][ct_inds,:],axis=0)
  flx_cv = np.interp(h_inds, ds.variables['zwv_0'+statdomain][:], flx_cv)
  cmp_flx_v = createTableRow(np.divide(flx_v,flx_cv)*100.,'wv (%)')

  flx = flx_u + flx_v
  flx_c = flx_cu + flx_cv
  cmp_flx = createTableRow(np.divide(flx,flx_c)*100., "u'w'+ v'w'")

  return [cmp_flx_u, cmp_flx_v, cmp_flx]

def compareProfilesTKE(statdomain, t_inds, ct_inds, h_inds, ds, cds):
  # Resolved perturbation energy and subgrid-scale TKE
  tke_e = np.mean(ds.variables['e*_0'+statdomain][t_inds,:],axis=0)
  tke_e = np.interp(h_inds, ds.variables['ze*_0'+statdomain][:], tke_e)
  tke_e_sg = np.mean(ds.variables['e_0'+statdomain][t_inds,:],axis=0)
  tke_e_sg = np.interp(h_inds, ds.variables['ze_0'+statdomain][:], tke_e_sg)
  tke_e = tke_e + tke_e_sg

  tke_ce = np.mean(cds.variables['e*_0'+statdomain][ct_inds,:],axis=0)
  tke_ce = np.interp(h_inds, cds.variables['ze*_0'+statdomain][:], tke_ce)
  tke_ce_sg = np.mean(cds.variables['e_0'+statdomain][ct_inds,:],axis=0)
  tke_ce_sg = np.interp(h_inds, cds.variables['ze_0'+statdomain][:], tke_ce_sg)
  tke_ce = tke_ce + tke_ce_sg

  cmp_tke_e = createTableRow(np.divide(tke_e,tke_ce)*100.,'TKE (%)')

  return [cmp_tke_e]

def compileDataListAverages(domain, t_inds, h_inds, ds, varlist):
  sep = [[''] * (len(h_inds) + 1)]  # empty line string
  dlist = []

  if ('ws' in varlist):
    pr_01 = averageProfilesWS(domain, t_inds, h_inds, ds)
    dlist = dlist + pr_01
  if ('var' in varlist):
    pr_02 = averageProfilesVariances(domain, t_inds, h_inds, ds)
    dlist = dlist + pr_02
  if ('flux' in varlist):
    pr_03 = averageProfilesMomentumFluxes(domain, t_inds, h_inds, ds)
    dlist = dlist + pr_03
  if ('tke' in varlist):
    pr_04 = averageProfilesTKE(domain, t_inds, h_inds, ds)
    dlist = dlist + pr_04
  if ('skew' in varlist):
    pr_05 = averageProfilesSkewness(domain, t_inds, h_inds, ds)
    dlist = dlist + pr_05
  return dlist

def compileDataListCompare(domain, t_inds, ct_inds, h_inds, ds, cds, varlist):
  sep = [[''] * (len(h_inds) + 1)]  # empty line string
  dlist = []

  if ('ws' in varlist):
    pr_01 = compareProfilesWS(domain, t_inds, ct_inds, h_inds, ds, cds)
    dlist = dlist + pr_01 + sep
  if ('var' in varlist):
    pr_02 = compareProfilesVariances(domain, t_inds, ct_inds, h_inds, ds, cds)
    dlist = dlist + pr_02 + sep
  if ('flux' in varlist):
    pr_03 = compareProfilesMomentumFluxes(domain, t_inds, ct_inds, h_inds, ds, cds)
    dlist = dlist + pr_03 + sep
  if ('tke' in varlist):
    pr_04 = compareProfilesTKE(domain, t_inds, ct_inds, h_inds, ds, cds)
    dlist = dlist + pr_04 + sep
  return dlist
