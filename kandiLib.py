import numpy as np
import sys

# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def createTableRow(qty,label):
  qty = qty.astype(str)
  qty = np.insert(qty, 0, label)
  return qty

def averageProfilesWS(statdomain, t_inds, h_inds, ds):
  # Mean values
  u = np.mean(ds.variables['u_'+statdomain][t_inds,h_inds],axis=0)
  str_u = createTableRow(u,'u')

  v = np.mean(ds.variables['v_'+statdomain][t_inds,h_inds],axis=0)
  str_v = createTableRow(v,'v')

  w = np.mean(ds.variables['w_'+statdomain][t_inds,h_inds],axis=0)
  str_w = createTableRow(w,'w')

  U = np.linalg.norm([u,v],axis=0)
  str_U = createTableRow(U,'U')
  return [str_u,str_v,str_w,str_U]

def averageProfilesVariances(statdomain, t_inds, h_inds, ds):
  # Variances
  var_u = np.mean(ds.variables['u*2_'+statdomain][t_inds,h_inds],axis=0)
  var_u = createTableRow(var_u,'u*2')

  var_v = np.mean(ds.variables['v*2_'+statdomain][t_inds,h_inds],axis=0)
  var_v = createTableRow(var_v,'v*2')

  var_w = np.mean(ds.variables['w*2_'+statdomain][t_inds,h_inds],axis=0)
  var_w = createTableRow(var_w,'w*2')

  return [var_u,var_v,var_w]

def averageProfilesMomentumFluxes(statdomain,t_inds, h_inds, ds):
  flx_u = np.mean(ds.variables['wu_'+statdomain][t_inds,h_inds],axis=0)
  flx_u = createTableRow(flx_u, "wu")

  flx_v = np.mean(ds.variables['wv_'+statdomain][t_inds,h_inds],axis=0)
  flx_v = createTableRow(flx_v, "wv")

  return [flx_u, flx_v]

def averageProfilesTKE(statdomain, t_inds, h_inds, ds):
  tke_e = np.mean(ds.variables['e_'+statdomain][t_inds,h_inds],axis=0)
  tke_e = createTableRow(tke_e, "TKE")

  return [tke_e]

def compareProfilesWS(statdomain, t_inds, ct_inds, h_inds, ds, cds):
  u = np.mean(ds.variables['u_'+statdomain][t_inds,h_inds],axis=0)
  cu = np.mean(cds.variables['u_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_u = createTableRow(np.divide(u,cu)*100.,'u (%)')

  v = np.mean(ds.variables['v_'+statdomain][t_inds,h_inds],axis=0)
  cv = np.mean(cds.variables['v_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_v = createTableRow(np.divide(v,cv)*100.,'v (%)')

  w = np.mean(ds.variables['w_'+statdomain][t_inds,h_inds],axis=0)
  cw = np.mean(cds.variables['w_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_w = createTableRow(np.divide(w,cw)*100.,'w (%)')

  U = np.linalg.norm([u,v],axis=0)
  cU = np.linalg.norm([cu,cv],axis=0)
  cmp_U = createTableRow(np.divide(U,cU)*100.,'U (%)')

  return [cmp_u, cmp_v, cmp_w, cmp_U]

def compareProfilesVariances(statdomain, t_inds, ct_inds, h_inds, ds, cds):
  var_u = np.mean(ds.variables['u*2_'+statdomain][t_inds,h_inds],axis=0)
  var_cu = np.mean(cds.variables['u*2_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_var_u = createTableRow(np.divide(var_u,var_cu)*100.,'u*2 (%)')

  var_v = np.mean(ds.variables['v*2_'+statdomain][t_inds,h_inds],axis=0)
  var_cv = np.mean(cds.variables['v*2_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_var_v = createTableRow(np.divide(var_v,var_cv)*100.,'v*2 (%)')

  var_w = np.mean(ds.variables['w*2_'+statdomain][t_inds,h_inds],axis=0)
  var_cw = np.mean(cds.variables['w*2_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_var_w = createTableRow(np.divide(var_w,var_cw)*100.,'w*2 (%)')

  return [cmp_var_u, cmp_var_v, cmp_var_w]

def compareProfilesMomentumFluxes(statdomain, t_inds, ct_inds, h_inds, ds, cds):
  flx_u = np.mean(ds.variables['wu_'+statdomain][t_inds,h_inds],axis=0)
  flx_cu = np.mean(cds.variables['wu_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_flx_u = createTableRow(np.divide(flx_u,flx_cu)*100.,'wu (%)')

  flx_v = np.mean(ds.variables['wv_'+statdomain][t_inds,h_inds],axis=0)
  flx_cv = np.mean(cds.variables['wv_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_flx_v = createTableRow(np.divide(flx_v,flx_cv)*100.,'wv (%)')

  return [cmp_flx_u, cmp_flx_v]

def compareProfilesTKE(statdomain, t_inds, ct_inds, h_inds, ds, cds):
  tke_e = np.mean(ds.variables['e_'+statdomain][t_inds,h_inds],axis=0)
  tke_ce = np.mean(cds.variables['e_'+statdomain][ct_inds,h_inds],axis=0)
  cmp_tke_e = createTableRow(np.divide(tke_e,tke_ce)*100.,'TKE (%)')

  return [cmp_tke_e]
