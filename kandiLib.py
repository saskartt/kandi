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
