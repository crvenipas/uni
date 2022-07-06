import pandas as pd
import numpy as np

A = np.deg2rad(192.86)
D = np.deg2rad(27.13)

df = pd.read_table("catalog_voitko.dat", sep="\s+", usecols=['RA','Dec','Dist','LogT','PMA','PMD','Vr'])

RA_rad = np.deg2rad(df.loc[:,'RA'])
Dec_rad = np.deg2rad(df.loc[:,'Dec'])

### 1 ###
b = np.arcsin(np.sin(Dec_rad) * np.sin(D) + np.cos(Dec_rad) * np.cos(D) * np.cos(RA_rad - A))
la =((np.cos(Dec_rad) * np.sin(RA_rad - A)) / np.cos(b))
la = np.arccos(la)

l_n = np.deg2rad(180) - la
l_n1 = np.deg2rad(360) + la

la = np.cos(la)
la[la < 0] = l_n
l_s = np.sin(np.arccos(la))
la[not la.empty > 0 and l_s < 0] = l_n1

la = np.rad2deg(la)
l = la + np.deg2rad(33)
l = np.deg2rad(l)
l_t = l - np.deg2rad(360)
l[l > np.deg2rad(360)] = l_t

b = np.rad2deg(b)
l = np.rad2deg(l)

### 2 ###
gc = pd.DataFrame({'l' : l , 'b' : b})
df['RA'] = gc['l']
df['Dec']= gc['b']
df.columns = ['l','b','Dist','LogT','PMA','PMD','Vr']
df = df.drop(df[(df.l < 10) & (df.l > -10)].index)
df = df.drop(df[(df.l < 190) & (df.l > 170)].index)

### 3 ###
b = np.deg2rad(b)
l = np.deg2rad(l)

sin_phi = (np.cos(D) * np.sin(RA_rad - A)) / np.cos(b)
cos_phi = (np.sin(D) - np.sin(b)*np.sin(Dec_rad)) / np.cos(b)*np.cos(Dec_rad)

pma = df.loc[:,'PMA']
pmd = df.loc[:,'PMD']

pml = pma * cos_phi + pmd * sin_phi
pmb = -(pma * sin_phi + pmd * cos_phi)

df.columns = ['l','b','Dist','LogT','pml','pmb','Vr']
### 4 ###
Vr = df.loc[:,'Vr']
r = (df.loc[:, 'Dist']) / 1000
Vl = 4.74 * r * pml
Vb = 4.74 * r * pmb

u = Vr*np.cos(b)*np.cos(l) - Vl*np.sin(l) - Vb*np.sin(b)*np.cos(l)
v = Vr*np.cos(b)*np.sin(l) - Vl*np.cos(l) - Vb*np.sin(b)*np.sin(l)
w = Vr*np.sin(b) + Vl*0 + Vb*np.cos(b)

### 5 ###
Vr_lsr = Vr + (u * np.cos(l)*np.cos(b) + v*np.sin(l)*np.cos(b) + w*np.sin(b))
Vl_lsr = Vl + -(u *np.sin(l) + v*np.cos(l))
Vb_lsr = Vb + -(u* np.cos(l)*np.sin(b) - v*np.sin(l)*np.sin(b) + w*np.cos(b))
del df['Vr']
del df['Dist']
#df['r']=r
df['Vr_lsr']= Vr_lsr
df['Vl_lsr']= Vl_lsr
df['Vb_lsr']= Vb_lsr

### 6 ###
R0 = 8.3
R = np.sqrt(R0**2 + r**2 * np.cos(b) * np.cos(b) - 2 * R0 * r * np.cos(l) *np.cos(b))
df['R']=R
df = df.drop(df[(df.R < 6.5)].index)
df = df.drop(df[(df.R > 15)].index)

### 7, 8 ### 
omega0 = (Vr_lsr * (R0 * np.cos(l) - r * np.cos(b)))/(R0 * r * np.sin(l)*np.cos(b)*np.cos(b)) - ((Vl_lsr) / (r * np.cos(b)))
df['omega0'] = omega0
'''df = df.drop(df[(df.omega0 < 0)].index)
df = df.drop(df[(df.omega0 > 100)].index)
'''
### 9, 10 ###


### 12 ###
omega = (Vr_lsr / (R0 * np.sin(l) * np.cos(b))) + omega0
df['omega'] = omega

### 11 ###
df['w-w0'] = omega - omega0

### 13 ###
df.to_csv('edited.dat', sep='\t')

print(df)
df.plot.scatter('R', 'w-w0')
