import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import math
plt.rcParams['figure.figsize']=16,12

def closest(lst, K):
    lst = np.asarray(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]

def spec(img_name):
    im = Image.open(img_name)
    pix = im.load()
    width, height=im.size
    x_all=[]
    y_all=[]
    wl_start=400
    wl_end=700
    wl_range=wl_end-wl_start
    wl_p_pix=wl_range/width

    for i in range(0,width):
        s = 0.0
        c = 0
        for j in range(0,height):
            RGB=pix[i,j]
            s+=0.21*int(RGB[0]) + 0.72*int(RGB[1]) + 0.07*int(RGB[2])#luminosity method
            c += 1.0

        x_all.append(i*wl_p_pix+wl_start)
        y_all.append(s/c)
        
    y_all=np.array(y_all)
    y_all=y_all[::-1]
    index_max=np.where(y_all==max(y_all))
    index_max=np.array(index_max,dtype='i')
    peak_start=(index_max[0][0]-1)*wl_p_pix+wl_start
    last_index=len(index_max)
    peak_end=(index_max[last_index-1][0])*wl_p_pix+wl_start
    max_x=(peak_end+peak_start)/2
    max_in_wl_err=(peak_end-peak_start)/2
    return[x_all,y_all,max_x]
    
x0,y0,max_x0=spec('CUVETTE.png')
x1,y1,max_x1=spec('REDSTOCK.png')
x2,y2,max_x2=spec('RED1.png')
x3,y3,max_x3=spec('RED2.png')
x4,y4,max_x4=spec('RED3.png')
x5,y5,max_x5=spec('RED4.png')
x6,y6,max_x6=spec('RED5.png')
x_unk,y_unk,max_unk=spec('REDSPORT.png')
x_unk1,y_unk1,max_unk1=spec('SLRED.png')



t1=y1/y0
t2=y2/y0
t3=y3/y0
t4=y4/y0
t5=y5/y0
t6=y6/y0
t_unk=y_unk/y0
t_unk1=y_unk1/y0



a1=-(np.log10(t1))
a2=-(np.log10(t2))
a3=-(np.log10(t3))
a4=-(np.log10(t4))
a5=-(np.log10(t5))
a6=-(np.log10(t6))
a_unk=-(np.log10(t_unk))
a_unk1=-(np.log10(t_unk1))




plt.plot(x0,a1,label='STOCK solution (4.53×[10]^(-6))')
plt.plot(x0,a2,label='3.29×[10]^(-6) M')
plt.plot(x0,a3,label='2.59×[10]^(-6) M')
plt.plot(x0,a4,label='2.13×[10]^(-6) M')
plt.plot(x0,a5,label='1.81×[10]^(-6) M')
plt.plot(x0,a6,label='1.58×[10]^(-6) M') 
plt.plot(x0,a_unk,label='YETI DRINK (RED)')
plt.plot(x0,a_unk,label='SL SPORT (RED)')


ab_wl=482 # selected wavelength
plt.axvline(x=ab_wl,linestyle='--',label='%.1f nm'%ab_wl)

c_wl=closest(x0, ab_wl)
print(c_wl)
index=np.where(x0==c_wl)
print(index)

ab1=a1[index[0]][0]
ab2=a2[index[0]][0]
ab3=a3[index[0]][0]
ab4=a4[index[0]][0]
ab5=a5[index[0]][0]
ab6=a6[index[0]][0]
ab_unk=a_unk[index[0]][0]
ab_unk1=a_unk1[index[0]][0]




print('absorbance 1 : '+str((ab1)))
print('absorbance 2 : '+str((ab2)))
print('absorbance 3 : '+str((ab3)))
print('absorbance 4 : '+str((ab4)))
print('absorbance 5 : '+str((ab5)))
print('absorbance 6 : '+str((ab6)))
print('absorbance (YETI drink) : '+str((ab_unk)))
print('absorbance (SL SPORT) : '+str((ab_unk1)))



con=[0.00000158,0.00000181,0.00000213,0.00000259,0.00000330,0.00000453] # concentrations of samples
abso=[ab6,ab5,ab4,ab3,ab2,ab1]
coefs=np.polyfit(con, abso, 1)
fitf = np.poly1d(coefs)
print(fitf)

l=0.025 # optical path length
m=coefs[0]
c=coefs[1]
print(m)
print(c)
print('gradient of graph: '+str(m))
e=m/l
print('absoption coefficient at %.1f nm: '%ab_wl+str(e))

unk_con=(ab_unk-c)/(m)
print('concentration in sportdrink (YETI): '+str(unk_con)+' moldm-3')


err_i=math.sqrt((0.21)**2+(0.72)**2+(0.07)**2)

def err_a(t,y):
    err_a=(np.log10(math.e))*(np.sqrt((1/y)**2+(1/y0)**2))*err_i
    return err_a

err_a1=err_a(t1,y1)
err_a2=err_a(t2,y2)
err_a3=err_a(t3,y3)
err_a4=err_a(t4,y4)
err_a5=err_a(t5,y5)
err_a6=err_a(t6,y6) 
err_a_unk=err_a(t_unk,y_unk)
err_a_unk1=err_a(t_unk1,y_unk1)


err_ab1=err_a1[index[0]][0]
err_ab2=err_a2[index[0]][0]
err_ab3=err_a3[index[0]][0]
err_ab4=err_a4[index[0]][0]
err_ab5=err_a5[index[0]][0]
err_ab6=err_a6[index[0]][0] 
err_ab_unk=err_a_unk[index[0]][0]
err_ab_unk1=err_a_unk1[index[0]][0]



err_abso=[err_ab1,err_ab2,err_ab3,err_ab4,err_ab5,err_ab6 ]
print('error of absorbance :'+str(err_abso))

n=len(err_abso)
s_x2=0
s_x=0
s=0
for i in range(0,n):
    s_x2=s_x2+((con[i]**2)/(err_abso[i]**2))
    s_x=s_x+(con[i]/(err_abso[i]**2))
    s=s+(1/(err_abso[i]**2))
 
d=s*s_x2-s_x**2
err_m=math.sqrt(s/d)
err_c=math.sqrt(s_x2/d)
print('error of gradient: '+str(err_m))
print('error of intersection: '+str(err_c))


err_l=0.05
err_e=e*math.sqrt((err_m/m)**2+(err_l/l)**2)
print('error of absorption coefficient: '+str(err_e))


err_unk_con=unk_con*(math.sqrt((err_ab_unk/ab_unk)**2+(err_e/e)**2+(err_l/l)**2))
print('error of YETI drink concentation: '+str(err_unk_con))

mol=[0.00000107,0.000000535,0.000000306,0.000000214,0.000000165,0.000000143] # mols
v=[0.1,0.1,0.1,0.1,0.1,0.1] #volumes
err_mol=0.00625
err_v=0.001
err_con=[]

for j in range(len(con)):
    err_con.append(con[j]*math.sqrt((err_mol/mol[j])**2+(err_v/v[j])**2))
print('errors of consentrations: '+str(err_con))

plt.legend()
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorbance')
plt.title('The graph of Absorbance vs Wavelength for E129 (Red)')
plt.show()
#pylab.savefig('Spectrum_REDabsorbance.png')
