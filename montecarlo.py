from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import odfinal as od
# np.random.seed(10052008)
np.random.seed(2)

table26 = fits.open('/home/mliang/Downloads/June26.fits')[1].data
table06 = fits.open('/home/mliang/Downloads/july6.fits')[1].data
table08 = fits.open('/home/mliang/Downloads/july8.fits')[1].data
fra26 = table26['field_ra']
ira26 = table26['index_ra']
fra06 = table06['field_ra']
ira06 = table06['index_ra']
fra08 = table08['field_ra']
ira08 = table08['index_ra']
# print(ra)
sum26 = 0
sum06 = 0
sum08 = 0
N26 = len(fra26)
N06 = len(fra06)
N08 = len(fra08)
for i in range(N26):
    sum26 += (fra26[i] - ira26[i])**2
for i in range(N06):
    sum06 += (fra06[i] - ira06[i])**2
for i in range(N08):
    sum08 += (fra08[i] - ira08[i])**2
sdra26 = (sum26/N26)**0.5
sdra06 = (sum06/N06)**0.5
sdra08 = (sum08/N08)**0.5


fdec26 = table26['field_dec']
idec26 = table26['index_dec']
fdec06 = table06['field_dec']
idec06 = table06['index_dec']
fdec08 = table08['field_dec']
idec08 = table08['index_dec']

sum26 = 0
sum06 = 0
sum08 = 0
N26 = len(fdec26)
N06 = len(fdec06)
N08 = len(fdec08)
for i in range(N26):
    sum26 += (fdec26[i] - idec26[i])**2
for i in range(N06):
    sum06 += (fdec06[i] - idec06[i])**2
for i in range(N08):
    sum08 += (fdec08[i] - idec08[i])**2
sddec26 = (sum26/N26)**0.5
sddec06 = (sum06/N06)**0.5
sddec08 = (sum08/N08)**0.5

alist = []
elist = []
ilist = []
Olist = []
wlist = []
Mlist = []
def montecarlo(k):
    n = 0
    for i in range(k):
        RA26 = np.random.normal(0,sdra26)
        Dec26 = np.random.normal(0,sddec26)
        RA06 = np.random.normal(0,sdra06)
        Dec06 = np.random.normal(0,sddec06)
        RA08 = np.random.normal(0,sdra08)
        Dec08 = np.random.normal(0,sddec08)
        # print(237.7970833+RA26,-12.9731667+Dec26)
        # print(236.9653417+RA06, -13.3891528+Dec06)
        # print(236.9172500+RA08, -13.5061889+Dec08)
        info = [f'2025 06 26 06:00:46 {237.7970833+RA26} {-12.9731667+Dec26} -8.068103668181402E-02   9.297295419704680E-01   4.029832670874380E-01', 
                f'2025 07 06 05:25:49 {236.9653417+RA06} {-13.3891528+Dec06} -2.470833223422502E-01   9.048284118859414E-01   3.921846503148367E-01', 
                f'2025 07 08 05:16:13 {236.9172500+RA08} {-13.5061889+Dec08} -2.796364694720073E-01   8.967895400402173E-01   3.886998032747111E-01']
        a,e,i,OMEGA,omega,M = od.finalmonte(info)
        i *= 180/np.pi
        OMEGA *= 180/np.pi
        omega *= 180/np.pi
        M *= 180/np.pi
        alist.append(a)
        elist.append(e)
        ilist.append(i)
        Olist.append(OMEGA)
        wlist.append(omega)
        Mlist.append(M)
        n+= 1
        print(n)
    amean = sum(alist)/k
    emean = sum(elist)/k
    imean = sum(ilist)/k
    Omean = sum(Olist)/k
    wmean = sum(wlist)/k
    Mmean = sum(Mlist)/k
    asd = np.std(alist)
    esd = np.std(elist)
    isd = np.std(ilist)
    Osd = np.std(Olist)
    wsd = np.std(wlist)
    Msd = np.std(Mlist)
    plt.figure(figsize=(13, 7)) 
    plt.hist(alist,color = 'powderblue', bins=int(k**0.5))
    plt.axvline(2.284084364220522,color = 'red', label = 'JPL')
    plt.axvline(amean + asd,color = 'black',linestyle='dotted',label = 'Standard Deviation')
    plt.axvline(amean - asd,color = 'black',linestyle='dotted')
    plt.axvline(amean,color = 'blue', label = 'Mean Value')
    plt.xlabel("Length (AU)")
    plt.ylabel("frequency")
    plt.title("Semi-Major Axis (a)")
    plt.legend()
    plt.figure(figsize=(13, 7)) 
    plt.hist(elist,color = 'powderblue', bins=int(k**0.5))
    plt.axvline(.2782618490944362,color = 'red', label = 'JPL')
    plt.axvline(emean + esd,color = 'black', linestyle='dotted',label = 'Standard Deviation')
    plt.axvline(emean - esd,color = 'black', linestyle='dotted')
    plt.axvline(emean,color = 'blue', label = 'Mean Value')
    plt.xlabel("Ratio")
    plt.ylabel("frequency")
    plt.title("Orbit Eccentricity (e)")
    plt.legend()
    plt.figure(figsize=(13, 7)) 
    plt.hist(ilist,color = 'powderblue', bins=int(k**0.5))
    plt.axvline(3.879772054206133,color = 'red', label = 'JPL')
    plt.axvline(imean + isd,color = 'black', linestyle='dotted',label = 'Standard Deviation')
    plt.axvline(imean - isd,color = 'black', linestyle='dotted')
    plt.axvline(imean,color = 'blue', label = 'Mean Value')
    plt.xlabel("Degrees (°)")
    plt.ylabel("frequency")
    plt.title("Inclination (i)")
    plt.legend()
    plt.figure(figsize=(13, 7)) 
    plt.hist(Olist,color = 'powderblue', bins=int(k**0.5))
    plt.axvline(143.95145801930183,color = 'red', label = 'JPL')
    plt.axvline(Omean + Osd,color = 'black', linestyle='dotted',label = 'Standard Deviation')
    plt.axvline(Omean - Osd,color = 'black',linestyle='dotted')
    plt.axvline(Omean,color = 'blue', label = 'Mean Value')
    plt.xlabel("Degrees (°)")
    plt.ylabel("frequency")
    plt.title("Longitude of Ascending Node (Ω)")
    plt.legend()
    plt.figure(figsize=(13, 7))  
    plt.hist(wlist,color = 'powderblue', bins=int(k**0.5))
    plt.axvline(175.7227984036659,color = 'red', label = 'JPL')
    plt.axvline(wmean + wsd,color = 'black', linestyle='dotted',label = 'Standard Deviation')
    plt.axvline(wmean - wsd,color = 'black', linestyle='dotted')
    plt.axvline(wmean,color = 'blue', label = 'Mean Value')
    plt.xlabel("Degrees (°)")
    plt.ylabel("frequency")
    plt.title("Argument of Periapsis (ω)")
    plt.legend()
    plt.figure(figsize=(13, 7)) 
    plt.hist(Mlist,color = 'powderblue', bins=int(k**0.5))
    plt.axvline(3.258787791429181E+02,color = 'red', label = 'JPL')
    plt.axvline(Mmean + Msd,color = 'black', linestyle='dotted',label = 'Standard Deviation')
    plt.axvline(Mmean - Msd,color = 'black', linestyle='dotted')
    plt.axvline(Mmean,color = 'blue', label = 'Mean Value')
    plt.xlabel("Degrees (°)")
    plt.ylabel("frequency")
    plt.title("Mean Anomaly (M)")
    plt.legend()
    plt.show()
    print(amean,emean,imean,Omean,wmean,Mmean)
    print(sddec06,sddec08,sddec26,sdra06,sdra08,sdra26)
    print(asd,esd,isd,Osd,wsd,Msd)
    print(100*(amean-2.284084364220522)/2.284084364220522, 
          100*(emean-.2782618490944362)/.2782618490944362,
          100*(imean-3.879772054206133)/3.879772054206133,
          100*(Omean-143.95145801930183)/143.95145801930183,
          100*(wmean-175.7227984036659)/175.7227984036659,
          100*(Mmean-3.258787791429181E+02)/3.258787791429181E+02)

    # return alist, elist, ilist, Olist, wlist, Mlist

print(montecarlo(7500000))
