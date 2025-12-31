import numpy as np
import math

# Constants
k = 0.01720209895
mu = k**2
cAU = 173.144643267
cAUGAU = cAU/k
eps = math.radians(23.4374)

# Helper Functions:
def angledet(sin,cos):
    y = math.asin(sin)
    x = math.acos(cos)
    if sin < 0:
        return -x
    elif cos < 0:
        return -y
    else:
        return y
def m_func(E,e,M):
    return E - e*math.sin(E) - M
def df(E,e):
    return 1 - e*math.cos(E)
def nr(e,M):
    xlist = [0]
    for i in range(10000):
        xnew = xlist[-1] - m_func(xlist[-1],e,M)/df(xlist[-1],e)
        xlist.append(xnew)
    return xlist[-1]
def rotation(v,a,b):
    #modified for OD
    x = v[0]
    y = v[1]
    z = v[2]
    xx = x*math.cos(a) - y*math.sin(a) 
    yy = x*math.sin(a)*math.cos(b) + y*math.cos(a)*math.cos(b) - z*math.sin(b)
    zz = x*math.sin(b)*math.sin(a) + y*math.sin(b)*math.cos(a) + z*math.cos(b)
    return [xx,yy,zz]
def counterrotate(v):
    # Convert from EQ coord to EC coord
    a = -23.43596*math.pi/180 
    x = v[0]
    y = v[1]
    z = v[2]    
    xx = x 
    yy = y*math.cos(a) - z*math.sin(a)
    zz = y*math.sin(a) + z*math.cos(a)
    return np.array([xx,yy,zz])
def spherical(ra,dec):
    theta = math.pi/180 * ra
    phi = math.pi/180 * dec
    x = math.cos(phi)*math.cos(theta)
    y = math.cos(phi)*math.sin(theta)
    z = math.sin(phi)
    return np.array([x,y,z])
def fg(t,r,v):
    r_mag = np.linalg.norm(r)
    f = 1-t**2/(2*r_mag**3)+np.dot(r,v)*t**3/(2*r_mag**5)+t**4/(24*r_mag**3)*(3*(np.dot(v,v)/r_mag**2-r_mag**-3)-15*(np.dot(r,v)/r_mag**2)**2+r_mag**-3)
    g = t-t**3/(6*r_mag**3)+np.dot(r,v)*t**4/(4*r_mag**5)
    return f,g
def fgend(t1,t2,r,v):
    f1,g1 = fg(t1,r,v)
    f3,g3 = fg(t2,r,v)
    return f1,g1,f3,g3
def triproduct(a,b,c):
    return np.dot(np.cross(a,b),c)
def rho_det(R1,R2,R3,h1rho,h2rho,h3rho,a1,a3):
    rho1 = float((a1*triproduct(R1,h2rho,h3rho)-triproduct(R2,h2rho,h3rho)+a3*triproduct(R3,h2rho,h3rho))/(a1*triproduct(h1rho,h2rho,h3rho)))
    rho2 = float((a1*triproduct(h1rho,R1,h3rho)-triproduct(h1rho,R2,h3rho)+a3*triproduct(h1rho,R3,h3rho))/(-triproduct(h1rho,h2rho,h3rho)))
    rho3 = float((a1*triproduct(h2rho,R1,h1rho)-triproduct(h2rho,R2,h1rho)+a3*triproduct(h2rho,R3,h1rho))/(a3*triproduct(h2rho,h3rho,h1rho)))
    
    return rho1,rho2,rho3
def hms_to_deg(a) -> float:
    '''Converts an angle in degrees, minutes, seconds to decimal degrees
    '''
    a = a.split(':')
    a = np.array([int(a[0]),int(a[1]),float(a[2])])
    hours = a[0]
    minutes = a[1]
    seconds = a[2]
    angle_hour = hours * 15
    angle_minute = minutes * 1/4
    angle_seconds = seconds * 1/240
    
    return angle_hour + angle_minute + angle_seconds
def dms_to_deg(a) -> float:
    '''Converts an angle in degrees, minutes, seconds to decimal degrees
    '''
    a = a.split(':')
    a = np.array([int(a[0]),int(a[1]),float(a[2])])
    degrees = a[0]
    minutes = a[1]
    seconds = a[2]
    if degrees < 0:
        is_neg = True
    else:
        is_neg = False
    angle_minute = minutes / 60
    angle_seconds = seconds/3600
    angle = np.abs(degrees) + angle_minute + angle_seconds
    if is_neg:
        return angle * -1
    else:
        return angle
def colon_to_list(a):
    a = a.split(':')
    a = np.array([int(a[0]),int(a[1]),float(a[2])])
    return a
def julian(time,date):
    [hour, minute, second] = time.split(':')
    [year, month, day] = [float(i) for i in date]
    time = (int(hour) + int(minute)/60 + float(second)/3600 + 6)/24
    j = 367 * year - (7*(year + ((month+9)//12)))//4 + (275*month)//9 + day + 1721013.5 + time
    return j
def T_Get(r,v,JD):
    a,M = orbital_components(r,v)[0], orbital_components(r,v)[5]
    T = JD - (M)/(mu/(a**3))**0.5
    return T
def orbital_components(r,v):
    # Setting up Values
    r_mag = np.linalg.norm(r)
    rx = r[0]
    ry = r[1]
    rz = r[2]
    v_mag = np.linalg.norm(v)
    h = np.cross(r,v)
    h_mag = np.linalg.norm(h)
    hx = h[0]
    hy = h[1]
    hz = h[2]

    # Find a
    a = 1/(2/r_mag-v_mag**2/mu)
    # Find e
    e = (1-h_mag**2/(mu*a))**0.5
    # Find i
    i = math.atan((hx**2+hy**2)**0.5/hz)
    if i < 0:
        i += 2*math.pi

    # Find OMEGA, using angle det function.
    OMEGA = math.atan2(hx/(h_mag*math.sin(i)),-hy/(h_mag*math.sin(i)))
    if OMEGA < 0:
        OMEGA += 2*math.pi

    # Find U for omega
    U = np.atan2(rz/(r_mag*math.sin(i)),(rx*math.cos(OMEGA) + ry*math.sin(OMEGA))/r_mag)

    # Find True Anomalies
    nu = np.atan2((a*(1-e**2)/(e*h_mag))*(np.dot(r,v)/r_mag), 1/e * (a*(1-e**2)/r_mag-1))

    # Find omega
    omega = U-nu
    if omega < 0:
        omega += 2*math.pi
    
    # Find M
    E = math.acos(1/(e)*(1-r_mag/(a)))
    if nu <= math.pi:
        E = -E
    M = E - e*math.sin(E)
    if M <0:
        M += 2*math.pi

    return float(a),float(e),float(i),float(OMEGA),float(omega),float(M)
def orbital_parameter(r,v):
    'Function used to output parameter with units'

    [a,e,i,OMEGA,omega,M] = orbital_components(r,v)
    # Getting ready for export
    a_return = f'{a} AU'
    e_return = f'{e}'
    i_return = f'{i*180/math.pi}°'
    O_return = f'{OMEGA*180/math.pi}°'
    w_return = f'{omega*180/math.pi}°'
    M_return = f'{M*180/math.pi}°'

    return [a_return,e_return,i_return,O_return,w_return,M_return]

# Input Function
def final(path):
    # Opening File
    with open(path, 'r') as file:
        raw = file.readlines()
    data = []

    # Turning into lists
    for i in raw:
        i = i.replace(',','')
        i = i.split()
        data.append(i)

    # Code for cases with more than three lines
    if len(data) < 3:
        return 'Error: Less than three lines of data'
    elif len(data) > 3:
        first = int(input('You have more than 3 data inputs, please enter the first index to choose your data set: '))
        second = int(input('You have more than 3 data inputs, please enter the second index to choose your data set: '))
        third = int(input('You have more than 3 data inputs, please enter the third index to choose your data set: '))
        data = [data[first],data[second],data[third]]
    
    # Define Variables
    RA1 = data[0][4]    
    RA2 = data[1][4]
    RA3 = data[2][4]

    Dec1 = data[0][5]
    Dec2 = data[1][5]
    Dec3 = data[2][5]

    R1 = np.array([float(i) for i in data[0][-3:]])
    R2 = np.array([float(i) for i in data[1][-3:]])
    R3 = np.array([float(i) for i in data[2][-3:]])

    # B/C
    RA1 = hms_to_deg(RA1)
    RA2 = hms_to_deg(RA2)
    RA3 = hms_to_deg(RA3)

    Dec1 = dms_to_deg(Dec1)
    Dec2 = dms_to_deg(Dec2)
    Dec3 = dms_to_deg(Dec3)

    # Figure out initial constants
    hrho1 = spherical(RA1,Dec1)
    hrho2 = spherical(RA2,Dec2)
    hrho3 = spherical(RA3,Dec3)

    # D
    t1 = julian(data[0][3],data[0][:3])
    t2 = julian(data[1][3],data[1][:3])
    t3 = julian(data[2][3],data[2][:3])

    tau1 = k*(t1-t2)
    tau0 = k*(t3-t1)
    tau3 = k*(t3-t2)

    # E
    a1 = tau3/tau0
    a3 = -tau1/tau0

    # F
    rho1,rho2,rho3 = rho_det(R1, R2, R3,hrho1,hrho2,hrho3,a1,a3)

    # G
    r1 = rho1 * hrho1 - R1
    r2 = rho2 * hrho2 - R2
    r3 = rho3 * hrho3 - R3

    v12 = (r2-r1)/(t2-t1)
    v23 = (r3-r2)/(t3-t2)

    rd2 = ((t3-t2)*v12 + (t2-t1)*v23)/(t3-t1)

    # H
    f1,g1,f3,g3 = fgend(tau1,tau3,r2,rd2)

    # I
    a1 = g3/(f1*g3-f3*g1)
    a3 = -g1/(f1*g3-f3*g1)

    # Loop Conditions:
    r2n = r2
    r2n1 = np.array([1,1,1])
    r2dn = rd2
    r2dn1 = np.array([1,1,1])
    
    # Begin Iterative Stage
    while (np.linalg.norm((r2n - r2n1)/r2n) > 10E-10) or (np.linalg.norm((r2dn - r2dn1)/r2dn) > 10E-10):
        # Light Correction
        ct1 = t1 - rho1/cAU
        ct2 = t2 - rho2/cAU
        ct3 = t3 - rho3/cAU
        ctau1 = k*(ct1-ct2)
        ctau3 = k*(ct3-ct2)
        # 1
        rho1,rho2,rho3 = rho_det(R1, R2, R3,hrho1,hrho2,hrho3,a1,a3)

        # 2
        r1 = rho1 * hrho1 - R1
        r2 = rho2 * hrho2 - R2
        r3 = rho3 * hrho3 - R3

        # 3
        r2 = (g3*r1-g1*r3)/(f1*g3-f3*g1)
        rd2 = k*(f3*r1-f1*r3)/(f3*g1-f1*g3)

        # 4
        f1,g1,f3,g3 = fgend(tau1,tau3,r2,rd2)

        # 5
        a1 = g3/(f1*g3-f3*g1)
        a3 = -g1/(f1*g3-f3*g1)
        
        # 6
        r2n1 = r2n
        r2n = r2
        r2dn1 = r2dn
        r2dn = rd2

        print(np.linalg.norm((r2n - r2n1)/r2n),np.linalg.norm((r2dn - r2dn1)/r2dn))
    # Return Values
    # Converting to Eccliptic 
    r2 = counterrotate(r2)
    rd2 = counterrotate(rd2)

    # Components & Mean Anomaly
    JD = julian('6:00:00.000',[2025, 7, 22])
    a,e,i,OMEGA,omega,M = orbital_components(r2,rd2)
    T = T_Get(r2,rd2,t2)
    M_date = (k/a**1.5*(JD-T))%(2*math.pi) * 180/math.pi
    return f'position vector = {r2}, velocity vector = {rd2}, range = {rho2}, orbital components(a,e,i,OMEGA,omega,M) = {orbital_parameter(r2,rd2)}, Mean Anomaly for July 22, 2025 @ 6:00:00 UTC: {M_date}'

def finalmonte(file):

    # Turning into lists
    data = []
    for i in file:
        # print((i))
        i = i.replace(',','')
        i = i.split()
        data.append(i)

    # Define Variables
    RA1 = float(data[0][4])
    RA2 = float(data[1][4])
    RA3 = float(data[2][4])

    Dec1 = float(data[0][5])
    Dec2 = float(data[1][5])
    Dec3 = float(data[2][5])

    R1 = np.array([float(i) for i in data[0][-3:]])
    R2 = np.array([float(i) for i in data[1][-3:]])
    R3 = np.array([float(i) for i in data[2][-3:]])


    # Figure out initial constants
    hrho1 = spherical(RA1,Dec1)
    hrho2 = spherical(RA2,Dec2)
    hrho3 = spherical(RA3,Dec3)

    # D
    t1 = julian(data[0][3],data[0][:3])
    t2 = julian(data[1][3],data[1][:3])
    t3 = julian(data[2][3],data[2][:3])

    tau1 = k*(t1-t2)
    tau0 = k*(t3-t1)
    tau3 = k*(t3-t2)

    # E
    a1 = tau3/tau0
    a3 = -tau1/tau0

    # F
    rho1,rho2,rho3 = rho_det(R1, R2, R3,hrho1,hrho2,hrho3,a1,a3)

    # G
    r1 = rho1 * hrho1 - R1
    r2 = rho2 * hrho2 - R2
    r3 = rho3 * hrho3 - R3

    v12 = (r2-r1)/(t2-t1)
    v23 = (r3-r2)/(t3-t2)

    rd2 = ((t3-t2)*v12 + (t2-t1)*v23)/(t3-t1)

    # H
    f1,g1,f3,g3 = fgend(tau1,tau3,r2,rd2)

    # I
    a1 = g3/(f1*g3-f3*g1)
    a3 = -g1/(f1*g3-f3*g1)

    # Loop Conditions:
    r2n = r2
    r2n1 = np.array([1,1,1])

    # Begin Iterative Stage
    while (np.linalg.norm((r2n - r2n1)/r2n) > 10E-7):
        # Light Correction
        ct1 = t1 - rho1/cAU
        ct2 = t2 - rho2/cAU
        ct3 = t3 - rho3/cAU
        ctau1 = k*(ct1-ct2)
        ctau3 = k*(ct3-ct2)
        # 1
        rho1,rho2,rho3 = rho_det(R1, R2, R3,hrho1,hrho2,hrho3,a1,a3)

        # 2
        r1 = rho1 * hrho1 - R1
        r3 = rho3 * hrho3 - R3

        # 3
        r2 = (g3*r1-g1*r3)/(f1*g3-f3*g1)
        rd2 = k*(f3*r1-f1*r3)/(f3*g1-f1*g3)

        # 4
        f1,g1,f3,g3 = fgend(ctau1,ctau3,r2,rd2)

        # 5
        a1 = g3/(f1*g3-f3*g1)
        a3 = -g1/(f1*g3-f3*g1)
        
        # 6
        r2n1 = r2n
        r2n = r2
        # print(np.linalg.norm((r2n - r2n1)/r2n))

    # Return Values
    # Converting to Eccliptic 
    r2 = counterrotate(r2)
    rd2 = counterrotate(rd2)

    return orbital_components(r2,rd2)
#     f'position vector = {r2n}, velocity vector = {rd2}, range = {rho2}, orbital components(a,e,i,OMEGA,omega,M) = {orbital_parameter(r2,rd2)}'

print(final('/home/mliang/Documents/ssp_code/od/odtext.txt'))

# jpllist = [2.284084364220522,.2782618490944362,3.879772054206133,143.9514580193018,175.7227984036659,3.258787791429181E+02]
# output = [2.3064024808689134, 0.28850863005937377, 3.876063894120386, 143.7352797642847, 175.78603584676117, 326.7591333548703]
# def error(jpl,intake):
#     result = []
#     for i in range(6):
#         result.append((intake[i]-jpl[i])/jpl[i] * 100)
#     return result
# print(error(jpllist,output))