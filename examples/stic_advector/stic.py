import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from math import sqrt, cos, sin, pi, pow, acos
from scipy.optimize import NonlinearConstraint, LinearConstraint, minimize
import imageio, os

VFTOL = 1.0e-12

def forward(x0,y0,hx,hy,anx,any,d):
    nx = anx* hx; ny = any * hy
    mag = abs(nx) + abs(ny)
    if (mag < 1.0e-15):
        mag = 1.0e-15
    alpha = (d - (anx*x0+any*y0)) / mag;
    nx /= mag
    ny /= mag
    if (nx < 0):
        nx = -nx
        alpha += nx
    if (ny < 0):
        ny = -ny
        alpha += ny
    if (alpha <= 0):
        return 0
    elif (alpha >= 1):
        return 1
    m = nx
    if (ny < nx): 
        m = ny
    alpha_0 = alpha
    if (alpha > 0.5):
        alpha_0 = 1.0 - alpha    
    tmp_vol = 0.0
    if (alpha_0 < m):
        tmp_vol = alpha_0**2/(2.0*m*(1.0-m))
    else:
        tmp_vol = alpha_0/(1.0-m)-m/(2.0*(1.0-m))
    if (alpha > 0.5):
        return (1.0 - tmp_vol)
    else:
        return tmp_vol

def forward3d(x0,y0,z0,hx,hy,hz,anx,any,anz,d):
    nx = anx* hx; ny = any * hy; nz = anz * hz
    mag = abs(nx) + abs(ny) + abs(nz)
    if (mag < 1.0e-15):
        mag = 1.0e-15
    alpha = (d - (anx*x0+any*y0+anz*z0)) / mag;
    nx /= mag
    ny /= mag
    nz /= mag
    if (nx < 0):
        nx = -nx
        alpha += nx
    if (ny < 0):
        ny = -ny
        alpha += ny
    if (nz < 0):
        nz = -nz
        alpha += nz
    if (alpha <= 0):
        return alpha
    elif (alpha >= 1):
        return alpha
    m_array = np.sort(np.array([nx,ny,nz]))
    m1 = m_array[0]
    m2 = m_array[1]
    m12 = m1+m2
    m3 = m_array[2]
    m = min(m12,m3)
    V1=m1**2/max(6.0*m2*m3,1.0e-15)
    alpha_0 = alpha
    if (alpha > 0.5):
        alpha_0 = 1.0 - alpha    
    tmp_vol = 0.0
    if (alpha_0 < m1):
        tmp_vol = alpha_0**3/(6.0*m1*m2*m3)
    elif (alpha_0 < m2):
        tmp_vol = alpha_0*(alpha_0-m1)/(2.0*(m2*m3))+V1
    elif (alpha_0 < m):
        tmp_vol = (alpha_0**2*(3.0*m12-alpha_0)+m1**2*(m1-3.0*alpha_0)+m2**2*(m2-3.0*alpha_0))/(6.0*m1*m2*m3)
    elif (m3 < m12):
        tmp_vol = (alpha_0**2*(3.0-2.0*alpha_0)+m1**2*(m1-3.0*alpha_0)+m2**2*(m2-3.0*alpha_0)+m3**2*(m3-3.0*alpha_0))/(6.0*m1*m2*m3)
    else:
        tmp_vol = (2.0*alpha_0-m12)/(2.0*m3)
    if (alpha > 0.5):
        return (1.0 - tmp_vol)
    else:
        return tmp_vol
    
def forward_hlvira(x0,y0,hx,hy,anx,any,d):
    nx = anx* hx; ny = any * hy
    mag = abs(nx) + abs(ny)
    if (mag < 1.0e-15):
        mag = 1.0e-15
    alpha = (d - (anx*x0+any*y0)) / mag;
    nx /= mag
    ny /= mag
    if (nx < 0):
        nx = -nx
        alpha += nx
    if (ny < 0):
        ny = -ny
        alpha += ny
    if (alpha <= 0):
        return alpha
    elif (alpha >= 1):
        return alpha
    m = nx
    if (ny < nx): 
        m = ny
    alpha_0 = alpha
    if (alpha > 0.5):
        alpha_0 = 1.0 - alpha    
    tmp_vol = 0.0
    if (alpha_0 < m):
        tmp_vol = alpha_0**2/(2.0*m*(1.0-m))
    else:
        tmp_vol = alpha_0/(1.0-m)-m/(2.0*(1.0-m))
    if (alpha > 0.5):
        return (1.0 - tmp_vol)
    else:
        return tmp_vol

def backward(x0,y0,h,nx,ny,vf):
    if (vf < VFTOL):
        return nx*x0+ny*y0
    elif (vf > 1.0-VFTOL):
        return nx*(x0+h)+ny*(y0+h)
    
    norm = abs(nx) + abs(ny)
    m = min(abs(nx/norm),abs(ny/norm))
    V1 = m/(2.0*(1.0-m))

    V0 = vf
    if (vf > 0.5):
        V0 = 1.0 - vf    
    tmp_alpha = 0.0
    if (V0 < V1):
        tmp_alpha = sqrt(2.0*m*(1.0-m)*V0)
    else:
        tmp_alpha = V0*(1.0-m)+0.5*m
    alpha = 0.0
    if (vf > 0.5):
        alpha = (1.0 - tmp_alpha)
    else:
        alpha = tmp_alpha
    if (nx < 0):
        alpha += nx/norm
    if (ny < 0):
        alpha += ny/norm
    return alpha * (h*norm) + (nx*x0+ny*y0)

def backward3d(x0,y0,z0,hx,hy,hz,anx,any,anz,vf):
    nx = anx* hx; ny = any * hy; nz = anz * hz
    if (vf < VFTOL):
        return anx*x0+any*y0+anz*z0
    elif (vf > 1.0-VFTOL):
        return anx*(x0+hx)+any*(y0+hy)+anz*(z0+hz)
    mag = abs(nx) + abs(ny) + abs(nz)
    if (mag < 1.0e-15):
        mag = 1.0e-15
    nx /= mag
    ny /= mag
    nz /= mag
    m_array = np.sort(np.array([abs(nx),abs(ny),abs(nz)]))
    m1 = m_array[0]
    m2 = m_array[1]
    m12 = m1+m2
    m3 = m_array[2]
    V1=m1**2/max(6.0*m2*m3,1.0e-15)
    V2=V1+0.5*(m2-m1)/m3
    V3=0
    if (m3<m12):
        V3=(m3**2*(3.0*m12-m3)+m1**2*(m1-3.0*m3)+m2**2*(m2-3.0*m3))/(6.0*m1*m2*m3)
    else:
        V3=0.5*m12/m3

    V0 = 0
    if (vf > 0.5):
        V0 = 1.0 - vf    
    else:
        V0 = vf    
    tmp_alpha = 0.0
    if (V0 < V1):
        tmp_alpha = pow(6.0*m1*m2*m3*V0,1.0/3.0)
    elif (V0 < V2):
        tmp_alpha = 0.5*(m1 + sqrt(m1**2+8.0*m2*m3*(V0-V1)))
    elif (V0 < V3):
        a0 = -(m1**3 + m2**3 - 6.0*m1*m2*m3 * V0)
        a1 = 3.0 * (m1**2 + m2**2)
        a2 = -3.0 * m12
        p0 = -(a1 / 3.0 - a2**2 / 9.0)
        q0 = (a1 * a2 - 3.0 * a0) / 6.0 - a2**3 / 27.0
        theta = acos(q0 / sqrt(p0**3)) / 3.0
        tmp_alpha = sqrt(p0) * (sqrt(3.0) * sin(theta) - cos(theta)) - a2 / 3.0
    elif (m3 < m12):
        a0 = -0.5 * (m1**3 + m2**3 + m3**3 - 6.0 * m1*m2*m3 * V0)
        a1 = 1.5 * (m1**2 + m2**2 + m3**2)
        a2 = -1.5;
        p0 = -(a1 / 3.0 - a2**2 / 9.0)
        q0 = (a1 * a2 - 3.0 * a0) / 6.0 - a2**3 / 27.0
        theta = acos(q0 / sqrt(p0**3)) / 3.0
        tmp_alpha = sqrt(p0) * (sqrt(3.0) * sin(theta) - cos(theta)) - a2 / 3.0
    else:
        tmp_alpha = m3*V0 + 0.5*m12
    if (vf > 0.5):
        tmp_alpha = (1.0 - tmp_alpha)
    if (nx < 0):
        tmp_alpha += nx
    if (ny < 0):
        tmp_alpha += ny
    if (nz < 0):
        tmp_alpha += nz
    alpha = tmp_alpha
    return alpha * mag + (anx*x0+any*y0+anz*z0)

def elvira(N,h,vf):
    norm = np.zeros((N, N, 2))
    dist = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            norm[i][j][0] = 1.0/sqrt(2.0)
            norm[i][j][1] = 1.0/sqrt(2.0)
            dist[i][j] = 0.0
    for i in range(1,N-1):
        for j in range(1,N-1):
            if (vf[i][j] >= 1.0-VFTOL):
                norm[i][j][0] = 1.0/sqrt(2.0)
                norm[i][j][1] = 1.0/sqrt(2.0)
                dist[i][j] = sqrt(2.0)
            elif (vf[i][j] > VFTOL):
                x0 = i / N; x1 = x0 + h
                y0 = j / N; y1 = y0 + h
                # Calculation differents normal options
                m = np.zeros((6,2))
                for k in range(0,3):
                    if (vf[i][j+1] > vf[i][j-1]):
                        m[k][1] = -1.0
                    else:
                        m[k][1] = 1.0
                    if (vf[i+1][j] > vf[i-1][j]):
                        m[3+k][0] = -1.0
                    else:
                        m[3+k][0] = 1.0
                for k in range(-1,2):
                    m[0][0] -= vf[i][j+k] - vf[i-1][j+k]
                    m[1][0] -= 0.5*(vf[i+1][j+k] - vf[i-1][j+k])
                    m[2][0] -= vf[i+1][j+k] - vf[i][j+k]
                    m[3][1] -= vf[i+k][j] - vf[i+k][j-1]
                    m[4][1] -= 0.5*(vf[i+k][j+1] - vf[i+k][j])
                    m[5][1] -= vf[i+k][j+1] - vf[i+k][j]
                # Normalize normals
                for k in range(6):
                    mag = sqrt(m[k][0]**2+m[k][1]**2)
                    m[k][0]/=mag; m[k][1]/=mag
                # Try normals
                E = np.zeros(6)
                maxE = 1000.0
                for k in range(6):
                    d = backward(x0,y0,h,m[k][0],m[k][1],vf[i][j])
                    for l in range(-1,2):
                        for n in range(-1,2):
                            x00 = x0+l*h; y00 = y0+n*h
                            tmp_vf = forward(x00,y00,h,h,m[k][0],m[k][1],d)
                            E[k] += (tmp_vf - vf[i+l][j+n])**2
                    E[k] = sqrt(E[k])
                    if (E[k] < maxE):
                        maxE = E[k]
                        norm[i][j][0] = m[k][0]
                        norm[i][j][1] = m[k][1]
                        dist[i][j] = d
                    
    return norm, dist

def cost_lvira(plane,N,h,vf,i,j):
    x0 = i / N; x1 = x0 + h
    y0 = j / N; y1 = y0 + h
    E = 0.0
    # if (vf[i][j] > 0 and vf[i][j] < 1):
    #     plane[2] = backward(x0,y0,h,plane[0],plane[1],vf[i][j])
    for l in range(-1,2):
        for n in range(-1,2):
            x00 = x0+l*h; y00 = y0+n*h
            tmp_vf = forward(x00,y00,h,h,plane[0],plane[1],plane[2])
            E += (tmp_vf - vf[i+l][j+n])**2
    return E

def lvira(N,h,vf,old_norm,old_dist):
    norm = old_norm
    dist = old_dist
    for i in range(1,N-1):
        for j in range(1,N-1):
            norm[i][j][0] = old_norm[i][j][0]
            norm[i][j][1] = old_norm[i][j][1]
            dist[i][j] = old_dist[i][j]
            if (vf[i][j] > VFTOL and vf[i][j] < 1.0-VFTOL):
                x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
                y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)
                plane_guess = np.zeros(3)
                plane_guess[0] = old_norm[i][j][0]
                plane_guess[1] = old_norm[i][j][1]
                plane_guess[2] = old_dist[i][j]

                def M0_constraint_lvira(plane):
                    tmp_vf = forward(x0,y0,h,h,plane[0],plane[1],plane[2])
                    if (tmp_vf == 0 and vf[i][j] > 0):
                        alpha = abs(plane[2] - (plane[0]*x+plane[1]*y))/h - 0.5;
                        return -alpha
                    elif (tmp_vf == 1  and vf[i][j] < 1):
                        alpha = abs(plane[2] - (plane[0]*x+plane[1]*y))/h - 0.5;
                        return alpha
                    else:
                        return tmp_vf - vf[i][j]
                    
                def callback_lvira(plane):
                    plane[2] = backward(x0,y0,h,plane[0],plane[1],vf[i][j])

                M0_constraint = NonlinearConstraint(M0_constraint_lvira, -VFTOL, +VFTOL)
                res = minimize(cost_lvira,plane_guess,args=(N,h,vf,i,j),callback=callback_lvira,constraints=[M0_constraint],options={'maxiter':100},tol=VFTOL)
                norm[i][j][0] = res.x[0]
                norm[i][j][1] = res.x[1]
                dist[i][j] = backward(x0,y0,h,norm[i][j][0],norm[i][j][1],vf[i][j])
                # if (not res.success):
                #     print(vf[i][j],M0_constraint_lvira(res.x))
                #     print(res)

    return norm, dist

def plot_plic(N,h,vf,norm,dist,ax):
    for i in range(1,N-1):
        for j in range(1,N-1):
            if (vf[i][j] > 1e-10 and vf[i][j] < 1-1e-10):
                x0 = i / N; x1 = x0 + h
                y0 = j / N; y1 = y0 + h
                ptx = []; pty = []
                npt = 0
                if (abs(norm[i][j][1]) > 0.0):
                    # X0
                    y = (dist[i][j] - norm[i][j][0]*x0)/norm[i][j][1]
                    if (y > y0 and y <= y1):
                        ptx.append(-0.5+x0/h)
                        pty.append(-0.5+y/h)
                        npt += 1
                    # X1
                    y = (dist[i][j] - norm[i][j][0]*x1)/norm[i][j][1]
                    if (y >= y0 and y < y1):
                        ptx.append(-0.5+x1/h)
                        pty.append(-0.5+y/h)
                        npt += 1
                if (abs(norm[i][j][0]) > 0.0):
                    # Y0
                    x = (dist[i][j] - norm[i][j][1]*y0)/norm[i][j][0]
                    if (x >= x0 and x < x1):
                        ptx.append(-0.5+x/h) 
                        pty.append(-0.5+y0/h)
                        npt += 1
                    # Y1
                    x = (dist[i][j] - norm[i][j][1]*y1)/norm[i][j][0]
                    if (x > x0 and x <= x1):
                        ptx.append(-0.5+x/h)
                        pty.append(-0.5+y1/h)
                        npt += 1
                if (npt == 2):
                    ax.plot(ptx,pty,'r')

def U(x,y,t):
    # return -2.0 * sin(pi * x)**2 * sin(pi * y) * cos(pi * y) * cos(pi * t / 8.0)
    return -2.0 * pi * (y-0.5)

def V(x,y,t):
    # return 2.0 * sin(pi * y)**2 * sin(pi * x) * cos(pi * x) * cos(pi * t / 8.0)
    return 2.0 * pi * (x-0.5);

def split_advect(N,h,old_vf,t,dt):
    vf = np.copy(old_vf)
    cc = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (vf[i][j] > 0.5):
                cc[i][j] = 1.0

    # ELVIRA
    norm, dist = elvira(N,h,vf)
    # # Store old normals and old distance
    # old_norm = norm
    # old_dist = dist
    # norm, dist = lvira(N,h,vf,old_norm,old_dist)
    fluxes = np.zeros(2*(N+1)*N)

    for i in range(1,N-1):
        for j in range(1,N-1):
            x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
            y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)
            ul = U(x0,y,t)
            ur = U(x1,y,t)

            # Sweep x
            Fl = 0; Fr = 0
            if (ul > 0):
                Fl += forward(x0-ul*dt,y0,ul*dt,h,norm[i-1][j][0],norm[i-1][j][1],dist[i-1][j])*abs(ul)*dt/h
            else:
                Fl -= forward(x0,y0,-ul*dt,h,norm[i][j][0],norm[i][j][1],dist[i][j])*abs(ul)*dt/h
            if (ur > 0):
                Fr -= forward(x1-ur*dt,y0,ur*dt,h,norm[i][j][0],norm[i][j][1],dist[i][j])*abs(ur)*dt/h
            else:
                Fr += forward(x1,y0,-ur*dt,h,norm[i+1][j][0],norm[i+1][j][1],dist[i+1][j])*abs(ur)*dt/h
            vf[i][j] += (Fl + Fr) + cc[i][j]*dt*(ur-ul)/h
            # if (vf[i][j] < 0 or vf[i][j] > 1):
            #     print(vf[i][j])
            fluxes[i+(N+1)*j] = Fl
            fluxes[(i+1)+(N+1)*j] = -Fr

    # ELVIRA
    norm, dist = elvira(N,h,vf)
    # # Store old normals and old distance
    # old_norm = norm
    # old_dist = dist
    # norm, dist = lvira(N,h,vf,old_norm,old_dist)

    for i in range(1,N-1):
        for j in range(1,N-1):
            x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
            y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)
            vb = V(x,y0,t)
            vt = V(x,y1,t)

            # Sweep y
            Fb = 0; Ft = 0
            if (vb > 0):
                Fb += forward(x0,y0-vb*dt,h,vb*dt,norm[i][j-1][0],norm[i][j-1][1],dist[i][j-1])*abs(vb)*dt/h
            else:
                Fb -= forward(x0,y0,h,-vb*dt,norm[i][j][0],norm[i][j][1],dist[i][j])*abs(vb)*dt/h
            if (vt > 0):
                Ft -= forward(x0,y1-vt*dt,h,vt*dt,norm[i][j][0],norm[i][j][1],dist[i][j])*abs(vt)*dt/h
            else:
                Ft += forward(x0,y1,h,-vt*dt,norm[i][j+1][0],norm[i][j+1][1],dist[i][j+1])*abs(vt)*dt/h
            vf[i][j] += (Fb + Ft) + cc[i][j]*dt*(vt-vb)/h
            # if (vf[i][j] < 0 or vf[i][j] > 1):
            #     print(vf[i][j])
            fluxes[N*(N+1)+j+(N+1)*i] = Fb
            fluxes[N*(N+1)+(j+1)+(N+1)*i] = -Ft

            # vf[i][j] += Fl + Fr + Fb + Ft
    
    def cost_boundedness(fluxes,N,old_fluxes):
        E = 0
        for i in range(2*(N+1)*N):
            E += (fluxes[i] - old_fluxes[i])**2
        return E
    
    old_fluxes = np.copy(fluxes)

    A = np.zeros((N*N+2*(N+1)*N,2*(N+1)*N))
    lb = np.zeros(N*N+2*(N+1)*N)
    ub = np.zeros(N*N+2*(N+1)*N)
    for i in range(N):
        for j in range(N):
            A[i+j*N][i+(N+1)*j] = 1
            A[i+j*N][i+1+(N+1)*j] = -1
            A[i+j*N][N*(N+1)+j+(N+1)*i] = 1
            A[i+j*N][N*(N+1)+j+1+(N+1)*i] = -1
            lb[i+j*N] = -vf[i][j]
            ub[i+j*N] = 1.0-vf[i][j]
    for i in range(2*(N+1)*N):
        A[N*N+i][i] = 1
        lb[N*N+i] = -1e15
        ub[N*N+i] = 1e15            
    for j in range(N):
        lb[N*N+0+(N+1)*j] = -1.0e-16
        ub[N*N+0+(N+1)*j] = 1.0e-16
        lb[N*N+N+(N+1)*j] = -1.0e-16
        ub[N*N+N+(N+1)*j] = 1.0e-16
    for i in range(N):
        lb[N*N+N*(N+1)+0+(N+1)*i] = -1.0e-16
        ub[N*N+N*(N+1)+0+(N+1)*i] = 1.0e-16
        lb[N*N+N*(N+1)+N+(N+1)*i] = -1.0e-16
        ub[N*N+N*(N+1)+N+(N+1)*i] = 1.0e-16

    def callback_flux_0_bdy(fluxes):
        for j in range(N):
            fluxes[0+(N+1)*j] = 0
            fluxes[N+(N+1)*j] = 0
        for i in range(N):
            fluxes[N*(N+1)+0+(N+1)*i] = 0
            fluxes[N*(N+1)+N+(N+1)*i] = 0
    
    boundedness_constraint = LinearConstraint(A,lb,ub)

    # res = minimize(cost_boundedness,fluxes,args=(N,old_fluxes),constraints=[boundedness_constraint],callback=callback_flux_0_bdy,options={'maxiter':100},tol=VFTOL)
    # fluxes = res.x
    # print(res.nit)

    # for i in range(N):
    #     for j in range(N):
    #         vf[i][j] += fluxes[i+(N+1)*j] - fluxes[(i+1)+(N+1)*j] + fluxes[N*(N+1)+j+(N+1)*i] - fluxes[N*(N+1)+(j+1)+(N+1)*i]

    return vf, norm, dist

def cost_hlvira(hyperplane,N,h,old_vf,guess_vf,i,j,t,dt):
    x0 = i / N; x1 = x0 + h
    y0 = j / N; y1 = y0 + h
    E = 0.0
    # Old time
    plane = np.zeros(3)
    plane[0] = hyperplane[0]
    plane[1] = hyperplane[1]
    plane[2] = hyperplane[3]-t*hyperplane[2]
    for l in range(-1,2):
        for n in range(-1,2):
            x00 = x0+l*h; y00 = y0+n*h
            tmp_vf = forward(x00,y00,h,h,plane[0],plane[1],plane[2])
            E += (tmp_vf - old_vf[i+l][j+n])**2
    # New time
    plane[0] = hyperplane[0]
    plane[1] = hyperplane[1]
    plane[2] = hyperplane[3]-(t+dt)*hyperplane[2]
    for l in range(-1,2):
        for n in range(-1,2):
            x00 = x0+l*h; y00 = y0+n*h
            tmp_vf = forward(x00,y00,h,h,plane[0],plane[1],plane[2])
            E += (tmp_vf - guess_vf[i+l][j+n])**2
    return E

def stic_advect(N,h,guess_norm,guess_dist,old_vf,guess_vf,t,dt):
    hyperplane = np.zeros((N,N,4))
    vf = np.copy(old_vf)
    for i in range(N):
        for j in range(N):
            hyperplane[i][j][0] = guess_norm[i][j][0]
            hyperplane[i][j][1] = guess_norm[i][j][1]
            hyperplane[i][j][2] = 0.0
            hyperplane[i][j][3] = guess_dist[i][j]

    for i in range(1,N-1):
        for j in range(1,N-1):
            x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
            y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)

            if (is_interface(old_vf[i][j]) or is_interface(guess_vf[i][j])):
                # We need to fit!
                hyperplane_guess = hyperplane[i][j]

                def M0_t0_constraint_hlvira(hyperplane):
                    plane = np.zeros(3)
                    plane[0] = hyperplane[0]
                    plane[1] = hyperplane[1]
                    plane[2] = hyperplane[3]-t*hyperplane[2]
                    tmp_vf = forward_hlvira(x0,y0,h,h,plane[0],plane[1],plane[2])
                    # if (tmp_vf == 0 and old_vf[i][j] > 0):
                    #     alpha = abs(plane[2] - (plane[0]*x+plane[1]*y))/h - 0.5;
                    #     return -alpha
                    # elif (tmp_vf == 1  and old_vf[i][j] < 1):
                    #     alpha = abs(plane[2] - (plane[0]*x+plane[1]*y))/h - 0.5;
                    #     return alpha
                    # else:
                    #     return tmp_vf - old_vf[i][j]
                    if (tmp_vf <= 0 and old_vf[i][j] < VFTOL):
                        return 0
                    elif (tmp_vf >= 1.0 and old_vf[i][j] > 1.0 - VFTOL):
                        return 0
                    else:
                        return tmp_vf - old_vf[i][j]
                def M0_t1_constraint_hlvira(hyperplane):
                    plane = np.zeros(3)
                    plane[0] = hyperplane[0]
                    plane[1] = hyperplane[1]
                    plane[2] = hyperplane[3]-(t+dt)*hyperplane[2]
                    tmp_vf = forward_hlvira(x0,y0,h,h,plane[0],plane[1],plane[2])
                    # if (tmp_vf == 0 and guess_vf[i][j] > 0):
                    #     alpha = abs(plane[2] - (plane[0]*x+plane[1]*y))/h - 0.5;
                    #     return -alpha
                    # elif (tmp_vf == 1  and guess_vf[i][j] < 1):
                    #     alpha = abs(plane[2] - (plane[0]*x+plane[1]*y))/h - 0.5;
                    #     return alpha
                    # else:
                    #     return tmp_vf - guess_vf[i][j]
                    if (tmp_vf <= 0 and guess_vf[i][j] < VFTOL):
                        return 0
                    elif (tmp_vf >= 1.0 and guess_vf[i][j] > 1.0 - VFTOL):
                        return 0
                    else:
                        return tmp_vf - guess_vf[i][j]
                def M0_3d_constraint_hlvira(hyperplane):
                    tmp_vf = forward3d(x0,y0,t,h,h,dt,hyperplane[0],hyperplane[1],hyperplane[2],hyperplane[3])
                    return tmp_vf - 0.5*(old_vf[i][j]+guess_vf[i][j])

                def callback_hlvira(hyperplane):
                    hyperplane[3] = backward3d(x0,y0,t,h,h,dt,hyperplane[0],hyperplane[1],hyperplane[2],0.5*(old_vf[i][j]+guess_vf[i][j]))

                M0_t0_constraint = NonlinearConstraint(M0_t0_constraint_hlvira, -VFTOL, +VFTOL)
                M0_t1_constraint = NonlinearConstraint(M0_t1_constraint_hlvira, -VFTOL, +VFTOL)
                M0_3d_constraint = NonlinearConstraint(M0_3d_constraint_hlvira, -VFTOL, +VFTOL)
                res = minimize(cost_hlvira,hyperplane_guess,args=(N,h,old_vf,guess_vf,i,j,t,dt),callback=callback_hlvira,options={'maxiter':100},tol=VFTOL)
                hyperplane[i][j] = res.x
                hyperplane[i][j][3] = backward3d(x0,y0,t,h,h,dt,hyperplane[i][j][0],hyperplane[i][j][1],hyperplane[i][j][2],0.5*(old_vf[i][j]+guess_vf[i][j]))
                # print(forward(x0,y0,h,h,hyperplane[i][j][0],hyperplane[i][j][1],hyperplane[i][j][3]-(t)*hyperplane[i][j][2]) - old_vf[i][j])
                # print(forward(x0,y0,h,h,hyperplane[i][j][0],hyperplane[i][j][1],hyperplane[i][j][3]-(t+dt)*hyperplane[i][j][2]) - guess_vf[i][j])
                # print(forward3d(x0,y0,t,h,h,dt,hyperplane[i][j][0],hyperplane[i][j][1],hyperplane[i][j][2],hyperplane[i][j][3]) - 0.5*(old_vf[i][j]+guess_vf[i][j]))

    fluxes = np.zeros(2*(N+1)*N)
    for i in range(1,N-1):
        for j in range(1,N-1):
            x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
            y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)
            ul = U(x0,y,t)
            ur = U(x1,y,t)
            vb = V(x,y0,t)
            vt = V(x,y1,t)
            Fl=0; Fr=0; Fb=0; Ft=0
            # UPWIND
            if (ul > 0):
                Fl += forward(y0,t,h,dt,hyperplane[i-1][j][1],hyperplane[i-1][j][2],hyperplane[i-1][j][3]-x0*hyperplane[i-1][j][0])*abs(ul)*dt/h
            else:
                Fl -= forward(y0,t,h,dt,hyperplane[i][j][1],hyperplane[i][j][2],hyperplane[i][j][3]-x0*hyperplane[i][j][0])*abs(ul)*dt/h
            if (ur > 0):
                Fr -= forward(y0,t,h,dt,hyperplane[i][j][1],hyperplane[i][j][2],hyperplane[i][j][3]-x1*hyperplane[i][j][0])*abs(ur)*dt/h
            else:
                Fr += forward(y0,t,h,dt,hyperplane[i+1][j][1],hyperplane[i+1][j][2],hyperplane[i+1][j][3]-x1*hyperplane[i+1][j][0])*abs(ur)*dt/h

            if (vb > 0):
                Fb += forward(t,x0,dt,h,hyperplane[i][j-1][2],hyperplane[i][j-1][0],hyperplane[i][j-1][3]-y0*hyperplane[i][j-1][1])*abs(vb)*dt/h
            else:
                Fb -= forward(t,x0,dt,h,hyperplane[i][j][2],hyperplane[i][j][0],hyperplane[i][j][3]-y0*hyperplane[i][j][1])*abs(vb)*dt/h
            if (vt > 0): 
                Ft -= forward(t,x0,dt,h,hyperplane[i][j][2],hyperplane[i][j][0],hyperplane[i][j][3]-y1*hyperplane[i][j][1])*abs(vt)*dt/h
            else:
                Ft += forward(t,x0,dt,h,hyperplane[i][j+1][2],hyperplane[i][j+1][0],hyperplane[i][j+1][3]-y1*hyperplane[i][j+1][1])*abs(vt)*dt/h

            # # CENTRAL
            # Fl += 0.5*forward(y0,t,h,dt,hyperplane[i-1][j][1],hyperplane[i-1][j][2],hyperplane[i-1][j][3]-x0*hyperplane[i-1][j][0])*ul*dt/h
            # Fl += 0.5*forward(y0,t,h,dt,hyperplane[i][j][1],hyperplane[i][j][2],hyperplane[i][j][3]-x0*hyperplane[i][j][0])*ul*dt/h
            # Fr -= 0.5*forward(y0,t,h,dt,hyperplane[i][j][1],hyperplane[i][j][2],hyperplane[i][j][3]-x1*hyperplane[i][j][0])*ur*dt/h
            # Fr -= 0.5*forward(y0,t,h,dt,hyperplane[i+1][j][1],hyperplane[i+1][j][2],hyperplane[i+1][j][3]-x1*hyperplane[i+1][j][0])*ur*dt/h
            # Fb += 0.5*forward(t,x0,dt,h,hyperplane[i][j-1][2],hyperplane[i][j-1][0],hyperplane[i][j-1][3]-y0*hyperplane[i][j-1][1])*vb*dt/h
            # Fb += 0.5*forward(t,x0,dt,h,hyperplane[i][j][2],hyperplane[i][j][0],hyperplane[i][j][3]-y0*hyperplane[i][j][1])*vb*dt/h
            # Ft -= 0.5*forward(t,x0,dt,h,hyperplane[i][j][2],hyperplane[i][j][0],hyperplane[i][j][3]-y1*hyperplane[i][j][1])*vt*dt/h
            # Fr -= 0.5*forward(t,x0,dt,h,hyperplane[i][j+1][2],hyperplane[i][j+1][0],hyperplane[i][j+1][3]-y1*hyperplane[i][j+1][1])*vt*dt/h
            
            fluxes[i+(N+1)*j] = Fl
            fluxes[(i+1)+(N+1)*j] = -Fr
            fluxes[N*(N+1)+j+(N+1)*i] = Fb
            fluxes[N*(N+1)+(j+1)+(N+1)*i] = -Ft

            # vf[i][j] += Fl + Fr + Fb + Ft
    
    def cost_boundedness(fluxes,N,old_fluxes):
        E = 0
        for i in range(2*(N+1)*N):
            E += (fluxes[i] - old_fluxes[i])**2
        return E
    
    old_fluxes = np.copy(fluxes)

    A = np.zeros((N*N+2*(N+1)*N,2*(N+1)*N))
    lb = np.zeros(N*N+2*(N+1)*N)
    ub = np.zeros(N*N+2*(N+1)*N)
    for i in range(N):
        for j in range(N):
            A[i+j*N][i+(N+1)*j] = 1
            A[i+j*N][i+1+(N+1)*j] = -1
            A[i+j*N][N*(N+1)+j+(N+1)*i] = 1
            A[i+j*N][N*(N+1)+j+1+(N+1)*i] = -1
            lb[i+j*N] = -vf[i][j]
            ub[i+j*N] = 1.0-vf[i][j]
    for i in range(2*(N+1)*N):
        A[N*N+i][i] = 1
        lb[N*N+i] = -1e15
        ub[N*N+i] = 1e15            
    for j in range(N):
        lb[N*N+0+(N+1)*j] = -1.0e-16
        ub[N*N+0+(N+1)*j] = 1.0e-16
        lb[N*N+N+(N+1)*j] = -1.0e-16
        ub[N*N+N+(N+1)*j] = 1.0e-16
    for i in range(N):
        lb[N*N+N*(N+1)+0+(N+1)*i] = -1.0e-16
        ub[N*N+N*(N+1)+0+(N+1)*i] = 1.0e-16
        lb[N*N+N*(N+1)+N+(N+1)*i] = -1.0e-16
        ub[N*N+N*(N+1)+N+(N+1)*i] = 1.0e-16

    def callback_flux_0_bdy(fluxes):
        for j in range(N):
            fluxes[0+(N+1)*j] = 0
            fluxes[N+(N+1)*j] = 0
        for i in range(N):
            fluxes[N*(N+1)+0+(N+1)*i] = 0
            fluxes[N*(N+1)+N+(N+1)*i] = 0
    
    boundedness_constraint = LinearConstraint(A,lb,ub)

    res = minimize(cost_boundedness,fluxes,args=(N,old_fluxes),constraints=[boundedness_constraint],callback=callback_flux_0_bdy,options={'maxiter':100},tol=VFTOL)
    fluxes = res.x
    print(res.nit)

    for i in range(N):
        for j in range(N):
            vf[i][j] += fluxes[i+(N+1)*j] - fluxes[(i+1)+(N+1)*j] + fluxes[N*(N+1)+j+(N+1)*i] - fluxes[N*(N+1)+(j+1)+(N+1)*i]

    return vf

def is_interface(vf):
    if (vf > VFTOL or vf < 1-VFTOL):
        return True
    else:
        return False

def cost_hlvira2(hyperplane,N,h,old_vf,guess_vf,i,j,t,dt,dir):
    x0 = i / N; x1 = x0 + h
    y0 = j / N; y1 = y0 + h
    E = 0.0
    # Old time
    plane = np.zeros(3)
    plane[0] = hyperplane[0]
    plane[1] = hyperplane[1]
    plane[2] = hyperplane[3]-t*hyperplane[2]
    for l in range(-1,1+dir):
        for n in range(-1,2-dir):
            x00 = x0+l*h; y00 = y0+n*h
            tmp_vf = forward(x00,y00,h,h,plane[0],plane[1],plane[2])
            E += (tmp_vf - old_vf[i+l][j+n])**2
    # New time
    plane[0] = hyperplane[0]
    plane[1] = hyperplane[1]
    plane[2] = hyperplane[3]-(t+dt)*hyperplane[2]
    for l in range(-1,1+dir):
        for n in range(-1,2-dir):
            x00 = x0+l*h; y00 = y0+n*h
            tmp_vf = forward(x00,y00,h,h,plane[0],plane[1],plane[2])
            E += (tmp_vf - guess_vf[i+l][j+n])**2
    return E

def stic_advect2(N,h,guess_norm,guess_dist,old_vf,guess_vf,t,dt):
    hyperplane_x = np.zeros((N+1,N,4))
    hyperplane_y = np.zeros((N,N+1,4))
    vf = np.copy(guess_vf)
    for i in range(1,N):
        for j in range(1,N-1):
            x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
            y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)
            hyperplane_x[i][j][0] = 0.5*(guess_norm[i-1][j][0]+guess_norm[i][j][0])
            hyperplane_x[i][j][1] = 0.5*(guess_norm[i-1][j][1]+guess_norm[i][j][1])
            hyperplane_x[i][j][2] = 0.0
            hyperplane_x[i][j][3] = 0.5*(guess_dist[i-1][j]+guess_dist[i][j])

            if (is_interface(old_vf[i-1][j]) or is_interface(old_vf[i][j]) or is_interface(guess_vf[i-1][j]) or is_interface(guess_vf[i][j])):
                # We need to fit!
                hyperplane_guess = hyperplane_x[i][j]
                res = minimize(cost_hlvira2,hyperplane_guess,args=(N,h,old_vf,guess_vf,i,j,t,dt,0),options={'maxiter':100},tol=VFTOL)
                hyperplane_x[i][j] = res.x
    for i in range(1,N-1):
        for j in range(1,N):
            x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
            y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)
            hyperplane_y[i][j][0] = 0.5*(guess_norm[i][j-1][0]+guess_norm[i][j][0])
            hyperplane_y[i][j][1] = 0.5*(guess_norm[i][j-1][1]+guess_norm[i][j][1])
            hyperplane_y[i][j][2] = 0.0
            hyperplane_y[i][j][3] = 0.5*(guess_dist[i][j-1]+guess_dist[i][j])

            if (is_interface(old_vf[i-1][j]) or is_interface(old_vf[i][j]) or is_interface(guess_vf[i-1][j]) or is_interface(guess_vf[i][j])):
                # We need to fit!
                hyperplane_guess = hyperplane_x[i][j]
                res = minimize(cost_hlvira2,hyperplane_guess,args=(N,h,old_vf,guess_vf,i,j,t,dt,1),options={'maxiter':100},tol=VFTOL)
                hyperplane_y[i][j] = res.x

    for i in range(1,N-1):
        for j in range(1,N-1):
            x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
            y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)
            ul = U(x0,y,t); ur = U(x1,y,t); vb = V(x,y0,t); vt = V(x,y1,t)
            Fl = forward(t,y0,dt,h,hyperplane_x[i][j][2],hyperplane_x[i][j][1],hyperplane_x[i][j][3]-x0*hyperplane_x[i][j][0])*ul*dt/h
            Fr = -forward(t,y0,dt,h,hyperplane_x[i+1][j][2],hyperplane_x[i+1][j][1],hyperplane_x[i+1][j][3]-x1*hyperplane_x[i+1][j][0])*ur*dt/h
            Fb = forward(t,x0,dt,h,hyperplane_y[i][j][2],hyperplane_y[i][j][1],hyperplane_y[i][j][3]-y0*hyperplane_y[i][j][0])*vb*dt/h
            Ft = -forward(t,x0,dt,h,hyperplane_y[i][j+1][2],hyperplane_y[i][j+1][1],hyperplane_y[i][j+1][3]-y1*hyperplane_y[i][j+1][0])*vt*dt/h
            vf[i][j] += Fl + Fr + Fb + Ft

    return vf

N = 25
h = 1/N
dt = 0.01
Nt = int(1/dt)

cx = 0.5; cy = 0.65; r = 0.15

vf = np.zeros((N, N))
norm = np.zeros((N, N, 2))
dist = np.zeros((N, N))

# initialisation of volume fractions
total_volume = 0.0
for i in range(N):
    for j in range(N):
        x0 = i / N; x1 = x0 + h; x = 0.5*(x0+x1)
        y0 = j / N; y1 = y0 + h; y = 0.5*(y0+y1)
        norm[i][j][0] = x-cx
        norm[i][j][1] = y-cy
        mag = sqrt(norm[i][j][0]**2 + norm[i][j][1]**2)
        norm[i][j][0] /= mag; norm[i][j][1] /= mag
        ptx = cx+norm[i][j][0]*r
        pty = cy+norm[i][j][1]*r
        dist[i][j] = norm[i][j][0]*ptx + norm[i][j][1]*pty
        vf[i][j] = forward(x0,y0,h,h,norm[i][j][0],norm[i][j][1],dist[i][j])
        dist[i][j] = backward(x0,y0,h,norm[i][j][0],norm[i][j][1],vf[i][j])
        total_volume += vf[i][j]


print("Time\t\t Mass error\tBound error")

filenames = []
for it in range(Nt):
    t = (it+0.5) * dt

    old_vf = np.copy(vf)
    guess_vf, norm, dist = split_advect(N,h,old_vf,t,dt)
    # vf = guess_vf
    vf = stic_advect(N,h,norm,dist,old_vf,guess_vf,t,dt)
    # guess_vf = np.copy(vf)
    # vf = stic_advect(N,h,norm,dist,old_vf,guess_vf,t,dt)
    # vf = stic_advect2(N,h,norm,dist,old_vf,guess_vf,t,dt)

    tmp_volume = 0.0
    bound_error = 0.0
    for i in range(N):
        for j in range(N):
            tmp_volume += vf[i][j]
            bound_error = max(bound_error, -vf[i][j])
            bound_error = max(bound_error, vf[i][j]-1.0)
            # vf[i][j] = max(0,vf[i][j])
            # vf[i][j] = min(1,vf[i][j])

    print("%2.6f \t %2.2e \t %2.2e" % ((it+1)*dt,abs(tmp_volume-total_volume),bound_error))

    figure, axes = plt.subplots()
    filename = f'frame_{it}.png'
    filenames.append(filename)
     # create discrete colormap
    cmap = mpl.colormaps['viridis']
    fig, ax = plt.subplots()
    plotvf = np.swapaxes(vf,0,1)
    ax.imshow(plotvf, cmap=cmap, origin='lower')
    plot_plic(N,h,vf,norm,dist,ax)
    # draw gridlines
    # ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
    plt.tick_params(labelleft=False, labelbottom=False)
    ax.set_xticks(np.arange(-.5, N, 1));
    ax.set_yticks(np.arange(-.5, N, 1));
    plt.savefig(filename)
    plt.rcParams.update({'figure.max_open_warning': 0})
    plt.close()

with imageio.get_writer('stic.gif', mode='I',duration=5) as writer:
    for filename in filenames:
        image = imageio.v2.imread(filename)
        writer.append_data(image)
