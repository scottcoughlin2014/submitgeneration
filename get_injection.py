#!/usr/bin/env python

import argparse
import numpy as np
from pylal import SimInspiralUtils
import lal
import lalsimulation as lalsim

def sph2cart(r,theta,phi):
    """
    Utiltiy function to convert r,theta,phi to cartesian co-ordinates.
    """
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x,y,z

def cart2sph(x,y,z):
    """
    Utility function to convert cartesian coords to r,theta,phi.
    """
    r = np.sqrt(x*x + y*y + z*z)
    theta = np.arccos(z/r)
    phi = np.fmod(2*np.pi + np.arctan2(y,x), 2*np.pi)

    return r,theta,phi

def array_dot(vec1, vec2):
    """
    Calculate dot products between vectors in rows of numpy arrays.
    """
    if vec1.ndim==1:
        product = (vec1*vec2).sum()
    else:
        product = (vec1*vec2).sum(axis=1).reshape(-1,1)
    return product



def array_ang_sep(vec1, vec2):
    """
    Find angles between vectors in rows of numpy arrays.
    """
    vec1_mag = np.sqrt(array_dot(vec1, vec1))
    vec2_mag = np.sqrt(array_dot(vec2, vec2))
    return np.arccos(array_dot(vec1, vec2)/(vec1_mag*vec2_mag))

def array_polar_ang(vec):
    """
    Find polar angles of vectors in rows of a numpy array.
    """
    if vec.ndim==1:
        z = vec[2]
    else:
        z = vec[:,2].reshape(-1,1)
    norm = np.sqrt(array_dot(vec,vec))
    return np.arccos(z/norm)

#}}}

def orbital_momentum(f_ref, mc, inclination, m1,m2,eta):
    #{{{
    Lmag = np.power(mc, 5.0/3.0) / np.power(np.pi * lal.MTSUN_SI * f_ref, 1.0/3.0)
    v0 = ((m1+m2)*lal.MTSUN_SI * np.pi *f_ref)**(1./3.)
    Lmag= Lmag*(1.0 + (v0**2) *  (2.5 -eta/6.) )

    Lx, Ly, Lz = sph2cart(Lmag, inclination, 0.0)
    return np.hstack((Lx,Ly,Lz))

def ROTATEZ(angle, vx, vy, vz):
    # This is the ROTATEZ in LALSimInspiral.c.
    tmp1 = vx*np.cos(angle) - vy*np.sin(angle);
    tmp2 = vx*np.sin(angle) + vy*np.cos(angle);
    return np.asarray([tmp1,tmp2,vz])

def ROTATEY(angle, vx, vy, vz):
    # This is the ROTATEY in LALSimInspiral.c
    tmp1 = vx*np.cos(angle) + vz*np.sin(angle);
    tmp2 = -1.0*vx*np.sin(angle) + vz*np.cos(angle);
    return np.asarray([tmp1,vy,tmp2])

def extract_inj_vals(sim_inspiral_event):
    a1, a2, spin1z, spin2z, theta_jn, phi_jl, tilt1, tilt2, phi12 = calculate_injected_sys_frame_params(sim_inspiral_event)
    injvals={
        'mc'          : sim_inspiral_event.mchirp,
	'q'           : sim_inspiral_event.mass2/sim_inspiral_event.mass1,
	'time'        : float(sim_inspiral_event.get_end()),
	'phi_orb'     : sim_inspiral_event.coa_phase,
        'dist'        : sim_inspiral_event.distance,
        'logdistance' : np.log(sim_inspiral_event.distance),
        'ra'          : sim_inspiral_event.longitude,
	'dec'         : sim_inspiral_event.latitude,
        'sindec'      : np.sin(sim_inspiral_event.latitude),
        'psi'         : np.mod(sim_inspiral_event.polarization, np.pi),
        'a1'          : a1,
        'a2'          : a2,
        'spin1'       : spin1z,
        'spin2'       : spin2z,
        'phi12'       : phi12,
        'tilt1'       : tilt1,
        'tilt2'       : tilt2,
        'costilt1'    : np.cos(tilt1),
        'costilt2'    : np.cos(tilt2),
        'theta_jn'    : theta_jn,
        'costheta_jn' : np.cos(theta_jn),
        'phi12'       : phi12,
        'phi_jl'      : phi_jl}
    return injvals

def calculate_injected_sys_frame_params(sim_inspiral_event, f_ref = 100.0):

    # To do this we should extract the parameters of
    # our binary system
    m1  = inj.mass1
    m2  = inj.mass2
    m1_sun = m1*lal.MSUN_SI
    m2_sun = m2*lal.MSUN_SI
    mc  = inj.mchirp
    eta = inj.eta
    s1x = inj.spin1x
    s1y = inj.spin1y
    s1z = inj.spin1z
    s2x = inj.spin2x
    s2y = inj.spin2y
    s2z = inj.spin2z
    inc = inj.inclination
    S1 = np.hstack((s1x, s1y, s1z))
    S2 = np.hstack((s2x, s2y, s2z))

    # Get Chi/a1
    a1, theta1,phi1 = cart2sph(s1x, s1y, s1z)
    a2, theta2,phi2 = cart2sph(s2x, s2y, s2z)

    # Spin Angles are retrieved through careful analysis of the following functions: SimInspiralTransformPrecessingNewInitialConditions

    # One thing we need first is to calculate Lmag the 
    # magnitude orbital Angular momentum which we only need 
    # eta m1 m2 mchirp to do.

    Lmag = orbital_momentum_leadord(f_ref, mc, m1,m2,eta)

    # Starting frame: LNhat is along the z-axis and the unit
    # spin vectors are defined from the angles relative to LNhat.
    # Note that we put s1hat in the x-z plane, and phi12
    # sets the azimuthal angle of s2hat measured from the x-axis.

    LNhat = np.array([0., 0., 1.])
    L = Lmag * LNhat
    Nhat = np.array([np.sin(inc), 0., np.cos(inc)])

    # Define S1, S2, J with proper magnitudes
    S1 *= m1**2
    S2 *= m2**2
    J = L + S1 + S2

    # Let us also calculate Jhat for the plot
    Jnorm = np.sqrt( J[0]**2 + J[1]**2 + J[2]**2)
    Jhat = J/Jnorm

    # alright let's plot this shit!
    #plot_M1M2_J_L_N(m1,m2,LNhat,Nhat,Jhat,inc)
    
    # Calculate some angles while in this frame before we need to rotate
    phi1 = np.arctan2(S1[1], S1[0])
    phi2 = np.arctan2(S2[1], S2[0])
    phi12 = phi2 - phi1

    tilt1 = angle(L, S1)
    tilt2 = angle(L, S2)
    theta_jn = angle(J, Nhat)
    beta  = angle(J, L)

    if (phi2 < phi1):
          phi12 = phi2 - phi1 + 2.*np.pi

    # we gotta get phi_jl
    # Rotation 1: Rotate about z-axis by -phi0 to put Jhat in x-z plane
    phi0 = np.arctan2(J[1], J[0])
    phi0 = np.pi - phi0
    theta0 = np.arccos(Jhat[2])

    J = ROTATEZ(phi0, J[0], J[1], J[2])
    L = ROTATEZ(phi0, L[0], L[1], L[2])
    Nhat = ROTATEZ(phi0, Nhat[0], Nhat[1], Nhat[2])

    # Rotation 2: Rotate about new y-axis by -theta to put Jhat along z-axis

    J = ROTATEY(theta0, J[0], J[1], J[2])
    L  = ROTATEY(theta0, L[0], L[1], L[2])
    Nhat = ROTATEY(theta0, Nhat[0], Nhat[1], Nhat[2])

    # Rotation 6: Now L is along z and we have to bring N in the x-z plane.
    phiN = np.arctan2(Nhat[1], Nhat[0])
    phiN = np.pi - phiN
    J = ROTATEZ(phiN, J[0], J[1], J[2])
    L = ROTATEZ(phiN, L[0], L[1], L[2])
    Nhat = ROTATEZ(phiN, Nhat[0], Nhat[1], Nhat[2])

    # Now calc phi_jl
    phi_jl = np.arctan2(L[1], L[0])
    if (phi_jl < 0):

       phi_jl[i] = phi_jl + 2.0*pi

    spins['phi_jl'] = phi_jl

    return a1, a2, spin1z, spin2z, theta_jn, phi_jl, tilt1, tilt2, phi12

###### DO STUFF

injection_table = SimInspiralUtils.ReadSimInspiralFromFiles([injname])
injection_object = injection_table[event]
injection_values = extract_inj_vals(injection_object)

for parameter in injection_values:
	print('injected %r: %r' % (parameter, injection_values[parameter]))

