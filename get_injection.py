#!/usr/bin/env python

import argparse
import numpy as np
from pylal import SimInspiralUtils
import lal
import lalsimulation as lalsim
from math import cos,ceil,floor,sqrt,pi as pi_constant

parser = argparse.ArgumentParser()
parser.add_argument("--inj", dest="inj", default=None, help="Injection XML Path")
parser.add_argument("--event", dest ="event", type=int, default=None, help="Event Number")
parser.add_argument("--frame", dest ="frame", default="OrbitalL", help="Spin Frame of injection")
args = parser.parse_args()

injname = args.inj
event = args.event
frame = args.frame

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

def orbital_momentum(fref, m1,m2,eta, inclination):
    """
    Calculate orbital angular momentum vector.
    Note: The units of Lmag are different than what used in lalsimulation.
    Mc must be called in units of Msun here.

    Note that if one wants to build J=L+S1+S2 with L returned by this function, S1 and S2
    must not get the Msun^2 factor.
    """
    Lmag = orbital_momentum_mag(fref, m1,m2,eta)
    Lx, Ly, Lz = sph2cart(Lmag, inclination, 0.0)
    return np.hstack((Lx,Ly,Lz))
#
#
def orbital_momentum_mag(fref, m1,m2,eta):
    v0 = np.power(pi_constant * lal.MTSUN_SI * fref, 1.0/3.0)
    #1 PN Mtot*Mtot*eta/v 
    PNFirst = (((m1+m2)**2)*eta)/v0
    PNSecond = 1+ (v0**2) * (3.0/2.0 -eta/6.0)
    Lmag= PNFirst*PNSecond
    return Lmag

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

def angle(x, y):

    norm1 = np.dot(x, x)
    norm2 = np.dot(y, y)
    ph = np.arccos(np.dot(x, y)/np.sqrt(norm1*norm2))

    return (ph)

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

def calculate_injected_sys_frame_params(sim_inspiral_event, f_ref = 20.0):

    inj = sim_inspiral_event

    axis = lalsim.SimInspiralGetFrameAxisFromString(frame)

    m1, m2 = inj.mass1, inj.mass2
    mc, eta = inj.mchirp, inj.eta

    # Convert to radiation frame
    iota, s1x, s1y, s1z, s2x, s2y, s2z = \
                lalsim.SimInspiralInitialConditionsPrecessingApproxs(inj.inclination,
                                                                     inj.spin1x, inj.spin1y, inj.spin1z,
                                                                     inj.spin2x, inj.spin2y, inj.spin2z,
                                                                     m1*lal.MSUN_SI, m2 * lal.MSUN_SI, f_ref, axis)

    # when we come out of SimInspiralInitialConditionsPrecessingApproxs 
    # (which is called before precessing waveforms such as 
    # SpinTaylorT* and IMRPv2 waveform drivers 
    # are called in XLALChooseTDWaveform) we are in the 
    # frame where N= (0,0,1) and L= (sin(inc),0,cos(inc))
    # In additon L mag is defined to the 1.5PN order

    a1, theta1, phi1 = cart2sph(s1x, s1y, s1z)
    a2, theta2, phi2 = cart2sph(s2x, s2y, s2z)

    L  = orbital_momentum(f_ref, m1,m2,eta,iota)
    S1 = np.hstack((s1x, s1y, s1z))
    S2 = np.hstack((s2x, s2y, s2z))

    S1 *= m1**2
    S2 *= m2**2
    J = L + S1 + S2

    tilt1 = array_ang_sep(L, S1)
    tilt2 = array_ang_sep(L, S2)
    beta  = array_ang_sep(J, L)

    # Need to do rotations of XLALSimInspiralTransformPrecessingInitialConditioin inverse order to go in the L frame
            # first rotation: bring J in the N-x plane, with negative x component
    phi0 = np.arctan2(J[1], J[0])
    phi0 = np.pi - phi0

    J = ROTATEZ(phi0, J[0], J[1], J[2])
    L = ROTATEZ(phi0, L[0], L[1], L[2])
    S1 = ROTATEZ(phi0, S1[0], S1[1], S1[2])
    S2 = ROTATEZ(phi0, S2[0], S2[1], S2[2])

    # now J in in the N-x plane and form an angle theta_jn with N, rotate by -theta_jn around y to have J along z
    theta_jn = array_polar_ang(J)

    J = ROTATEY(theta_jn, J[0], J[1], J[2])
    L = ROTATEY(theta_jn, L[0], L[1], L[2])
    S1 = ROTATEY(theta_jn, S1[0], S1[1], S1[2])
    S2 = ROTATEY(theta_jn, S2[0], S2[1], S2[2])

    # J should now be || z and L should have a azimuthal angle phi_jl
    phi_jl = np.arctan2(L[1], L[0])
    phi_jl = np.pi - phi_jl

    # bring L in the Z-X plane, with negative x
    J = ROTATEZ(phi_jl, J[0], J[1], J[2])
    L = ROTATEZ(phi_jl, L[0], L[1], L[2])
    S1 = ROTATEZ(phi_jl, S1[0], S1[1], S1[2])
    S2 = ROTATEZ(phi_jl, S2[0], S2[1], S2[2])

    theta0 = array_polar_ang(L)
    J = ROTATEY(theta0, J[0], J[1], J[2])
    L = ROTATEY(theta0, L[0], L[1], L[2])
    S1 = ROTATEY(theta0, S1[0], S1[1], S1[2])
    S2 = ROTATEY(theta0, S2[0], S2[1], S2[2])

    # The last rotation is useless as it won't change the differenze in spins' azimuthal angles
    phi1 = np.arctan2(S1[1], S1[0])
    phi2 = np.arctan2(S2[1], S2[0])
    if phi2 < phi1:
        phi12 = phi2 - phi1 + 2.*np.pi
    else:
        phi12 = phi2 - phi1

    return a1, a2, s1z, s2z, theta_jn, phi_jl, tilt1, tilt2, phi12

###### DO STUFF

injection_table = SimInspiralUtils.ReadSimInspiralFromFiles([injname])
injection_object = injection_table[event]
injection_values = extract_inj_vals(injection_object)

for parameter in injection_values:
	print('injected %r: %r' % (parameter, injection_values[parameter]))

