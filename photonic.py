# from __future__ import division, print_function, absolute_import
# import numpy as np
# import itertools
# from numpy import sqrt, pi, cos, sin, log, exp, sinh
# from scipy.special import iv as besseli
# from scipy.optimize import fmin, fminbound
# from scipy import integrate
# # from scipy.interpolate import interp1d

# import gdspy
# from phidl.device_layout import Device, Port
# from phidl.device_layout import _parse_layer, DeviceReference
# import phidl.routing as pr
# import copy as python_copy
# from collections import OrderedDict
# import pickle

# from skimage import draw, morphology

from phidl.geometry import *

# waveguide: 1st para is length, second is width
def waveguide(width = 1, length = 10, layer = 0):
    WG = Device('waveguide')
    WG.add_polygon( [(0, -width/2), (0, width/2), (length, width/2), (length, -width/2)], layer = layer)
    WG.add_port(name = 1, midpoint = [0,0], width = width, orientation = 180)
    WG.add_port(name = 2, midpoint = [length,0], width = width, orientation = 0)
    return WG

# new ring as an alternative of ring
def new_ring(radius = 10, width = 0.5, angle_resolution = 1, layer = 0):
    D = Device(name = 'new_ring')
    inner_radius = radius - width/2
    outer_radius = radius + width/2
    angle = np.linspace(0, 360, np.ceil(360/angle_resolution))
    angle.sort()
    t=angle*np.pi/180
    inner_points_x = (inner_radius*cos(t)).tolist()
    inner_points_y = (inner_radius*sin(t)).tolist()
    outer_points_x = (outer_radius*cos(t)).tolist()
    outer_points_y = (outer_radius*sin(t)).tolist()
    xpts = inner_points_x + outer_points_x[::-1]
    ypts = inner_points_y + outer_points_y[::-1]
    D.add_polygon(points = (xpts,ypts), layer = layer)
    return D

def extend(port, length = 20,layer = 0):
	'''
	Extend the port with a straight waveguide.
	'''
	D = waveguide(port.width, length, layer=layer)
	D.rotate(port.orientation).move(origin = D.ports[1], destination = port)
	return D

# racetrack
def racetrack(radius = 10, width = 0.5, lc = 5, angle_resolution = 2.5, layer = 0):
    D = Device(name = 'racetrack')
    inner_radius = radius - width/2
    outer_radius = radius + width/2
    angle = np.append(np.linspace(0, 360, np.ceil(360/angle_resolution)),[90,90,270,270])
    angle.sort()
    t=angle*np.pi/180
    inner_points_x = [i+np.sign(i)*lc/2 for i in inner_radius*cos(t)]
    inner_points_y = (inner_radius*sin(t)).tolist()
    outer_points_x = [i+np.sign(i)*lc/2 for i in outer_radius*cos(t)]
    outer_points_y = (outer_radius*sin(t)).tolist()
    xpts = inner_points_x + outer_points_x[::-1]
    ypts = inner_points_y + outer_points_y[::-1]
    D.add_polygon(points = (xpts,ypts), layer = layer)
    return D

# ring with triple slot
def slot_ring(radius=50,CGS=[1.0,0.1,0.3], angle_resolution = 2.5, layer = 0):
    D = Device(name = 'slot_ring')
    cR=radius
    iR=radius-CGS[0]*0.5-CGS[1]-CGS[2]*0.5
    oR=radius+CGS[0]*0.5+CGS[1]+CGS[2]*0.5
    D << new_ring(radius=cR,width=CGS[0], layer=layer)
    D << new_ring(radius=iR,width=CGS[2], layer=layer)
    D << new_ring(radius=oR,width=CGS[2], layer=layer)
    return D

# arc shape grating coupler
@device_lru_cache
def arc_grating(num_periods = 20, period = 0.75, fill_factor = 0.5, angle = 45, length_taper = 5, width = 0.5, layer = 0):
    #returns a fiber grating
    G = Device('grating')

    # make the grating teeth
    cradius = length_taper + period*(1-0.5*fill_factor)
    for i in range(num_periods):
        cgrating = G.add_ref(arc(radius=cradius, start_angle=-angle/2, theta=angle, width=period*fill_factor, layer = layer))
        cradius += period

    # make the taper
    out_len = width*0.5/np.tan(angle/360*np.pi)
    A = bbox([(0,-width/2),(out_len,width/2)])
    B = arc(radius=length_taper/2,start_angle=-angle/2,theta=angle,width=length_taper,layer=layer)
    G.add_ref(boolean(A, B, operation = 'a+b'))
    p = G.add_port(name = 1, midpoint = (0,0), width = width, orientation = 180)

    G.flatten()
    return G

# from numpy import cos, sin, arctan
r2d = lambda x : x*180/np.pi
d2r = lambda x : x*np.pi/180
arg = lambda x, y: np.arctan((y[1]-y[0])/(x[1]-x[0]))
def archimedes(bent = 20, width = 0.5, n = 1, distance = 10, angle_resolution = 1, layer = 1):
    shift = distance/np.pi
    D = Device('archimedes')
    sp_t = np.linspace(0, 360*n, np.ceil(360*n/angle_resolution))*np.pi/180

    rho = 2*bent + sp_t*shift
    inner_rho = rho-width*.5
    outer_rho = rho+width*.5
    inner_sp_x = (inner_rho*cos(sp_t)).tolist()
    inner_sp_y = (inner_rho*sin(sp_t)).tolist()
    outer_sp_x = (outer_rho*cos(sp_t)).tolist()
    outer_sp_y = (outer_rho*sin(sp_t)).tolist()
    # arc section
    inner_t = arg(inner_sp_x,inner_sp_y)
    arc_radius = bent/sin(inner_t)
    arc_t = np.linspace(1.5*np.pi-inner_t, 1.5*np.pi+inner_t,100)[:-1]
    inner_arc_x = ( arc_radius*sin(inner_t) + (arc_radius-width*0.5)*cos(arc_t) ).tolist()
    inner_arc_y = ( arc_radius*cos(inner_t) + (arc_radius-width*0.5)*sin(arc_t) ).tolist()
    outer_arc_x = ( arc_radius*sin(inner_t) + (arc_radius+width*0.5)*cos(arc_t) ).tolist()
    outer_arc_y = ( arc_radius*cos(inner_t) + (arc_radius+width*0.5)*sin(arc_t) ).tolist()
    
    outer_t = arg(outer_sp_x[::-1], outer_sp_y[::-1]) 
    ext_radius = (inner_sp_x[-1]+outer_sp_x[-1])*0.5/sin(outer_t)
    ext_t = np.linspace(outer_t-np.pi*0.5, np.pi*.5,100)[0:]
    inner_ext_x = ( (ext_radius-width*0.5)*cos(ext_t) ).tolist()
    inner_ext_y = ( ext_radius*cos(outer_t) + (ext_radius-width*0.5)*sin(ext_t) ).tolist()
    outer_ext_x = ( (ext_radius+width*0.5)*cos(ext_t) ).tolist()
    outer_ext_y = ( ext_radius*cos(outer_t) + (ext_radius+width*0.5)*sin(ext_t) ).tolist()
    
    xpts = inner_arc_x + inner_sp_x + inner_ext_x + outer_ext_x[::-1] + outer_sp_x[::-1] + outer_arc_x[::-1]
    ypts = inner_arc_y + inner_sp_y + inner_ext_y + outer_ext_y[::-1] + outer_sp_y[::-1] + outer_arc_y[::-1]
    D.add_polygon(points = (xpts,ypts), layer = layer)
    D.add_polygon(points = (xpts,ypts), layer = layer).rotate(180)
    
    D.add_port(name = 1, midpoint = [0,D.ymax-width/2], width = width, orientation = 180)
    D.add_port(name = 2, midpoint = [0,D.ymin+width/2], width = width, orientation = 0)
    return D

# def semi_spiral(bend=20,shift=10,width=1,layer=1, n=4, angle_resolution=1):
#     D = Device('semi_spiral')
#     inn = arc(radius=(2*bend-shift)/4,start_angle=0,theta=180,width=width).rotate(180).movex(-(2*bend-shift)/4)
#     D << copy(inn).rotate(180)
#     D << inn
#     for i in range(n):
# #         radius = bend + i*shift
#         out = arc(radius=bend+shift*i,start_angle=0,theta=180,width=width).movex(shift/2)
#         D << out
#     for i in range(n):
#         out = arc(radius=bend+shift*i,start_angle=180,theta=180,width=width).movex(-shift/2)
#         D << out
#     WG = waveguide(length=bend+shift*n,width=width).rotate(90).movex(-bend-shift*n+shift/2)
#     D << WG
#     D << copy(WG).rotate(180)
#     D.add_port(name = 1, midpoint = [D.xmin+width/2,D.ymax], width = width, orientation = 90)
#     D.add_port(name = 2, midpoint = [-D.xmin-width/2,-D.ymax], width = width, orientation = 270)
#     return D.rotate(90)

# ring resonator and coupling bus, all pass
def ALLPASS(width_rg=1, width_bus=1, radius=30, gap=0.1, layer=0):
    D = Device('MRRBUS')
    RG = D << ring(radius=radius,width=width_rg, layer=layer) 
    BUS = D << waveguide(length=2*radius+width_rg,width=width_bus, layer=layer).movex(-radius-width_rg/2)
    # align ring & bus, rotate
    BUS.ymax = RG.ymin - gap
    D.add_port(name = 1, midpoint = [D.xmax,D.ymin+width_bus/2], width = width_bus, orientation = 0)
    D.add_port(name = 2, midpoint = [D.xmin,D.ymin+width_bus/2], width = width_bus, orientation = 180)
    return D

# ring resonator and coupling bus, four part
def FOURPORT(width_rg=1, width_bus=1, radius=30, gap=0.1, layer=0):
    D = Device('MRRBUS')
    RG = D << ring(radius=radius,width=width_rg, layer=layer) 
    BUS1 = D << waveguide(length=2*radius+width_rg,width=width_bus, layer=layer).movex(-radius-width_rg/2)
    BUS2 = D << waveguide(length=2*radius+width_rg,width=width_bus, layer=layer).movex(-radius-width_rg/2)
    # align ring & bus, rotate
    BUS1.ymax = RG.ymin - gap
    BUS2.ymin = RG.ymax + gap    
    D.add_port(name = 1, midpoint = [D.xmax,D.ymin+width_bus/2], width = width_bus, orientation = 0)
    D.add_port(name = 2, midpoint = [D.xmin,D.ymin+width_bus/2], width = width_bus, orientation = 180)
    D.add_port(name = 3, midpoint = [D.xmax,D.ymax-width_bus/2], width = width_bus, orientation = 0)
    D.add_port(name = 4, midpoint = [D.xmin,D.ymax-width_bus/2], width = width_bus, orientation = 180)
    return D