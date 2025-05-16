#!/usr/bin/env python

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import argparse
import os
import pickle

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def plot_hotspots(m, base_path = '.', lon360 = False, **kwargs):
    hotspots = pickle.load(open('%s/data/hotspots.pkl' % (base_path), 'rb'))
    if lon360:
        hotspots[:,0] = (hotspots[:,0] + 360) % 360.0
    x, y = m(hotspots[:,0], hotspots[:,1])
    if kwargs:
        m.scatter(x, y, **kwargs)
    else:
        m.scatter(x, y)

def plot_plates(m, base_path = '.', lon360 = False, **kwargs):
    for bound in ['ridge', 'transform', 'trench']:
        name, segs = pickle.load(open('%s/data/%s.pkl' % (base_path,bound), 'rb'))
        ind_nan, = np.nonzero(np.isnan(segs[:,0]))
        segs[ind_nan,0] = 0
        segs[ind_nan,1] = 0
        if lon360:
            segs[:,0] = (segs[:,0] + 360) % 360.0
        x, y = m(segs[:,0], segs[:,1])
        x[ind_nan] = np.nan
        y[ind_nan] = np.nan
        dx = np.abs(x[1:] - x[:-1])
        ind_jump, = np.nonzero(dx > 1000000)
        x[ind_jump] = np.nan
        if kwargs:
            m.plot(x, y, '-', **kwargs)
        else:
            m.plot(x, y, '-')

def plot_model(fname, vmin, vmax):
    # load the model data; extract the lon, lat, and dlnVs columns
    tmp = np.loadtxt(fname)
    lon = tmp[:,1].reshape((181,361))
    lat = tmp[:,2].reshape((181,361))
    dvs = tmp[:,3].reshape((181,361))
    # initialize figure and axes
    fig = plt.figure(figsize=(8,4))
    ax_map = fig.add_axes([0,0,0.9,1.0])
    ax_cbr = fig.add_axes([0.91,0.3,0.01,0.4])
    # set up map
    m = Basemap(projection='robin', lon_0=150, resolution='c', ax=ax_map)
    clip_path = m.drawmapboundary()
    m.drawcoastlines()
    m.drawparallels(np.arange(-90,90,30))
    m.drawmeridians(np.arange(-180,180,30))
    # plot the model
    s = m.transform_scalar(dvs, lon[0,:], lat[:,0], 1000, 500)
    im = m.imshow(s, cmap=plt.cm.jet_r, clip_path=clip_path, vmin=vmin, vmax=vmax)
    # add plates and hotspots
    base_path = os.path.dirname(os.path.abspath(__file__))
    plot_plates(m, base_path=base_path, color='k')
    plot_hotspots(m, base_path=base_path, s=30, color='m', edgecolor='k')
    # add a colorbar
    cb = plt.colorbar(im, cax=ax_cbr)
    cb.set_label('dlnVs (%)')
    # done
    plt.show()
    
def main():
    parser = argparse.ArgumentParser(description='plot exported model samples')
    parser.add_argument('-u', '--upper-bound', type=float, default=4.0,
        help='Upper bound for color scale saturation level (percent)')
    parser.add_argument('-l', '--lower-bound', type=float, default=-4.0,
        help='Lower bound for color scale saturation level (percent)')
    parser.add_argument('-f', '--file', type=str, default='model-samples.out',
        help='File name from which to read')
    arg = parser.parse_args()
    plot_model(arg.file, arg.lower_bound, arg.upper_bound)

if __name__ == "__main__":
    main()
