import peakutils
from peakutils.plot import plot as pplot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
import math
from scipy.sparse.linalg import spsolve

PATH = "PDha_1.csv"  # csv or txt format
df = pd.read_csv( PATH, sep="\t", skiprows=[0, 1, 2, 3, 4],delimiter="," )  # in skiprows write the rows which contain text in order to elliminate them
df.to_numpy()
class baselinecorrection:
    def __init__(self,path):
        self.path=path
        self.read_data()
    def read_data(self):
        df = pd.read_csv( self.path, sep="\t", skiprows=[0, 1, 2, 3, 4],delimiter="," )  # in skiprows write the rows which contain text in order to eliminate them
        df.to_numpy()
        data = df.values.T
        y = data
        x = data[0]  # wavelength
        self.x = np.nan_to_num( x )
        y = np.nan_to_num( y )
        self.y = np.delete( y, 0, 0 ) #spectral data


    def baseline_subtraction(self,rawspc):
        b = len( rawspc )
        t = int( input( "plz enter polynomial degree " ) )
        c = np.zeros( rawspc.shape )
        for i in range( b ):
            base = peakutils.baseline(rawspc[i], t)  # choose order of polynomial here
            c[i] = rawspc[i] - base
            c[c < 0] = 0
        return c, base

    def baseline_als(self,rawspc, lam, p, niter):
        """

        :param lam:  λ for smoothness default 10**5
        :param p:  p for asymmetry default: 0.1
        :param niter: number of iteration default: 10
        :return: least square baseline error
        """
        y = rawspc[0, :]
        l = len( y )
        eye = np.eye( l )
        D = sparse.csc_matrix( np.diff( eye, 2 ) )
        w = np.ones( l )
        for i in range( niter ):
            W = sparse.spdiags( w, 0, l, l )
            Z = W + lam * D.dot( D.transpose() )
            z = spsolve( Z, w * y )
            w = p * (y > z) + (1 - p) * (y < z)
        return(z)
    def find_index(self, wavelenght_arr, number):
        for index, value in enumerate(wavelenght_arr):
            if abs(value - number) < 0.0000001:
                index = float(index)
                return index
        return index

    def region_select(self):
        """

        :return: new x and y due to selected region
        """
        wu = float( input( "plz enter upper boundry wavelenght " ) )
        ou = self.find_index(self.x, wu )
        wl = float( input( "plz enter lower boundry wavelenght " ) )
        ol = self.find_index(self.x, wl )
        ou=int(ou)
        ol=int(ol)
        xi = self.x[ol:ou, ]
        yi = self.y[:, ol:ou]
        return xi, yi
    def plot_sub(self,raw,wavelebgth):
        fig, ax = plt.subplots(3)
        ax[0].plot(wavelebgth, raw.T)
        ax[0].set_title( "Raw spectra", fontsize=15 )
        ax[1].plot(wavelebgth, c.T)
        ax[1].set_title( "corrected spectra", fontsize=15 )
        ax[2].plot(wavelebgth, base)
        ax[2].set_title( "noise", fontsize=15 )
        plt.savefig( 'FTIR subtraction poly.png' )
        plt.tight_layout()
        plt.show()
        plt.plot(wavelebgth, c.T)
        plt.title("corrected graph")
        plt.savefig('corrected graph.png')
        plt.show()
        return
    def plot_als(self,raw,wavelebgth):
        c = self.baseline_als(raw,lam=10**5, p=0.1, niter=10)
        fig, ax = plt.subplots(3)
        ax[0].plot(wavelebgth, raw.T)
        ax[0].set_title("Raw spectra", fontsize=15)
        ax[1].plot(wavelebgth, (raw - c).T)
        ax[1].set_title("corrected spectra", fontsize=15)
        ax[2].plot(wavelebgth, c)
        ax[2].set_title("noise", fontsize=15)
        plt.savefig('corrected graph.png')
        plt.tight_layout()
        plt.show()
        plt.plot(wavelebgth, (raw - c).T)
        plt.title("FTIR corrected")
        plt.savefig('FTIR corrected.png')
        plt.show()
        return




c1 = baselinecorrection(PATH)
wl = c1.x
raw = c1.y
plt.plot(wl, raw.T)
plt.title("raw spectra")
plt.show()
r=np.zeros(df.shape)


l = input( "Do you need to select desire regions? yes =y no =n" ).lower()
if l == "n":  # 1. whole spectrum fitting
    l1 = input( "least square or subtraction method? least square =a subtraction=s" ).lower()
    if l1 == "s":
        c, base = c1.baseline_subtraction( raw )
        r[:, 0] = wl.transpose()
        r[:, 1:] = c.transpose()
        np.savetxt( 'baseline corrected spectra.txt', r, delimiter='\t' )  # change plot title base on your poly degree
        c1.plot_sub(raw,wl)
    elif l1 == "a":
        G = c1.baseline_als( lam=10 ** 5, p=0.1, niter=10 )
        G1 = raw - G
        r[:, 0] = wl.transpose()
        r[:, 1:] = G1.transpose()
        np.savetxt( 'baseline corrected spectra.txt', r, delimiter='\t' )  # change plot title base on your poly degree
        c1.plot_als(raw,wl)
    else:
        raise ValueError( "this value is not valid for arg method(method must be s or a)" )
elif l == "y":
    l1 = input( "least square or subtraction method? least square =a subtraction=s" ).lower()
    if l1 == "s":
        wli, rawi = c1.region_select()
        c, base = c1.baseline_subtraction( rawi )
        c1.plot_sub(rawi,wli)
    elif l1 == "a":
        wli, rawi = c1.region_select()
        c = c1.baseline_als(rawi,lam=10**5,p=0.1,niter=10)
        c1.plot_als(rawi,wli)
    else:
        raise ValueError("this value is not valid for arg method(method must be s or a)")
else:
    raise ValueError("this value is not valid for arg method(method must be y or n)")
