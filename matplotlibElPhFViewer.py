## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "matplotlibElPhFViewer.py"
 #                                     created: 11/1/06 {4:34:01 PM}
 #                                 last update: 7/24/07 {8:57:32 AM}
 # Author: Jonathan Guyer
 # E-mail: <guyer@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are 
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 # 
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 2006-11- 1 JEG 1.0 original
 # 
 # ########################################################################
 ##

import pylab
 
from fipy.viewers.matplotlibViewer.matplotlibViewer import AbstractMatplotlibViewer

class MatplotlibElPhFViewer(AbstractMatplotlibViewer):
    def __init__(self, phase, potential, components, charge, potentials, limits=None):
        AbstractMatplotlibViewer.__init__(self, 
                                  vars=[phase] + components + [potential, charge],
                                  limits=limits)
        self.phase = phase
        self.potential = potential
        self.components = components
        self.charge = charge
        
        pylab.ioff()
        
        pylab.subplot(511)
        self.lines = [pylab.plot(*datum) for datum in self._getData([phase])]
        pylab.ylim(ymin=0., ymax=1.)
        pylab.xlim(xmin = self._getLimit('xmin'),
                   xmax = self._getLimit('xmax'))
        pylab.ylabel(r'$\xi$')

        pylab.subplot(512)
#         self.lines += [pylab.plot(*datum) for datum in self._getData(self.components)]
        self.lines += [pylab.semilogy(*datum) for datum in self._getData(self.components)]
        pylab.legend([var.name for var in self.components])
        pylab.xlim(xmin = self._getLimit('xmin'),
                   xmax = self._getLimit('xmax'))
        pylab.ylabel(r'$\bar{V}_S C_j$')

        pylab.subplot(513)
        self.lines += [pylab.plot(*datum) for datum in self._getData([potential])]
        pylab.ylabel(r'$\phi\cal{F} / R T$')
        pylab.xlim(xmin = self._getLimit('xmin'),
                   xmax = self._getLimit('xmax'))

        pylab.subplot(514)
        self.lines += [pylab.plot(*datum) for datum in self._getData([charge])]
        pylab.xlim(xmin = self._getLimit('xmin'),
                   xmax = self._getLimit('xmax'))
        pylab.ylabel(r'$\bar{V}_S  \sum_j z_j C_j$')

        pylab.subplot(515)
        self.chemicalPotentials = potentials # [var.muBar for var in potentials] #  + [self.phase]
        self.lines += [pylab.plot(*datum) for datum in self._getData(self.chemicalPotentials)]
        pylab.legend([var.name for var in potentials]) #  + [self.phase]
        pylab.xlim(xmin = self._getLimit('xmin'),
                   xmax = self._getLimit('xmax'))
        pylab.ylabel(r'$\bar{\mu}_j / R T$')

        pylab.xlabel(r'$x/L$')
        
        pylab.draw()
        pylab.ion()

    def _getData(self, vars):
        from fipy.tools.numerix import array
        return [[array(var.mesh.cellCenters[0]), array(var)] for var in vars]

    def _plot(self):
        from fipy.tools import numerix
        
        pylab.subplot(511)
        ymin, ymax = self._autoscale([self.phase], 
                                     datamin=self._getLimit('phasemin'), 
                                     datamax=self._getLimit('phasemax'))
        pylab.ylim(ymin=ymin, ymax=ymax)

        pylab.subplot(512)
        ymin, ymax = self._autoscale(self.components, 
                                     datamin=self._getLimit('Cmin'), 
                                     datamax=self._getLimit('Cmax'))
        pylab.ylim(ymin=ymin, ymax=ymax)

        pylab.subplot(513)
        ymin, ymax = self._autoscale([self.potential], 
                                     datamin=self._getLimit('Vmin'), 
                                     datamax=self._getLimit('Vmax'))
        pylab.ylim(ymin=ymin, ymax=ymax)

        pylab.subplot(514)
        ymin, ymax = self._autoscale([self.charge], 
                                     datamin=self._getLimit('rhomin'), 
                                     datamax=self._getLimit('rhomax'))
        pylab.ylim(ymin=ymin, ymax=ymax)

        pylab.subplot(515)
        ymin, ymax = self._autoscale(self.chemicalPotentials, 
                                     datamin=self._getLimit('mumin'), 
                                     datamax=self._getLimit('mumax'))
        pylab.ylim(ymin=ymin, ymax=ymax)


        for line, datum in zip(self.lines, 
                               self._getData([self.phase] 
                                             + self.components 
                                             + [self.potential, self.charge]
                                             + self.chemicalPotentials)):
            line[0].set_xdata(datum[0])
            line[0].set_ydata(datum[1])
