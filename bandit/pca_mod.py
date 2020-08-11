# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

r"""Principal Component Analysis (PCA) --- :mod:`MDAnalysis.analysis.pca`
=====================================================================

:Authors: John Detlefs
:Year: 2016
:Copyright: GNU Public License v3

.. versionadded:: 0.16.0

This module contains the linear dimensions reduction method Principal Component
Analysis (PCA). PCA sorts a simulation into 3N directions of descending
variance, with N being the number of atoms. These directions are called
the principal components. The dimensions to be analyzed are reduced by only
looking at a few projections of the first principal components. To learn how to
run a Principal Component Analysis, please refer to the :ref:`PCA-tutorial`.

The PCA problem is solved by solving the eigenvalue problem of the covariance
matrix, a :math:`3N \times 3N` matrix where the element :math:`(i, j)` is the
covariance between coordinates :math:`i` and :math:`j`. The principal
components are the eigenvectors of this matrix.

For each eigenvector, its eigenvalue is the variance that the eigenvector
explains. Stored in :attr:`PCA.cumulated_variance`, a ratio for each number of
eigenvectors up to index :math:`i` is provided to quickly find out how many
principal components are needed to explain the amount of variance reflected by
those :math:`i` eigenvectors. For most data, :attr:`PCA.cumulated_variance`
will be approximately equal to one for some :math:`n` that is significantly
smaller than the total number of components, these are the components of
interest given by Principal Component Analysis.

From here, we can project a trajectory onto these principal components and
attempt to retrieve some structure from our high dimensional data.

For a basic introduction to the module, the :ref:`PCA-tutorial` shows how
to perform Principal Component Analysis.

.. _PCA-tutorial:

PCA Tutorial
------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF` and
:data:`~MDAnalysis.tests.datafiles.DCD`). This tutorial shows how to use the
PCA class.

First load all modules and test data

    >>> import MDAnalysis as mda
    >>> import MDAnalysis.analysis.pca as pca
    >>> from MDAnalysis.tests.datafiles import PSF, DCD

Given a universe containing trajectory data we can perform Principal Component
Analyis by using the class :class:`PCA` and retrieving the principal
components.

    >>> u = mda.Universe(PSF, DCD)
    >>> PSF_pca = pca.PCA(u, select='backbone')
    >>> PSF_pca.run()

Inspect the components to determine the principal components you would like
to retain. The choice is arbitrary, but I will stop when 95 percent of the
variance is explained by the components. This cumulated variance by the
components is conveniently stored in the one-dimensional array attribute
``cumulated_variance``. The value at the ith index of `cumulated_variance`
is the sum of the variances from 0 to i.

    >>> n_pcs = np.where(PSF_pca.cumulated_variance > 0.95)[0][0]
    >>> atomgroup = u.select_atoms('backbone')
    >>> pca_space = PSF_pca.transform(atomgroup, n_components=n_pcs)

From here, inspection of the ``pca_space`` and conclusions to be drawn from the
data are left to the user.

Classes and Functions
---------------------

.. autoclass:: PCA
.. autofunction:: cosine_content

"""
from __future__ import division, absolute_import
from six.moves import range
import warnings

import numpy as np
import scipy.integrate

from MDAnalysis import Universe
from MDAnalysis.analysis.align import _fit_to
from MDAnalysis.lib.log import ProgressMeter

from .base import AnalysisBase
from .pca import PCA

class PCA_MOD(PCA):
    def transform(self, atomgroup, n_components=None, start=None, stop=None,
                  step=None,eigen_vec=None):
        """Apply the dimensionality reduction on a trajectory

        Parameters
        ----------
        atomgroup: MDAnalysis atomgroup/ Universe
            The atomgroup or universe containing atoms to be PCA transformed.
        n_components: int, optional
            The number of components to be projected onto, Default none: maps
            onto all components.
        start: int, optional
            The frame to start on for the PCA transform. Default: None becomes
            0, the first frame index.
        stop: int, optional
            Frame index to stop PCA transform. Default: None becomes n_frames.
            Iteration stops *before* this frame number, which means that the
            trajectory would be read until the end.
        step: int, optional
            Number of frames to skip over for PCA transform. Default: None
            becomes 1.

        Returns
        -------
        pca_space : array, shape (number of frames, number of components)

        .. versionchanged:: 0.19.0
           Transform now requires that :meth:`run` has been called before,
           otherwise a :exc:`ValueError` is raised.
        """
        if np.any(eigen_vec) != None:
            self.p_components = eigen_vec
            print(eigen_vec - self.p_components)
        if not self._calculated:
            raise ValueError('Call run() on the PCA before using transform')

        if isinstance(atomgroup, Universe):
            atomgroup = atomgroup.atoms

        if(self._n_atoms != atomgroup.n_atoms):
            raise ValueError('PCA has been fit for'
                             '{} atoms. Your atomgroup'
                             'has {} atoms'.format(self._n_atoms,
                                                   atomgroup.n_atoms))
        if not (self._atoms.types == atomgroup.types).all():
            warnings.warn('Atom types do not match with types used to fit PCA')

        traj = atomgroup.universe.trajectory
        start, stop, step = traj.check_slice_indices(start, stop, step)
        n_frames = len(range(start, stop, step))

        dim = (n_components if n_components is not None else
               self.p_components.shape[1])

        dot = np.zeros((n_frames, dim))

        for i, ts in enumerate(traj[start:stop:step]):
            xyz = atomgroup.positions.ravel() - self.mean
            dot[i] = np.dot(xyz, self.p_components[:, :n_components])

        return dot


