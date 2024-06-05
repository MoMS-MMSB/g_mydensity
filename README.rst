==============================================
g_mydensity: Calculate local partial densities
==============================================

``g_mydensity`` is a modified version of the ``g_density`` GROMACS tool. It
adds the possibility to calculate partial density landscape on a plane and
partial density profile as a function of the distance to a group of atoms.

Only number and mass densities are supported yet.

Even if ``g_density`` regular features are available in ``g_mydensity`` you
should use the original tool to use them.

Installation
============

``g_mydensity`` depends on the GROMACS software package, which needs to be
installed.  Only versions 4.5.x have been tested, but ``g_mydensity`` might be
compatible with other versions of GROMACS.

To install ``g_mydensity``, GROMACS needs to be loaded. You can load
it using:

    source /path_to_gromacs/bin/GMXRC

Go into the source directory of the program, then run ``make``. The
``g_mydensity`` executable should be created.  Make sure that
this executable is in the research path of your shell.

Usage
=====
Here we assume that ``g_mydensity`` is in the research path of your shell. To
get some help just run ``g_mydensity -h``. All available options will be
listed.

A classical use would be:

    g_mydensity -f traj.xtc -s topol.tpr -n index.ndx -o order.xvg

GROMACS needs to be loaded for ``g_mydensity`` to work.

Basic arguments
----------------

So as most GROMACS tools, the basic arguments are:

* ``-f``: the path to the trajectory to read;
* ``-s``: the path to the topology (tpr file);
* ``-n``: the path to an index file that describe the group of atoms you are
  interested in;

The ``-dens`` argument can be use to switch from mass density (when set to
"mass", default) to number density (when set to "number"). The "charge" and
"electron" options are not implemented.

Output control
--------------

The following arguments control the output. You can get either one or both of
the possible outputs but you need to select at least one of them.

* ``-og``: produce the partial density landscape; a file path can be given as
  argument. The produced file can be converted into a picture. See the
  `Generate picture from landscapes`_ section to know more
  about that.
* ``-od``: produce the partial dentity profile as a function of distance to a
  group.  Distance is calculated as a function of the center of mass of a
  reference group. The distance is calculated in 2D by default, the normal axis
  is ignored in the calculation. To calculate distances in 3D, use the ``-3d``
  option.

Generate pictures from landscapes
---------------------------------

The landscape output is a text file describing the order parameter values on a
grid. The file format is not XPM like most grid outputs produced by GROMACS
tools so the ``xpm2ps`` utility can not be used to produce usable pictures. The
``dispgrid`` python script aims to exploit the data and to produce pictures from
them.

You need python 2.7, with the numpy and matplotlib modules to run dispgrid.

Basic usage of dispgrid is:

    dispgrid input.dat output.png

See the help available by typing ``dispgrid -h`` for more features.
