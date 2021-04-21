Installing the Cubit plugin
==================================

There is a common location that stores all the plugins that are available for
use. Go to the `DAGMC Plugins <DAGMC_plugins_>`_ page and download the plugin
appropriate for your operating system.

Linux install
~~~~~~~~~~~~~

Obtain the Linux plugin from the link above. These instructions assume that you
have installed Cubit in ``/opt/Coreform-Cubit-2020.2``. If you installed it somewhere
else, modify these instructions appropriately.

These operations may need to be performed as the root user (sudo).

1.  Unpack the tarball in ``/opt/Coreform-Cubit-2020.2/bin/plugins``.
2.  Change directory to ``/opt/Coreform-Cubit-2020.2/bin/plugins/dagmc``.
3.  Run the install script: ``./install.sh``.

OS/X install
~~~~~~~~~~~~

Obtain the OS/X plugin from the link above. These instructions assume that you
have installed Cubit in ``/Applications/Coreform-Cubit-2020.2.app``. If you installed it
somewhere else, modify these instructions appropriately.

These operations may need to be performed as the root user (sudo).

1.  If it does not already exist, create a folder called ``plugins`` in
    ``/Applications/Trelis-16.0/Contents/MacOS``.
2.  Unpack the tarball in ``/Applications/Trelis-16.0/Contents/MacOS/plugins``.
3.  Change directory to
    ``/Applications/Trelis-16.0/Contents/MacOS/plugins/dagmc``.
4.  Run the install script: ``./install.sh``.

Windows install
~~~~~~~~~~~~~~~

A plugin for Windows is currently being developed; it will be posted to the
DAGMC plugins page when it is ready.

..  _DAGMC_plugins: https://go.wisc.edu/svalinn-cubit-plugin
