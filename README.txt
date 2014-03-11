 ___       ____  __  __ 
/ __) ___ (  _ \(  \/  )  
\__ \( _ ) )___/ )    (   Statistical nonParametric Mapping toolbox
(___/(_)_)(__)  (_/\/\_)  SnPM -  http://www.fil.ion.ucl.ac.uk/spm/snpm

                              R E A D M E

________________________________________________________________________

This README gives a brief introduction to the installation and use of
the SnPM toolbox for SPM. Full details can be found on the SnPMweb site:
                 http://www.fil.ion.ucl.ac.uk/spm/snpm

A manifest for this release is contained in the file Contents.m
The release is described in the file snpm.man
 
________________________________________________________________________
                                                                    SnPM

   Statistical nonParametric Mapping refers to the enterprise
   of making statistical inferences on volumetric statistic
   images with minimal assumptions using non-parametric statistical
   techniques. SnPM refers to an implementation of Statistical
   nonParametric Mapping by Andrew Holmes and Thomas Nichols.

________________________________________________________________________
                                                            Installation

The SnPM8 software is a suite of MatLab functions and scripts, which
utilise the SPM API to implement Statistical nonParametric Mapping.

SnPM8 therefore *must* be installed alongside SPM8.

Installation is relatively straightforward. Download the gzipped tar
file and unzip and untar it; it will create a directory named snpm99.
Add the full pathname of this directory to your MATLABPATH,
ensuring that the first SPM distribution on the MATLABPATH is
SPM8. (You can check this by typing which spm in MatLab and examining
the path, or try spm ver).  Note that none of the names of SnPM3
routines clash with those of SPM8 so the positioning of SnPM on the
MATLABPATH is not important

For step by step instructions on downloading and installing SPM
distributions, refer to the SPM distribution page:

	http://www.fil.ion.ucl.ac.uk/spm/distrib.html
	

________________________________________________________________________
                                                         Getting started
                                                         
SnPM is invoked with the command `snpm` at the MatLab prompt. SnPM8
requires SPM8, and will launch the SPM environment if not already
started.

We recommend you start by reviewing the help system, by selecting "About
SPM" from the splash screen. This initially displays the "snpm.man"
topic, detailing the toolbox. Press the "Menu" button to display a
representation of the SPM menu window, with buttons linked to
appropriate help pages.

Before attempting to analyze data using SnPM, we recommend you spend
some time reading. It is essential to understand the concepts of
Statistical nonParametric Mapping in order to effectively use the
software as a research tool. You should begin with the SnPMweb pages,
which give a theoretical primer, an overview of the toolbox, a manual,
and a worked example. You are also urged to read the appropriate literature.

________________________________________________________________________
                                                               Resources

The SnPMweb site is the central repository for all SnPM resources:
                 http://www.fil.ion.ucl.ac.uk/spm/snpm
Introductory material, installation details, documentation, examples
and patches are published on the site.

SnPMweb is part of SPMweb, the central repository for all SPM resources:
                 http://www.fil.ion.ucl.ac.uk/spm

There is an SPM eMail discussion list, hosted at <spm@mailbase.ac.uk>.
The list is monitored by the authors, and discusses theoretical,
methodological and practical issues of Statistical Parametric Mapping,
SPM and SnPM. See the SPMweb site for details on subscribing:

                 http://www.fil.ion.ucl.ac.uk/spm/help

Please report bugs to the authors at <snpm-authors@fil.ion.ucl.ac.uk>.
Peculiarities may actually be features, and should be raised on the SPM
eMail discussion list, <spm@mailbase.ac.uk>.

________________________________________________________________________
                                                              Disclaimer

SnPM (being the collection of files given in the manifest in the
Contents.m file) is free but copyrighted software, (c) 2009 Thomas
Nichols & Andrew Holmes.  
SnPM is distributed under the terms of the GNU General
Public Licence as published by the Free Software Foundation (either
version 2, as given in file spm_LICENCE.man, or at your option, any
later version). Further details on "copyleft" can be found at
http://www.gnu.org/copyleft/. In particular, SnPM is supplied as is.
No formal support or maintenance is provided or implied.

________________________________________________________________________
SnPM is developed by 
   Andrew Holmes
   Thomas Nichols, Department of Statistics & Warwick Manufacturing
   Group, Univeristy of Warwick

% $Id $








