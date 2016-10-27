pyFlamesTk

The development of this package was developed initially by Eugene J. Ragasa
at the University of Florida starting in 2013 and currently continuing to the
present.  

The name pyFlamesTk comes from 

"python" 
+ "Florida Laboratory for Advanced Materials Simulation" 
+ "ToolKit"

The Florida Laboratory for Advanced Materials Simulation(FLAMES) is a 
computational materials reserach group in the Materials Science and Engineering
Department at the University of Florida.  While there I was mentored by
Dr Simon Philpot, Dr. Susan Sinnot, Dr. Richad Hennig, Dr. Alex Chernatysnkiy, 
and Dr. Tao Liang.

This package as written to require few dependencies other than an Anaconda
distribution of python which includes some key scientific libraries to prevent
problems of using this package on supercomputer clusters which lack the ability
to install packages onto the computers.

Some of the code within this package is pedagogical in nature which provided 
the ability for me to learn material science.  As a programmer-at-heart, I
detest redoing the same task twice, and like to develop reusable tools.

Major Functionalities:

Raman Spectrum Calculator

Interatomic Potential Fitter.  This portion of the software was developed for
a project for SANDIA which formed a major portion of my dissertation.  It 
uses DAKOTA for conjugate gradient minimization, pareto set deterimination.
The molecular dynamics calculations are provided by an interface to LAMMPS.  
Dr. Steven Foiles (SANDIA), Dr. Stephen O'Brien (SANDIA), Dr. Simon Philpot (UF),
and Dr. Richard Hennig (UF) were my collaborators on the project.  Special thanks 
to Dr. Jacky Martinez(UF) whose work on posMAT and whose insight into fitting
potentials were invaluable.

Authors:

Eugene J. Ragasa (eragsa@ufl.edu)

References: