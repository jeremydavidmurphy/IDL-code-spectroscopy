# IDL-code-spectroscopy

This is code I wrote while in graduate school. Most of it is designed to reduce spectroscoptic data from the Mitchell Spectrograph (formally VIRUS-P) currently in use on the Harlan J. Smith telescope at McDonald Observatory. My observational work was focused on measuring the velocity and velocity dispersion of massive elliptical galaxies, and then using that data to constrain the amount of dark matter present in the galaxy. This was done running large number of n-body simulations using the TACC supercomputers.

The primary code in this directory is the not-so-cleverly-named pipe1.pro and pipe2.pro. These routines take the partly-reduced data from the telescope and complete the reduction. The versions pipe1C and pipe2C are modifications of those files to work on spectra that has been collapsed in other routines.

Much of the rest of the code you'll find here are designed to carry out specfic tasks, many of which were focused on data visualization.
