# OCAPE Implementation
This project serves as a generalized implementation of Ocean Convective Available Potential Energy (OCAPE) that will work with given values for potential temperature, salinity (in psu), and depth.
OCAPE (described in detail at https://authors.library.caltech.edu/66427/1/jpo-d-14-0155%252E1.pdf) is the maximum energy of a water column that can be converted into kinetic energy. OCAPE is the basis of thermobaric convection which comes from the nonlinear equation of state (EOS) of water. 
In this implementation, we used the Jackett and McDougall (1995) water equation of state. However, using another equation of state would be fairly trivial.
In order to calulate the value for OCAPE, we need the current enthalpy and reference value for enthalpy (the minimum). In order to do so, we take advantage of the Hungarian algorithm to find the minimum of the matrix composed of all of the different rearrangements of column parcels. 
