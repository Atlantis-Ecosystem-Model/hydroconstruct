# hydroconstruct

HydroConstruct creates the final format of the temperature, salinity and transport files required by Atlantis. It will read in files in a number of formats but its generally transport across each face as well as temperature, salinity and vertical flux per box.

The majority of the data read in has the level in the water column in the opposite order to how Atlantis reads them in, so HydroConstruct reverses this order.

It mainly changes the transport data into data per layer per box and handles hyperdiffusion.


