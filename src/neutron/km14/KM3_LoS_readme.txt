Dear Michal,

thanks for letting me know about your progress.

In our modelling of the neutron spectrum components, we take into the account the actual line of sight.


_KM3.los is a txt file. Each line describes a cell of the line of sight. For each line you will find 8 numbers,



x,y,z,C,V,u,v,w



with the following meaning:



x,y,z are the x,y,z coordinates of the cell, in meters

V is the volume of the cell, in meter**3

C is related to the emission solid angle omega = (4*pi*C/V). Omega is the solid angle associated to each cell. Neutrons emitted from the cell and within Omega will reach the detector, the others will not.

u,v,w are the x,y,z components of the versor that describes the emission direction that neutrons produced from a cell must have to reach the detector.

In our Monte Carlo calculation with GENESIS, we randomly extract a cell from the line of sight. We then find the TRANSP distribution that is associated to that cell and use that to calculate the neutron emission from the cell, using the densities, the solid angle etc. as weighting factors. The total spectrum that impinges on the detector will be the weighted sum of the spectrum from all the cells belonging to the LoS. 

Regards 

                      Massimo
					  
Thank you Massimo, 
this will help me very much. 
I am trying to map the coordinates with PSIR and PSIZ coordinates from EFIT. 
From the values I infer that
x corresponds to local coordinate in the tangential toroidal direction, 
y correponds PSIR in meters,
z corresponds to PSIZ in meters.

Is that correct? 

Thank you,
Michał 

Hi Michal,

this is a set of cartesian right-handed orthogonal coordinates. The poloidal plane is in the x-z plane. The toroidal plane is in the x-y plane.
 
The x,y,z definition of the LoS is NOT discharge, nor time dependent, as it depends on the geometry of the LoS, which does not depend on the plasma. 
The (rho,theta,phi) to (x,y,z) conversion is discharge and time dependent.



             Massimo