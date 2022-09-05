# DynamicSoapfilmsWithEvolvingThickness

This is an implementation of the paper **“A Model for Soap Film Dynamics with Evolving Thickness
”** Ishida\*, Synak\*, Narita, Hachisuka, Wojtan, Transactions on Graphics (SIGGRAPH 2020).

[[Project site]][P]  
<a href="https://sadashigeishida.bitbucket.io/soapfilm_with_thickness">  <img src="https://sadashigeishida.bitbucket.io/soapfilm_with_thickness/teaser_white_background.jpeg" height="150px"> </a>  
[[Youtube video]][Y]  
<a href="https://www.youtube.com/watch?v=Pr1zibwxAKU"><img src="https://i.ytimg.com/vi/Pr1zibwxAKU/0.jpg" width="200px"></a>

[Y]:https://www.youtube.com/watch?v=Pr1zibwxAKU
[P]:https://sadashigeishida.bitbucket.io/soapfilm_with_thickness  
Author: Sadashige Ishida, Fumiya Narita, and Peter Synak
Lisence: MPL-2.0

## Basic Usage
Running the executables without command line arguments will display the usage. Data files are located in the fes_assets folder.

[KEY SETTING]  
Space: Turn on/off the clock.  
s: Proceed one time step.   
i: Change simulation mode 0. Eulerian flow only 1. Lagrangian deformation only 2. Both of them  

C: Add random velocities using curl noise  
P: Add random thickness using Perlin noise

b: Burst the bubble with the thinnest film among.    
B: Toggle thickness-dependent bubble burst.  
e: Toggle evaporation.
 
\/: Switch the equation between our soap film model and the incompressible Euler equation  
D: Turn on/off the scene specific update.   

m: Change rendering mode.  
v: Toggle visualization of Eulerian velocities on vertices  
t: Toggle visualization of Eulerian velocities on triangles   

o: Save the state as files containing information of mesh, labels, and constrained vertices.  
O: Save the state as above, but with ghost vertices and faces.  
and etc.

## Dependencies
This program is built by standard procedures using CMAKE (http://www.cmake.org).
The following external libraries are required:   
Eigen (http://eigen.tuxfamily.org)  
LAPACK (http://www.netlib.org/lapack/)  
libigl (http://libigl.github.io/libigl/)  
OpenGL (https://www.opengl.org/)  
GLEW (http://glew.sourceforge.net/) for non-mac OS
