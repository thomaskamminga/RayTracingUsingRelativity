# RayTracingUsingRelativity
 C++ project that generates images of black holes from space-time curved by gravity in the Schwarzschild and Kerr metric. This code was for my bachelor thesis "Imaging a Kerr Black Hole" available at repository.tudelft.nl . The thesis explains the working of the raytracing algorithm and the mathematics of light around a black hole. The project is an extention of Peter Shirley's course "Ray Tracing in One Weekend".
 
 ## Building
  The project relies on the boost header library to build. This library was to big to include. Download the latest version from boost.org and put the "boost" folder in the "RayTracingUsingRelativity" folder. Use Visual Studio to open the project and build to release. If the Kerr metric is selected in geodesicODE function, the building proces will take several minutes.

## Running
To run build to relase and run "RayTracingUsingRelativity.exe" in RayTracingUsingRelativity/release. After running the image will appear at RayTracingUsingRelativity/release/images/image.png . The image name can be changed in the source code.

# Resulting images
![Schwarzschild black hole with milkyway background](https://github.com/thomaskamminga/RayTracingUsingRelativity/blob/main/Release/images/front_page.png)
![Schwarzschild black hole with accretion disk](https://github.com/thomaskamminga/RayTracingUsingRelativity)
![Kerr black hole with accretion disk](https://github.com/thomaskamminga/RayTracingUsingRelativity/blob/main/Release/images/Kerr_correct_ring.png)
