# PYFURCATION: Create 3D Branching Respiratory Tree Structures

## Getting Started
I actually don't know exactly how you would start using this on your local device. Copy this repository to your device?

## Dependencies 
The PYFURCATION Library utilizes the following python packages:
* Open3d (version >= 0.18.0)
* NumPy (version >= 1.26.4)

See Open3d [GitHub](https://github.com/isl-org/Open3D) for Python installation guide and [Documentation](https://www.open3d.org/docs/release/) for more information on the library.

## About the Model: Defining the Geometry

### Bifurcation Units
<p align="center">
  <img src="/supporting_documents/pointcloud_unit.png" width = 500>
</p>

The PYFURCATION Library generates a series of symmetric Bifurcation Units (Single Inlet, Double Outlet unit trees) using the methodology set out in [[1]](#1) and Appendix A of [[2]](#2). 
Each Bifurcation Unit is defined by the following parameters:
* R<sub>p</sub> - Parent Radius; the radius of the parent section (inlet) of the bifurcation unit
* L<sub>p</sub> - Parent Length; the length of the parent section of the bifurcation unit
* R<sub>d</sub> - Daughter Radius; the radius of the daughter sections (each outlet) of the bifurcation unit
* L<sub>d</sub> - Daughter Length; the length of the daughter sections of the bifurcation unit
* R<sub>o</sub> - Outer Radius; the external radius of the bending/merging pipe section
* &Phi;<sub>&Beta;</sub> - Branching Angle; the angle between the parent and daughter pipe sections
* &Delta;&gamma; - Carina Angle; the angle of the smooth carina connecting the daughter pipe sections

See Figures Below. Note: all other parameters shown in the images are calculated from the above.

![](/supporting_documents/figure1.jpg)  |  ![](/supporting_documents/figure2.jpg)
:-------------------------:|:-------------------------:
**Figure 1**               |  **Figure 2**

If you are feeling adventurous, there are bounds on these values for which a valid geometry can be generated, based on the equations found in Appendix A of [[2]](#2). It's a pretty interesting problem,
but would require *a lot* of algebra to fully determine. There might be more information in [[1]](#1) too, that I have yet to fully dive into.

See code documentation in the [bifurcation unit generation directory](/bifurcating_unit_modules) for code workflow.

In prac

### Respiratory Tree


## How To: Creating a Bifurcating Respiratory Tree Surface Mesh

## References
