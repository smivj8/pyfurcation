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
  <img src="/supporting_documents/pointcloud_unit2.png" width = 500>
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

To generate a complete mesh, an additional parameter, n<sub>Cont. Outlets</sub>, is necessary. This parameter defines the number of outlets (0, 1, or 2) that will be left "open," in that they are not closed off by triangles. This parameter is important for connecting successive bifurcation units in a complete tree, and has no utility in the overall geometry generation of a single unit (and in which case is should be set to 0).

Additional parameters, such as pointcloud density or Poisson Mesh characteristic parameters, can be manipulated within the code (they are left as optional arguments in their corresponding functions) if necessary, although it is not recommended. I have found that the provided values generate adequate meshes in a short time. Any push/pull requests that alter these values, without good justification, will be rejected. Probably. I'm not sure how much I will be keeping track of this in the future. 

#### Example Parameters

At the start of this section, a sample bifurcation unit, with the pointcloud, was shown. The parameters used were: 
* R<sub>p</sub> = 420
* L<sub>p</sub> = 2010
* R<sub>d</sub> = 335
* L<sub>d</sub> = 2010
* R<sub>o</sub> = 1675
* &Phi;<sub>&Beta;</sub> = 35&deg;
* &Delta;&gamma; = 3.6&deg;
* n<sub>Cont. Outlets</sub> = 0

The radius and length parameters have units &mu;m, although the unit generation should work so long as they are within the necessary relative bounds (of which I don't know). 

### Respiratory Tree


## How To: Creating a Bifurcating Respiratory Tree Surface Mesh

## References
