# MultiscaleTopOpt
A 3D multsicale topology optimization code using surrogate models of lattice microscale response

Background

This Matlab code is contained in Appendix D of "Simple, accurate surrogate models of the elastic response of three-dimensional open truss micro-architectures with applications to multiscale topology design" by Watts et al., appearing in 2019 in Structural and Multidisciplinary Optimization (full reference below). This article derives polynomial surrogate models of the homogenized elasticity tensor for each of the Isotruss, octet truss, and ORC truss architectures as a function of the relative density and material properties of their constituent solid. 

This code shows how these surrogate models can be incorporated into a standard topology optimization code, converting it into a multiscale topology design tool at negligible additional runtime cost (although the microstructure is restricted to one for which a surrogate model has been generated). Ole Sigmund’s two-dimensional 99-line Matlab topology optimization code  (Structural and Multidisciplinary Optimization, vol. 21, pp. 120-7, 2001) is used as a starting point, and it is modified as described in Section 5 of the 2019 article. Quoting from that section, “[w]e echo Sigmund’s comments that this code is written for its pedagogical attributes, not its speed. Our goal in providing the modified code is to make plain what changes are needed to existing topology optimization codes to design with surrogate models of the micro-architecture, not to streamline or optimize the code. Many techniques, including those mentioned by Sigmund, exist towards these ends.”
 
Use

This code is used to design light, stiff structures via topology optimization using a density parameterization. A uniform-density structure is initialized, and its performance under a prescribed loading is analyzed using finite element analysis. The design is updated iteratively using an optimality criterion method with filtered sensitivities, such that the design maximizes its stiffness subject to a resource constraint. Additional details on the topology optimization problem and this solution approach are contained in the two articles cited in the previous section, and in their references.

The code is executed in Matlab by calling

    top(nelx, nely, nelz, volfrac, rmin, truss, Es, vs, minVF, maxVF, maxit) 
 
    where: nelx, nely, and nelz are the number (nel* >= 1) of unit-size square or cube finite elements in the mesh of the design domain
           volfrac is the upper bound on the volume fraction (0.0 < volfrac < 1.0) of the domain that may be solid
           rmin is the radius of the sensitivity filter, in units of elements (rmin >= 1)
           truss is the type of microstructure (‘simp’, ‘bound’, ‘iso’, ‘octet’, ‘orc’)
           Es is the constituent material’s Young’s modulus (Es > 0.0)
           vs is the constituent material’s Poisson’s ratio (0.0 < vs < 0.5 typically)
           minVF and maxVF are respectively lower and upper bounds on the local element density (0.0 < minVF < maxVF < 1.0)
           maxit is the maximum number of optimization iterations allowed for convergence (maxit > 1)

The five designs (a)-(e) of Example 6.1 of the 2019 article can be obtained respectively by calling

    top(1, 60, 20, 0.50, 1.5, ‘simp’,  1.0, 0.3, 0.001, 1.0, 150)
    top(1, 60, 20, 0.50, 1.5, ‘bound’, 1.0, 0.3, 0.001, 1.0, 150) 
    top(1, 60, 20, 0.50, 1.5, ‘iso’,   1.0, 0.3, 0.001, 1.0, 150) 
    top(1, 60, 20, 0.50, 1.5, ‘octet’, 1.0, 0.3, 0.001, 1.0, 150) 
    top(1, 60, 20, 0.50, 1.5, ‘orc’,   1.0, 0.3, 0.001, 1.0, 150) 

Modification of the boundary conditions imposed on the design domain, e.g. to obtain the designs of Example 6.2 or 6.3 of the 2019 article, requires changing the assigned loaded nodes (for non-homogeneous natural boundary conditions) and fixed nodes (for homogeneous essential boundary conditions). See lines 79 – 94. 

Future Development & Contributing

This code is being provided as a convenience to readers of the 2019 article who wish to experiment with the code of its Appendix D. Since this code matches a published version, we intend to keep it static and do not anticipate development of this code, except for bug fixes. Additionally, we are not soliciting contributions, although we encourage others to use this code as a basis for their own projects. 

If you find this code useful, we would appreciate a citation of our work. Bibliographical details are below and will be updated when a volume and page numbers are assigned.

Seth Watts, William Arrighi, Jun Kudo, Daniel A. Tortorelli, and Daniel A. White, “Simple, accurate surrogate models of the elastic response of three-dimensional open truss micro-architectures with applications to multiscale topology design," Structural and Multidisciplinary Optimization, 2019, DOI 10.1007/s00158-019-02297-5

LLNL-CODE-757968
