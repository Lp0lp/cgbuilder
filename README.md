CGBuilder — Visual tool to build CG molecule models
===================================================

This project is a fork of [Original CGBuilder](https://github.com/jbarnoud/cgbuilder), with additional
functionality for complex bead placement, and Shaker format export.

Double-bond rendering notes
---------------------------

- For consistent display of double and triple bonds prefer input formats that include bond-order metadata such as SDF (.sdf, .sd) or MOL2 (.mol2). PDB and GRO files often lack bond-order information and may not show multiple bonds correctly.
