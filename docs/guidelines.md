# Guidelines & Resources

This page collects resources that are commonly consulted when mapping, bead typing, and parameterising molecules for Martini 3. They are gathered here for convenience — no need to hunt through papers and supplementary files mid-workflow.

<details>
<summary>Martini 3 interaction matrix</summary>

The Martini 3 non-bonded interaction matrix. Values indicate the interaction level between bead types (0 = hyper-attractive → 21 = super-repulsive).

![Martini 3 interaction matrix](docs/img/martini3_matrix.png)

<p class="info-footnote">Adapted from Souza, P.C.T. et al. (2021). Martini 3: a general purpose force field for coarse-grained molecular dynamics. <em>Nature Methods</em>, 18, 382–388. Supplementary Table 1. <a href="https://doi.org/10.1038/s41592-021-01098-3" target="_blank">doi:10.1038/s41592-021-01098-3</a></p>

<details class="nested">
<summary>Label corrections</summary>

Interaction level corrections applied on top of the base matrix when beads carry special labels. Adapted from Souza et al. (2021) *Nature Methods* SI, Tables 11–16. [doi:10.1038/s41592-021-01098-3](https://doi.org/10.1038/s41592-021-01098-3)

**Hydrogen bonding (d/a)**

| Label pair | RR/RS | SS/ST/RT | TT |
|---|---|---|---|
| d-d / a-a † | +2 | +1 | +1 |
| d-a | 0 | −1 | −3 |
| W-d / W-a | +1 | +1 | +1 |

† New level cannot be higher than 12.

**Electron polarizability (v/e)**

| Label pair | RR/RS | SS/ST/RT | TT |
|---|---|---|---|
| v-v / e-e † | +1 | +1 | +1 |
| e-v | −1 | −2 | −3 |
| e/v – others | −1 | −1 | −1 |

† New level cannot be higher than 14.

**Positive and negative (p/n)**

| Label pair | RR/RS | SS/ST/RT | TT |
|---|---|---|---|
| p-p / n-n † | +3 | +2 | +1 |
| p-n | 0 | 0 | 0 |
| p–C/X beads | −1 | −1 | −1 |
| n–P/N beads | −1 | −1 | −1 |

† New level cannot be higher than 12.

**Cross-interactions between labels**

| Label pair | RR/RS | SS/ST/RT | TT |
|---|---|---|---|
| d-e / a-v | −1 | −1 | −1 |
| d-v / a-e | +2 | +1 | +1 |
| p-d / n-a / n-e / p-v | +2 | +1 | +1 |
| p-a / n-d / p-e / n-v | −1 | −2 | −3 |

**Self-interaction (h/r)**

| Label pair | RR/RS | SS/ST/RT | TT |
|---|---|---|---|
| h-h | −1 | −1 | −1 |
| r-r | +1 | +1 | +1 |
| h-r | 0 | 0 | 0 |
| h/r – others | 0 | 0 | 0 |

**Partial charge (q)**

| Label pair | RR/RS | SS/ST/RT | TT |
|---|---|---|---|
| q-q | −3 | −2 | −1 |
| q–W/Q | −3 | −3 | −3 |
| q–others ‡ | −2 | −2 | −2 |

‡ D beads are not included.

</details>

</details>

<details>
<summary>General mapping rules</summary>

The flowchart below summarises the recommended mapping procedure for Martini 3 molecules.

![Martini 3 mapping rules flowchart](docs/img/mapping_rules.png)

<p class="info-footnote">Adapted from Souza, P.C.T. et al. (2021). Martini 3: a general purpose force field for coarse-grained molecular dynamics. <em>Nature Methods</em>, 18, 382–388. Supplementary Information.</p>

</details>

<details>
<summary>Bead type assignments</summary>

<details class="nested">
<summary>Neutral molecules/fragments</summary>

| Bead | Chemical Group | Bead | Chemical Group |
|------|---------------|------|---------------|
| C1 | linear alkane | N6d | primary amine |
| C2 | branched alkane | N4 | secondary amine |
| C3 | cyclic alkane | N3a | tertiary amine |
| C4 | alkene | P1 | alcohol |
| C5 | dienes | P2 | enol |
| C5 | aromatics | P2 | carboxylic acid |
| C5e | aromatic (no hydrogens) | P2 | hemiacetal/hemiketal |
| C6r | alkynes | P4 | diol |
| C6 | thiol | P5 | primary amide |
| C6 | sulfide | P3 | secondary amide |
| N1 | methyl pyrrole | P3a | tertiary amide |
| N2 | primary imine | P6 | sulfoxide |
| N1a | secondary imine | P6 | amino acid |
| N2ah | acyl chloride | W | water |
| N2a | nitrile | X1 | primary organoiodine |
| N3r/N3a† | ether | X2 | primary organobromide |
| N4a | acetal/ketal | X3 | primary organochloride |
| N4a | ester | X4e | primary organofluorine |
| N5a | ketone | X3h | dichloroethane |
| N5a | aldehyde | X2 | trichloromethane |
| N6 | phenol | X1 | tetrachloromethane |

† Both options could work for ethers depending on the fragments attached.

<p class="info-footnote">Adapted from Souza, P.C.T. et al. (2021). Martini 3. <em>Nature Methods</em>, 18, 382–388. Supplementary Table 24.</p>

</details>

<details class="nested">
<summary>Charged fragments / ions</summary>

| Bead (charge) | Cation | Bead (charge) | Anion |
|---------------|--------|---------------|-------|
| TQ5 (+1) | sodium | TQ5 (−1) | chloride |
| SQ4 (+1) | potassium | SQ4 (−1) | bromide |
| TQ5p (+1) | ammonium | SQ2 (−1) | iodine |
| Q4p/SQ4p (+1) | alkyl ammonium | Q2 (−1) | tetrafluoroborate |
| Q3p/SQ3p (+1) | alkyl methyl ammonium | Q1 (−1) | hexafluorophosphate |
| Q2p/SQ2p (+1) | alkyl dimethyl ammonium | SQ1 (−1) | thiocyanate |
| Q2 (+1) | tetramethyl ammonium | SQ3 (−1) | nitrate |
| Q1 (+1) | choline (lipid head) | Q2 (−1) | perchlorate |
| Q1 (+1) | tetramethyl phosphonium | Q5n/SQ5n (−1) | carboxylate |
| SQ3p (+1) | guanidinium | Q4n/SQ4n (−1) | sulfonate |
| TD (+2) | magnesium | Q5 (−1) | phosphate (lipid head) |
| SD/TD‡ (+2) | calcium | D (−2) | phosphate |

‡ Depending on the application, TD could be useful.

<p class="info-footnote">Adapted from Souza, P.C.T. et al. (2021). Martini 3. <em>Nature Methods</em>, 18, 382–388. Supplementary Information.</p>

</details>

<details class="nested">
<summary>Protein sidechains</summary>

Martini 3 bead type assignments for the 20 standard amino acid sidechains.

![Martini 3 protein sidechain bead types](docs/img/protein_sidechains.png)

<p class="info-footnote">Reprinted from Souza, P.C.T., Borges-Araújo, L., Brasnett, C. et al. GōMartini 3: From large conformational changes in proteins to environmental bias corrections. <em>Nat Commun</em> 16, 4051 (2025). <a href="https://doi.org/10.1038/s41467-025-58719-0" target="_blank">doi:10.1038/s41467-025-58719-0</a>. Licensed under <a href="https://creativecommons.org/licenses/by-nc-nd/4.0/" target="_blank">CC BY-NC-ND 4.0</a>.</p>

</details>

<details class="nested">
<summary>Small molecule building blocks</summary>

Bead type assignments for a wide range of small molecule fragments ([Alessandri et al., 2022](https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/adts.202100391), Apache 2.0). Source: [M3-Bible on GitHub](https://github.com/Martini-Force-Field-Initiative/M3-Bible).

<iframe src="docs/assets/building_block_table.pdf" width="100%" height="700" style="border:1px solid var(--border-input);border-radius:6px;display:block;"></iframe>

</details>

<details class="nested">
<summary>The Martini Bible</summary>

A comprehensive community-maintained reference for Martini 3 bead type assignments and parameterisation guidelines, covering a wide range of chemical building blocks. Source: [M3-Bible on GitHub](https://github.com/Martini-Force-Field-Initiative/M3-Bible).

<iframe src="docs/assets/MartiniBible.pdf" width="100%" height="700" style="border:1px solid var(--border-input);border-radius:6px;display:block;"></iframe>

<p class="info-footnote">Maintained by Riccardo Alessandri and Fabian Grünewald. Last accessed: 18 June 2026.</p>

</details>

</details>

<details>
<summary>Key papers</summary>

The Martini 3 force field and the key papers that have parameterised large collections of molecules.

<ul class="ref-list">
<li>
Souza, P.C.T. et al. <em>Nat. Methods</em> 2021, 18, 382–388.
<a href="https://doi.org/10.1038/s41592-021-01098-3" target="_blank">Read →</a>
</li>
<li>
Alessandri, R. et al. <em>Adv. Theory Simul.</em> 2022, 5, 2100391.
<a href="https://doi.org/10.1002/adts.202100391" target="_blank">Read →</a>
&nbsp;·&nbsp;
<a href="https://github.com/Martini-Force-Field-Initiative/M3-Small-Molecules/tree/main" target="_blank">GitHub →</a>
</li>
<li>
Grünewald, F. et al. <em>J. Chem. Theory Comput.</em> 2022, 18, 7555.
<a href="https://pubs.acs.org/doi/10.1021/acs.jctc.2c00757" target="_blank">Read →</a>
</li>
<li>
Vainikka, P. et al. <em>ACS Sustainable Chem. Eng.</em> 2021, 9, 17338.
<a href="https://pubs.acs.org/doi/10.1021/acssuschemeng.1c06521" target="_blank">Read →</a>
</li>
<li>
Kjølbye, L.R. et al. <em>ChemRxiv</em> 2024.
<a href="https://doi.org/10.26434/chemrxiv-2024-bf4n8" target="_blank">Read →</a>
</li>
<li>
Pedersen, K.P. et al. <em>ACS Cent. Sci.</em> 2025, 11, 1598–1610.
<a href="https://pubs.acs.org/doi/10.1021/acscentsci.5c00755" target="_blank">Read →</a>
</li>
<li>
Martini Bible — community-maintained parameterisation reference.
<a href="https://github.com/Martini-Force-Field-Initiative/M3-Bible/tree/main" target="_blank">GitHub →</a>
</li>
</ul>

</details>

