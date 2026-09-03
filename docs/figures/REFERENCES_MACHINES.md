# References for the machine figure (svg/01_all_systems_complete.svg)

One section per machine, one entry per drawn component, in the order the parts
are listed in `scripts/machine_config.py`. Each machine is a **composite of
separately deposited structures**; no entry here is a structure of a whole
secretion system, because for most of these no such structure exists.

**Provenance note.** Every citation below was retrieved on 2026-09-03 from the
RCSB PDB entry record for the exact accession the figure draws, and the author
list, journal, year and DOI are as deposited. Nothing here was written from
memory. Companion file: `REFERENCES.md`, which covers the substrate figure.

Components marked **[AF]** are AlphaFold models, drawn with a dashed outline in
the figure. For those the *fold* is predicted but the *copy number* is evidence;
the evidence is named in the entry and recorded per part in
`scripts/build_af_parts.py`.

---

## T1SS — Haemolysin A secretion system (*E. coli*)

**TolC** — PDB **1EK9** (X-ray, 2.1 Å)  
outer-membrane exit duct; the OPM entry sets the OM plane.

- Koronakis V, Sharff A, Koronakis E, Luisi B, Hughes C (2000). Crystal structure of the bacterial membrane protein TolC central to multidrug efflux and protein export. *Nature* 405:914. doi:10.1038/35016007

**HlyD hexamer** — **[AF]** AlphaFold model of HlyD, UniProt **P09986** (*E. coli*), 6 copies  
periplasmic adaptor. Its alpha-hairpin and lipoyl body are disordered in 7SGR, so the fold is AlphaFold and only the copy number (6) comes from the structure.

- Cite the AlphaFold method and database entries listed under *Software and databases* below, plus UniProt P09986.

**HlyB-HlyD inner-membrane complex** — PDB **7SGR** (cryo-EM, 2.9 Å)  
the ABC transporter and the resolved fragments of its adaptor.

- Zhao H, Lee J, Chen J (2022). The hemolysin A secretion system is a multi-engine pump containing three ABC transporters. *Cell* 185:3329. doi:10.1016/j.cell.2022.07.017

---

## T2SS — Pullulanase secretion system (*K. pneumoniae* / ETEC)

**GspK-GspI-GspJ tip** — PDB **3CI0** (X-ray, 2.2 Å)  
pseudopilus tip heterotrimer.

- Korotkov KV, Hol WG (2008). Structure of the GspK-GspI-GspJ complex from the enterotoxigenic Escherichia coli type 2 secretion system. *Nat Struct Mol Biol* 15:462. doi:10.1038/nsmb.1426

**Secretin GspD + PulS** — PDB **6HCG** (cryo-EM, 4.3 Å)  
C15 outer-membrane channel; sets the OM plane.

- Chernyatina AA, Low HH (2019). Core architecture of a bacterial type II secretion system. *Nat Commun* 10:5437. doi:10.1038/s41467-019-13301-3

**Pseudopilus GspG/PulG** — PDB **5WDA** (cryo-EM, 5.0 Å)  
helical filament, truncated to 15 of 25 subunits to fit the periplasm.

- Lopez-Castilla A, Thomassin JL, Bardiaux B, Zheng W, Nivaskumar M, Yu X et al. (2017). Structure of the calcium-dependent type 2 secretion pseudopilus. *Nat Microbiol* 2:1686. doi:10.1038/s41564-017-0041-2

**GspL platform (C6)** — **[AF]** AlphaFold model of GspL/PulL, UniProt **P15751** (*K. pneumoniae*), 6 copies  
inner-membrane platform; TM helices unresolved in every deposit.

- Cite the AlphaFold method and database entries listed under *Software and databases* below, plus UniProt P15751.

**GspM platform (C6)** — **[AF]** AlphaFold model of GspM/PulM, UniProt **P15752** (*K. pneumoniae*), 6 copies  
fold rescued by the GspL-GspM multimer job (mean pLDDT 63.4 to 75.1); the predicted interface scored ipTM 0.43 and is NOT used.

- Cite the AlphaFold method and database entries listed under *Software and databases* below, plus UniProt P15752.

**GspF dimer** — **[AF]** AlphaFold model of GspF/PulF, UniProt **P15745** (*K. pneumoniae*), 2 copies  
no membrane-domain structure exists anywhere in the PDB.

- Cite the AlphaFold method and database entries listed under *Software and databases* below, plus UniProt P15745.

**GspE ATPase + GspL cytoplasmic domain** — PDB **4PHT** (X-ray, 2.83 Å)  
3 copies in the deposit, not a closed hexamer.

- Lu C, Korotkov KV, Hol WG (2014). Crystal structure of the full-length ATPase GspE from the Vibrio vulnificus type II secretion system in complex with the cytoplasmic domain of GspL. *J Struct Biol* 187:223. doi:10.1016/j.jsb.2014.07.006

---

## T3SS — SPI-1 injectisome (*S.* Typhimurium)

**SipD needle tip** — PDB **7RYE** (cryo-EM, 3.9 Å)  
needle filament-tip complex.

- Guo EZ, Galan JE (2021). Cryo-EM structure of the needle filament tip complex of the Salmonella type III secretion injectisome. *Proc Natl Acad Sci U S A* 118. doi:10.1073/pnas.2114552118

**Needle complex** — PDB **7AH9** (cryo-EM, 3.3 Å)  
substrate-engaged basal body spanning both membranes; the source of the 220 Å periplasm width.

- Miletic S, Fahrenkamp D, Goessweiner-Mohr N, Wald J, Pantel M, Vesper O et al. (2021). Substrate-engaged type III secretion system structures reveal gating mechanism for unfolded protein translocation. *Nat Commun* 12:1546. doi:10.1038/s41467-021-21143-1

**Secretin InvG (membrane frame)** — PDB **6DV3** (cryo-EM, 4.1 Å)  
open-gate secretin, used as the OPM reference that places 7AH9's OM plane.

- Hu J, Worrall LJ, Hong C, Vuckovic M, Atkinson CE, Caveney N et al. (2018). Cryo-EM analysis of the T3S injectisome reveals the structure of the needle and open secretin. *Nat Commun* 9:3840. doi:10.1038/s41467-018-06298-8

**SctV/InvA export apparatus (C9)** — **[AF]** AlphaFold model of InvA/SctV, UniProt **P0A1I3** (*S.* Typhimurium), 9 copies  
AlphaFold supplies the 8 TM helices no experimental entry resolves.

- Cite the AlphaFold method and database entries listed under *Software and databases* below, plus UniProt P0A1I3.

**Sorting platform pods** — **[AF]** AlphaFold model of SpaO, UniProt **P40699** (*S.* Typhimurium), 6 copies  
drawn as six bare SpaO monomers; the assembly has no atomic structure.

- Cite the AlphaFold method and database entries listed under *Software and databases* below, plus UniProt P40699.

**ATPase SctN/EscN** — PDB **6NJP** (cryo-EM, 3.29 Å)  
isolated EPEC hexamer with its central stalk.

- Majewski DD, Worrall LJ, Hong C, Atkinson CE, Vuckovic M, Watanabe N et al. (2019). Cryo-EM structure of the homohexameric T3SS ATPase-central stalk complex reveals rotary ATPase-like asymmetry. *Nat Commun* 10:626. doi:10.1038/s41467-019-08477-7

---

## T4SS — R388 conjugative system (*E. coli*)

**Pilus TrwL** — PDB **8S6H** (cryo-EM, 3.39 Å)  
conjugative pilus helical polymer.

- Vadakkepat AK, Xue S, Redzej A, Smith TK, Ho BT, Waksman G (2024). Cryo-EM structure of the R388 plasmid conjugative pilus reveals a helical polymer characterized by an unusual pilin/phospholipid binary complex. *Structure* 32:1335. doi:10.1016/j.str.2024.06.009

**OM core complex** — PDB **8RT6** (cryo-EM, 3.18 Å)  
full-length outer-membrane core complex.

- Mace K, Waksman G (2024). Cryo-EM structure of a conjugative type IV secretion system suggests a molecular switch regulating pilus biogenesis. *EMBO J* 43:3287. doi:10.1038/s44318-024-00135-z

**O-layer (membrane frame)** — PDB **7O3J** (cryo-EM, 2.6 Å)  
OPM reference that places 8RT6's OM plane.

- Mace K, Vadakkepat AK, Redzej A, Lukoyanova N, Oomen C, Braun N et al. (2022). Cryo-EM structure of a type IV secretion system. *Nature* 607:191. doi:10.1038/s41586-022-04859-y

**VirB11 ATPase** — PDB **1G6O** (X-ray, 2.5 Å)  
*H. pylori* HP0525 hexamer; absent from every near-complete T4SS model, so taken from a separate entry. Biological assembly, not the asymmetric unit.

- Yeo HJ, Savvides SN, Herr AB, Lanka E, Waksman G (2000). Crystal structure of the hexameric traffic ATPase of the Helicobacter pylori type IV secretion system. *Mol Cell* 6:1461. doi:10.1016/S1097-2765(00)00142-8

**VirD4 coupling protein** — PDB **1E9R** (X-ray, 2.4 Å)  
TrwB hexamer, likewise from a separate entry.

- Gomis-Rueth FX, Moncalian G, Perez-Luque R, Gonzalez A, Cabezon E, De La Cruz F et al. (2001). The Bacterial Conjugation Protein Trwb Resembles Ring Helicases and F1-ATPase. *Nature* 409:637. doi:10.1038/35054586

**Stalk + arches + IMC** — PDB **8RTD** (cryo-EM, 4.33 Å)  
inner-membrane half of the assembled machine.

- Mace K, Waksman G (2024). Cryo-EM structure of a conjugative type IV secretion system suggests a molecular switch regulating pilus biogenesis. *EMBO J* 43:3287. doi:10.1038/s44318-024-00135-z

**IMC protomer (membrane frame)** — PDB **7OIU** (cryo-EM, 3.7 Å)  
OPM reference that places 8RTD's IM plane.

- Mace K, Vadakkepat AK, Redzej A, Lukoyanova N, Oomen C, Braun N et al. (2022). Cryo-EM structure of a type IV secretion system. *Nature* 607:191. doi:10.1038/s41586-022-04859-y

---

## T5SS — EspP autotransporter (*E. coli*)

**BAM + EspP barrel** — PDB **8BO2** (cryo-EM, 3.1 Å)  
engineered BamA-S425C/EspP-S1299C disulfide trap caught mid-insertion; a legend note is warranted.

- Shen C, Chang S, Luo Q, Chan KC, Zhang Z, Luo B et al. (2023). Structural basis of BAM-mediated outer membrane beta-barrel protein assembly. *Nature* 617:185. doi:10.1038/s41586-023-05988-8

**EspP passenger** — PDB **3SZE** (X-ray, 2.5 Å)  
the secreted beta-helical passenger domain.

- Khan S, Mian HS, Sandercock LE, Chirgadze NY, Pai EF (2011). Crystal Structure of the Passenger Domain of the Escherichia coli Autotransporter EspP. *J Mol Biol* 413:985. doi:10.1016/j.jmb.2011.09.028

**SurA chaperone** — PDB **1M5Y** (X-ray, 3.0 Å)  
fills the real 117 Å periplasmic gap with what occupies it; NOT a physical conduit.

- Bitto E, McKay DB (2002). Crystallographic Structure of SurA, a Molecular Chaperone that  
Facilitates Folding of Outer Membrane Porins. *Structure* 10:1489. doi:10.1016/S0969-2126(02)00877-8

**SecYEG translocon** — PDB **5AWW** (X-ray, 2.724 Å)  
*T. thermophilus*, chosen because no *E. coli* SecYEG is cleanly OPM-aligned.

- Tanaka Y, Sugano Y, Takemoto M, Mori T, Furukawa A, Kusakizako T et al. (2015). Crystal Structures of SecYEG in Lipidic Cubic Phase Elucidate a Precise Resting and a Peptide-Bound State. *Cell Rep* 13:1561. doi:10.1016/j.celrep.2015.10.025

---

## T6SS — *E. coli* / *V. cholerae* / *P. aeruginosa* composite

**TssJM membrane core** — PDB **6HS7** (cryo-EM, 4.6 Å)  
TssJ x15 + TssM x10; TssM modelled only from residue 569, so its TM helices are absent.

- Rapisarda C, Cherrak Y, Kooger R, Schmidt V, Pellarin R, Logger L et al. (2019). In situand high-resolution cryo-EM structure of a bacterial type VI secretion system membrane complex. *EMBO J* 38. doi:10.15252/embj.2018100886

**TssM N-terminal half (C10)** — **[AF]** AlphaFold model of TssM/VasK, UniProt **B7LFU0** (*E. coli* 55989), 10 copies  
the weakest element in the figure: TM helices predict at ~65 pLDDT, drawn at a lowered cutoff of 55.

- Cite the AlphaFold method and database entries listed under *Software and databases* below, plus UniProt B7LFU0.

**TssL (C10)** — **[AF]** AlphaFold model of TssL, UniProt **D3GU40** (*E. coli* 042), 10 copies  
absent from every PDB entry; 10 copies inferred from the 3:2:2 TssJ:TssL:TssM stoichiometry on the C5 map.

- Cite the AlphaFold method and database entries listed under *Software and databases* below, plus UniProt D3GU40.

**VgrG spike** — PDB **4UHV** (X-ray, 2.0 Å)  
deposits 2 chains rather than the C3 trimer, so the drawn spike is thinner than the real one.

- Spinola-Amilibia M, Davo-Siguero I, Ruiz FM, Santillana E, Medrano FJ, Romero A (2016). The Structure of Vgrg1 from Pseudomonas Aeruginosa, the Needle Tip of the Bacterial Type Vi Secretion System. *Acta Crystallogr D Biol Crystallogr* 72:34. doi:10.1107/S205979831502149X

**Baseplate TssKFG** — PDB **6GIY** (cryo-EM, 4.3 Å)  
TssEFGK is split across 6GIY (KFG) and 6GJ1 (EFG); only 6GIY is drawn.

- Cherrak Y, Rapisarda C, Pellarin R, Bouvier G, Bardiaux B, Allain F et al. (2018). Biogenesis and structure of a type VI secretion baseplate. *Nat Microbiol* 3:1404. doi:10.1038/s41564-018-0260-1

**Sheath + Hcp tube** — PDB **5OJQ** (cryo-EM, 3.7 Å)  
extended (non-contractile) sheath-tube complex.

- Wang J, Brackmann M, Castano-Diez D, Kudryashev M, Goldie KN, Maier T et al. (2017). Cryo-EM structure of the extended type VI secretion system sheath-tube complex. *Nat Microbiol* 2:1507. doi:10.1038/s41564-017-0020-7

---

## Supporting citations (not drawn, but relied on)

These do not contribute coordinates to the figure. They supply a copy number, a
membrane position, or an architectural claim the figure makes.

**7ALW** — nonameric SctV cytoplasmic ring; sets the 9-fold copy number and ring width for the predicted SctV.

- Kuhlen L, Johnson S, Cao J, Deme JC, Lea SM (2021). Nonameric structures of the cytoplasmic domain of FlhA and SctV in the context of the full-length protein. *PLoS One* 16:e0252800. doi:10.1371/journal.pone.0252800

**7OSL** — second nonameric SctV-C structure, same purpose.

- Yuan B, Portaliou AG, Parakra R, Smit JH, Wald J, Li Y et al. (2021). Structural Dynamics of the Functional Nonameric Type III Translocase Export Gate. *J Mol Biol* 433:167188. doi:10.1016/j.jmb.2021.167188

**5NIK** — MacAB-TolC, the closest assembled analogue to the T1SS pump and the only entry in this set that OPM places with both bilayers, 275 Å apart.

- Fitzpatrick AWP, Llabres S, Neuberger A, Blaza JN, Bai XC, Okada U et al. (2017). Structure of the MacAB-TolC ABC-type tripartite multidrug efflux pump. *Nat Microbiol* 2:17070. doi:10.1038/nmicrobiol.2017.70

**5WQ7** — T2SS secretin GspD; the independent check that the secretin reaches ~203 Å across the periplasm.

- Yan Z, Yin M, Xu D, Zhu Y, Li X (2017). Structural insights into the secretin translocation channel in the type II secretion system. *Nat Struct Mol Biol* 24:177. doi:10.1038/nsmb.3350

**6GJ1** — the TssEFG half of the T6SS baseplate, not drawn but named because the baseplate is split across two entries.

- Cherrak Y, Rapisarda C, Pellarin R, Bouvier G, Bardiaux B, Allain F et al. (2018). Biogenesis and structure of a type VI secretion baseplate. *Nat Microbiol* 3:1404. doi:10.1038/s41564-018-0260-1

**EMD-8544** — in-situ Salmonella injectisome; the only evidence for where the
ATPase sits under the machine, and the source of the six-pod sorting-platform cage.

- Hu B, Lara-Tejero M, Kong Q, Galan JE, Liu J (2017). In situ molecular architecture of the Salmonella type III secretion machine. *Cell* 168:1065-1074.e10. doi:10.1016/j.cell.2017.02.022

**Sorting-platform stoichiometry** (1 OrgA + 4 SpaO + 2 OrgB per pod, six pods) —
quantitative superresolution counting in live cells.

- Zhang Y, Lara-Tejero M, Bewersdorf J, Galan JE (2017). Visualization and characterization of individual type III protein secretion machines in live bacteria. *Proc Natl Acad Sci USA* 114:6098-6103. doi:10.1073/pnas.1705823114

- Soto JE, Galan JE, Lara-Tejero M (2022). Assembly and architecture of the type III secretion sorting platform. *Proc Natl Acad Sci USA* 119:e2218010119. doi:10.1073/pnas.2218010119 — the assembly map. Deposited no coordinates, which is why the pods are drawn as bare SpaO monomers.

**EMD-4562** — in-situ subtomogram average of assembled TssJLM at 25 Å, the only map
of the complete T6SS membrane complex. Same paper as 6HS7.

- Rapisarda C, Cherrak Y, Kooger R, Schmidt V, Pellarin R, Logger L et al. (2019). In situand high-resolution cryo-EM structure of a bacterial type VI secretion system membrane complex. *EMBO J* 38. doi:10.15252/embj.2018100886

**T6SS 3:2:2 TssJ:TssL:TssM stoichiometry** — the basis for drawing 10 copies each of
TssL and the TssM N-terminal half on a C5 ring.

- Durand E, Nguyen VS, Zoued A, Logger L, Pehau-Arnaudet G, Aschtgen MS et al. (2015). Biogenesis and structure of a type VI secretion membrane core complex. *Nature* 523:555-60. doi:10.1038/nature14667

---

## Software and databases

- Lomize MA, Pogozheva ID, Joo H, Mosberg HI, Lomize AL (2012). OPM database and PPM web server: resources for positioning of proteins in membranes. *Nucleic Acids Res* 40:D370-6. doi:10.1093/nar/gkr703 — every membrane plane in the figure is an OPM bilayer position, not a hand-drawn line.
- Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583-589. doi:10.1038/s41586-021-03819-2
- Varadi M, Bertoni D, Magana P, Paramval U, Pidruchna I, Radhakrishnan M et al. (2024). AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. *Nucleic Acids Res* 52:D368-D375. doi:10.1093/nar/gkad1011 — source of the monomer models (v6).
- Abramson J, Adler J, Dunger J, Evans R, Green T, Pritzel A et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493-500. doi:10.1038/s41586-024-07487-w — the AlphaFold Server model behind the GspL-GspM job. Non-commercial use only; model seed fixed at 1.

- cellscape (github.com/jordisr/cellscape), MIT licence, last commit 2022-06-22.
  No accompanying paper; cite as software. The figure was produced with seven
  local fixes in `scripts/cellscape-shapely2-and-fixes.patch`, two of which are
  genuine upstream bugs affecting output geometry (an OPM view matrix with
  determinant -1, i.e. a reflection that reversed depth order and helix
  handedness, and a depth-slice loop that dropped the back-most 3 Å slab). A
  methods section should say the patched version was used.

Structure coordinates: RCSB PDB (rcsb.org). Predicted models: AlphaFold DB
(alphafold.ebi.ac.uk). Membrane placements: OPM (opm.phar.umich.edu).
