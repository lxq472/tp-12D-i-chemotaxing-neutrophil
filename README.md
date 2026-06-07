# Week 3 – Chemotaxing Neutrophil  

**Course:** Physics of Molecular Diseases  
**Subunit:** Receptor Adaptation in Chemotaxis of White Blood Cells

---

## Overview

This script implements the **receptor-mediated cell orientation model** of
Lin & Butcher (J. Immunol. 2008) for a chemotaxing neutrophil, including
homologous receptor desensitization, internalization and recycling.

It reproduces the single-gradient results of the paper (Fig. 1) and adds a
migration simulation showing how a cell with normal (desensitizable) vs
nondesensitizable receptors moves through a single ligand gradient.

---

## Biological Background

G-protein-coupled chemoattractant receptors signal transiently on ligand
binding, then are rapidly desensitized. The importance of desensitization
was long unclear, because nondesensitizable receptor mutants chemotax
efficiently in a single gradient. The model tests the hypothesis that
desensitization matters for navigation in COMPETING gradients (the next
question); this script establishes the single-gradient building block.

---

## Model

The cell is two receptor units (front / back) separated by `l = 10 um`.
Each unit obeys the receptor kinetics:

```
R + L  <kf, kr1>  LR*  --kdes-->  LRd  --ki-->  (LRi -> Ri instantly)
Ri  --kup-->  R                                  (recycling / up-regulation)
```

**ODEs (Eqs 1-3), per unit:**

```
dLR*/dt = kf*L*R - (kr1 + kdes)*LR*
dLRd/dt = kdes*LR* - ki*LRd
dR/dt   = kr1*LR* - kf*L*R + kup*Ri
```

**Conservation (Eq 4):**  `Rtot = R + Ri + LR* + LRd`  →  `Ri = Rtot - R - LR* - LRd`

**Orientation signal (Eq 10):**  `ΔLR* = LR*_front - LR*_back`

The cell orients up-gradient if `ΔLR* >= +threshold`, down-gradient if
`<= -threshold`, otherwise it migrates randomly (`threshold = 10` receptors).

---

## Parameters

| Symbol | Value | Meaning |
|--------|-------|---------|
| `kf` | 8.4e7 M⁻¹s⁻¹ | ligand-receptor association |
| `kr1` | 0.37 s⁻¹ | low-affinity dissociation |
| `kdes` | 0.065 s⁻¹ (0 if nondesensitizable) | desensitization |
| `ki` | 0.0033 s⁻¹ | internalization |
| `kup` | 0.004 s⁻¹ | up-regulation / recycling |
| `Rtot` | 25,000 per unit | total receptors |
| `Lmax` | 17.6 nM | peak ligand (mean over field = kd) |
| `W` | 1000 um | gradient width |
| `n` | 3 | gradient nonlinearity |
| `l` | 10 um | cell length (front-back separation) |
| `kd = kr1/kf` | 4.4 nM | dissociation constant |
| threshold | 10 receptors | orientation threshold on \|ΔLR*\| |

**Gradient profile (Eq 14, L0 = 0):**  `L(x) = Lmax * (x/W)^n`, source at `x = W`.

---

## Script Structure

```
week3_chemotaxis_neutrophil.py
│
├── PARAMETERS + gradient profile Lgrad(x)
│
├── RECEPTOR DYNAMICS (one unit)
│     rhs()              — the 3-ODE system
│     trajectory_unit()  — time course of [LR*, LRd, R]
│     steady_LRstar()    — equilibrium active complexes
│     orientation_signal()— ΔLR* = front - back at steady state
│
├── MIGRATION
│     chemotax()         — quasi-static walk using the orientation rule
│
├── FIGURE (4 panels)
│     A — receptor dynamics at x = 0.2 W (Fig 1A)
│     B — equilibrium ΔLR* vs position, desens vs nondesens (Fig 1B)
│     C — migration trajectories of both receptor types
│     D — final-position distribution over 60 cells
│
└── CONSOLE SUMMARY
      ΔLR* table vs position + interpretation
```

---

## Expected Output

| Panel | What it shows |
|-------|---------------|
| A | total LR climbs to ~7000-7600; active LR* stays low (~350-400) because desensitization converts complexes to LRd |
| B | desensitizable ΔLR* peaks ~28 near x=150 µm then falls below threshold at high [L]; nondesensitizable grows monotonically to ~333 |
| C | nondesensitizable cell drives straight to the source (1000 µm); desensitizable cell stalls in the ~350-550 µm zone (signal below threshold) |
| D | nondesensitizable cells cluster at the source; desensitizable cells spread across mid-field |

**Console:** ΔLR* vs position table and a short interpretation.

---

## Important Interpretation

In a SINGLE gradient the **nondesensitizable** cell actually performs
"better" — it drives straight to the source while the desensitizable cell
stalls in the mid-field where its signal drops below threshold. This is
**not** a bug: it is exactly the paradox the paper opens with. Nondesensitizable
receptor mutants chemotax efficiently (even more strongly) in single
gradients, which is why the role of desensitization was unclear.

The functional ADVANTAGE of desensitization appears only in **competing
gradients**, where the desensitizable cell can integrate signals and orient
to a distant or balanced source, while the nondesensitizable cell is
trapped by the nearest local source. That is the subject of the next (final)
question of the week — and Panels C/D here set up the contrast.

Why the desensitizable signal collapses at high [L]: rapid desensitization
plus slow recycling keeps the active pool LR* small and saturated, so the
front-back difference shrinks once both sides are far above kd. The
nondesensitizable receptor has no such ceiling, so its difference keeps
growing — useful in one gradient, but unable to be "switched off" to let a
distant competing source win.

---

## Requirements

```
numpy
matplotlib
scipy
```

---

## Connection to the Broader Week 3 Analysis

| Concept | Link |
|---------|------|
| Receptor desensitization (GPCR) | Chemoattractant receptors are GPCRs; all desensitize |
| Adaptation | The transient LR* spike then decay is adaptation at the receptor level |
| Fold-change / relative sensing | Desensitization lets cells respond to relative differences, not absolute [L] |
| Competing gradients (next question) | Where desensitization becomes essential for navigation |
| Disease relevance | Leukocyte trafficking in inflammation, cancer metastasis, lymphocyte homing |

The take-home: receptor desensitization is not needed for simple
single-gradient chemotaxis — it is the mechanism that lets cells integrate
competing chemoattractant signals and navigate complex fields step by step.

---

## References

- Ala Trusina, Lecture slides from the course *Physics of Molecular Diseases* (Niels Bohr Institute, 2020). The broader week covers:
- Receptor desensitization (GPCR, insulin signalling, WNT pathway)
- Adaptation and fold-change detection (Weber's law)
- Feedback algebra and conditions for perfect adaptation
- Receptor adaptation in neutrophil chemotaxis (Lin & Butcher, *Modeling the Role of Homologous Receptor Desensitization in Cell Gradient Sensing*,
  J Immunol. 2008 December 15; 181(12): 8335–8343.)
