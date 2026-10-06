# Dynamic Leading- and Trailing-Edge Morphing Aerofoil While Periodically Pitching

> **Status:** Under construction

This repository demonstrates the use of a compiled User-Defined Function (UDF) in **Ansys Fluent** to simulate an unsteady computational fluid dynamics (CFD) case involving a dynamically morphing aerofoil. The aerofoil undergoes both:

- leading-edge (LE) and trailing-edge (TE) morphing, and
- periodic pitching motion.

The UDF is intended to be compiled in Fluent launched through a **Visual Studio environment**.

This work is inspired by the framework provided in [Dynamic-Morphing-Wing](https://github.com/chawkiabd/Dynamic-Morphing-Wing), and extends the approach to combined LE-TE morphing with periodic pitching.

For a detailed explanation of the UDF inputs and mathematical formulation, see [README_UDF_BREAKDOWN.md](README_UDF_BREAKDOWN.md)

## Examples

1- 

2- Madi, F. (2026). *Numerical study of a pitching aerofoil with a morphing leading- and trailing-edge section for improved aerodynamic performance* [Conference presentation video]. ECCOMAS 2026, STS425B: Disruptive Aircraft’s Wing Design through Innovative Electroactive Morphing towards Sustainable Aviation II.


## Overview

The repository uses a second-order polynomial deformation function, modified from Wang et al. [1], to deform the aerofoil surface. A harmonic function is used to control the time-dependent morphing amplitude.

The UDF includes control parameters for:

- morphing frequency,
- morphing amplitude,
- leading-edge and trailing-edge morphing locations,
- periodic pitching motion,
- upper, lower, and blunt aerofoil surface motion.

In Fluent, two types of motion are used:

- **Grid motion** controls the morphing aerofoil surface sections.
- **Zone motion** controls the periodic pitching of the inner fluid zone.

## Geometry Requirements

The Fluent case should contain two fluid zones separated by an interface:

1. **Inner moving zone**  
   This zone should be centred around the pitching location, for example at `0.25C`, where `C` is the aerofoil chord length. This zone is hooked to the pitching UDF.

2. **Outer static zone**  
   This larger surrounding zone remains stationary.

The aerofoil should be a 4-digit NACA aerofoil with separated boundary names so that each surface can be hooked individually after compiling the UDF. For example:

- `wing_upper`
- `wing_lower`
- `blunt`

The exact zone and boundary names in your case file must match the names expected by the UDF, or the UDF must be edited accordingly.

## Repository Contents

This repository includes example Fluent case/data files and the UDF source code.

Suggested files:

- `<case-file.cas.h5>` - Fluent case file
- `<data-file.dat.h5>` - Fluent data file
- `<udf-source-file.c>` - UDF source code for morphing and pitching motion

Replace the placeholder names above with the actual files included in this repository.

## Running the Example Case

### 1. Launch Fluent

Launch Ansys Fluent from a Visual Studio environment, preferably as administrator if required by your local installation.

Then read the example case and data files:

```text
File > Read > Case & Data
```

Load:

```text
<case-file.cas.h5>
<data-file.dat.h5>
```

### 2. Allocate User-Defined Memory

Before compiling the UDF, allocate memory for node-based motion data.

In Fluent:

```text
User-Defined > Memory
```

Then:

1. Set **Node Memory Locations** to `7`.
2. Enable **Zone-Based Memory Allocation**.
3. Click **Edit**.
4. Select the aerofoil surface zones, for example:

```text
blunt
wing_lower
wing_upper
```

5. Click **OK**.

### 3. Compile the UDF

Go to:

```text
User-Defined > Functions > Compiled
```

Then:

1. Click **Add**.
2. Select the UDF source file, for example:

```text
<udf-source-file.c>
```

3. Click **Build**.
4. Click **Load**.

After loading successfully, Fluent should list four compiled UDF functions in the TUI:

```text
morphing_upper
morphing_lower
morphing_blunt
pitching_zone_motion
```

### 4. Hook the Morphing Grid Motion

The morphing motion should be hooked before switching to the transient pitching case.

Go to:

```text
Dynamic Mesh
```

Hook each surface to the corresponding grid motion UDF:

| Aerofoil surface | UDF function |
|---|---|
| `wing_upper` | `morphing_upper` |
| `wing_lower` | `morphing_lower` |
| `blunt` | `morphing_blunt` |

### 5. Run the Initial Steady-State Solution

Initialize the case and run a steady-state solution until convergence is reached.

A typical starting point is:

```text
2000 iterations
```

This steady-state solution is used as the initial condition for the transient morphing and pitching simulation.

### 6. Switch to a Transient Simulation

After the steady-state solution has converged, change the solver settings from steady to transient.

The pitching zone motion can only be hooked after the case is set to transient.

### 7. Hook the Pitching Zone Motion

Open:

```text
Cell Zone Conditions > Fluid > Contact_region_2-trg
```

For the example case, this corresponds to:

```text
ID = 8
```

Then:

1. Enable **Mesh Motion**.
2. Go to the mesh motion settings.
3. Under **Zone Motion Function**, select:

```text
pitching_zone_motion
```

4. Click **Apply**.

Note that the cell zone name and ID may differ in your own case file.

### 8. Preview the Motion

Before running the transient simulation, it is recommended to preview the motion.

Fluent provides two useful options:

- **Zone Motion Preview**  
  Allows you to view the aerofoil pitching and morphing behaviour over time.

- **Preview Mesh Motion**  
  Actually deforms the mesh during preview. Use this carefully, as the previewed deformation may not be reversible without reloading the case.

### 9. Run the Transient Simulation

For the attached example UDF, the recommended transient settings are:

```text
Time-step size: 1e-3 s
Inner iterations per time step: 150
Total time steps: 1200
```

These settings produce approximately 10 cycles of motion for the default frequency used in the UDF.

## Related Publication

The example UDF is linked to the open-access paper:

**A computational study on the aerodynamic properties of a pitching aerofoil with a morphing leading edge section at high Reynolds numbers**  
International Journal of Heat and Fluid Flow, 2026.  
https://doi.org/10.1016/j.ijheatfluidflow.2026.110531

## Reference

[1] Wang, R., Ma, X., Zhang, G., Ying, P. and Wang, X., 2023.  
Numerical simulation of continuous morphing wing with leading edge and trailing edge parabolic flaps.  
*Journal of Aerospace Engineering*, 36(5), 04023051.

## Author

Fakhreddine Madi  
Fakhreddine.Madi@uwe.ac.uk
