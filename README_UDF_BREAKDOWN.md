``markdown
# UDF Code Explanation: Dynamic LE-TE Morphing and Periodic Pitching

This document explains how the example Fluent UDF works, how the user inputs control the simulation, and the mathematical approach used to deform the aerofoil.

Example UDF file:

``text
LE_12DEG_PHI15_TE_13P731DEG_PHI165_TEXC_75.c
``

The UDF applies three independently controlled motions:

1. Leading-edge morphing
2. Trailing-edge morphing
3. Whole-aerofoil pitching motion

The morphing is applied through Fluent grid-motion hooks, while the pitching can be applied either through the grid-motion hooks, the zone-motion hook, or both, depending on the selected switches.

## Compiled UDF Outputs

After the UDF is compiled and loaded in Fluent, the following functions should be available:

```text
morphing_upper
morphing_lower
morphing_blunt
pitching_zone_motion
```

The functions should be hooked as follows:

| UDF function | Fluent use |
|---|---|
| `morphing_upper` | Grid motion for the upper aerofoil surface |
| `morphing_lower` | Grid motion for the lower aerofoil surface |
| `morphing_blunt` | Grid motion for the blunt trailing-edge face |
| `pitching_zone_motion` | Zone motion for rigid-body pitching of the inner fluid zone |

## User-Defined Node Memory Requirement

The UDF requires at least **7 User-Defined Node Memory locations**.

These are used internally as:

| UDNM index | Stored quantity |
|---|---|
| `0` | Original local x-coordinate |
| `1` | Original local y-coordinate |
| `2` | Original local z-coordinate |
| `3` | Node initialization marker |
| `4` | Original NACA thickness |
| `5` | Leading-edge spanwise weight |
| `6` | Trailing-edge spanwise weight |

The original coordinates are stored the first time the mesh-motion hook visits each node. All later node positions are rebuilt from these original coordinates, which prevents the mesh from accumulating deformation error over time.

## 1. Motion Mode Inputs

The UDF defines several motion modes for the LE and TE morphing:

```c
#define MOTION_SINE                  0
#define MOTION_PULSE                 1
#define MOTION_UP_ONLY               2
#define MOTION_DOWN_ONLY             3
#define MOTION_RAMP_HOLD             4
#define MOTION_RAMP_HOLD_WITH_PITCH  5
```

### `MOTION_SINE`

Applies continuous sinusoidal morphing.

The general form is:

```text
waveform(t) = sin(2 pi f t + phase)
```

This mode allows both positive and negative morphing.

### `MOTION_PULSE`

Applies a periodic pulse motion.

Each period contains:

```text
ramp up -> hold -> ramp down -> zero motion
```

The pulse shape is controlled by:

```c
#define PULSE_RAMP_FRAC    0.11
#define PULSE_DUR_FRAC     0.50
```

`PULSE_RAMP_FRAC` controls the duration of each smooth ramp as a fraction of the period.

`PULSE_DUR_FRAC` controls the total active pulse duration as a fraction of the period.

### `MOTION_UP_ONLY`

Uses only the positive part of the sine wave:

```text
waveform(t) = max(sin(2 pi f t + phase), 0)
```

### `MOTION_DOWN_ONLY`

Uses only the negative part of the sine wave:

```text
waveform(t) = min(sin(2 pi f t + phase), 0)
```

### `MOTION_RAMP_HOLD`

Smoothly ramps the morphing from zero to the target morphing angle, then holds the final morphed shape.

### `MOTION_RAMP_HOLD_WITH_PITCH`

Legacy mode. It behaves like `MOTION_RAMP_HOLD`, but can also request pitch ramp-and-hold when `PITCH_MOTION_MODE` is set to `PITCH_MODE_AUTO`.

## 2. Geometry Inputs

The main geometry inputs are:

```c
#define GEOMETRY_MODE_3D 0
#define Thick            0.12
#define chord            0.535
#define AoA_deg          -9.97
#define x_center          0.13375
#define y_center          0.0
```

### `GEOMETRY_MODE_3D`

Controls whether the UDF behaves as a 2D or 3D case.

```c
#define GEOMETRY_MODE_3D 0
```

means 2D mode.

```c
#define GEOMETRY_MODE_3D 1
```

enables 3D spanwise controls.

### `Thick`

This is the NACA thickness ratio.

For example:

```c
#define Thick 0.12
```

corresponds to a NACA 0012-type thickness distribution.

### `chord`

The aerofoil chord length in metres.

```c
#define chord 0.535
```

### `AoA_deg`

The fixed reference or mounting angle of attack of the clean aerofoil.

This angle is included in the final coordinate rotation.

### `x_center` and `y_center`

These define the pitching centre.

For example:

```c
#define x_center 0.13375
#define y_center 0.0
```

For a chord of `0.535 m`, this corresponds to:

```text
x_center = 0.25C
```

## 3. Leading-Edge Morphing Inputs

The LE morphing inputs are:

```c
#define LE_ENABLE         1
#define LE_MOTION_MODE    MOTION_SINE
#define freq_le           8.459
#define W_le_deg          12.0
#define phase_le_deg      15.0
#define LE_SIGN           (1.0)
#define x_s_le            0.13375
#define yLE               0.0
```

### `LE_ENABLE`

Turns LE morphing on or off.

```c
#define LE_ENABLE 1
```

enables LE morphing.

```c
#define LE_ENABLE 0
```

disables LE morphing.

### `LE_MOTION_MODE`

Selects the LE waveform.

Example:

```c
#define LE_MOTION_MODE MOTION_SINE
```

### `freq_le`

The LE morphing frequency in Hz.

```c
#define freq_le 8.459
```

### `W_le_deg`

The LE morphing amplitude or target angle in degrees.

```c
#define W_le_deg 12.0
```

This value is converted internally to radians:

```text
W_le_amp = W_le_deg pi / 180
```

### `phase_le_deg`

The phase angle of the LE morphing motion in degrees.

```c
#define phase_le_deg 15.0
```

This is converted internally to radians:

```text
phase_le = phase_le_deg pi / 180
```

### `LE_SIGN`

Controls the LE morphing direction.

```c
#define LE_SIGN (1.0)
```

uses the default direction.

```c
#define LE_SIGN (-1.0)
```

reverses the LE morphing direction.

### `x_s_le`

The x-location where the LE morphing section joins the undeformed aerofoil.

```c
#define x_s_le 0.13375
```

For the supplied chord, this is approximately `0.25C`.

Nodes with:

```text
x <= x_s_le
```

are treated as part of the LE morphing region.

### `yLE`

The LE camber reference value at the join location.

For a symmetric aerofoil this is usually:

```c
#define yLE 0.0
```

## 4. Trailing-Edge Morphing Inputs

The TE morphing inputs are:

```c
#define TE_ENABLE         1
#define TE_MOTION_MODE    MOTION_SINE
#define freq_te           8.459
#define W_te_deg          -13.73125
#define phase_te_deg      165.0
#define TE_SIGN           (1.0)
#define x_s               0.40125
#define yTE               0.0
#define TE_ANGLE_EPS      1.0e-6
```

### `TE_ENABLE`

Turns TE morphing on or off.

```c
#define TE_ENABLE 1
```

enables TE morphing.

```c
#define TE_ENABLE 0
```

disables TE morphing.

### `TE_MOTION_MODE`

Selects the TE waveform.

Example:

```c
#define TE_MOTION_MODE MOTION_SINE
```

### `freq_te`

The TE morphing frequency in Hz.

```c
#define freq_te 8.459
```

### `W_te_deg`

The TE morphing amplitude or target angle in degrees.

```c
#define W_te_deg -13.73125
```

The sign controls the direction of the TE deflection. The code also applies a sign convention internally when sending the TE angle into the Wang trailing-edge mapping.

### `phase_te_deg`

The TE phase angle in degrees.

```c
#define phase_te_deg 165.0
```

### `TE_SIGN`

Controls the TE morphing direction.

```c
#define TE_SIGN (1.0)
```

uses the default direction.

```c
#define TE_SIGN (-1.0)
```

reverses the TE morphing direction.

### `x_s`

The x-location where TE morphing starts.

```c
#define x_s 0.40125
```

For the supplied chord, this is approximately:

```text
x_s = 0.75C
```

Nodes with:

```text
x >= x_s
```

are treated as part of the TE morphing region.

### `yTE`

The TE camber reference value at the join location.

For a symmetric aerofoil this is usually:

```c
#define yTE 0.0
```

### `TE_ANGLE_EPS`

Small threshold used to skip the TE mapping when the TE angle is very close to zero.

```c
#define TE_ANGLE_EPS 1.0e-6
```

This avoids unnecessary calculations and protects against numerical noise.

## 5. Pitching Inputs

The pitching inputs are:

```c
#define PITCH_MOTION_MODE  PITCH_MODE_AUTO
#define Freq_Pitching      8.459
#define PitchAmpDeg        12.0
#define phase_pitching_deg 0.0
#define PITCH_RAMP_PERIOD  0.0
```

The available pitch modes are:

```c
#define PITCH_MODE_AUTO              0
#define PITCH_MODE_HARMONIC          1
#define PITCH_MODE_RAMP_HOLD         2
```

### `PITCH_MODE_AUTO`

Preserves the legacy behaviour.

If either LE or TE uses:

```c
MOTION_RAMP_HOLD_WITH_PITCH
```

then the pitch motion becomes ramp-and-hold.

Otherwise, the pitch remains harmonic.

### `PITCH_MODE_HARMONIC`

Forces continuous sinusoidal pitching.

The pitch law is:

```text
theta(t) = -PitchAmpRad sin(2 pi Freq_Pitching tau + phase_pitching)
```

where:

```text
tau = time - pitch_start
```

and:

```text
pitch_start = FTT_MOTION + FTT_PITCHING
```

### `PITCH_MODE_RAMP_HOLD`

Smoothly ramps the whole aerofoil from zero pitch to the target pitch angle, then holds it.

The held value is:

```text
theta = -PitchAmpRad
```

The negative sign preserves the sign convention used by the original harmonic pitching equation.

### `Freq_Pitching`

Pitching frequency in Hz.

```c
#define Freq_Pitching 8.459
```

### `PitchAmpDeg`

Pitching amplitude in degrees.

```c
#define PitchAmpDeg 12.0
```

Converted internally to radians:

```text
PitchAmpRad = PitchAmpDeg pi / 180
```

### `phase_pitching_deg`

Pitching phase angle in degrees.

```c
#define phase_pitching_deg 0.0
```

### `PITCH_RAMP_PERIOD`

Optional smoothing time for the start of harmonic pitching.

```c
#define PITCH_RAMP_PERIOD 0.0
```

A value of zero means no additional smoothing is applied to the harmonic pitch start.

## 6. Grid Pitch and Zone Pitch Switches

The code contains two important switches:

```c
#define GRID_PITCH_ENABLE  1
#define ZONE_PITCH_ENABLE  1
```

### `GRID_PITCH_ENABLE`

If set to `1`, the surface-node grid-motion hooks rotate the aerofoil surface coordinates.

```c
#define GRID_PITCH_ENABLE 1
```

If set to `0`, the surface grid-motion hooks apply morphing only, without applying pitch rotation to the surface nodes.

### `ZONE_PITCH_ENABLE`

If set to `1`, the `pitching_zone_motion` function outputs angular velocity for the moving fluid zone.

```c
#define ZONE_PITCH_ENABLE 1
```

If set to `0`, the zone-motion function returns zero angular velocity.

## 7. Motion Start and Timing Inputs

The timing inputs are:

```c
#define FTT                0.0
#define PRE_MORPH_DT       1.0e-3
#define PRE_MORPH_PERIOD   (1.0 / freq_le)
#define FTT_PITCHING       0.0
```

### `FTT`

The global time at which morphing is allowed to begin.

```c
#define FTT 0.0
```

Before this time, morphing and pitch are zero.

### `PRE_MORPH_PERIOD`

The duration of the smooth initial morphing ramp.

```c
#define PRE_MORPH_PERIOD (1.0 / freq_le)
```

By default, this is one LE morphing period.

The ramp uses a smoothstep function:

```text
alpha = 3s^2 - 2s^3
```

where:

```text
s = (time - FTT) / PRE_MORPH_PERIOD
```

This gives a smooth start with zero slope at the beginning and end of the ramp.

### `FTT_MOTION`

This is calculated internally as:

```text
FTT_MOTION = FTT + PRE_MORPH_PERIOD
```

After `FTT_MOTION`, the selected LE and TE waveform functions are used directly.

### `FTT_PITCHING`

Additional delay before harmonic pitching begins.

```c
#define FTT_PITCHING 0.0
```

The harmonic pitch start time is:

```text
pitch_start = FTT_MOTION + FTT_PITCHING
```

## 8. 3D Spanwise Controls

The code can also be used for 3D morphing cases.

The main switches are:

```c
#define GEOMETRY_MODE_3D 0
#define SPANWISE_MORPH_ENABLE       1
#define PITCH_ONLY_ON_MORPH_SPANS   0
```

When `GEOMETRY_MODE_3D = 0`, spanwise controls are ignored and the spanwise weight is always `1`.

When `GEOMETRY_MODE_3D = 1`, the code assumes:

```text
x-y plane = aerofoil section
z direction = spanwise direction
```

### Spanwise Morphing Bands

The LE and TE can each have up to three spanwise morphing bands.

Example LE band:

```c
#define LE_SPAN_1_ENABLE 1
#define LE_Z_1_MIN       0.16
#define LE_Z_1_MAX       1.34
#define LE_Z_1_RAMP      0.10
```

Example TE band:

```c
#define TE_SPAN_1_ENABLE 1
#define TE_Z_1_MIN       0.16
#define TE_Z_1_MAX       1.34
#define TE_Z_1_RAMP      0.10
```

The spanwise weighting works as:

```text
z < Z_MIN                  weight = 0
Z_MIN to Z_MIN + RAMP      smooth ramp from 0 to 1
Z_MIN + RAMP to Z_MAX      weight = 1
Z_MAX to Z_MAX + RAMP      smooth ramp from 1 to 0
z > Z_MAX + RAMP           weight = 0
```

This avoids a sharp discontinuity at the edges of the morphing span.

## 9. NACA Thickness Reconstruction

The UDF reconstructs the upper and lower surface positions using the original NACA thickness distribution.

The thickness equation is:

```text
y_t = 5 c t [
  0.2969 sqrt(x/c)
  - 0.1260 (x/c)
  - 0.3516 (x/c)^2
  + 0.2843 (x/c)^3
  - 0.1015 (x/c)^4
]
```

where:

```text
c = chord
t = Thick
x = local x-coordinate
```

For this case:

```text
c = 0.535 m
t = 0.12
```

The thickness is stored once in User-Defined Node Memory and then reused during the morphing calculation.

## 10. Leading-Edge Morphing Method

The LE deformation is based on the Wang leading-edge mapping.

For nodes in the LE morphing region:

```text
x <= x_s_le
```

the code calculates a morphed camber-line position using the instantaneous LE morphing angle.

The LE morphing angle is calculated as:

```text
W_LE(t) = W_le_amp * LE_SIGN * waveform_LE(t)
```

where:

```text
W_le_amp = W_le_deg pi / 180
```

The Wang mapping gives a new camber-line point:

```text
(xc, yc)
```

and a local camber-line slope:

```text
dyc/dx
```

The surface point is then reconstructed normal to the morphed camber line.

For a surface thickness offset `yt`, the mapped point is:

```text
x_morph = xc - yt / sqrt(1 + (dyc/dx)^2) * (dyc/dx)

y_morph = yc + yt / sqrt(1 + (dyc/dx)^2)
```

For the lower surface, the thickness offset is negative.

At shared leading-edge points, the offset is set to zero to avoid hook-order dependence between the upper and lower surface hooks.

## 11. Trailing-Edge Morphing Method

For nodes in the TE morphing region:

```text
x >= x_s
```

the code applies the Wang trailing-edge mapping.

The TE angle is calculated from:

```text
W_TE_raw(t) = W_te_amp * TE_SIGN * waveform_TE(t)
```

where:

```text
W_te_amp = W_te_deg pi / 180
```

Inside the grid-motion hooks, the mapping uses:

```text
W_TE(t) = -W_TE_raw(t)
```

This sign convention is used to preserve the intended physical TE deflection direction.

The Wang TE mapping gives:

```text
(xc, yc)
```

and:

```text
dyc/dx
```

The surface point is reconstructed using the signed offset from the original TE camber reference:

```text
signed_offset = y_original - yTE
```

Then:

```text
x_morph = xc - signed_offset / sqrt(1 + (dyc/dx)^2) * (dyc/dx)

y_morph = yc + signed_offset / sqrt(1 + (dyc/dx)^2)
```

The same TE mapping is used for the upper surface, lower surface, and blunt trailing-edge face. This helps keep shared trailing-edge corner nodes consistent.

## 12. Pitching Rotation

After LE and TE morphing are applied, the code rotates the resulting local coordinates about the pitching centre.

The rotation centre is:

```text
(x_center, y_center)
```

The final rotation angle is:

```text
theta_total = theta_pitch(t) + AoA_rad
```

The rotation is:

```text
x_rot = cos(theta_total)(x - x_center)
        - sin(theta_total)(y - y_center)
        + x_center

y_rot = sin(theta_total)(x - x_center)
        + cos(theta_total)(y - y_center)
        + y_center
```

The final node coordinates written back to Fluent are:

```text
NODE_X = x_rot
NODE_Y = y_rot
NODE_Z = z_original
```

In 2D, the z-coordinate is ignored.

## 13. Zone-Motion Pitching

The `pitching_zone_motion` function does not return an absolute pitch angle. Fluent zone motion requires angular velocity.

The code calculates angular velocity using a backward finite difference:

```text
omega(t) = [theta(t) - theta(t - dt)] / dt
```

This ensures that the zone-motion output remains consistent with the same pitch law used by the grid-motion rotation.

For 2D cases, the pitching origin is:

```text
origin = (x_center, y_center)
```

For 3D cases, the pitching axis is the global z-axis:

```text
axis = (0, 0, 1)
```

## 14. Example Input Values in This UDF

The attached example uses:

```c
#define chord       0.535
#define Thick       0.12
#define AoA_deg     -9.97
#define x_center    0.13375
#define y_center    0.0
```

Leading-edge morphing:

```c
#define LE_ENABLE      1
#define LE_MOTION_MODE MOTION_SINE
#define freq_le        8.459
#define W_le_deg       12.0
#define phase_le_deg   15.0
#define x_s_le         0.13375
```

Trailing-edge morphing:

```c
#define TE_ENABLE      1
#define TE_MOTION_MODE MOTION_SINE
#define freq_te        8.459
#define W_te_deg       -13.73125
#define phase_te_deg   165.0
#define x_s            0.40125
```

Pitching:

```c
#define PITCH_MOTION_MODE  PITCH_MODE_AUTO
#define Freq_Pitching      8.459
#define PitchAmpDeg        12.0
#define phase_pitching_deg 0.0
```

For this example:

```text
x_center = 0.25C
x_s_le   = 0.25C
x_s      = 0.75C
```

Therefore, the LE morphing region extends from the leading edge to `0.25C`, and the TE morphing region extends from `0.75C` to the trailing edge.

## 15. Typical Workflow for Editing the UDF

To create a new case, most users only need to edit the section labelled:

```c
USER INPUTS - EDIT ONLY THIS SECTION
```

A typical modification workflow is:

1. Set the chord and NACA thickness:

```c
#define chord 0.535
#define Thick 0.12
```

2. Set the pitching centre:

```c
#define x_center 0.13375
#define y_center 0.0
```

3. Choose whether LE and TE morphing are enabled:

```c
#define LE_ENABLE 1
#define TE_ENABLE 1
```

4. Choose the LE and TE motion modes:

```c
#define LE_MOTION_MODE MOTION_SINE
#define TE_MOTION_MODE MOTION_SINE
```

5. Set LE and TE frequency, amplitude, and phase:

```c
#define freq_le      8.459
#define W_le_deg     12.0
#define phase_le_deg 15.0

#define freq_te      8.459
#define W_te_deg     -13.73125
#define phase_te_deg 165.0
```

6. Set where morphing starts or ends:

```c
#define x_s_le 0.13375
#define x_s    0.40125
```

7. Set the pitching motion:

```c
#define Freq_Pitching      8.459
#define PitchAmpDeg        12.0
#define phase_pitching_deg 0.0
```

8. Compile and load the UDF in Fluent.

9. Hook the grid-motion and zone-motion functions.

10. Preview the motion before running the transient simulation.

## 16. Important Notes

Hook the UDF while the mesh is still in its original, undeformed shape.

The first mesh-motion call stores the original coordinates in User-Defined Node Memory.

If the mesh has already been previewed or deformed and you want to restart cleanly, reload the original case/data before re-hooking and previewing the motion again.

For production simulations, it is recommended to disable detailed profiling:

```c
#define UDF_PROFILE_ENABLE 0
```

The profiler is useful for debugging and performance checks, but it adds overhead because timing calls are made during the grid-motion loops.

## 17. Citation

The morphing approach is based on the leading-edge and trailing-edge mapping method described in:

Wang, R., Ma, X., Zhang, G., Ying, P. and Wang, X., 2023.  
Numerical simulation of continuous morphing wing with leading edge and trailing edge parabolic flaps.  
*Journal of Aerospace Engineering*, 36(5), 04023051.
````
