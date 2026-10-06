/*********************************************************************/
/*UDF for 2D Harmonic Morphing Trailing- and Leading Edge Morphing Aerofoil by
   Fakhreddine Madi, UWE Bristol (Supervised by Yufeng Y. , Shine W.N. , Matt O'd.)

   Inspired by the work of Chawki Abdessemed (BASE code and structure):
   https://github.com/chawkiabd/Dynamic-Morphing-Wing

   This routine deforms the grid on the upper/lower/blunt surface of the aerofoil
   (blunt Trailing-Edge only), by updating X and Y node from the

   Equation below taken from:
   "Numerical Simulation of Continuous Morphing Wing with Leading Edge and Trailing Edge Parabolic Flaps"
   (Ruochen Wang; Xiaoping Ma; Guoxin Zhang; Pei Ying; and Xiangyu Wang)
   2D/3D MORPHING AEROFOIL WITH INDEPENDENT LE, TE AND PITCH CONTROL
 * -------------------------------------------------------------------
 *
 * PURPOSE
 * -------
 * This UDF deforms a NACA-type aerofoil using the Wang leading-edge and
 * trailing-edge mappings, and can also apply rigid-body pitching about a
 * user-defined rotation centre.
 *
 * The important design feature of this version is that the three motions are
 * controlled independently:
 *
 *     1) LE_MOTION_MODE     -> controls leading-edge morphing.
 *     2) TE_MOTION_MODE     -> controls trailing-edge morphing.
 *     3) PITCH_MOTION_MODE  -> controls whole-wing pitching.
 *
 * This allows combinations such as:
 *
 *     LE  : continuously sinusoidal morphing
 *     TE  : pulse morphing
 *     PITCH: smoothly move to PitchAmpDeg once and then hold there
 *
 * without forcing the morphing to stop when pitching becomes stationary.
 *
 * -------------------------------------------------------------------
 * HOW THE CODE WORKS - HIGH LEVEL FLOW
 * -------------------------------------------------------------------
 *
 * A. USER INPUTS
 *    Every value you should normally edit is grouped immediately below this
 *    header in the section labelled "USER INPUTS - EDIT ONLY THIS SECTION".
 *
 * B. ORIGINAL COORDINATE STORAGE
 *    On the first visit to each surface node, the UDF stores the clean/local
 *    x, y and z coordinates in User-Defined Node Memory (UDNM).  All later
 *    mesh positions are rebuilt from those original coordinates.  This avoids
 *    accumulating geometric error from repeatedly deforming an already moved
 *    node.
 *
 * C. MORPH ANGLES
 *    get_le_angle() and get_te_angle() calculate the instantaneous LE and TE
 *    morph angles from the selected waveform and the current Fluent time.
 *
 * D. MORPH GEOMETRY
 *    wang_le()/wang_te() calculate the deformed camber-line geometry.  The
 *    upper/lower surface thickness offset is then reconstructed normal to the
 *    local morphed camber line.
 *
 * E. SPANWISE WEIGHTING (3D ONLY)
 *    le_span_weight()/te_span_weight() apply smooth 0->1->0 windows in z so
 *    only selected spanwise regions morph.  In 2D these weights are always 1.
 *
 * F. PITCHING
 *    get_pitch_angle() independently calculates either harmonic pitching or
 *    the one-way ramp-and-hold pitch law.  Pitch is applied AFTER LE/TE
 *    morphing, so the already-morphed local coordinates rotate rigidly about
 *    (x_center,y_center).
 *
 * G. GRID-MOTION HOOKS
 *    morphing_upper, morphing_lower and morphing_blunt repeat the same basic
 *    sequence for each node:
 *
 *       clean stored coordinate
 *              -> LE/TE deformation
 *              -> optional span weighting
 *              -> rigid pitch rotation
 *              -> write new Fluent node coordinate
 *
 * -------------------------------------------------------------------
 * FLUENT SETUP REQUIREMENTS
 * -------------------------------------------------------------------
 *
 *  - Allocate AT LEAST 7 User-Defined Node Memory locations.
 *  - Hook morphing_upper to the upper surface dynamic-mesh zone.
 *  - Hook morphing_lower to the lower surface dynamic-mesh zone.
 *  - Hook morphing_blunt to the complete blunt trailing-edge face zone.
 *  - If using the zone-motion output, hook pitching_zone_motion as required.
 *  - Hook the functions while the mesh is still in its original clean shape.
 *  - The first mesh-motion/preview call stores the original coordinates.
 *  - In 3D, this code assumes x-y is the aerofoil plane and global z is span.
 *
 *********************************************************************/

#include "udf.h"   /* Fluent UDF macros, Node/Thread definitions and motion hooks. */
#include "math.h"  /* sin, cos, tan, atan, sqrt, fabs, fmod, asinh and pow. */
#include "time.h"  /* clock() used only by the optional diagnostic profiler. */

#ifndef M_PI
/* Some C compilers do not expose M_PI by default, so define it if necessary. */
#define M_PI 3.14159265358979323846
#endif

/*****************************************************************************
 * USER INPUTS - EDIT ONLY THIS SECTION
 * ---------------------------------------------------------------------------
 * Everything that should normally be changed between Fluent cases is here.
 * Internal/derived constants, UDNM locations and helper macros are deliberately
 * kept OUTSIDE this section so they are less likely to be edited accidentally.
 *****************************************************************************/

/*===========================================================================*/
/* 1. AVAILABLE MODE NAMES                                                   */
/*===========================================================================*/
/* These names make the selections below readable.  Normally do not edit the
 * numeric values themselves; choose one of the names in the actual selectors.
 */
#define MOTION_SINE                  0  /* Full positive/negative sine morphing. */
#define MOTION_PULSE                 1  /* Periodic ramp-up -> hold -> ramp-down pulse. */
#define MOTION_UP_ONLY               2  /* Keep only the positive half of the sine. */
#define MOTION_DOWN_ONLY             3  /* Keep only the negative half of the sine. */
#define MOTION_RAMP_HOLD             4  /* Morph once to the requested edge angle, then hold. */
#define MOTION_RAMP_HOLD_WITH_PITCH  5  /* Legacy Mode 5: Mode 4 + AUTO pitch hold request. */

#define PITCH_MODE_AUTO              0  /* Legacy-compatible automatic pitch choice. */
#define PITCH_MODE_HARMONIC          1  /* Continuous sinusoidal pitching. */
#define PITCH_MODE_RAMP_HOLD         2  /* Pitch once to target angle, then hold. */

/*===========================================================================*/
/* 2. DIMENSION / GENERAL GEOMETRY                                           */
/*===========================================================================*/
#define GEOMETRY_MODE_3D 0   /* 0 = 2D behaviour; 1 = enable 3D spanwise controls. */
#define Thick            0.12 /* NACA thickness ratio t/c; 0.12 corresponds to NACA0012. */
#define chord            0.535 /* Aerofoil chord length [m]. */

/* The clean/reference aerofoil can already be mounted at a fixed AoA.  AoA_deg
 * is included in every final coordinate rotation in addition to dynamic pitch.
 */
#define AoA_deg          -9.97  /* Fixed reference/mounting angle of attack [deg]. */
#define x_center          0.13375 /* x-coordinate of pitch centre in local coordinates [m]. */
#define y_center          0.0    /* y-coordinate of pitch centre in local coordinates [m]. */

/*===========================================================================*/
/* 3. LEADING-EDGE MORPHING INPUTS                                           */
/*===========================================================================*/
#define LE_ENABLE         1      /* 1 = enable LE morphing; 0 = keep LE clean. */
#define LE_MOTION_MODE    MOTION_SINE /* Select one of MOTION_* above. */
#define freq_le           8.459    /* LE morphing frequency [Hz] for periodic modes. */
#define W_le_deg          12.0   /* Full LE morph amplitude/target angle [deg]. */
#define phase_le_deg      15.0  /* LE phase angle [deg] for periodic modes. */
#define LE_SIGN           (1.0)  /* +1 = normal direction; -1 = reverse LE morph direction. */
#define x_s_le            0.13375 /* x-position where LE morphing joins clean aerofoil [m]. */
#define yLE               0.0    /* y-coordinate of LE morphing join/camber reference [m]. */

/*===========================================================================*/
/* 4. TRAILING-EDGE MORPHING INPUTS                                          */
/*===========================================================================*/
#define TE_ENABLE         1      /* 1 = enable TE morphing; 0 = keep TE clean. */
#define TE_MOTION_MODE    MOTION_SINE /* Select one of MOTION_* above. */
#define freq_te           8.459    /* TE morphing frequency [Hz] for periodic modes. */
#define W_te_deg          -13.73125  /* Full TE morph amplitude/target angle [deg]. */
#define phase_te_deg      165.0  /* TE phase angle [deg] for periodic modes. */
#define TE_SIGN           (1.0)  /* +1 = normal direction; -1 = reverse TE morph direction. */
#define x_s               0.40125 /* x-position where TE morphing starts [m]. */
#define yTE               0.0    /* y-coordinate of TE morphing join/camber reference [m]. */
#define TE_ANGLE_EPS      1.0e-6 /* Skip TE mapping when |angle| is smaller than this [rad]. */

/*===========================================================================*/
/* 5. INDEPENDENT WHOLE-WING PITCHING INPUTS                                 */
/*===========================================================================*/
/*
 * PITCH_MODE_AUTO:
 *   Preserves the previous Mode-5 behaviour.  If LE or TE uses
 *   MOTION_RAMP_HOLD_WITH_PITCH, pitch ramps once and holds. Otherwise pitch
 *   remains harmonic.
 *
 * PITCH_MODE_HARMONIC:
 *   Forces continuous sinusoidal pitch irrespective of LE/TE modes.
 *
 * PITCH_MODE_RAMP_HOLD:
 *   Smoothly changes whole-wing pitch from 0 to -PitchAmpDeg during
 *   FTT -> FTT_MOTION and then holds it there. LE and TE keep executing their
 *   own selected morph laws independently.
 */
#define PITCH_MOTION_MODE  PITCH_MODE_AUTO /* Choose AUTO, HARMONIC or RAMP_HOLD. */
#define Freq_Pitching      8.459    /* Harmonic pitch frequency [Hz]; ignored by pitch hold mode. */
#define PitchAmpDeg        12.0   /* Pitch amplitude/held target magnitude [deg]. */
#define phase_pitching_deg 0.0    /* Harmonic pitch phase [deg]; ignored by pitch hold mode. */
#define PITCH_RAMP_PERIOD  0.0    /* Optional harmonic-start smoothing time [s]; 0 = none. */

/* These switches control where the calculated pitch is output.  Normally both
 * remain 1 to preserve the supplied behaviour.
 */
#define GRID_PITCH_ENABLE  1 /* 1 = grid-motion hooks rotate surface-node coordinates. */
#define ZONE_PITCH_ENABLE  1 /* 1 = pitching_zone_motion returns angular velocity. */

/*===========================================================================*/
/* 6. MOTION START / PRE-MORPH TIMING                                        */
/*===========================================================================*/
#define FTT                0.0    /* Global time at which morphing is allowed to begin [s]. */
#define PRE_MORPH_DT       1.0e-3 /* Retained legacy preview/time-scale input; not used directly. */
#define PRE_MORPH_PERIOD   (1.0 / freq_le) /* Smooth initial morph duration [s]; default = one LE cycle. */
#define FTT_PITCHING       0.0    /* Extra delay after FTT_MOTION before harmonic pitch begins [s]. */

/*===========================================================================*/
/* 7. PULSE-WAVEFORM SHAPE                                                   */
/*===========================================================================*/
/* These fractions apply only when an edge uses MOTION_PULSE.  They are
 * fractions of that edge's own period T = 1/frequency.
 */
#define PULSE_RAMP_FRAC    0.11 /* Fraction of T used for EACH cosine ramp. */
#define PULSE_DUR_FRAC     0.50 /* Total fraction of T for which the pulse is active. */

/*===========================================================================*/
/* 8. 3D SPANWISE MORPH / PITCH CONTROLS                                     */
/*===========================================================================*/
/* These settings are ignored when GEOMETRY_MODE_3D = 0. */
#define SPANWISE_MORPH_ENABLE       1 /* 1 = use LE/TE span bands in 3D; 0 = morph entire span. */
#define PITCH_ONLY_ON_MORPH_SPANS   0 /* 1 = pitch only where LE/TE span weight exists; 0 = pitch whole span. */

/* Each enabled band uses:
 *   Z_MIN = start of smooth ramp from weight 0 -> 1
 *   Z_MAX = start of smooth ramp from weight 1 -> 0
 *   RAMP  = axial/spanwise length of each transition
 *
 * Therefore a band with MIN=0.16, MAX=1.34 and RAMP=0.10 gives:
 *   z < 0.16       : weight 0
 *   0.16 -> 0.26   : smooth 0 -> 1
 *   0.26 -> 1.34   : weight 1
 *   1.34 -> 1.44   : smooth 1 -> 0
 *   z > 1.44       : weight 0
 */
#define LE_SPAN_1_ENABLE 1    /* Enable LE span band 1. */
#define LE_Z_1_MIN       0.16 /* LE band 1 ramp-in start [m]. */
#define LE_Z_1_MAX       1.34 /* LE band 1 ramp-out start [m]. */
#define LE_Z_1_RAMP      0.10 /* LE band 1 transition length [m]. */

#define LE_SPAN_2_ENABLE 0    /* Enable LE span band 2. */
#define LE_Z_2_MIN       0.00 /* LE band 2 ramp-in start [m]. */
#define LE_Z_2_MAX       0.00 /* LE band 2 ramp-out start [m]. */
#define LE_Z_2_RAMP      0.00 /* LE band 2 transition length [m]. */

#define LE_SPAN_3_ENABLE 0    /* Enable LE span band 3. */
#define LE_Z_3_MIN       0.00 /* LE band 3 ramp-in start [m]. */
#define LE_Z_3_MAX       0.00 /* LE band 3 ramp-out start [m]. */
#define LE_Z_3_RAMP      0.00 /* LE band 3 transition length [m]. */

#define TE_SPAN_1_ENABLE 1    /* Enable TE span band 1. */
#define TE_Z_1_MIN       0.16 /* TE band 1 ramp-in start [m]. */
#define TE_Z_1_MAX       1.34 /* TE band 1 ramp-out start [m]. */
#define TE_Z_1_RAMP      0.10 /* TE band 1 transition length [m]. */

#define TE_SPAN_2_ENABLE 0    /* Enable TE span band 2. */
#define TE_Z_2_MIN       0.00 /* TE band 2 ramp-in start [m]. */
#define TE_Z_2_MAX       0.00 /* TE band 2 ramp-out start [m]. */
#define TE_Z_2_RAMP      0.00 /* TE band 2 transition length [m]. */

#define TE_SPAN_3_ENABLE 0    /* Enable TE span band 3. */
#define TE_Z_3_MIN       0.00 /* TE band 3 ramp-in start [m]. */
#define TE_Z_3_MAX       0.00 /* TE band 3 ramp-out start [m]. */
#define TE_Z_3_RAMP      0.00 /* TE band 3 transition length [m]. */

/*===========================================================================*/
/* 9. OPTIONAL UDF PERFORMANCE PROFILER                                      */
/*===========================================================================*/
/* Profiling is useful while checking UDF cost, but the detailed per-section
 * timers call clock() inside node loops and therefore add measurable overhead.
 * For production simulations use UDF_PROFILE_ENABLE=0.
 */
#define UDF_PROFILE_ENABLE      1 /* 1 = print timing diagnostics; 0 = compile profiler out. */
#define UDF_PROFILE_SECTIONS    1 /* 1 = detailed section timing; 0 = totals/counts only. */
#define UDF_PROFILE_PRINT_EVERY 1 /* Print one profiler line every N hook calls. */

/*****************************************************************************
 * END OF USER INPUTS
 * ---------------------------------------------------------------------------
 * Normally do not edit below this line unless changing the UDF implementation.
 *****************************************************************************/

/*===========================================================================*/
/* INPUT VALIDATION AND COMPILE-TIME CONTROL                                 */
/*===========================================================================*/
/* Stop compilation immediately when a selector contains an invalid value. */
#if (GEOMETRY_MODE_3D != 0) && (GEOMETRY_MODE_3D != 1)
#error "GEOMETRY_MODE_3D must be 0 (2D) or 1 (3D)."
#endif

#if GEOMETRY_MODE_3D && !RP_3D
#error "GEOMETRY_MODE_3D=1 requires a 3D Fluent case. Use GEOMETRY_MODE_3D=0 for 2D."
#endif

#if (LE_MOTION_MODE < MOTION_SINE) || (LE_MOTION_MODE > MOTION_RAMP_HOLD_WITH_PITCH)
#error "LE_MOTION_MODE must be one of MOTION_SINE through MOTION_RAMP_HOLD_WITH_PITCH."
#endif

#if (TE_MOTION_MODE < MOTION_SINE) || (TE_MOTION_MODE > MOTION_RAMP_HOLD_WITH_PITCH)
#error "TE_MOTION_MODE must be one of MOTION_SINE through MOTION_RAMP_HOLD_WITH_PITCH."
#endif

#if (PITCH_MOTION_MODE < PITCH_MODE_AUTO) || (PITCH_MOTION_MODE > PITCH_MODE_RAMP_HOLD)
#error "PITCH_MOTION_MODE must be PITCH_MODE_AUTO, PITCH_MODE_HARMONIC or PITCH_MODE_RAMP_HOLD."
#endif

/* If profiling is completely disabled, forcibly disable intrusive section
 * timers as well. This removes their code at preprocessing time.
 */
#if !UDF_PROFILE_ENABLE
#undef UDF_PROFILE_SECTIONS
#define UDF_PROFILE_SECTIONS 0
#endif

/* Resolve the requested pitch law ONCE at compile time.  This avoids checking
 * pitch-mode selectors repeatedly inside every node update.
 */
#if (PITCH_MOTION_MODE == PITCH_MODE_RAMP_HOLD)
#define PITCH_RAMP_HOLD_ACTIVE 1
#elif (PITCH_MOTION_MODE == PITCH_MODE_HARMONIC)
#define PITCH_RAMP_HOLD_ACTIVE 0
#else
/* AUTO preserves legacy Mode-5 coupling: either edge using Mode 5 requests
 * pitch ramp-and-hold; otherwise pitch remains harmonic.
 */
#if (LE_MOTION_MODE == MOTION_RAMP_HOLD_WITH_PITCH) || \
    (TE_MOTION_MODE == MOTION_RAMP_HOLD_WITH_PITCH)
#define PITCH_RAMP_HOLD_ACTIVE 1
#else
#define PITCH_RAMP_HOLD_ACTIVE 0
#endif
#endif

/*===========================================================================*/
/* DERIVED CONSTANTS - CALCULATED FROM THE USER INPUTS                       */
/*===========================================================================*/
#define FTT_MOTION       (FTT + PRE_MORPH_PERIOD) /* Time at which initial smooth morph ramp is complete. */
#define AoA_rad          (AoA_deg * M_PI / 180.0) /* Fixed AoA converted to radians. */
#define PitchAmpRad      (PitchAmpDeg * M_PI / 180.0) /* Pitch amplitude/target converted to radians. */
#define phase_le         (phase_le_deg * M_PI / 180.0) /* LE phase converted to radians. */
#define phase_te         (phase_te_deg * M_PI / 180.0) /* TE phase converted to radians. */
#define phase_pitching   (phase_pitching_deg * M_PI / 180.0) /* Pitch phase converted to radians. */
#define W_le_amp         (W_le_deg * M_PI / 180.0) /* LE target/amplitude converted to radians. */
#define W_te_amp         (W_te_deg * M_PI / 180.0) /* TE target/amplitude converted to radians. */

/*===========================================================================*/
/* OPTIONAL PROFILER IMPLEMENTATION                                          */
/*===========================================================================*/
#if UDF_PROFILE_ENABLE
/* Convert processor-clock ticks returned by clock() into milliseconds. */
static double udf_profile_ms(clock_t ticks)
{
    return 1000.0 * ((double)ticks / (double)CLOCKS_PER_SEC);
}

/* Print one compact diagnostic line for a grid-motion hook.  The routine is
 * intentionally outside the node loops so formatting/Message0 overhead occurs
 * only once per reported hook call.
 */
static void udf_profile_report(const char *hook_name, real flow_time, long call_id,
                               clock_t total_ticks, clock_t init_ticks,
                               clock_t le_ticks, clock_t te_ticks,
                               clock_t rotate_ticks, long face_node_visits,
                               long updated_nodes, long le_nodes, long te_nodes)
{
    double total_ms = udf_profile_ms(total_ticks);  /* Total hook CPU time. */
    double init_ms = udf_profile_ms(init_ticks);    /* Original-coordinate initialization time. */
    double le_ms = udf_profile_ms(le_ticks);        /* LE mapping time. */
    double te_ms = udf_profile_ms(te_ticks);        /* TE mapping time. */
    double rotate_ms = udf_profile_ms(rotate_ticks);/* Pitch rotation + node write time. */
    double accounted_ms = init_ms + le_ms + te_ms + rotate_ms;
    double traversal_ms = total_ms - accounted_ms;  /* Face/node traversal and miscellaneous work. */

    if (traversal_ms < 0.0)
        traversal_ms = 0.0; /* Timer resolution can otherwise create tiny negative residuals. */

    Message0("[UDF-PROFILE] %s call=%ld flow-time=%.9g total=%.6f ms "
             "traversal/other=%.6f ms init=%.6f ms LE=%.6f ms "
             "TE=%.6f ms rotate/write=%.6f ms visits=%ld updated=%ld "
             "LE-nodes=%ld TE-nodes=%ld\n",
             hook_name, call_id, flow_time, total_ms, traversal_ms, init_ms,
             le_ms, te_ms, rotate_ms, face_node_visits, updated_nodes,
             le_nodes, te_nodes);
}
#endif

/*------------------------------------------------------------
 * SHARED NODE MEMORY LOCATIONS
 *-----------------------------------------------------------
 * Allocate at least seven User-Defined Node Memory locations:
 *   0 = original local x
 *   1 = original local y
 *   2 = original local z
 *   3 = initialization marker
 *   4 = original NACA thickness
 *   5 = leading-edge span weight (legacy/unused by motion hooks)
 *   6 = trailing-edge span weight (legacy/unused by motion hooks)
 *
 * Span weights are recalculated from ORIG_Z_UDMI during every motion call,
 * so changing span MIN/MAX/RAMP values does not depend on cached weights.
 */

#define ORIG_X_UDMI     0
#define ORIG_Y_UDMI     1
#define ORIG_Z_UDMI     2
#define INIT_FLAG_UDMI  3
#define THICKNESS_UDMI  4
#define LE_WEIGHT_UDMI  5
#define TE_WEIGHT_UDMI  6

/* Exact, representable marker used to distinguish initialized nodes. */
#define INIT_MAGIC      12345.0
#define INIT_MAGIC_TOL  0.25

/*------------------------------------------------------------
 * 2D/3D NODE ACCESS HELPERS
 *-----------------------------------------------------------*/

#if RP_3D
#define GET_NODE_Z(v)       NODE_Z(v)
#define SET_NODE_Z(v, zv)   (NODE_Z(v) = (zv))
#else
#define GET_NODE_Z(v)       0.0
#define SET_NODE_Z(v, zv)   ((void)0)
#endif

/*
 * In 2D mode, spanwise effects are deliberately disabled even if the
 * same source is compiled in a 3D Fluent session.
 */
#if GEOMETRY_MODE_3D && RP_3D
#define GET_NODE_SPAN_COORD(v) GET_NODE_Z(v)
#else
#define GET_NODE_SPAN_COORD(v) 0.0
#endif

/*------------------------------------------------------------
 * GLOBAL TO LOCAL COORDINATE TRANSFORMATION
 *-----------------------------------------------------------*/

#define GLOBAL_TO_LOCAL_3D(xg, yg, zg, xl, yl, zl)                 \
    do                                                              \
    {                                                               \
        (xl) = ((xg) - x_center) * cos(AoA_rad) +                  \
               ((yg) - y_center) * sin(AoA_rad) + x_center;        \
        (yl) = -((xg) - x_center) * sin(AoA_rad) +                 \
               ((yg) - y_center) * cos(AoA_rad) + y_center;        \
        (zl) = (zg);                                                \
    } while (0)

#define MY_ISNAN(a) ((a) != (a))

/*------------------------------------------------------------
 * UTILITY FUNCTIONS
 *-----------------------------------------------------------*/

/* Clamp any scalar to [0,1].  Used before smoothstep so transition weights
 * cannot overshoot when the current time/z lies outside a ramp interval. */
static real clamp01(real v)
{
    if (v < 0.0) return 0.0;
    if (v > 1.0) return 1.0;
    return v;
}

/* Cubic smoothstep s = 3v^2-2v^3.  It gives zero slope at both ends, which
 * avoids an abrupt velocity change at the start/end of a morphing ramp. */
static real smoothstep01(real v)
{
    v = clamp01(v);
    return v * v * (3.0 - 2.0 * v);
}

/* Small helper because Fluent's supported compiler environment should not
 * depend on C++/library max() behaviour. */
static real max_real(real a, real b)
{
    return (a > b) ? a : b;
}

static real span_band_weight(real z, real z_min_in, real z_max_in, real ramp_in)
{
    real z_min = z_min_in;
    real z_max = z_max_in;
    real ramp = fabs(ramp_in);

    if (z_max < z_min)
    {
        real tmp = z_min;
        z_min = z_max;
        z_max = tmp;
    }

    if (z_max <= z_min)
        return 0.0;

    if (ramp <= 0.0)
        return (z >= z_min && z <= z_max) ? 1.0 : 0.0;

    /*
     * Interpret z_min and z_max as the START locations of the two
     * transitions:
     *
     *   z_min ........ z_min+ramp ........ z_max ........ z_max+ramp
     *     0    ramp-in      1      full       1    ramp-out      0
     *
     * This makes, for example:
     *   MIN=0.16, MAX=1.44, RAMP=0.10
     * ramp in over 0.16-0.26 and ramp out over 1.44-1.54.
     */
    return smoothstep01((z - z_min) / ramp) *
           smoothstep01((z_max + ramp - z) / ramp);
}

static real le_span_weight(real z)
{
#if GEOMETRY_MODE_3D && RP_3D && SPANWISE_MORPH_ENABLE
    real w = 0.0;

    if (LE_SPAN_1_ENABLE)
        w = max_real(w, span_band_weight(z, LE_Z_1_MIN, LE_Z_1_MAX, LE_Z_1_RAMP));

    if (LE_SPAN_2_ENABLE)
        w = max_real(w, span_band_weight(z, LE_Z_2_MIN, LE_Z_2_MAX, LE_Z_2_RAMP));

    if (LE_SPAN_3_ENABLE)
        w = max_real(w, span_band_weight(z, LE_Z_3_MIN, LE_Z_3_MAX, LE_Z_3_RAMP));

    return clamp01(w);
#else
    (void)z;
    return 1.0;
#endif
}

static real te_span_weight(real z)
{
#if GEOMETRY_MODE_3D && RP_3D && SPANWISE_MORPH_ENABLE
    real w = 0.0;

    if (TE_SPAN_1_ENABLE)
        w = max_real(w, span_band_weight(z, TE_Z_1_MIN, TE_Z_1_MAX, TE_Z_1_RAMP));

    if (TE_SPAN_2_ENABLE)
        w = max_real(w, span_band_weight(z, TE_Z_2_MIN, TE_Z_2_MAX, TE_Z_2_RAMP));

    if (TE_SPAN_3_ENABLE)
        w = max_real(w, span_band_weight(z, TE_Z_3_MIN, TE_Z_3_MAX, TE_Z_3_RAMP));

    return clamp01(w);
#else
    (void)z;
    return 1.0;
#endif
}

static real pitch_span_weight(real z)
{
#if GEOMETRY_MODE_3D && RP_3D && PITCH_ONLY_ON_MORPH_SPANS
    return clamp01(max_real(le_span_weight(z), te_span_weight(z)));
#else
    (void)z;
    return 1.0;
#endif
}

/* Return the clean NACA thickness magnitude at x_local.  This stored value
 * is later placed normal to the morphed LE camber line to reconstruct the
 * upper/lower surface without cumulatively deforming the old mesh. */
static real naca_thickness(real x_local)
{
    real x_norm = clamp01(x_local / chord);
    real polynomial;

    polynomial = x_norm *
                 (-0.1260 + x_norm *
                 (-0.3516 + x_norm *
                 ( 0.2843 - 0.1015 * x_norm)));

    return 5.0 * chord * Thick *
           (0.2969 * sqrt(x_norm) + polynomial);
}

/* Evaluate the Wang mapping denominator while protecting the W -> 0 limit
 * from division by tan(W).  The analytic limiting value is 2. */
static real safe_denom(real W)
{
    real t = tan(W);
    real eps = 1.0e-12;

    if (fabs(t) < eps)
        return 2.0;

    return sqrt(4.0 * t * t + 1.0) + asinh(2.0 * t) / (2.0 * t);
}

/*------------------------------------------------------------
 * AUTOMATIC ORIGINAL-COORDINATE STORAGE
 *-----------------------------------------------------------
 * This function is called only for nodes encountered by the hooked surface
 * zones. It stores the clean coordinates before the node is moved for the
 * first time. The initialization marker is shared by all three surfaces.
 */

static int initial_coords_stored(Node *v)
{
    real flag = N_UDMI(v, INIT_FLAG_UDMI);

    if (MY_ISNAN(flag))
        return 0;

    return (fabs(flag - INIT_MAGIC) <= INIT_MAGIC_TOL);
}

static void ensure_initial_coords(Node *v)
{
    if (!initial_coords_stored(v))
    {
        real xg = NODE_X(v);
        real yg = NODE_Y(v);
        real zg = GET_NODE_Z(v);
        real xl, yl, zl;

        GLOBAL_TO_LOCAL_3D(xg, yg, zg, xl, yl, zl);

        N_UDMI(v, ORIG_X_UDMI) = xl;
        N_UDMI(v, ORIG_Y_UDMI) = yl;

#if GEOMETRY_MODE_3D && RP_3D
        N_UDMI(v, ORIG_Z_UDMI) = zl;
#else
        N_UDMI(v, ORIG_Z_UDMI) = 0.0;
#endif

        /* These values depend only on the clean, original coordinates. */
        N_UDMI(v, THICKNESS_UDMI) = naca_thickness(xl);
        N_UDMI(v, LE_WEIGHT_UDMI) = le_span_weight(zl);
        N_UDMI(v, TE_WEIGHT_UDMI) = te_span_weight(zl);

        /* Write the marker last so partially initialized nodes are rejected. */
        N_UDMI(v, INIT_FLAG_UDMI) = INIT_MAGIC;
    }
}

/*------------------------------------------------------------
 * MORPHING WAVEFORM FUNCTIONS
 *-----------------------------------------------------------*/

/*
 * Periodic cosine-ramp pulse.
 *
 * For each period T:
 *   0 -> ramp up -> hold -> ramp down -> 0
 *
 * td = total ON duration
 * R1 = ramp-up duration
 * R2 = ramp-down duration
 * t_shift = phase expressed as a time shift
 */
static real cosine_pulse_periodic(real t, real t0, real T,
                                  real td, real R1, real R2,
                                  real t_shift)
{
    real tau, tss, eps;

    if (t < t0)
        return 0.0;

    if (T <= 0.0)
        return 0.0;

    tau = fmod(t - t0 - t_shift, T);

    if (tau < 0.0)
        tau += T;

    if (td <= 0.0)
        return 0.0;

    if (td > T)
        td = T;

    if (R1 < 0.0)
        R1 = 0.0;

    if (R2 < 0.0)
        R2 = 0.0;

    /* Prevent the two ramps from exceeding the available ON time. */
    if (R1 + R2 > td)
    {
        real scale = td / (R1 + R2 + 1.0e-30);
        R1 *= scale;
        R2 *= scale;
    }

    tss = td - R1 - R2;

    if (tss < 0.0)
        tss = 0.0;

    if (tau > td)
        return 0.0;

    eps = 1.0e-30;

    /* Half-cosine ramp from 0 to 1. */
    if (tau <= R1)
    {
        if (R1 <= eps)
            return 1.0;

        return 0.5 * (1.0 - cos(M_PI * tau / R1));
    }

    /* Constant maximum during the hold section. */
    if (tau <= (R1 + tss))
        return 1.0;

    /* Half-cosine ramp from 1 back to 0. */
    {
        real xi = tau - (R1 + tss);

        if (R2 <= eps)
            return 0.0;

        return 0.5 * (1.0 + cos(M_PI * xi / R2));
    }
}

/*
 * Return the selected dimensionless morphing waveform.
 *
 * SINE:
 *     sin(phase)
 *
 * PULSE:
 *     positive cosine-ramp pulse from 0 to 1
 *
 * UP_ONLY:
 *     max(sin(phase),0)
 *
 * DOWN_ONLY:
 *     min(sin(phase),0)
 *
 * RAMP_HOLD and RAMP_HOLD_WITH_PITCH:
 *     1.0 after their start time. The actual one-way 0 -> 1 edge ramp is
 *     already supplied by get_le_angle()/get_te_angle() during
 *     FTT -> FTT_MOTION, so returning 1.0 here makes the final morph angle
 *     remain fixed forever. Mode 5 additionally changes get_pitch_angle()
 *     to the one-way global pitch ramp-and-hold law.
 *
 * The returned value is later multiplied by the user-specified
 * amplitude and LE_SIGN/TE_SIGN.
 */
/*------------------------------------------------------------
 * PHASE -> TIME SHIFT
 * Converts phase in radians to an equivalent time shift in seconds.
 *-----------------------------------------------------------*/
static real phase_to_shift_seconds(real phase_rad, real T)
{
    if (T <= 0.0)
        return 0.0;

    return (phase_rad / (2.0 * M_PI)) * T;
}

static real selected_morph_waveform(real time, real t0, real freq,
                                    real phase_rad, int mode)
{
    real T;
    real arg;
    real sine_value;

    if (time < t0 || freq <= 1.0e-12)
        return 0.0;

    T = 1.0 / freq;

    arg = 2.0 * M_PI * freq * (time - t0) + phase_rad;
    sine_value = sin(arg);

    if (mode == MOTION_SINE)
    {
        return sine_value;
    }
    else if (mode == MOTION_PULSE)
    {
        real td = PULSE_DUR_FRAC * T;
        real R1 = PULSE_RAMP_FRAC * T;
        real R2 = PULSE_RAMP_FRAC * T;
        real t_shift = phase_to_shift_seconds(phase_rad, T);

        return cosine_pulse_periodic(time, t0, T, td, R1, R2, t_shift);
    }
    else if (mode == MOTION_UP_ONLY)
    {
        return (sine_value > 0.0) ? sine_value : 0.0;
    }
    else if (mode == MOTION_DOWN_ONLY)
    {
        /* Keep only the negative half of the sinusoid. */
        return (sine_value < 0.0) ? sine_value : 0.0;
    }
    else if (mode == MOTION_RAMP_HOLD ||
             mode == MOTION_RAMP_HOLD_WITH_PITCH)
    {
        /*
         * Both one-way hold modes deliberately return a constant edge target
         * of 1.0.
         *
         * Before FTT_MOTION, get_le_angle()/get_te_angle() multiply this
         * target by smoothstep01((time-FTT)/PRE_MORPH_PERIOD), producing a
         * smooth one-way morph from zero to the complete requested edge angle.
         *
         * At and after FTT_MOTION this function continues returning 1.0, so
         * the final edge geometry is held indefinitely.
         *
         * Mode 4 stops there and leaves pitching unchanged.
         * Mode 5 remains a legacy morph-and-hold selector. In PITCH_MODE_AUTO
         * it also requests pitch ramp-and-hold, preserving previous behaviour.
         * With the independent pitch controller, pitch can now be overridden.
         */
        return 1.0;
    }

    /* Invalid mode: fail safely to zero motion. Compile-time checks above
     * should normally prevent execution from ever reaching this fallback. */
    return 0.0;
}

/*------------------------------------------------------------
 * TIME-DEPENDENT ANGLES
 *-----------------------------------------------------------*/

static real get_le_angle(real time)
{
    real waveform;

    if (!LE_ENABLE)
        return 0.0;

    if (time <= FTT)
        return 0.0;

    /*
     * Preserve the original pre-morphing start behaviour.
     * During FTT -> FTT_MOTION the code smoothly introduces the LE angle.
     *
     * For MOTION_RAMP_HOLD, selected_morph_waveform() evaluated at
     * FTT_MOTION returns exactly 1.0. Consequently alpha below performs a
     * smooth 0 -> full W_le_amp ramp. After FTT_MOTION the selected waveform
     * remains 1.0, so the LE stays at that final morphed shape.
     *
     * All other motion modes follow their original behaviour unchanged.
     */
    if (time < FTT_MOTION)
    {
        real alpha = smoothstep01((time - FTT) / PRE_MORPH_PERIOD);
        real start_wave = selected_morph_waveform(
            FTT_MOTION, FTT_MOTION, freq_le, phase_le, LE_MOTION_MODE);

        return alpha * W_le_amp * LE_SIGN * start_wave;
    }

    waveform = selected_morph_waveform(
        time, FTT_MOTION, freq_le, phase_le, LE_MOTION_MODE);

    return W_le_amp * LE_SIGN * waveform;
}

static real get_te_angle(real time)
{
    real waveform;

    if (!TE_ENABLE)
        return 0.0;

    if (time <= FTT)
        return 0.0;

    /*
     * Preserve the original pre-morphing start behaviour for TE as well.
     * With MOTION_RAMP_HOLD this becomes a smooth 0 -> full W_te_amp ramp,
     * followed by a permanent hold at that final TE morph angle.
     * Other TE motion modes are unchanged.
     */
    if (time < FTT_MOTION)
    {
        real alpha = smoothstep01((time - FTT) / PRE_MORPH_PERIOD);
        real start_wave = selected_morph_waveform(
            FTT_MOTION, FTT_MOTION, freq_te, phase_te, TE_MOTION_MODE);

        return alpha * W_te_amp * TE_SIGN * start_wave;
    }

    waveform = selected_morph_waveform(
        time, FTT_MOTION, freq_te, phase_te, TE_MOTION_MODE);

    return W_te_amp * TE_SIGN * waveform;
}

static real get_pitch_start_time(void)
{
    return FTT_MOTION + FTT_PITCHING;
}

static real get_pitch_angle(real time)
{
#if PITCH_RAMP_HOLD_ACTIVE
    /*
     * INDEPENDENT ONE-WAY PITCH RAMP + PERMANENT HOLD
     * ------------------------------------------------------------
     * This branch is selected when:
     *
     *   PITCH_MOTION_MODE = PITCH_MODE_RAMP_HOLD
     *
     * or when PITCH_MOTION_MODE = PITCH_MODE_AUTO and one of the edges uses
     * legacy morph Mode 5.
     *
     * The pitch motion is intentionally independent of LE/TE morphing:
     *
     * Before FTT:
     *     theta = 0
     *
     * FTT -> FTT_MOTION:
     *     theta smoothly ramps from 0 to -PitchAmpRad
     *
     * At/after FTT_MOTION:
     *     theta = -PitchAmpRad permanently
     *
     * LE and TE continue to use get_le_angle() and get_te_angle(), so selecting
     * PITCH_MODE_RAMP_HOLD does NOT freeze, alter, phase-shift or otherwise
     * interfere with their selected morphing waveforms.
     *
     * The minus sign preserves the physical sign convention of the original
     * harmonic equation: -PitchAmpRad*sin(...).
     */
    if (time <= FTT)
        return 0.0;

    if (time < FTT_MOTION)
    {
        real alpha = smoothstep01((time - FTT) / PRE_MORPH_PERIOD);
        return -PitchAmpRad * alpha;
    }

    return -PitchAmpRad;
#else
    /*
     * ORIGINAL CONTINUOUS HARMONIC PITCHING LAW
     * ------------------------------------------------------------
     * Used when PITCH_MOTION_MODE = PITCH_MODE_HARMONIC, or when AUTO resolves
     * to harmonic because neither edge uses legacy Mode 5.
     *
     * FTT_PITCHING, Freq_Pitching, phase_pitching_deg and PITCH_RAMP_PERIOD
     * retain exactly their original meanings here.
     */
    real pitch_start = get_pitch_start_time();
    real tau = time - pitch_start;
    real ramp = 1.0;

    if (tau <= 0.0)
        return 0.0;

    if (PITCH_RAMP_PERIOD > 0.0)
        ramp = smoothstep01(tau / PITCH_RAMP_PERIOD);

    return -PitchAmpRad * ramp *
           sin(2.0 * M_PI * Freq_Pitching * tau + phase_pitching);
#endif
}

/* Fluent zone motion requires angular velocity rather than absolute angle.
 * A backward finite difference of the SAME get_pitch_angle() law guarantees
 * that harmonic and ramp-hold pitch remain consistent with the grid motion. */
static real get_pitch_angular_velocity(real time, real dtime)
{
    real theta_now, theta_old;

    if (dtime <= 0.0)
        return 0.0;

    theta_now = get_pitch_angle(time);
    theta_old = get_pitch_angle(time - dtime);

    return (theta_now - theta_old) / dtime;
}

/* Return the pitch angle actually applied to a surface node.  In the normal
 * configuration every span station receives the same pitch.  When
 * PITCH_ONLY_ON_MORPH_SPANS=1 in 3D, pitch is multiplied by the span weight. */
static real get_grid_pitch_angle(real time, real z)
{
#if GRID_PITCH_ENABLE
    return get_pitch_angle(time) * pitch_span_weight(z);
#else
    (void)time;
    (void)z;
    return 0.0;
#endif
}

/*------------------------------------------------------------
 * ROTATION
 *-----------------------------------------------------------*/

static void rotate_about_center_3d(real x, real y, real z, real theta,
                                   real *xr, real *yr, real *zr)
{
    real cos_t = cos(theta + AoA_rad);
    real sin_t = sin(theta + AoA_rad);

    *xr = cos_t * (x - x_center) - sin_t * (y - y_center) + x_center;
    *yr = sin_t * (x - x_center) + cos_t * (y - y_center) + y_center;
    *zr = z;
}

/* Use this when every node in the hooked zone has the same pitch angle.
 * The sine and cosine are calculated once before the face/node loops.
 */
static void rotate_about_center_precomputed(real x, real y, real z,
                                            real cos_t, real sin_t,
                                            real *xr, real *yr, real *zr)
{
    *xr = cos_t * (x - x_center) - sin_t * (y - y_center) + x_center;
    *yr = sin_t * (x - x_center) + cos_t * (y - y_center) + y_center;
    *zr = z;
}

/*------------------------------------------------------------
 * RIGID-BODY PITCHING: ZONE MOTION
 *-----------------------------------------------------------*/

DEFINE_ZONE_MOTION(pitching_zone_motion, omega, axis, origin, velocity, time, dtime)
{
#if ZONE_PITCH_ENABLE
    *omega = get_pitch_angular_velocity(time, dtime);
#else
    *omega = 0.0;
#endif

#if RP_3D
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;

    origin[0] = x_center;
    origin[1] = y_center;
    origin[2] = 0.0;

    axis[0] = 0.0;
    axis[1] = 0.0;
    axis[2] = 1.0;
#else
    velocity[0] = 0.0;
    velocity[1] = 0.0;

    origin[0] = x_center;
    origin[1] = y_center;

    axis[0] = 0.0;
    axis[1] = 1.0;
#endif
}

/*------------------------------------------------------------
 * WANG LEADING-EDGE MAPPING
 *-----------------------------------------------------------*/

static void wang_le(real x_org, real W, real *xc, real *yc, real *dyc_dx_out)
{
    real denom = safe_denom(W);
    real cofp = 2.0 * sqrt(x_s_le * x_s_le + yLE * yLE) / denom;
    real xi = cofp * ((x_s_le - x_org) / x_s_le);
    real eta = -xi * ((x_s_le - x_org) / x_s_le) * tan(W);
    real phi = -atan(yLE / x_s_le);
    real delta_yc = -yLE * (1.0 - (x_s_le - x_org) / x_s_le);
    real ang = atan(2.0 * (x_s_le - x_org) / x_s_le * tan(W));

    *xc = x_s_le - xi * cos(phi) + eta * sin(phi) - delta_yc * sin(ang);
    *yc = yLE + xi * sin(phi) + eta * cos(phi) + delta_yc * cos(ang);

    *dyc_dx_out = cofp * (tan(W) / (x_s_le * x_s_le)) *
                  2.0 * (x_s_le - x_org);
}

/*------------------------------------------------------------
 * WANG TRAILING-EDGE MAPPING
 *-----------------------------------------------------------*/

/* Common leading-edge mapping for both upper and lower hooks.
 * The previous version let the upper hook own the shared LE node while
 * the lower hook skipped it. That can create a visible glitch because
 * the result depends on hook/node visitation order.
 *
 * This version uses the original node side to assign the signed NACA
 * thickness. A shared node with y_org == yLE gets zero offset, so both
 * hooks calculate exactly the same LE position.
 */
static void map_le_point(real x_org, real y_org, real W, real thickness,
                         real *x_morph, real *y_morph)
{
    real xc, yc, dyc_dx;
    real signed_offset;
    real s;
    const real join_tol = 1.0e-10;

    wang_le(x_org, W, &xc, &yc, &dyc_dx);
    s = sqrt(1.0 + dyc_dx * dyc_dx);

    if (fabs(y_org - yLE) <= join_tol || thickness <= 0.0)
        signed_offset = 0.0;
    else if (y_org > yLE)
        signed_offset = fabs(thickness);
    else
        signed_offset = -fabs(thickness);

    *x_morph = xc - (signed_offset / s) * dyc_dx;
    *y_morph = yc + (signed_offset / s);
}

static void wang_te(real x_org, real W, real *xc, real *yc, real *dyc_dx_out)
{
    real L = chord - x_s;
    real denom = safe_denom(W);
    real R = sqrt(L * L + yTE * yTE);
    real cofp = 2.0 * R / denom;
    real xi = cofp * ((x_org - x_s) / L);
    real eta = -cofp * pow((x_org - x_s) / L, 2.0) * tan(W);
    real phi = -atan(yTE / L);
    real delta_yc = -yTE * (1.0 - (x_org - x_s) / L);
    real ang = atan(2.0 * (x_org - x_s) / L * tan(W));

    *xc = x_s + xi * cos(phi) - eta * sin(phi) + delta_yc * sin(ang);
    *yc = yTE + xi * sin(phi) + eta * cos(phi) + delta_yc * cos(ang);

    *dyc_dx_out = -cofp * (tan(W) / (L * L)) *
                   2.0 * (x_org - x_s);
}

/* Map a point using its actual signed offset from the original TE camber
 * line. This keeps the upper, lower and blunt mappings identical at shared
 * trailing-edge corner nodes.
 */
static void map_te_point(real x_org, real y_org, real W,
                         real *x_morph, real *y_morph)
{
    real xc, yc, dyc_dx;
    real signed_offset = y_org - yTE;
    real s;

    wang_te(x_org, W, &xc, &yc, &dyc_dx);
    s = sqrt(1.0 + dyc_dx * dyc_dx);

    *x_morph = xc - (signed_offset / s) * dyc_dx;
    *y_morph = yc + (signed_offset / s);
}

/*------------------------------------------------------------
 * UPPER SURFACE MOTION
 *-----------------------------------------------------------*/

DEFINE_GRID_MOTION(morphing_upper, domain, dt, time, dtime)
{
#if UDF_PROFILE_ENABLE
    static long profile_call_id = 0;
    clock_t profile_total_start = clock();
    clock_t profile_init_ticks = 0;
    clock_t profile_le_ticks = 0;
    clock_t profile_te_ticks = 0;
    clock_t profile_rotate_ticks = 0;
    long profile_face_node_visits = 0;
    long profile_updated_nodes = 0;
    long profile_le_nodes = 0;
    long profile_te_nodes = 0;
#endif
    Thread *tf = DT_THREAD(dt);
    face_t f;
    int n;
    const real tol_x = 1.0e-6;
    const int motion_active = (time > FTT);
    const real W_le_now = motion_active ? get_le_angle(time) : 0.0;
    const real W_te_now = motion_active ? -get_te_angle(time) : 0.0;
#if GRID_PITCH_ENABLE
    const real pitch_now = get_pitch_angle(time);
#else
    const real pitch_now = 0.0;
#endif
#if !PITCH_ONLY_ON_MORPH_SPANS
    const real pitch_cos = cos(pitch_now + AoA_rad);
    const real pitch_sin = sin(pitch_now + AoA_rad);
#endif

    (void)domain;
    (void)dtime;

    SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));

    begin_f_loop(f, tf)
    {
        f_node_loop(f, tf, n)
        {
            Node *node_p = F_NODE(f, tf, n);
#if UDF_PROFILE_ENABLE
            profile_face_node_visits++;
#endif

            if (NODE_POS_NEED_UPDATE(node_p))
            {
                real x_org, y_org, z_org;
                real x_local, y_local, z_local;

#if UDF_PROFILE_ENABLE
                profile_updated_nodes++;
#endif
#if UDF_PROFILE_SECTIONS
                {
                    clock_t profile_t0 = clock();
                    ensure_initial_coords(node_p);
                    profile_init_ticks += clock() - profile_t0;
                }
#else
                ensure_initial_coords(node_p);
#endif

                x_org = N_UDMI(node_p, ORIG_X_UDMI);
                y_org = N_UDMI(node_p, ORIG_Y_UDMI);
#if GEOMETRY_MODE_3D && RP_3D
                z_org = N_UDMI(node_p, ORIG_Z_UDMI);
#else
                z_org = 0.0;
#endif

                NODE_POS_UPDATED(node_p);

                x_local = x_org;
                y_local = y_org;
                z_local = z_org;

                if (LE_ENABLE && x_org <= (x_s_le + tol_x) && motion_active)
                {
#if UDF_PROFILE_ENABLE
                    profile_le_nodes++;
#endif
#if UDF_PROFILE_SECTIONS
                    clock_t profile_t0 = clock();
#endif
                    real w_span = le_span_weight(z_org);
                    real t = N_UDMI(node_p, THICKNESS_UDMI);
                    real W_le = W_le_now;
                    real x_morph, y_morph;

                    map_le_point(x_org, y_org, W_le, t, &x_morph, &y_morph);

                    x_local = x_org + w_span * (x_morph - x_org);
                    y_local = y_org + w_span * (y_morph - y_org);
#if UDF_PROFILE_SECTIONS
                    profile_le_ticks += clock() - profile_t0;
#endif
                }

                if (TE_ENABLE && x_org >= (x_s - tol_x) && motion_active)
                {
                    real W_te = W_te_now;

                    if (fabs(W_te) > TE_ANGLE_EPS)
                    {
#if UDF_PROFILE_ENABLE
                        profile_te_nodes++;
#endif
#if UDF_PROFILE_SECTIONS
                        clock_t profile_t0 = clock();
#endif
                        real w_span = te_span_weight(z_org);
                        real x_morph, y_morph;

                        map_te_point(x_org, y_org, W_te, &x_morph, &y_morph);

                        x_local = x_org + w_span * (x_morph - x_org);
                        y_local = y_org + w_span * (y_morph - y_org);
#if UDF_PROFILE_SECTIONS
                        profile_te_ticks += clock() - profile_t0;
#endif
                    }
                }

                if (MY_ISNAN(x_local) || MY_ISNAN(y_local) || MY_ISNAN(z_local))
                {
                    x_local = x_org;
                    y_local = y_org;
                    z_local = z_org;
                }

                {
#if UDF_PROFILE_SECTIONS
                    clock_t profile_t0 = clock();
#endif
                    real x_rot, y_rot, z_rot;

#if PITCH_ONLY_ON_MORPH_SPANS
                    real theta = pitch_now * pitch_span_weight(z_org);

                    rotate_about_center_3d(x_local, y_local, z_local, theta,
                                           &x_rot, &y_rot, &z_rot);
#else
                    rotate_about_center_precomputed(x_local, y_local, z_local,
                                                    pitch_cos, pitch_sin,
                                                    &x_rot, &y_rot, &z_rot);
#endif

                    NODE_X(node_p) = x_rot;
                    NODE_Y(node_p) = y_rot;
                    SET_NODE_Z(node_p, z_rot);
#if UDF_PROFILE_SECTIONS
                    profile_rotate_ticks += clock() - profile_t0;
#endif
                }
            }
        }
    }
    end_f_loop(f, tf);
#if UDF_PROFILE_ENABLE
    profile_call_id++;
    if ((profile_call_id % UDF_PROFILE_PRINT_EVERY) == 0)
    {
        udf_profile_report("UPPER", time, profile_call_id,
                           clock() - profile_total_start,
                           profile_init_ticks, profile_le_ticks,
                           profile_te_ticks, profile_rotate_ticks,
                           profile_face_node_visits, profile_updated_nodes,
                           profile_le_nodes, profile_te_nodes);
    }
#endif
}

/*------------------------------------------------------------
 * LOWER SURFACE MOTION
 *-----------------------------------------------------------*/

DEFINE_GRID_MOTION(morphing_lower, domain, dt, time, dtime)
{
#if UDF_PROFILE_ENABLE
    static long profile_call_id = 0;
    clock_t profile_total_start = clock();
    clock_t profile_init_ticks = 0;
    clock_t profile_le_ticks = 0;
    clock_t profile_te_ticks = 0;
    clock_t profile_rotate_ticks = 0;
    long profile_face_node_visits = 0;
    long profile_updated_nodes = 0;
    long profile_le_nodes = 0;
    long profile_te_nodes = 0;
#endif
    Thread *tf = DT_THREAD(dt);
    face_t f;
    int n;
    const real tol_x = 1.0e-6;
    const int motion_active = (time > FTT);
    const real W_le_now = motion_active ? get_le_angle(time) : 0.0;
    const real W_te_now = motion_active ? -get_te_angle(time) : 0.0;
#if GRID_PITCH_ENABLE
    const real pitch_now = get_pitch_angle(time);
#else
    const real pitch_now = 0.0;
#endif
#if !PITCH_ONLY_ON_MORPH_SPANS
    const real pitch_cos = cos(pitch_now + AoA_rad);
    const real pitch_sin = sin(pitch_now + AoA_rad);
#endif

    (void)domain;
    (void)dtime;

    SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));

    begin_f_loop(f, tf)
    {
        f_node_loop(f, tf, n)
        {
            Node *node_p = F_NODE(f, tf, n);
#if UDF_PROFILE_ENABLE
            profile_face_node_visits++;
#endif

            if (NODE_POS_NEED_UPDATE(node_p))
            {
                real x_org, y_org, z_org;
                real x_local, y_local, z_local;

#if UDF_PROFILE_ENABLE
                profile_updated_nodes++;
#endif
#if UDF_PROFILE_SECTIONS
                {
                    clock_t profile_t0 = clock();
                    ensure_initial_coords(node_p);
                    profile_init_ticks += clock() - profile_t0;
                }
#else
                ensure_initial_coords(node_p);
#endif

                x_org = N_UDMI(node_p, ORIG_X_UDMI);
                y_org = N_UDMI(node_p, ORIG_Y_UDMI);
#if GEOMETRY_MODE_3D && RP_3D
                z_org = N_UDMI(node_p, ORIG_Z_UDMI);
#else
                z_org = 0.0;
#endif

                /* Shared LE nodes are handled by the same LE mapping as
                 * upper/lower surface nodes. This removes hook-order dependence. */
                NODE_POS_UPDATED(node_p);

                x_local = x_org;
                y_local = y_org;
                z_local = z_org;

                if (LE_ENABLE && x_org <= (x_s_le + tol_x) && motion_active)
                {
#if UDF_PROFILE_ENABLE
                    profile_le_nodes++;
#endif
#if UDF_PROFILE_SECTIONS
                    clock_t profile_t0 = clock();
#endif
                    real w_span = le_span_weight(z_org);
                    real t = N_UDMI(node_p, THICKNESS_UDMI);
                    real W_le = W_le_now;
                    real x_morph, y_morph;

                    map_le_point(x_org, y_org, W_le, t, &x_morph, &y_morph);

                    x_local = x_org + w_span * (x_morph - x_org);
                    y_local = y_org + w_span * (y_morph - y_org);
#if UDF_PROFILE_SECTIONS
                    profile_le_ticks += clock() - profile_t0;
#endif
                }

                if (TE_ENABLE && x_org >= (x_s - tol_x) && motion_active)
                {
                    real W_te = W_te_now;

                    if (fabs(W_te) > TE_ANGLE_EPS)
                    {
#if UDF_PROFILE_ENABLE
                        profile_te_nodes++;
#endif
#if UDF_PROFILE_SECTIONS
                        clock_t profile_t0 = clock();
#endif
                        real w_span = te_span_weight(z_org);
                        real x_morph, y_morph;

                        map_te_point(x_org, y_org, W_te, &x_morph, &y_morph);

                        x_local = x_org + w_span * (x_morph - x_org);
                        y_local = y_org + w_span * (y_morph - y_org);
#if UDF_PROFILE_SECTIONS
                        profile_te_ticks += clock() - profile_t0;
#endif
                    }
                }

                if (MY_ISNAN(x_local) || MY_ISNAN(y_local) || MY_ISNAN(z_local))
                {
                    x_local = x_org;
                    y_local = y_org;
                    z_local = z_org;
                }

                {
#if UDF_PROFILE_SECTIONS
                    clock_t profile_t0 = clock();
#endif
                    real x_rot, y_rot, z_rot;

#if PITCH_ONLY_ON_MORPH_SPANS
                    real theta = pitch_now * pitch_span_weight(z_org);

                    rotate_about_center_3d(x_local, y_local, z_local, theta,
                                           &x_rot, &y_rot, &z_rot);
#else
                    rotate_about_center_precomputed(x_local, y_local, z_local,
                                                    pitch_cos, pitch_sin,
                                                    &x_rot, &y_rot, &z_rot);
#endif

                    NODE_X(node_p) = x_rot;
                    NODE_Y(node_p) = y_rot;
                    SET_NODE_Z(node_p, z_rot);
#if UDF_PROFILE_SECTIONS
                    profile_rotate_ticks += clock() - profile_t0;
#endif
                }
            }
        }
    }
    end_f_loop(f, tf);
#if UDF_PROFILE_ENABLE
    profile_call_id++;
    if ((profile_call_id % UDF_PROFILE_PRINT_EVERY) == 0)
    {
        udf_profile_report("LOWER", time, profile_call_id,
                           clock() - profile_total_start,
                           profile_init_ticks, profile_le_ticks,
                           profile_te_ticks, profile_rotate_ticks,
                           profile_face_node_visits, profile_updated_nodes,
                           profile_le_nodes, profile_te_nodes);
    }
#endif
}

/*------------------------------------------------------------
 * BLUNT TRAILING EDGE: SINGLE HOOK
 *-----------------------------------------------------------
 * Hook this function to the complete blunt trailing-edge face zone.
 * Every node in that hooked face zone is handled by this motion.
 */

DEFINE_GRID_MOTION(morphing_blunt, domain, dt, time, dtime)
{
#if UDF_PROFILE_ENABLE
    static long profile_call_id = 0;
    clock_t profile_total_start = clock();
    clock_t profile_init_ticks = 0;
    clock_t profile_le_ticks = 0;
    clock_t profile_te_ticks = 0;
    clock_t profile_rotate_ticks = 0;
    long profile_face_node_visits = 0;
    long profile_updated_nodes = 0;
    long profile_le_nodes = 0;
    long profile_te_nodes = 0;
#endif
    Thread *tf = DT_THREAD(dt);
    face_t f;
    int n;
    const int motion_active = (time > FTT);
    const real W_te_now = motion_active ? -get_te_angle(time) : 0.0;
#if GRID_PITCH_ENABLE
    const real pitch_now = get_pitch_angle(time);
#else
    const real pitch_now = 0.0;
#endif
#if !PITCH_ONLY_ON_MORPH_SPANS
    const real pitch_cos = cos(pitch_now + AoA_rad);
    const real pitch_sin = sin(pitch_now + AoA_rad);
#endif

    (void)domain;
    (void)dtime;

    SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));

    begin_f_loop(f, tf)
    {
        f_node_loop(f, tf, n)
        {
            Node *node_p = F_NODE(f, tf, n);
#if UDF_PROFILE_ENABLE
            profile_face_node_visits++;
#endif

            if (NODE_POS_NEED_UPDATE(node_p))
            {
                real x_org, y_org, z_org;
                real x_local, y_local, z_local;

#if UDF_PROFILE_ENABLE
                profile_updated_nodes++;
#endif
#if UDF_PROFILE_SECTIONS
                {
                    clock_t profile_t0 = clock();
                    ensure_initial_coords(node_p);
                    profile_init_ticks += clock() - profile_t0;
                }
#else
                ensure_initial_coords(node_p);
#endif

                x_org = N_UDMI(node_p, ORIG_X_UDMI);
                y_org = N_UDMI(node_p, ORIG_Y_UDMI);
#if GEOMETRY_MODE_3D && RP_3D
                z_org = N_UDMI(node_p, ORIG_Z_UDMI);
#else
                z_org = 0.0;
#endif

                NODE_POS_UPDATED(node_p);

                x_local = x_org;
                y_local = y_org;
                z_local = z_org;

                if (TE_ENABLE && motion_active)
                {
                    real W_te = W_te_now;

                    if (fabs(W_te) > TE_ANGLE_EPS)
                    {
#if UDF_PROFILE_ENABLE
                        profile_te_nodes++;
#endif
#if UDF_PROFILE_SECTIONS
                        clock_t profile_t0 = clock();
#endif
                        real w_span = te_span_weight(z_org);
                        real x_morph, y_morph;

                        map_te_point(x_org, y_org, W_te, &x_morph, &y_morph);

                        x_local = x_org + w_span * (x_morph - x_org);
                        y_local = y_org + w_span * (y_morph - y_org);
#if UDF_PROFILE_SECTIONS
                        profile_te_ticks += clock() - profile_t0;
#endif
                    }
                }

                if (MY_ISNAN(x_local) || MY_ISNAN(y_local) || MY_ISNAN(z_local))
                {
                    x_local = x_org;
                    y_local = y_org;
                    z_local = z_org;
                }

                {
#if UDF_PROFILE_SECTIONS
                    clock_t profile_t0 = clock();
#endif
                    real x_rot, y_rot, z_rot;

#if PITCH_ONLY_ON_MORPH_SPANS
                    real theta = pitch_now * pitch_span_weight(z_org);

                    rotate_about_center_3d(x_local, y_local, z_local, theta,
                                           &x_rot, &y_rot, &z_rot);
#else
                    rotate_about_center_precomputed(x_local, y_local, z_local,
                                                    pitch_cos, pitch_sin,
                                                    &x_rot, &y_rot, &z_rot);
#endif

                    NODE_X(node_p) = x_rot;
                    NODE_Y(node_p) = y_rot;
                    SET_NODE_Z(node_p, z_rot);
#if UDF_PROFILE_SECTIONS
                    profile_rotate_ticks += clock() - profile_t0;
#endif
                }
            }
        }
    }
    end_f_loop(f, tf);
#if UDF_PROFILE_ENABLE
    profile_call_id++;
    if ((profile_call_id % UDF_PROFILE_PRINT_EVERY) == 0)
    {
        udf_profile_report("BLUNT", time, profile_call_id,
                           clock() - profile_total_start,
                           profile_init_ticks, profile_le_ticks,
                           profile_te_ticks, profile_rotate_ticks,
                           profile_face_node_visits, profile_updated_nodes,
                           profile_le_nodes, profile_te_nodes);
    }
#endif
}
