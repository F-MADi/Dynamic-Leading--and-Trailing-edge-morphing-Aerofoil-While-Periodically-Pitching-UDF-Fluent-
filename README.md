# Dynamic-Leading--and-Trailing-edge-morphing-Aerofoil-While-Periodically-Pitching-UDF-Fluent-

[UNDER CONSTRUCTION]

This repository describes and demonstrates the use of a User-Defined Function (UDF) compiled in Fluent Ansys (launched in a Visual Studio environment) to perform an unsteady computational fluid dynamics (CFD) study on the dynamically morphing leading edge (LE) and trailing edge (TE) of an aerofoil while undergoing periodic pitching.

This repository is inspired by https://github.com/chawkiabd/Dynamic-Morphing-Wing, building on this framework.

This Repository uses a second-order polynomial (modified from [1]) to deform the aerofoil shape by introducing a harmonic function to control the deflection amplitude of the morphing LE-TE sections. This includes parameters such as frequency, amplitude, and location where morphing is implemented. A separate UDF is needed to control the periodic pitching motion using a zone motion approach, while the main UDF is needed to control morphing sections.

Geometry Requirements: Two zones, separated by an interface where the inner zone is centred at the pitching motion location (e.g. 0.25C). This will allow you to perfectly hook the zone to the UDF “XXX”. While the larger zone remains static. The aerofoil needs to be a 4-digit NACAXXXX with separated surface names (e.g., Upper, Lower, Blunt) to allow you to hook them after code completion.

Process: Launch Fluent using a Visual Studio environment. Allocate user-defined memory (7) to save grid node locations. Copy and save "Morphing-LE-TE-Pitching-Aerofoil" and "ZONE_MOTION" into Notepad and save as a .c file. Compile the downloaded .c file “XXXX” (not using the built-in compiler), and press Load. Initialise storing the node coordinates by hooking “XXX” UDF output run-on-demand and active (confirm by reading UI output the number of stored nodes) Hooking zone UDF “XXX” to your rotating zone domain and hooking your upper, lower, and blunt trailing-edge surfaces to their corresponding names.

Published open-access paper linked to example UDF seen in: A computational study on the aerodynamic properties of a pitching aerofoil with a morphing leading edge section at high Reynolds numbers. https://doi.org/10.1016/j.ijheatfluidflow.2026.110531

[1] Wang, R., Ma, X., Zhang, G., Ying, P. and Wang, X., 2023. Numerical simulation of continuous morphing wing with leading edge and trailing edge parabolic flaps. Journal of Aerospace Engineering, 36(5), p.04023051.

Author: Fakhreddine Madi (Fakhreddine.Madi@uwe.ac.uk)
