# Absorbing Boundaries
AxiSEM3D features a combination of an absorbing boundary condition (Clayton & Engquist) and a sponge boundary (Kosloff & Kosloff). These boundary conditions can be applied to any desired subset of medium boundaries (top, bottom and right). Note that the left edge does not need absorption as this is the axial boundary which is always placed in the center of the simulated medium.

The Clayton & Engquist boundary is merely an on-off switch and does not affect the model in any way. 

The Kosloff & Kosloff sponge boundary requires a bit of thought: The boundary consists of an attenuating region which extends into the modeled medium. If this affects critical regions of the model, mesh and model have to be expanded as necessary. There are two user-defined input parameters required for the sponge boundary:
1. Boundary width. The effectiveness of the Kosloff & Kosloff boundary scales directly with the boundary width relative to the maximum wavelength within the absorbing region. A noticeable increase of absorption (compared to a model using Clayton & Engquist boundaries alone) can be expected from a sponge boundary width of 2 wavelengths, although the performance of the absorbing boundary is hard to predict and always depends on the model. 99% absorption is reliably achieved at a width of 10 wavelengths. 
2. Maximum attenuating factor. The ideal value for absorption in the sponge has been shown to depend on 3D medium properties, sponge boundary width and other factors. Hence, there is an option to use parsed equations for the maximum attenuation inside the sponge in solid and fluid media. Extensive testing has been conducted to find an optimum attenuation value. The default expressions provided in *inparam.model.yaml* have been shown to work well if the sponge boundary and the Clayton & Engquist boundary are used simultaneously.

> **Note on Efficiency.**  
>The ABC has poor stability properties under some circumstances. Hence the automatically determined time step (explained in the Sources section) is made smaller if this ABC is detected.
