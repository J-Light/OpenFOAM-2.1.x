    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    mrfZones.addCoriolis(rhok, UEqn());

    UEqn().relax();

    pZones.addResistance(UEqn());

    if (simple.momentumPredictor())
    {
        solve
        (
            UEqn()
            ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
    }
