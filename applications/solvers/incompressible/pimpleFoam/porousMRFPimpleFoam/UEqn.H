// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  ==
    sources(U)
);

mrfZones.addCoriolis(UEqn());

UEqn().relax();
        
sources.constrain(UEqn());

pZones.addResistance(UEqn());

volScalarField rAU(1.0/UEqn().A());

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));
}
