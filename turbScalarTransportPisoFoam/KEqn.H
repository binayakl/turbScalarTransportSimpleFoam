volScalarField DKK ("DKK", DK + turbulence->nut()/Sct) ;

while (piso.correctNonOrthogonal())
        {
            fvScalarMatrix KEqn
            (
                fvm::ddt(K)
              + fvm::div(phi, K)
              - fvm::laplacian(DKK, K)
             ==
                fvOptions(K)
            );

            KEqn.relax();
            fvOptions.constrain(KEqn);
            KEqn.solve();
            fvOptions.correct(K);
        }
