// Solve the Momentum equation

fvVectorMatrix UEqn
(
    fvm::ddt(U) 
  + fvm::div(phi, U)
  - fvm::laplacian(nu, U)
);

UEqn.relax();

//volVectorField gradP = IBM.pressureGradientField();

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p) );//+ gradP);
}

// 假设从 controlDict 中读取关键字，可根据实际情况修改
IOdictionary IBMDict
(
    IOobject
    (
        "IBMDict",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
const word ibForceMethod = IBMDict.subDict("IBModel").lookupOrDefault<word>("ibForceMethod", "ibForce");

if (ibForceMethod == "ibForce")
{
    ibForce = IBM.ibForce(U);
}
else if (ibForceMethod == "ibForceInt")
{
    ibForce = IBM.ibForceInt();
}
else
{
    FatalErrorInFunction
        << "Invalid ibForceMethod specified in dictionary. "
        << "Valid options are 'ibForce' and 'ibForceInt'."
        << abort(FatalError);
}
//volVectorField ibForce = IBM.ibForce(U);
//ibForce = IBM.ibForceInt(); // 从 IBM 对象获取 ibForce 数据

solve(UEqn == -fvc::grad(p) + ibForce);

//IBM.multiDirectForcing(U, ibForce);
