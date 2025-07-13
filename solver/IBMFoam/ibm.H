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
else
{
    ibForce = IBM.ibForceInt();
}

solve(UEqn == -fvc::grad(p) + ibForce);

IBM.multiDirectForcing(U, ibForce);


IBM.update();

alphaSolid = IBM.calculateAlphaSolid();

mesh.update();
