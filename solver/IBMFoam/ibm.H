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
    ibForce = IBM.ibForce(U);  //计算ibforce
}
else
{
    ibForce = IBM.ibForceInt();
}

solve(UEqn == -fvc::grad(p) + ibForce);

IBM.multiDirectForcing(U, ibForce);

IBM.write();

runTime.write();

meshRefine = IBM.calculateAlphaSolid();  //计算自适应区域

mesh.update(); //更新网格

IBM.update();  //更新固体单元，网格info等等