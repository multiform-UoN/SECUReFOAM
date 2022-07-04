
vector cellSize(0,0,0);

const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

forAll ( mesh.C(), celli)
{
    const cell & cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi)
           pLocal[pointi] = pp[pLabels[pointi]];

    cellSize[0] = max(cellSize[0], Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0)) );
    cellSize[1] = max(cellSize[1], Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0)) );
    cellSize[2] = max(cellSize[2], Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1)) );
}