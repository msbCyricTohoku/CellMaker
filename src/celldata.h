#ifndef CELLDATA_H
#define CELLDATA_H

struct CompleteCell {
    double x, y, z;        // cytoplasm Center
    double rx, rz;         // cytoplasm Radii
    double majorX, majorY; // scaling factors

    double nx, ny, nz;     // nucleus Center
    double nrx, nrz;       // nucleus Radii

    int cellSurfId;        // id for PHITS surface
    int nucSurfId;         // id for PHITS nucleus
};

#endif // CELLDATA_H
