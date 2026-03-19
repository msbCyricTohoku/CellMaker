#ifndef CELLDATA_H
#define CELLDATA_H

struct CompleteCell {
    double x, y, z;        // cytoplasm center
    double rx, rz;         // cytoplasm radii
    double majorX, majorY; // scaling factors

    double nx, ny, nz;     // nucleus center
    double nrx, nrz;       // nucleus radii

    int cellSurfId;        // id for PHITS surface
    int nucSurfId;         // id for PHITS nucleus
};

#endif // CELLDATA_H
