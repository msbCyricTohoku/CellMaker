#include "manualarrangedialog.h"
#include <QVBoxLayout>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QRandomGenerator>
#include <QMessageBox>
#include <cmath>

ManualArrangeDialog::ManualArrangeDialog(QWidget *parent) : QDialog(parent) {
    setWindowTitle("3D In-Vivo Cell Array Setup");
    resize(420, 380);

    QVBoxLayout *mainLayout = new QVBoxLayout(this);
    QGridLayout *gridLayout = new QGridLayout();

    gridLayout->addWidget(new QLabel("Distribution Type:"), 0, 0);
    comboDistType = new QComboBox(this);
    comboDistType->addItem("Uniform");
    comboDistType->addItem("Random");
    gridLayout->addWidget(comboDistType, 0, 1);

    gridLayout->addWidget(new QLabel("Random Radius (um):"), 1, 0);
    spinRandomRadius = new QDoubleSpinBox(this);
    spinRandomRadius->setRange(0, 1000);
    spinRandomRadius->setEnabled(false);
    gridLayout->addWidget(spinRandomRadius, 1, 1);

    gridLayout->addWidget(new QLabel("Depth In-Vivo (um):"), 2, 0);
    spinDepth = new QDoubleSpinBox(this);
    spinDepth->setRange(0, 50000);
    spinDepth->setValue(10.0);
    gridLayout->addWidget(spinDepth, 2, 1);

    gridLayout->addWidget(new QLabel("Cells in X-axis:"), 3, 0);
    spinX = new QSpinBox(this); spinX->setRange(1, 100);
    gridLayout->addWidget(spinX, 3, 1);

    gridLayout->addWidget(new QLabel("Cells in Y-axis:"), 4, 0);
    spinY = new QSpinBox(this); spinY->setRange(1, 100);
    gridLayout->addWidget(spinY, 4, 1);

    gridLayout->addWidget(new QLabel("Cells in Z-axis (Layers):"), 5, 0);
    spinZ = new QSpinBox(this); spinZ->setRange(1, 100); spinZ->setValue(3);
    gridLayout->addWidget(spinZ, 5, 1);

    gridLayout->addWidget(new QLabel("Z-Axis Pitch (um):"), 6, 0);
    spinPitchZ = new QDoubleSpinBox(this); spinPitchZ->setRange(0, 1000);
    gridLayout->addWidget(spinPitchZ, 6, 1);

    // ratios
    gridLayout->addWidget(new QLabel("Cell Flatness Ratio:"), 7, 0);
    spinCellZRatio = new QDoubleSpinBox(this);
    spinCellZRatio->setRange(0.01, 2.0);
    spinCellZRatio->setSingleStep(0.1);
    gridLayout->addWidget(spinCellZRatio, 7, 1);

    gridLayout->addWidget(new QLabel("Nucleus Flatness Ratio:"), 8, 0);
    spinNucZRatio = new QDoubleSpinBox(this);
    spinNucZRatio->setRange(0.01, 2.0);
    spinNucZRatio->setSingleStep(0.1);
    gridLayout->addWidget(spinNucZRatio, 8, 1);

    mainLayout->addLayout(gridLayout);

    QPushButton *btnOk = new QPushButton("Generate 3D Array", this);
    btnOk->setMinimumHeight(35);
    btnOk->setStyleSheet("font-weight: bold; background-color: #0078D7; color: white;");

    connect(comboDistType, &QComboBox::currentTextChanged, [=](const QString &text) {
        spinRandomRadius->setEnabled(text == "Random");
    });

    connect(btnOk, &QPushButton::clicked, [=]() {
        accept();
    });

    mainLayout->addStretch();
    mainLayout->addWidget(btnOk);
}

void ManualArrangeDialog::setDefaultParams(double cSize, double nSize,
                                           int defaultXYCount, double defaultPitch,
                                           double randomRadius, double defaultCellRatio, double defaultNucRatio) {
    defCellSize = cSize;
    defNucSize = nSize;
    defPitch = defaultPitch;

    spinX->setValue(defaultXYCount);
    spinY->setValue(defaultXYCount);
    spinPitchZ->setValue(defaultPitch);
    spinRandomRadius->setValue(randomRadius);
    spinCellZRatio->setValue(defaultCellRatio);
    spinNucZRatio->setValue(defaultNucRatio);
}

QList<CompleteCell> ManualArrangeDialog::getFinalCells() {
    QList<CompleteCell> cells;
    QRandomGenerator *gen = QRandomGenerator::global();

    int nxCount = spinX->value();
    int nyCount = spinY->value();
    int nzCount = spinZ->value();
    double depth = spinDepth->value();

    bool isRandom = (comboDistType->currentText() == "Random");
    double userRadius = spinRandomRadius->value();
    const int maxAttempts = 2000;

    // calc Z dimens. based on the user ratio
    double actualCellRz = (defCellSize / 2.0) * spinCellZRatio->value();
    double actualNucRz = (defNucSize / 2.0) * spinNucZRatio->value();

    double safeCollisionDist = std::max(defCellSize, actualCellRz * 2.0);

    double spacingXY = defCellSize + defPitch;
    double spacingZ = (actualCellRz * 2.0) + spinPitchZ->value();

    double halfSpanX = (nxCount - 1) * spacingXY / 2.0;
    double halfSpanY = (nyCount - 1) * spacingXY / 2.0;

    double zOffset = isRandom ? userRadius : 0.0;

    int idCounter = 1;
    bool placementOk = true;

    for (int k = 0; k < nzCount && placementOk; ++k) {
        for (int i = 0; i < nxCount && placementOk; ++i) {
            for (int j = 0; j < nyCount && placementOk; ++j) {

                double baseCx = (i * spacingXY) - halfSpanX;
                double baseCy = (j * spacingXY) - halfSpanY;
                double baseCz = -depth - actualCellRz - zOffset - (k * spacingZ);

                double cx = baseCx;
                double cy = baseCy;
                double cz = baseCz;

                if (isRandom) {
                    int attempts = 0;
                    bool placed = false;
                    while (attempts < maxAttempts && !placed) {
                        double u, v, w, distSq;
                        do {
                            u = (gen->generateDouble() * 2.0) - 1.0;
                            v = (gen->generateDouble() * 2.0) - 1.0;
                            w = (gen->generateDouble() * 2.0) - 1.0;
                            distSq = (u*u) + (v*v) + (w*w);
                        } while (distSq > 1.0 || distSq == 0.0);

                        double candCx = baseCx + (u * userRadius);
                        double candCy = baseCy + (v * userRadius);
                        double candCz = baseCz + (w * userRadius);

                        bool validPlacement = true;
                        for (const CompleteCell &placedCell : cells) {
                            double dx = candCx - placedCell.x;
                            double dy = candCy - placedCell.y;
                            double dz = candCz - placedCell.z;
                            if (std::sqrt(dx*dx + dy*dy + dz*dz) < safeCollisionDist) {
                                validPlacement = false;
                                break;
                            }
                        }

                        if (validPlacement) {
                            cx = candCx;
                            cy = candCy;
                            cz = candCz;
                            placed = true;
                        }
                        ++attempts;
                    }

                    if (!placed) {
                        QMessageBox::warning(this, "Placement Error",
                                             "Unable to place cells without overlap. Try changing the Random Radius or the Cell Pitch.");
                        placementOk = false;
                        cells.clear();
                        break;
                    }
                }

                CompleteCell cc;
                cc.x = cx;
                cc.y = cy;
                cc.z = cz;

                cc.rx = defCellSize / 2.0;
                cc.rz = actualCellRz;

                cc.majorX = 0.0;
                cc.majorY = 0.0;

                cc.nrx = defNucSize / 2.0;
                cc.nrz = actualNucRz;

                //normalized scale factor to ensure nucleus does not go out of bound
                double ratioX = cc.nrx / cc.rx;
                double ratioZ = cc.nrz / cc.rz;
                double S = std::max(ratioX, ratioZ); // find the tightest constraint

                double safeRx = cc.rx * (1.0 - S) - 0.000001;
                double safeRz = cc.rz * (1.0 - S) - 0.000001;

                if (safeRx < 0) safeRx = 0;
                if (safeRz < 0) safeRz = 0;

                double u, v, w, distSq;
                do {
                    u = (gen->generateDouble() * 2.0) - 1.0;
                    v = (gen->generateDouble() * 2.0) - 1.0;
                    w = (gen->generateDouble() * 2.0) - 1.0;
                    distSq = (u*u) + (v*v) + (w*w);
                } while (distSq > 1.0);

                cc.nx = cc.x + (u * safeRx);
                cc.ny = cc.y + (v * safeRx);
                cc.nz = cc.z + (w * safeRz);

                cc.cellSurfId = idCounter;
                cc.nucSurfId = idCounter + (nxCount * nyCount * nzCount);
                cells.append(cc);

                idCounter++;
            }
        }
    }
    return cells;
}
