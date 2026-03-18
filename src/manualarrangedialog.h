#ifndef MANUALARRANGEDIALOG_H
#define MANUALARRANGEDIALOG_H
#include <QDialog>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QPushButton>
#include <QDialog>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QList>
#include <QComboBox>
#include "celldata.h"

class ManualArrangeDialog : public QDialog {
    Q_OBJECT

public:
    explicit ManualArrangeDialog(QWidget *parent = nullptr);

    void setDefaultParams(double cSize, double nSize,
                          int defaultXYCount, double defaultPitch,
                          double randomRadius, double defaultCellRatio, double defaultNucRatio);

    QList<CompleteCell> getFinalCells();

private:
    QSpinBox *spinX;
    QSpinBox *spinY;
    QSpinBox *spinZ;
    QDoubleSpinBox *spinPitchZ;
    QComboBox *comboDistType;
    QDoubleSpinBox *spinRandomRadius;
    QDoubleSpinBox *spinDepth;

    // here the ratio controls
    QDoubleSpinBox *spinCellZRatio;
    QDoubleSpinBox *spinNucZRatio;

    double defCellSize;
    double defNucSize;
    double defPitch;
};

#endif // MANUALARRANGEDIALOG_H
