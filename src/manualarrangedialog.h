#ifndef MANUALARRANGEDIALOG_H
#define MANUALARRANGEDIALOG_H
#include <QDialog>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QPushButton>
#include "celldata.h"

class ManualArrangeDialog : public QDialog {
    Q_OBJECT
public:
    explicit ManualArrangeDialog(QWidget *parent = nullptr);
    void setDefaultParams(double cSize, double nSize, double mZ);
    QList<CompleteCell> getFinalCells();

private slots:
    void addCell();
    void removeCell();

private:
    QTableWidget *table;
    double defCellSize, defNucSize, defMajorZ;
};

#endif // MANUALARRANGEDIALOG_H
