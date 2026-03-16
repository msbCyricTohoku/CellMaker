#include "manualarrangedialog.h"
#include <QHeaderView>
#include <QMessageBox>

ManualArrangeDialog::ManualArrangeDialog(QWidget *parent) : QDialog(parent) {
    setWindowTitle("Manual 3D Cell Arrangement");
    resize(700, 450); // slightly wider to accommodate columns

    QVBoxLayout *layout = new QVBoxLayout(this);
    table = new QTableWidget(0, 5, this);
    table->setHorizontalHeaderLabels({"X (um)", "Y (um)", "Z (um)", "Cell Dia.", "Nuc Dia."});

    // make table columns stretch evenly
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    // setup buttons
    QHBoxLayout *btnLayout1 = new QHBoxLayout();
    QPushButton *btnAdd = new QPushButton("Add New Cell", this);
    QPushButton *btnDup = new QPushButton("Duplicate Selected", this);
    QPushButton *btnRemove = new QPushButton("Remove Selected", this);
    QPushButton *btnClear = new QPushButton("Clear All", this);

    btnLayout1->addWidget(btnAdd);
    btnLayout1->addWidget(btnDup);
    btnLayout1->addWidget(btnRemove);
    btnLayout1->addWidget(btnClear);

    QPushButton *btnOk = new QPushButton("Confirm Arrangement", this);
    btnOk->setMinimumHeight(35);
    btnOk->setStyleSheet("font-weight: bold; background-color: #0078D7; color: white;");

    layout->addWidget(table);
    layout->addLayout(btnLayout1);
    layout->addWidget(btnOk);

    // original connections
    connect(btnAdd, &QPushButton::clicked, this, &ManualArrangeDialog::addCell);
    connect(btnRemove, &QPushButton::clicked, this, &ManualArrangeDialog::removeCell);

    // duplicate selected row
    connect(btnDup, &QPushButton::clicked, [=]() {
        int currentRow = table->currentRow();
        if (currentRow >= 0) {
            int newRow = table->rowCount();
            table->insertRow(newRow);
            // copy all 5 columns to new row
            for(int col = 0; col < 5; ++col) {
                QString existingText = table->item(currentRow, col) ? table->item(currentRow, col)->text() : "0.0";
                table->setItem(newRow, col, new QTableWidgetItem(existingText));
            }
            table->selectRow(newRow); // auto-select duplicated row
        }
    });

    // clear all rows
    connect(btnClear, &QPushButton::clicked, [=]() {
        table->setRowCount(0);
    });

    // validation before accepting
    connect(btnOk, &QPushButton::clicked, [=]() {
        if (table->rowCount() == 0) {
            QMessageBox::warning(this, "Empty Arrangement", "Please add at least one cell before confirming.");
            return;
        }
        accept();
    });
}

// set default params
void ManualArrangeDialog::setDefaultParams(double cSize, double nSize, double mZ) {
    defCellSize = cSize;
    defNucSize = nSize;
    defMajorZ = mZ;
}

// add new cell
void ManualArrangeDialog::addCell() {
    int row = table->rowCount();
    table->insertRow(row);
    table->setItem(row, 0, new QTableWidgetItem("0.0"));
    table->setItem(row, 1, new QTableWidgetItem("0.0"));
    table->setItem(row, 2, new QTableWidgetItem("0.0"));
    table->setItem(row, 3, new QTableWidgetItem(QString::number(defCellSize)));
    table->setItem(row, 4, new QTableWidgetItem(QString::number(defNucSize)));

    table->selectRow(row);
}

// remove selected cell
void ManualArrangeDialog::removeCell() {
    int currentRow = table->currentRow();
    if (currentRow >= 0) {
        table->removeRow(currentRow);
    } else if (table->rowCount() > 0) {
        table->removeRow(table->rowCount() - 1);
    }
}

// convert table data back to cell list
QList<CompleteCell> ManualArrangeDialog::getFinalCells() {
    QList<CompleteCell> cells;
    for(int i = 0; i < table->rowCount(); ++i) {
        CompleteCell cc;

        // safe extraction to prevent blank cell crash
        QString xStr = table->item(i, 0) ? table->item(i, 0)->text() : "0.0";
        QString yStr = table->item(i, 1) ? table->item(i, 1)->text() : "0.0";
        QString zStr = table->item(i, 2) ? table->item(i, 2)->text() : "0.0";
        QString cDiaStr = table->item(i, 3) ? table->item(i, 3)->text() : QString::number(defCellSize);
        QString nDiaStr = table->item(i, 4) ? table->item(i, 4)->text() : QString::number(defNucSize);

        // basic coordinates
        cc.x = xStr.toDouble();
        cc.y = yStr.toDouble();
        cc.z = zStr.toDouble();

        // radii calculations
        cc.rx = cDiaStr.toDouble() / 2.0;
        cc.rz = defMajorZ;
        cc.nrx = nDiaStr.toDouble() / 2.0;
        cc.nrz = cc.nrx / 2.0;

        // set nucleus center
        cc.nx = cc.x;
        cc.ny = cc.y;
        cc.nz = cc.z + cc.nrz + 1.0; // adjusted for custom z placement

        // assign temp ids
        cc.cellSurfId = i + 1;
        cc.nucSurfId = i + 1000;

        cells.append(cc);
    }
    return cells;
}
