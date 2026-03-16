#include "cellmaker.h"
#include "ui_cellmaker.h"
#include <QDebug>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QGraphicsPathItem>
#include <QGroupBox>
#include <QMessageBox>
#include <QOpenGLFunctions>
#include <QOpenGLWidget>
#include <QPainter>
#include <QPainterPath>
#include <QPointF>
#include <QProcess>
#include <QRandomGenerator>
#include <QSizePolicy>
#include <QStringList>
#include <QTextStream>
#include <QVector>
#include <QWheelEvent>
#include <QtMath>
#include <QIcon>
#include <QFileDialog>
#include <QDialog>
#include <QVBoxLayout>
#include <QTextBrowser>
#include <QPushButton>
#include <QPdfWriter>
#include "manualarrangedialog.h"

cellmaker::cellmaker(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::cellmaker) {
  ui->setupUi(this);

    // Enable click-and-drag panning
    ui->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
    // Tell the view to listen to our mouse wheel event filter
    ui->graphicsView->viewport()->installEventFilter(this);

  ui->lineEdit_6->setText("3");

  ui->lineEdit->setText("30");

  ui->lineEdit_10->setText("10");

  ui->lineEdit_11->setText("0");

  ui->lineEdit_2->setText("100");

  ui->lineEdit_5->setText("20");

  ui->checkBox->setChecked(true);
  ui->checkBox_2->setChecked(false);

  ui->checkBox_3->setChecked(false);

  ui->lineEdit_3->setText("0");

  ui->lineEdit_4->setText("0");

  ui->comboBox->addItem("Uniform");
  ui->comboBox->addItem("Random");
  ui->comboBox->addItem("Manual");

  ui->comboBox_2->addItem("proton");
  ui->comboBox_2->addItem("neutron");
  ui->comboBox_2->addItem("photon");
  ui->comboBox_2->addItem("alpha");

  ui->lineEdit_7->setText("70");

  ui->comboBox_3->addItem("Disk");
  ui->comboBox_3->addItem("Point");

  ui->lineEdit_8->setText("1");

  ui->lineEdit_9->setText("100");

  ui->comboBox_4->addItem("TISSUE-SOFT(ICRU-44)");
  ui->comboBox_4->addItem("WATER");
  ui->comboBox_4->addItem("BONE-CORTICAL(ICRU-44)");
  ui->comboBox_4->addItem("ADIPOSE-TISSUE");
  ui->comboBox_4->addItem("A-150-PLASTIC");
  ui->comboBox_4->addItem("B-100-PLASTIC");
  ui->comboBox_4->addItem("MUSCLE-SKELETAL");

  ui->comboBox_5->addItem("TISSUE-SOFT(ICRU-44)");
  ui->comboBox_5->addItem("WATER");
  ui->comboBox_5->addItem("BONE-CORTICAL(ICRU-44)");
  ui->comboBox_5->addItem("ADIPOSE-TISSUE");
  ui->comboBox_5->addItem("A-150-PLASTIC");
  ui->comboBox_5->addItem("B-100-PLASTIC");
  ui->comboBox_5->addItem("MUSCLE-SKELETAL");

  ui->comboBox_6->addItem("TISSUE-SOFT(ICRU-44)");
  ui->comboBox_6->addItem("WATER");
  ui->comboBox_6->addItem("BONE-CORTICAL(ICRU-44)");
  ui->comboBox_6->addItem("ADIPOSE-TISSUE");
  ui->comboBox_6->addItem("A-150-PLASTIC");
  ui->comboBox_6->addItem("B-100-PLASTIC");
  ui->comboBox_6->addItem("MUSCLE-SKELETAL");

  ui->comboBox_7->addItem("AIR-DRY-NIST");

  ui->lineEdit_12->setText("10000");
  ui->lineEdit_13->setText("10");

    ui->checkBox_4->setChecked(true);
}

cellmaker::~cellmaker() { delete ui; }


bool cellmaker::eventFilter(QObject *obj, QEvent *event) {
    if (obj == ui->graphicsView->viewport() && event->type() == QEvent::Wheel) {
        QWheelEvent *wheelEvent = static_cast<QWheelEvent*>(event);

        // determine zoom direction
        double scaleFactor = 1.15; // 15% zoom per scroll click
        if (wheelEvent->angleDelta().y() < 0) {
            scaleFactor = 1.0 / scaleFactor; // Zoom out
        }

        ui->graphicsView->scale(scaleFactor, scaleFactor);
        return true; // handled the event
    }
    return QMainWindow::eventFilter(obj, event);
}


void cellmaker::phitsScriptGen(const QString &path, const QString &maxcas,
                               const QString &maxbch, const QString sourceType,
                               const QString proj, const QString r0,
                               const QString z0, const QString e0,
                               QList<CompleteCell> cells, double bufH,
                               double majorZ, int cytoMatNo, int nucMatNo,
                               int buffMatNo) {
    const double micro_factor = 1; // 1E-4; // factor to scale down to micron

    for (int i = 0; i < cells.size(); ++i) {

        cells[i].x *= micro_factor;
        cells[i].y *= micro_factor;
        cells[i].z *= micro_factor;
        cells[i].rx *= micro_factor;
        cells[i].rz *= micro_factor;
        cells[i].majorX *= micro_factor;
        cells[i].majorY *= micro_factor;

        cells[i].nx *= micro_factor;
        cells[i].ny *= micro_factor;
        cells[i].nz *= micro_factor;
        cells[i].nrx *= micro_factor;
        cells[i].nrz *= micro_factor;
    }

    // pre-calculate dimensions for parameters
    bufH *= micro_factor;
    majorZ *= micro_factor;

    // below fixing issue for randomized cells, the last() function cannot work
    // for randomized find the furthest point any cell reaches (center+rad)
    double maxExtent = 0.0;
    for (const auto &cell : cells) {
        double reachX = std::abs(cell.x) + cell.rx;
        double reachY = std::abs(cell.y) + cell.rx; // assuming rx is used for Y plane
        maxExtent = std::max({maxExtent, reachX, reachY});
    }
    double halfGrid = maxExtent + (20.0 * micro_factor);

    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate)) {
        qWarning() << "Failed to open" << path << ":" << f.errorString();
        return;
    }
    QTextStream out(&f);
    out.setRealNumberNotation(QTextStream::FixedNotation);
    out.setRealNumberPrecision(6);

    // header of inp
    // out << "$OMP=40" << Qt::endl;
    out << "[ T i t l e ]" << Qt::endl;
    out << "Generated by Cell Maker Program" << "\n\n";

    out << "[ P a r a m e t e r s ]" << Qt::endl;
    out << "icntl  =  0" << Qt::endl;
    out << "maxcas =  " << maxcas << Qt::endl;
    out << "maxbch =  " << maxbch << "\n\n";

    // parameter geomtry constants in form of c*
    out << "$--- Cell Geometry Parameters (cm) ---" << Qt::endl;
    // c10: cell radius, c11: cell height
    out << QString("set: c10[%1] $ Cell Radius").arg(cells[0].rx, 8, 'f', 6) << Qt::endl;
    out << QString("set: c11[%1] $ Cell Height (Major Z)").arg(majorZ, 8, 'f', 6) << Qt::endl;

    // c20: nucleus radius, c21: nucleus height
    out << QString("set: c20[%1] $ Nucleus Radius").arg(cells[0].nrx, 8, 'f', 6) << Qt::endl;
    out << QString("set: c21[%1] $ Nucleus Height").arg(cells[0].nrz, 8, 'f', 6) << Qt::endl;

    // c30: buffer medium height, c31: buffer half-width
    out << QString("set: c30[%1] $ Total Buffer Height").arg(majorZ + bufH, 8, 'f', 6) << Qt::endl;
    out << QString("set: c31[%1] $ Buffer Half-Width").arg(halfGrid, 8, 'f', 6) << Qt::endl;
    out << "\n";

    // source card here, it is dynamic for disk and point
    out << "[ S o u r c e ]" << Qt::endl;

    if (sourceType.trimmed().compare(QLatin1String("Point"),
                                     Qt::CaseInsensitive) == 0) {
        // point source (s-type = 9)
        out << "s-type = 9" << Qt::endl;
        out << "proj = " << proj << Qt::endl;
        out << "dir = all" << Qt::endl;
        out << "r1 = " << r0 << Qt::endl;
        out << "r2 = " << r0 << Qt::endl;
        out << "x0 = 0.00" << Qt::endl;
        out << "y0 = 0.00" << Qt::endl;
        out << "z0 = " << z0 << Qt::endl;
    } else {
        //s-type = 1
        out << "s-type = 1" << Qt::endl;
        out << "proj = " << proj << Qt::endl;
        out << "dir = -1.00" << Qt::endl;
        out << "r0 = " << r0 << Qt::endl;
        out << "x0 = 0.00" << Qt::endl;
        out << "y0 = 0.00" << Qt::endl;
        out << "z0 = " << z0 << Qt::endl;
        out << "z1 = " << z0 << Qt::endl;
    }
    out << "e0 = " << e0 << "\n\n";

    out << "[ M a t e r i a l ]" << Qt::endl;
    out << "mat[1] $TISSUE-SOFT(ICRU-44)" << Qt::endl;
    out << "1000 -10.5 " << Qt::endl;
    out << "6000 -25.6" << Qt::endl;
    out << "7000 -2.7 " << Qt::endl;
    out << "8000 -60.2" << Qt::endl;
    out << "11000 -0.1" << Qt::endl;
    out << "15000 -0.2" << Qt::endl;
    out << "16000 -0.3" << Qt::endl;
    out << "17000 -0.2" << Qt::endl;
    out << "19000 -0.2\n" << Qt::endl;

    out << "mat[2] $WATER" << Qt::endl;
    out << "1000 -11.1" << Qt::endl;
    out << "8000 -88.9\n" << Qt::endl;

    out << "mat[3] $BONE-CORTICAL(ICRU-44)" << Qt::endl;
    out << "1000 -3.4" << Qt::endl;
    out << "6000 -15.5" << Qt::endl;
    out << "7000 -4.2" << Qt::endl;
    out << "8000 -43.5" << Qt::endl;
    out << "11000 -0.1" << Qt::endl;
    out << "12000 -0.2" << Qt::endl;
    out << "15000 -10.3" << Qt::endl;
    out << "16000 -0.3" << Qt::endl;
    out << "20000 -22.5\n" << Qt::endl;

    out << "mat[4] $ADIPOSE-TISSUE" << Qt::endl;
    out << "1000 -11.4" << Qt::endl;
    out << "6000 -59.8" << Qt::endl;
    out << "7000 -0.7" << Qt::endl;
    out << "8000 -27.8" << Qt::endl;
    out << "11000 -0.1" << Qt::endl;
    out << "16000 -0.1" << Qt::endl;
    out << "17000 -0.1\n" << Qt::endl;

    out << "mat[5] $A-150-PLASTIC" << Qt::endl;
    out << "1000 -10.1" << Qt::endl;
    out << "6000 -77.6" << Qt::endl;
    out << "7000 -3.5" << Qt::endl;
    out << "8000 -5.2" << Qt::endl;
    out << "9000 -1.7" << Qt::endl;
    out << "20000 -1.9\n" << Qt::endl;

    out << "mat[6] $B-100-PLASTIC" << Qt::endl;
    out << "1000 -6.6" << Qt::endl;
    out << "6000 -53.7" << Qt::endl;
    out << "7000 -2.1" << Qt::endl;
    out << "8000 -3.2" << Qt::endl;
    out << "9000 -16.7" << Qt::endl;
    out << "20000 -17.7\n" << Qt::endl;

    out << "mat[7] $MUSCLE-SKELETAL" << Qt::endl;
    out << "1000 -10.2" << Qt::endl;
    out << "6000 -14.3" << Qt::endl;
    out << "7000 -3.4" << Qt::endl;
    out << "8000 -71.0" << Qt::endl;
    out << "11000 -0.1" << Qt::endl;
    out << "15000 -0.2" << Qt::endl;
    out << "16000 -0.3" << Qt::endl;
    out << "17000 -0.1" << Qt::endl;
    out << "19000 -0.4\n" << Qt::endl;

    out << "mat[8] $AIR-DRY-NIST" << Qt::endl;
    out << "6000 -0.0124" << Qt::endl;
    out << "7000 -75.5267" << Qt::endl;
    out << "8000 -23.1781" << Qt::endl;
    out << "18000 -1.2827\n\n" << Qt::endl;

    out << "[ S u r f a c e ]" << Qt::endl;

    // nucleus surf section ehere
    int totalCells = cells.size();
    for (int i = 0; i < totalCells; ++i) {
        const auto &cell = cells[i];

        out << QString("%1  ELL  %2 %3 %4  0 0 c21  -c20")
                   .arg(cell.nucSurfId, -3)
                   .arg(cell.nx, 8, 'f', 6)
                   .arg(cell.ny, 8, 'f', 6)
                   .arg(cell.nz, 8, 'f', 6)
            << " $nucleus" << Qt::endl;
    }

    // cell cyto surfaces
    for (int i = 0; i < totalCells; ++i) {
        const auto &cell = cells[i];

        out << QString("%1  ELL  %2 %3 %4  %5 %6 c11  -c10")
                   .arg(cell.cellSurfId, -3)
                   .arg(cell.x, 8, 'f', 6)
                   .arg(cell.y, 8, 'f', 6)
                   .arg(cell.z, 8, 'f', 6)
                   .arg(cell.majorX, 8, 'f', 6)
                   .arg(cell.majorY, 8, 'f', 6)
            << " $cytoplasm" << Qt::endl;
    }

    // next id
    int nextId = (totalCells * 2) + 1;

    int containerId = nextId++;
    // c31 for the lateral bounds and c30 for the top
    out << QString("%1  RPP  -c31 c31 -c31 c31 -1.0e-5 c30")
               .arg(containerId)
        << " $buffer" << Qt::endl;

    // the planes
    int pz_bottom = nextId++; // equivalent to 11 in old code
    int pz_top = nextId++;    // equivalent to 12
    out << QString("%1  PZ  -2.0e-5").arg(pz_bottom) << Qt::endl;
    out << QString("%1  PZ  0.0").arg(pz_top) << Qt::endl;

    int py_p = nextId++;
    out << QString("%1  PY  0.50").arg(py_p) << Qt::endl;
    int py_m = nextId++;
    out << QString("%1  PY -0.50").arg(py_m) << Qt::endl;
    int px_p = nextId++;
    out << QString("%1  PX  0.50").arg(px_p) << Qt::endl;
    int px_m = nextId++;
    out << QString("%1  PX -0.50").arg(px_m) << Qt::endl;

    // le void
    out << "4000   SO   500.0 $outer boundary" << Qt::endl;

    out << "\n[ C e l l ]" << Qt::endl;

    // instead of a string of negative # operators, string of
    // positive surface IDs representing the outsides of the cytoplasms
    QString positiveCytoSurfaces = "";
    QString domainsREG = "";
    int cellCounter = 0;

    double Cytodensity = 0.0;
    double Nucdensity = 0.0;
    double Bufdensity = 0.0;
    double Airdensity = -0.00129;
    static const double densities[] = {1.06, 1.00, 1.92, 0.95, 1.127, 1.45, 1.05};
    const int num_materials = sizeof(densities) / sizeof(densities[0]);

    // qDebug() << num_materials << Qt::endl;

    if (cytoMatNo >= 0 && cytoMatNo < num_materials) {
        Cytodensity = densities[cytoMatNo];
    }

    if (nucMatNo >= 0 && nucMatNo < num_materials) {
        Nucdensity = densities[nucMatNo];
    }

    if (buffMatNo >= 0 && buffMatNo < num_materials) {
        Bufdensity = densities[buffMatNo];
    }

    // nuc and cyto pairs
    for (const auto &cell : cells) {

        out << QString("%1  %2 -%3  -%4")
                   .arg(cell.nucSurfId, -3)
                   .arg(nucMatNo + 1)
                   .arg(Nucdensity)
                   .arg(cell.nucSurfId)
            << " $nucleus" << Qt::endl;

        out << QString("%1  %2 -%3 (-%4 %5) #%6")
                   .arg(cell.cellSurfId, -3)
                   .arg(cytoMatNo + 1)
                   .arg(Cytodensity)
                   .arg(cell.cellSurfId)
                   .arg(pz_top)
                   .arg(cell.nucSurfId)
            << " $cytoplasm" << Qt::endl;

        // line break with 5 space to tackle phits limit
        if (cellCounter > 0 && cellCounter % 8 == 0) {
            positiveCytoSurfaces += "\n     ";
            domainsREG += "\n     ";
        }

        //appending the positive cellSurfId (outside cytoplasm)
        positiveCytoSurfaces += QString(" %1").arg(cell.cellSurfId);
        domainsREG += QString(" %1 %2").arg(cell.cellSurfId).arg(cell.nucSurfId);
        cellCounter++;
    }

    // optimized buffer cell (3000)
    out << QString("3000  %1 -%2  (-%3 %4)")
               .arg(buffMatNo + 1)
               .arg(Bufdensity)
               .arg(containerId)
               .arg(pz_top)
        << positiveCytoSurfaces << " $optimized buffer" << Qt::endl;

    // le buffer
    out << QString("3001  %1 -%2 %3 -%4 -%5 %6 -%7 %8")
               .arg(buffMatNo + 1)
               .arg(Bufdensity)
               .arg(pz_bottom)
               .arg(pz_top)
               .arg(py_p)
               .arg(py_m)
               .arg(px_p)
               .arg(px_m)
        << Qt::endl;

    // optimized outer world cell (4001)
    // complement operator (:) to define everything outside the RPP container
    out << QString("4001 %1 %2 -4000  ( %3 : -%4 ) #3001")
               .arg(num_materials + 1)
               .arg(Airdensity)
               .arg(containerId)
               .arg(pz_top)
        << " $optimized outer boundary" << Qt::endl;

    // the void
    out << "4000  -1  4000" << Qt::endl;

    if (ui->checkBox_3->isChecked()) {

        out << "\n[ T - Gshow ]" << Qt::endl;
        out << "title = generated cell array" << Qt::endl;
        out << "mesh =  xyz" << Qt::endl;
        out << "x-type = 2" << Qt::endl;
        out << "xmin = -c31" << Qt::endl;
        out << "xmax = c31" << Qt::endl;
        out << "nx = 200" << Qt::endl;
        out << "y-type = 2" << Qt::endl;
        out << "ymin = -c31" << Qt::endl;
        out << "ymax = c31" << Qt::endl;
        out << "ny = 200" << Qt::endl;
        out << " z-type = 1" << Qt::endl;
        out << "nz = 1" << Qt::endl;
        out << "0.000 c30/40.0" << Qt::endl;
        out << "axis = xy" << Qt::endl;
        out << "file = geometry_topview.out" << Qt::endl;
        out << "output = 6" << Qt::endl;
        out << "epsout = 1" << Qt::endl;
    }

    if (ui->checkBox_2->isChecked()) {
        out << "\n[ T - Track ]" << Qt::endl;
        out << "title = track on cell array" << Qt::endl;
        out << "mesh =  xyz" << Qt::endl;
        out << "x-type = 2" << Qt::endl;
        out << "xmin = -c31" << Qt::endl;
        out << "xmax = c31" << Qt::endl;
        out << "nx = 200" << Qt::endl;
        out << "y-type = 2" << Qt::endl;
        out << "ymin = -c31" << Qt::endl;
        out << "ymax = c31" << Qt::endl;
        out << "ny = 200" << Qt::endl;
        out << " z-type = 1" << Qt::endl;
        out << "nz = 1" << Qt::endl;
        out << "0.000 c30" << Qt::endl;
        out << "axis = xy" << Qt::endl;
        out << "file = tracks.out" << Qt::endl;
        out << "e-type = 3" << Qt::endl;
        out << "emin = 1.0E-3" << Qt::endl;
        out << "emax = " << e0 << Qt::endl;
        out << "ne = 1" << Qt::endl;
        out << "part = all" << Qt::endl;
        out << "material = all" << Qt::endl;
        out << "epsout = 1" << Qt::endl;
    }

    if (ui->checkBox->isChecked()) {

        out << "\n[ T - Deposit ]" << Qt::endl;
        out << "title = dose on cell array" << Qt::endl;
        out << "mesh =  xyz" << Qt::endl;
        out << "x-type = 2" << Qt::endl;
        out << "xmin = -c31" << Qt::endl;
        out << "xmax = c31" << Qt::endl;
        out << "nx = 200" << Qt::endl;
        out << "y-type = 2" << Qt::endl;
        out << "ymin = -c31" << Qt::endl;
        out << "ymax = c31" << Qt::endl;
        out << "ny = 200" << Qt::endl;
        out << " z-type = 1" << Qt::endl;
        out << "nz = 1" << Qt::endl;
        out << "0.000 c30" << Qt::endl;
        out << "axis = xy" << Qt::endl;
        out << "file = top_view_deposit.out" << Qt::endl;
        out << "part = all" << Qt::endl;
        out << "material = all" << Qt::endl;
        out << "output = dose" << Qt::endl;
        out << "unit = 2" << Qt::endl;
        out << "epsout = 1" << Qt::endl;

        out << Qt::scientific;
        out.setRealNumberPrecision(
            3); // good enough for cells in micron range since will be written as cm

        out << "\n[ T - Deposit ]" << Qt::endl;
        out << "title = dose in cell constituents" << Qt::endl;
        out << "mesh = reg" << Qt::endl;
        out << "reg = " << domainsREG << Qt::endl;
        out << "volume" << Qt::endl;
        out << "reg      vol" << Qt::endl;

        for (const auto &cell : cells) {
            double nucVol = M_PI * (4.0 / 3.0) * (cell.nrx * cell.nrx * cell.nrz);
            double cellVol =
                M_PI * (2.0 / 3.0) * (cell.rx * cell.rx * cell.rz) -
                nucVol; // to be precise remove out volume of nucleus from cytoplasm

            out << cell.cellSurfId << "   " << cellVol << Qt::endl;
            out << cell.nucSurfId << "   " << nucVol << Qt::endl;
        }

        out << "axis = reg" << Qt::endl;
        out << "file = dose_cells.out" << Qt::endl;
        out << "unit = 0" << Qt::endl;
        out << "part = photon neutron proton alpha" << Qt::endl;
        out << "material = all" << Qt::endl;
        out << "output = dose" << Qt::endl;
        out << "epsout = 0" << Qt::endl;
    }

    out << "\n[ E n d ]";
}


void cellmaker::on_pushButton_clicked() {
    ui->textBrowser->setText("");
    QRandomGenerator *generator = QRandomGenerator::global();

    // Read UI inputs
    int cellNo = ui->lineEdit_6->text().toInt();
    double cellSize = ui->lineEdit->text().toDouble();
    double nucleusSize = ui->lineEdit_10->text().toDouble();
    double CellPitch = ui->lineEdit_11->text().toDouble();

    const double majorZ = ui->lineEdit_5->text().toDouble();
    const double majorX = ui->lineEdit_3->text().toDouble();
    const double majorY = ui->lineEdit_4->text().toDouble();
    const double bufH = ui->lineEdit_2->text().toDouble();

    double spacing = cellSize + CellPitch;
    double totalSpan = (cellNo - 1) * spacing;
    double halfSpan = totalSpan / 2.0;

    QString cytoplasmType = ui->comboBox->currentText();
    double userRadius = ui->spinBoxRandomRadius->value();
    const int maxAttempts = 1000;

    QVector<QPointF> placedCenters;
    bool placementOk = true;
    currentCellList.clear();


    //maual cell array deifition
    if (cytoplasmType == "Manual") {
        ManualArrangeDialog dlg(this);
        dlg.setDefaultParams(ui->lineEdit->text().toDouble(),
                             ui->lineEdit_10->text().toDouble(),
                             ui->lineEdit_5->text().toDouble());

        if (dlg.exec() == QDialog::Accepted) {
            currentCellList = dlg.getFinalCells();
            renderManualCells(); // Drawing happens here
        }
        return;
    }


    //cell array data generation
    for (int i = 0; i < cellNo && placementOk; ++i) {
        for (int j = 0; j < cellNo && placementOk; ++j) {
            double baseCx = (i * spacing) - halfSpan;
            double baseCy = (j * spacing) - halfSpan;
            double cz = 0.0;
            double cx = baseCx;
            double cy = baseCy;

            if (cytoplasmType == "Random") {
                int attempts = 0;
                bool placed = false;
                while (attempts < maxAttempts && !placed) {
                    double thetaPos = 2.0 * M_PI * generator->generateDouble();
                    double rPos = userRadius * std::sqrt(generator->generateDouble());

                    double candCx = baseCx + rPos * std::cos(thetaPos);
                    double candCy = baseCy + rPos * std::sin(thetaPos);

                    bool overlap = false;
                    for (const QPointF &p : placedCenters) {
                        if (std::hypot(candCx - p.x(), candCy - p.y()) < cellSize) {
                            overlap = true;
                            break;
                        }
                    }

                    if (!overlap) {
                        cx = candCx;
                        cy = candCy;
                        placed = true;
                        placedCenters.append(QPointF(cx, cy));
                    }
                    ++attempts;
                }

                if (!placed) {
                    ui->textBrowser->setText(tr("<font color='red'>Error: Unable to place cells without overlap.</font>"));
                    placementOk = false;
                    break;
                }
            } else {
                placedCenters.append(QPointF(cx, cy));
            }

            CompleteCell cc;
            cc.x = cx; cc.y = cy; cc.z = cz;
            cc.rx = cellSize / 2.0;
            cc.rz = majorZ;
            cc.majorX = majorX; cc.majorY = majorY;
            cc.nrx = nucleusSize / 2.0;
            cc.nrz = nucleusSize / 4.0;

            double minNz = cc.nrz + 0.000001;
            double maxNz = cc.rz * std::sqrt(1.0 - std::pow(cc.nrx / cc.rx, 2)) - cc.nrz - 0.000001;

            if (maxNz < minNz || cc.nrx >= cc.rx) {
                ui->textBrowser->setText(tr("<font color='red'>Error: Nucleus is too large.</font>"));
                placementOk = false;
                continue;
            }

            cc.nz = minNz + (maxNz - minNz) * generator->generateDouble();
            double safeCellRx = cc.rx * std::sqrt(1.0 - std::pow((cc.nz + cc.nrz) / cc.rz, 2));

            double maxLateralOffset = safeCellRx - cc.nrx - 0.000001;
            if (maxLateralOffset < 0) maxLateralOffset = 0;

            double theta = 2.0 * M_PI * generator->generateDouble();
            double r = maxLateralOffset * std::sqrt(generator->generateDouble());
            cc.nx = cc.x + r * std::cos(theta);
            cc.ny = cc.y + r * std::sin(theta);

            cc.cellSurfId = currentCellList.size() + 1;
            cc.nucSurfId = cc.cellSurfId + (cellNo * cellNo);
            currentCellList.append(cc);
        }
    }

    if (!placementOk) return;

    renderManualCells();


    ui->textBrowser->append("\n<font color='green'>2D Models generated successfully! Scroll to zoom.</font>");

}

//view toggles
void cellmaker::on_pushButton_5_clicked() {
    // XY (top)
    if (sceneXY) {
        ui->graphicsView->setScene(sceneXY);
        ui->graphicsView->resetTransform();
        ui->graphicsView->fitInView(sceneXY->itemsBoundingRect(), Qt::KeepAspectRatio);
    }
}

void cellmaker::on_pushButton_6_clicked() {
    // XZ (side)
    if (sceneXZ) {
        ui->graphicsView->setScene(sceneXZ);
        ui->graphicsView->resetTransform();
        ui->graphicsView->fitInView(sceneXZ->itemsBoundingRect(), Qt::KeepAspectRatio);
    }
}



void cellmaker::on_pushButton_4_clicked() { QApplication::quit(); }

void cellmaker::on_pushButton_3_clicked() {
  QMessageBox::about(
      this, tr("CellMaker Program"),
      tr("Generate Cell Arrays for Monte Carlo Radiation Transport Studies.\n"
         "Developed by: Mehrdad S. Beni and Hiroshi Watabe, RARiS, Tohoku "
         "University, JAPAN -- 2026\n"
         "Version 1.0"));
}

void cellmaker::on_actionQuit_triggered() { on_pushButton_4_clicked(); }

void cellmaker::on_actionAbout_triggered() { on_pushButton_3_clicked(); }

void cellmaker::on_actionGenerate_Cell_Array_triggered() {
  on_pushButton_clicked();
}

void cellmaker::on_actionZoom_triggered() { on_pushButton_5_clicked(); }

void cellmaker::on_actionZoom_2_triggered() { on_pushButton_6_clicked(); }


void cellmaker::on_pushButton_2_clicked() {

    QString defaultPath = qApp->applicationDirPath();

    QString runDir = QFileDialog::getExistingDirectory(
        this,
        tr("Select Directory Containing cell.inp (Run Directory)"),
        defaultPath,
        QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
        );

    if (runDir.isEmpty()) {
        return;
    }

    QProcess *process = new QProcess(this);


    process->setWorkingDirectory(runDir);

    connect(process, &QProcess::readyRead,
            [=]() { ui->textBrowser->append(process->readAll()); });

    connect(process, qOverload<int, QProcess::ExitStatus>(&QProcess::finished),
            [=]() {
                ui->textBrowser->append("\n--- RUN DONE ---");
                process->deleteLater();
            });

    process->start("phits.sh", {"cell.inp"});
}

void cellmaker::on_pushButton_7_clicked() {

    if (currentCellList.isEmpty()) {
        ui->textBrowser->setText("<font color='red'>Error: Please generate a model first!</font>");
        return;
    }


    QString maxcas = ui->lineEdit_12->text();
    QString maxbch = ui->lineEdit_13->text();
    QString sourceType = ui->comboBox_3->currentText();
    QString proj = ui->comboBox_2->currentText();
    QString r0 = ui->lineEdit_8->text();
    QString z0 = ui->lineEdit_9->text();
    QString e0 = ui->lineEdit_7->text();

    QString bufHString = ui->lineEdit_2->text();
    double bufH = bufHString.toDouble();
    QString majorZString = ui->lineEdit_5->text();
    double majorZ = majorZString.toDouble();

    int cytoMatNo = ui->comboBox_4->currentIndex();
    int nucMatNo = ui->comboBox_5->currentIndex();
    int buffMatNo = ui->comboBox_6->currentIndex();


    QString defaultPath = qApp->applicationDirPath() + "/program_output/cell.inp";
    QString saveFilePath = QFileDialog::getSaveFileName(
        this,
        tr("Save PHITS Input File"),
        defaultPath,
        tr("PHITS Input Files (*.inp);;All Files (*)")
        );

    if (saveFilePath.isEmpty()) {
        return; // cancelled
    }

    phitsScriptGen(saveFilePath, maxcas, maxbch, sourceType, proj, r0,
                   z0, e0, currentCellList, bufH, majorZ, cytoMatNo, nucMatNo,
                   buffMatNo);

    ui->textBrowser->append("\n<font color='blue'>File saved to: " + saveFilePath + "</font>");
}

void cellmaker::on_checkBox_4_stateChanged(int arg1)
{
    bool showBuffer = (arg1 != 0);

    if (sceneXY) {
        for (QGraphicsItem *item : sceneXY->items()) {
            if (item->zValue() == 0) {
                item->setVisible(showBuffer);
            }
        }
    }

    if (sceneXZ) {
        for (QGraphicsItem *item : sceneXZ->items()) {
            if (item->zValue() == 0) {
                item->setVisible(showBuffer);
            }
        }
    }
}


void cellmaker::on_pushButton_8_clicked() {
    QDialog *tutorialDialog = new QDialog(this);
    tutorialDialog->setWindowTitle(tr("PHITS CellMaker User Manual"));
    tutorialDialog->resize(800, 700);

    QVBoxLayout *layout = new QVBoxLayout(tutorialDialog);
    QHBoxLayout *topButtons = new QHBoxLayout();
    QTextBrowser *textBrowser = new QTextBrowser(tutorialDialog);

    struct ManualText {
        QString content;
    };

    // --- ENGLISH MANUAL ---
    ManualText english = {
        "<h2 style='color: #0078D7;'>PHITS CellMaker Overview</h2>"
        "<p>This utility automates the creation of complex 3D semi-ellipsoid cell geometries for PHITS. "
        "The tool ensures all components are mathematically contained within their respective boundaries to prevent geometry errors.</p>"

        "<h3>1. Physical Geometry and Shape</h3>"
        "<ul>"
        "<li><b>Cell Shape:</b> Each cell is modeled as a <b>semi-ellipsoid (dome)</b> sitting on the Z=0 plane. "
        "The <b>'Major Z'</b> specifically defines the <b>peak vertical height</b> of the cell dome.</li>"
        "<li><b>Nucleus Shape & Placement:</b> The nucleus is a 3D <b>ellipsoid</b>. By default, it is slightly flattened (oblate) "
        "to mimic biological reality. The software randomly positions it within the cytoplasm while ensuring it never "
        "protrudes through the cell membrane (dome surface).</li>"
        "<li><b>Medium / Buffer (Extracellular Region):</b> The 'Buffer' represents the <b>Culture Medium</b> (for in-vitro studies) "
        "or the <b>Interstitium/Extracellular fluid</b> (for in-vivo modeling). It is the volume of liquid that "
        "submerges the cell array, extending from the floor up to the specified buffer height.</li>"
        "</ul>"

        "<h3>2. Array Distribution and Spacing</h3>"
        "<ul>"
        "<li><b>Cell Array Size:</b> Creates an N x N grid of cells.</li>"
        "<li><b>Cell Pitch:</b> Specifies the clear distance between the outer boundaries of adjacent cells.</li>"
        "<li><b>Uniform vs. Random:</b> 'Uniform' creates a perfect lattice. 'Random' allows cells to shift within a defined radius "
        "to simulate realistic, non-idealized biological samples.</li>"
        "<li><b>Manual Arrangement:</b> Allows you to bypass the grid and input exact (X, Y, Z) coordinates for individual cells using a custom data table.</li>"
        "</ul>"

        "<h3>3. Simulation Parameters & Scoring</h3>"
        "<ul>"
        "<li><b>Materials:</b> Assign materials to each region. The software handles material IDs for the PHITS <code>[ Material ]</code> section.</li>"
        "<li><b>Source Configuration:</b> Define the particle beam as either a <b>Point</b> source or a broad parallel beam. <code>maxcas</code> and <code>maxbch</code> control the statistical history generation.</li>"
        "<li><b>Tallies (Outputs):</b> Use the checkboxes to automatically generate cards for geometry visualization (<code>[ T - Gshow ]</code>), particle tracking (<code>[ T - Track ]</code>), and highly precise microscopic dose deposition inside the distinct cellular compartments (<code>[ T - Deposit ]</code>).</li>"
        "</ul>"

        "<h3>4. Workflow: Preview to Execution</h3>"
        "<ol>"
        "<li><b>Generate:</b> Process the model and visualize it in the 2D viewer.</li>"
        "<li><b>Verify:</b> Use <b>XY (Top View)</b> for lateral spacing and <b>XZ (Side View)</b> for vertical profile and nucleus depth.</li>"
        "<li><b>Export:</b> Generate the <code>.inp</code> file with automatic <code>[ Surface ]</code> and <code>[ Cell ]</code> cards. <i>Note: The exported file is completely parameterized using <code>c*</code> variables, meaning you can easily tweak cell radii and heights directly in the text file later!</i></li>"
        "<li><b>Execute:</b> Run the local PHITS executable directly from the tool.</li>"
        "</ol>"

        "<h3>Important Technical Notes</h3>"
        "<ul>"
        "<li><b>Geometry Validation:</b> If the 'Nucleus Size' exceeds the 'Major Z' or 'Cell Size', generation will fail to prevent overlaps.</li>"
        "<li><b>Units:</b> Inputs are in <b>micrometers (µm)</b>, but are automatically exported to PHITS in <b>centimeters (cm)</b>.</li>"
        "</ul>"
    };

    // --- JAPANESE MANUAL ---
    ManualText japanese = {
        "<h2 style='color: #0078D7;'>PHITS CellMaker 概要</h2>"
        "<p>このユーティリティは、モンテカルロ放射線輸送シミュレーション用の複雑な3D半楕円体細胞形状の作成を自動化します。"
        "PHITSのジオメトリエラーを防止するため、すべての構成要素がそれぞれの境界内に数学的に収まるように設計されています。</p>"

        "<h3>1. 物理的なジオメトリと形状</h3>"
        "<ul>"
        "<li><b>細胞の形状:</b> 各細胞は、Z=0平面上に配置された<b>半楕円体（ドーム型）</b>としてモデル化されます。"
        "<b>「Major Z」</b>は、この細胞ドームの底面から頂点までの<b>垂直方向の高さ</b>を定義します。</li>"
        "<li><b>核の形状と配置:</b> 核は3次元の<b>楕円体</b>としてモデル化されています。デフォルトでは、実際の生物学的細胞に"
        "近い、やや平らな形状をしています。核が細胞壁（ドーム表面）を突き抜けないよう、安全な範囲内でランダムに配置されます。</li>"
        "<li><b>培養液・細胞外領域 (Medium / Buffer):</b> 「バッファー」は、細胞を浸している液体環境を表します。"
        "試験管内（in-vitro）実験における<b>培養液（Medium）</b>、または生体内（in-vivo）モデルにおける<b>間質液や細胞外領域</b>に相当します。"
        "これは床面から指定された高さまで細胞配列全体を覆う液体の体積を定義します。</li>"
        "</ul>"

        "<h3>2. 配列の分布と間隔</h3>"
        "<ul>"
        "<li><b>細胞配列サイズ:</b> N x N の細胞グリッドを作成します。</li>"
        "<li><b>セルピッチ:</b> 隣接する細胞の外境界間の距離を指定します。</li>"
        "<li><b>均一 vs ランダム:</b> 「Uniform」は一様な格子形状を作成します。「Random」は位置をランダムにシフトさせ、実際の生体サンプルをシミュレートします。</li>"
        "<li><b>手動配置 (Manual):</b> グリッドを使用せず、専用のデータテーブルから個々の細胞の正確な (X, Y, Z) 座標を直接入力して配置できます。</li>"
        "</ul>"

        "<h3>3. シミュレーションパラメータとスコアリング</h3>"
        "<ul>"
        "<li><b>材料:</b> 各領域に材料を割り当てる必要があります。材料識別番号（Material ID）は自動で処理されます。</li>"
        "<li><b>線源設定:</b> 粒子ビームを<b>点線源 (Point)</b> または平行ビームとして定義します。統計精度は <code>maxcas</code> と <code>maxbch</code> で制御します。</li>"
        "<li><b>タリー (出力):</b> チェックボックスを使用して、形状確認 (<code>[ T - Gshow ]</code>)、飛跡の可視化 (<code>[ T - Track ]</code>)、および各細胞コンパートメント内の高精度な微視的付与線量計算 (<code>[ T - Deposit ]</code>) のカードを自動生成します。</li>"
        "</ul>"

        "<h3>4. ワークフロー: プレビューから実行まで</h3>"
        "<ol>"
        "<li><b>生成 (Generate):</b> 数学モデルを処理し、2Dビューアで視覚化します。</li>"
        "<li><b>検証 (Verify):</b> <b>XY (平面図)</b> で横方向の間隔を、<b>XZ (側面図)</b> で垂直プロファイルと核の深さを確認します。</li>"
        "<li><b>保存 (Export):</b> <code>.inp</code> ファイルを生成し、<code>[ Surface ]</code> および <code>[ Cell ]</code> カードを自動作成します。<i>注：エクスポートされたファイルは <code>c*</code> 変数を使用して完全にパラメータ化されているため、後でテキストファイル上で細胞の半径や高さを簡単に微調整できます！</i></li>"
        "<li><b>実行 (Execute):</b> PHITS実行ファイルを呼び出し、輸送計算を開始します。</li>"
        "</ol>"

        "<h3>重要な技術的注意点</h3>"
        "<ul>"
        "<li><b>ジオメトリの検証:</b> 核のサイズが「Major Z」や「細胞サイズ」を超える場合、形状の不整合を防ぐため生成に失敗します。</li>"
        "<li><b>単位:</b> 入力は<b>マイクロメートル (µm)</b> ですが、PHITS出力は<b>センチメートル (cm)</b> に自動変換されます。</li>"
        "</ul>"
    };

    textBrowser->setHtml(english.content);

    //language toggle
    QPushButton *btnLang = new QPushButton("日本語に切替 (Switch to Japanese)", tutorialDialog);
    btnLang->setMinimumHeight(30);
    connect(btnLang, &QPushButton::clicked, [=]() mutable {
        static bool isEng = true;
        isEng = !isEng;
        textBrowser->setHtml(isEng ? english.content : japanese.content);
        btnLang->setText(isEng ? "日本語に切替 (Switch to Japanese)" : "Englishに切替 (Switch to English)");
    });

    //print to pdf
    QPushButton *btnPrint = new QPushButton("Save as PDF", tutorialDialog);
    btnPrint->setMinimumHeight(30);
    connect(btnPrint, &QPushButton::clicked, [=]() {
        QString fileName = QFileDialog::getSaveFileName(this, "Export Manual to PDF", "PHITS_CellMaker_Manual.pdf", "*.pdf");
        if (!fileName.isEmpty()) {
            QPdfWriter writer(fileName);
            writer.setPageSize(QPageSize(QPageSize::A4));
            writer.setPageMargins(QMarginsF(15, 15, 15, 15));
            textBrowser->document()->print(&writer);
            QMessageBox::information(this, "Export Success", "The user manual has been saved as a PDF.");
        }
    });


    topButtons->addWidget(btnLang);

    topButtons->addWidget(btnPrint);


    layout->addLayout(topButtons);

    layout->addWidget(textBrowser);

    QPushButton *closeButton = new QPushButton(tr("Close Manual"), tutorialDialog);
    closeButton->setMinimumHeight(40);
    connect(closeButton, &QPushButton::clicked, tutorialDialog, &QDialog::accept);
    layout->addWidget(closeButton);

    tutorialDialog->exec();
    tutorialDialog->deleteLater();
}


void cellmaker::renderManualCells() {
    if (currentCellList.isEmpty()) return;


    if (sceneXY) { delete sceneXY; sceneXY = nullptr; }
    if (sceneXZ) { delete sceneXZ; sceneXZ = nullptr; }

    sceneXY = new QGraphicsScene(this);
    sceneXZ = new QGraphicsScene(this);


    sceneXY->setBackgroundBrush(QColor(20, 25, 30));
    sceneXZ->setBackgroundBrush(QColor(20, 25, 30));

    QPen cellOutline(QColor(0, 150, 255, 150));
    cellOutline.setWidthF(0.5);
    QPen nucleusOutline(Qt::NoPen);
    QPen bufferOutline(QColor(255, 200, 0, 150));
    bufferOutline.setStyle(Qt::DashLine);
    QBrush bufferFill(QColor(255, 200, 0, 15));


    double minX = currentCellList[0].x - currentCellList[0].rx;
    double maxX = currentCellList[0].x + currentCellList[0].rx;
    double minY = currentCellList[0].y - currentCellList[0].rx;
    double maxY = currentCellList[0].y + currentCellList[0].rx;
    double maxZ = currentCellList[0].rz;

    for(const auto& cc : currentCellList) {
        minX = qMin(minX, cc.x - cc.rx);
        maxX = qMax(maxX, cc.x + cc.rx);
        minY = qMin(minY, cc.y - cc.rx);
        maxY = qMax(maxY, cc.y + cc.rx);
        maxZ = qMax(maxZ, cc.rz);
    }

    double bufH = ui->lineEdit_2->text().toDouble();
    double totalBufZ = maxZ + bufH;

    //XY view
    sceneXY->addRect(minX-20, minY-20, (maxX-minX)+40, (maxY-minY)+40, bufferOutline, bufferFill)->setZValue(0);

    for (const auto& cc : currentCellList) {
        QRadialGradient cytoGradXY(cc.x, cc.y, cc.rx);
        cytoGradXY.setColorAt(0.0, QColor(0, 255, 255, 30));
        cytoGradXY.setColorAt(1.0, QColor(0, 120, 220, 120));
        sceneXY->addEllipse(cc.x - cc.rx, cc.y - cc.rx, cc.rx * 2, cc.rx * 2, cellOutline, QBrush(cytoGradXY))->setZValue(1);

        QRadialGradient nucGradXY(cc.nx, cc.ny, cc.nrx);
        nucGradXY.setFocalPoint(cc.nx - cc.nrx * 0.3, cc.ny - cc.nrx * 0.3);
        nucGradXY.setColorAt(0.0, QColor(255, 120, 120));
        nucGradXY.setColorAt(1.0, QColor(180, 0, 0));
        sceneXY->addEllipse(cc.nx - cc.nrx, cc.ny - cc.nrx, cc.nrx * 2, cc.nrx * 2, nucleusOutline, QBrush(nucGradXY))->setZValue(2);
    }

    //XZ view
    sceneXZ->addLine(minX-40, 0, maxX+40, 0, QPen(QColor(150, 150, 150), 1.0))->setZValue(0);
    sceneXZ->addRect(minX-20, -totalBufZ, (maxX-minX)+40, totalBufZ, bufferOutline, bufferFill)->setZValue(0);

    for (const auto& cc : currentCellList) {
        QLinearGradient cytoGradXZ(cc.x, -cc.rz, cc.x, 0);
        cytoGradXZ.setColorAt(0.0, QColor(0, 255, 255, 100));
        cytoGradXZ.setColorAt(1.0, QColor(0, 120, 220, 160));

        QPainterPath dome;
        dome.moveTo(cc.x + cc.rx, 0);
        dome.arcTo(cc.x - cc.rx, -cc.rz, cc.rx * 2, cc.rz * 2, 0, 180);
        dome.lineTo(cc.x - cc.rx, 0);
        sceneXZ->addPath(dome, cellOutline, QBrush(cytoGradXZ))->setZValue(1);

        double nucCenterZ = -(cc.nz + cc.nrz);
        QRadialGradient nucGradXZ(cc.nx, nucCenterZ, cc.nrx);
        nucGradXZ.setFocalPoint(cc.nx - cc.nrx * 0.3, nucCenterZ - cc.nrz * 0.3);
        nucGradXZ.setColorAt(0.0, QColor(255, 120, 120));
        nucGradXZ.setColorAt(1.0, QColor(180, 0, 0));
        sceneXZ->addEllipse(cc.nx - cc.nrx, nucCenterZ - cc.nrz, cc.nrx * 2, cc.nrz * 2, nucleusOutline, QBrush(nucGradXZ))->setZValue(2);
    }

    ui->graphicsView->setScene(sceneXY);
    ui->graphicsView->resetTransform();
    ui->graphicsView->fitInView(sceneXY->itemsBoundingRect(), Qt::KeepAspectRatio);
    ui->graphicsView->scale(0.9, 0.9);
}
