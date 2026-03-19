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

    //panning with drag
    ui->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);

    ui->graphicsView->viewport()->installEventFilter(this);

    ui->lineEdit_6->setText("3");

    ui->lineEdit->setText("30");

    ui->lineEdit_10->setText("10");

    ui->lineEdit_11->setText("0");

    ui->lineEdit_2->setText("10");

    //now removed ui->lineEdit_5->setText("20");

    ui->checkBox->setChecked(true);
    ui->checkBox_2->setChecked(false);

    ui->checkBox_3->setChecked(false);

    ui->lineEdit_3->setText("0.9"); // cell z to xy ratio
    ui->lineEdit_4->setText("0.3"); // nucleus z to xy ratio

    ui->comboBox->addItem("Uniform");
    ui->comboBox->addItem("Random");
    ui->comboBox->addItem("3D in-vivo");

    ui->comboBox_2->addItem("proton");
    ui->comboBox_2->addItem("neutron");
    ui->comboBox_2->addItem("photon");
    ui->comboBox_2->addItem("alpha");

    ui->lineEdit_7->setText("10");

    ui->comboBox_3->addItem("Disk");
    ui->comboBox_3->addItem("Point");

    ui->lineEdit_8->setText("0.1");

    ui->lineEdit_9->setText("1");

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

    ui->comboBox_6->setCurrentIndex(1);

    ui->comboBox_7->addItem("AIR-DRY-NIST");

    ui->lineEdit_12->setText("100000");
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
            scaleFactor = 1.0 / scaleFactor; // zoom out
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
                               int cytoMatNo, int nucMatNo,
                               int buffMatNo, bool is3DMode) {


    const double micro_factor = 1E-4; // factor to scale down to micron

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

    // scaled used here
    double maxExtent = 0.0;
    double minZ = 0.0;

    for (const auto &cell : cells) {

        double reachX = std::abs(cell.x) + cell.rx;

        double reachY = std::abs(cell.y) + cell.rx; // assuming rx is used for Y plane

        maxExtent = std::max({maxExtent, reachX, reachY});

        minZ = std::min(minZ, cell.z - cell.rz);

    }

    double halfGrid = maxExtent + (20.0 * micro_factor);

    double bufferBottom = minZ - (20.0 * micro_factor);

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
    out << "Generated by CellMaker Program" << "\n\n";

    out << "[ P a r a m e t e r s ]" << Qt::endl;
    out << "icntl  =  0" << Qt::endl;
    out << "maxcas =  " << maxcas << Qt::endl;
    out << "maxbch =  " << maxbch << "\n\n";

    // parameter geomtry constants in form of c*
    out << "$--- Cell Geometry Parameters (cm) ---" << Qt::endl;
    // c10: cell radius, c11: cell height
    // c10: cell radius, c11: cell height (Now ratio-scaled)
    out << QString("set: c10[%1] $ Cell XY Radius").arg(cells[0].rx, 8, 'f', 6) << Qt::endl;
    out << QString("set: c11[%1] $ Cell Z Height").arg(cells[0].rz, 8, 'f', 6) << Qt::endl;

    // c20: nucleus radius, c21: nucleus height (Now ratio-scaled)
    out << QString("set: c20[%1] $ Nucleus XY Radius").arg(cells[0].nrx, 8, 'f', 6) << Qt::endl;
    out << QString("set: c21[%1] $ Nucleus Z Height").arg(cells[0].nrz, 8, 'f', 6) << Qt::endl;

    // c30: buffer medium height, c31: buffer half-width
    /*
    out << QString("set: c30[%1] $ Total Buffer Height").arg(majorZ + bufH, 8, 'f', 6) << Qt::endl;
    out << QString("set: c31[%1] $ Buffer Half-Width").arg(halfGrid, 8, 'f', 6) << Qt::endl;
    out << "\n";
*/

    if (is3DMode) {

        out << QString("set: c30[%1] $ Buffer Bottom Boundary").arg(bufferBottom, 8, 'f', 6) << Qt::endl;

    } else {

        double topBufferBoundary = cells[0].rz + bufH;
        out << QString("set: c30[%1] $ Total Buffer Height").arg(topBufferBoundary, 8, 'f', 6) << Qt::endl;

    }

    out << QString("set: c31[%1] $ Buffer Half Width").arg(halfGrid, 8, 'f', 6) << Qt::endl;

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

    double ref_rz = cells[0].rz;
    if (ref_rz == 0.0) ref_rz = 1.0; // fallback to prevent divide zero

    // nucleus surf section ehere
    int totalCells = cells.size();

    /*******************************************************/
    // nucleus surfaces
    for (int i = 0; i < totalCells; ++i) {
        const auto &cell = cells[i];

        QString nucZStr;

        if (is3DMode) {

            double bottomGap = (cell.z - cell.rz) - bufferBottom;

            QString gapStr = QString(bottomGap >= 0 ? "+%1" : "%1").arg(bottomGap, 0, 'f', 6);

            // proportional Z offset of nucleus relative to cell center
            double dnz = cell.nz - cell.z;

            double nucRatio = dnz / ref_rz;

            QString ratioStr = QString(nucRatio >= 0 ? "+%1*c11" : "%1*c11").arg(nucRatio, 0, 'f', 6);

            // Nucleus Z = BufferBottom + Gap + CellRadius + NucleusRelativeOffset
            nucZStr = QString("(c30%1+c11%2)").arg(gapStr).arg(ratioStr);
        } else {

            double nucRatio = cell.nz / ref_rz;

            nucZStr = QString("(%1*c11)").arg(nucRatio, 0, 'f', 6);
        }

        out << QString("%1  ELL  %2 %3 %4  0 0 c21  -c20")
                    .arg(cell.nucSurfId, -3)
                    .arg(cell.nx, 8, 'f', 6)
                    .arg(cell.ny, 8, 'f', 6)
                    .arg(nucZStr)
            << " $nucleus" << Qt::endl;
    }

    // cytoplasm surfaces
    for (int i = 0; i < totalCells; ++i) {
        const auto &cell = cells[i];

        QString cellZStr;
        if (is3DMode) {

            double bottomGap = (cell.z - cell.rz) - bufferBottom;

            QString gapStr = QString(bottomGap >= 0 ? "+%1" : "%1").arg(bottomGap, 0, 'f', 6);

            // Cell Z = BufferBottom + Gap + CellRadius
            cellZStr = QString("(c30%1+c11)").arg(gapStr);

        } else {

            // 2d mode
            cellZStr = "0.0";
        }

        // cell cytoplasm
        out << QString("%1  ELL  %2 %3 %4  0 0 c11  -c10")
                    .arg(cell.cellSurfId, -3)
                    .arg(cell.x, 8, 'f', 6)
                    .arg(cell.y, 8, 'f', 6)
                    .arg(cellZStr)
            << " $cytoplasm" << Qt::endl;
    }
    /*******************************************************/

    /*
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

        // ellipsoid of revolution (0 0 c11)
        out << QString("%1  ELL  %2 %3 %4  0 0 c11  -c10")
                   .arg(cell.cellSurfId, -3)
                   .arg(cell.x, 8, 'f', 6)
                   .arg(cell.y, 8, 'f', 6)
                   .arg(cell.z, 8, 'f', 6)
            << " $cytoplasm" << Qt::endl;
    }
    */

    // next id
    int nextId = (totalCells * 2) + 1;

    int containerId = nextId++;

    /*
    // c31 for the lateral bounds and c30 for the top
    out << QString("%1  RPP  -c31 c31 -c31 c31 -1.0e-5 c30")
               .arg(containerId)
        << " $buffer" << Qt::endl;
*/
    if (is3DMode) {

        // buffer goes from Z=0 down to the bottom boundary
        out << QString("%1  RPP  -c31 c31 -c31 c31 c30 0.0")
                   .arg(containerId)
            << " $tissue buffer" << Qt::endl;

    } else {

        // dome cells resting at Z=0, buffer goes up
        out << QString("%1  RPP  -c31 c31 -c31 c31 -1.0e-5 c30")
                   .arg(containerId)
            << " $in-vitro medium buffer" << Qt::endl;
    }

    /*

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
    */

    // bug fix here, no cutting plane needed when cells are full shape
    //int pz_bottom = 0, pz_top = 0, py_p = 0, py_m = 0, px_p = 0, px_m = 0;

    int pz_top = 0;

    if (!is3DMode) {
        //pz_bottom = nextId++; // equivalent to 11 in old code
        pz_top = nextId++;    // equivalent to 12
        //out << QString("%1  PZ  -2.0e-5").arg(pz_bottom) << Qt::endl;
        out << QString("%1  PZ  0.0").arg(pz_top) << Qt::endl;
        //py_p = nextId++;
        //out << QString("%1  PY  0.50").arg(py_p) << Qt::endl;
        //py_m = nextId++;
        //out << QString("%1  PY -0.50").arg(py_m) << Qt::endl;
        //px_p = nextId++;
        //out << QString("%1  PX  0.50").arg(px_p) << Qt::endl;
        //px_m = nextId++;
        //out << QString("%1  PX -0.50").arg(px_m) << Qt::endl;
    }

    // le void
    out << "4000   SO   500.0 $outer boundary" << Qt::endl;

    out << "\n[ C e l l ]" << Qt::endl;


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


        if (is3DMode) {

            // full ellipsoid when 3D mode, no cutting plane (-pz_top is removed)
            out << QString("%1  %2 -%3 -%4 #%5")
                       .arg(cell.cellSurfId, -3)
                       .arg(cytoMatNo + 1)
                       .arg(Cytodensity)
                       .arg(cell.cellSurfId)
                       .arg(cell.nucSurfId)
                << " $full cytoplasm" << Qt::endl;

        } else {

            // 2D dome cut by pz_top
            out << QString("%1  %2 -%3 (-%4 %5) #%6")
                       .arg(cell.cellSurfId, -3)
                       .arg(cytoMatNo + 1)
                       .arg(Cytodensity)
                       .arg(cell.cellSurfId)
                       .arg(pz_top)
                       .arg(cell.nucSurfId)
                << " $cytoplasm dome" << Qt::endl;
        }

        /*
        out << QString("%1  %2 -%3 (-%4 %5) #%6")
                   .arg(cell.cellSurfId, -3)
                   .arg(cytoMatNo + 1)
                   .arg(Cytodensity)
                   .arg(cell.cellSurfId)
                   .arg(pz_top)
                   .arg(cell.nucSurfId)
            << " $cytoplasm" << Qt::endl;
        */


        // line break in phits script
        if (cellCounter > 0 && cellCounter % 8 == 0) {
            positiveCytoSurfaces += "\n     ";
            domainsREG += "\n     ";
        }

        // appending the positive cellSurfId (outside cytoplasm)
        positiveCytoSurfaces += QString(" %1").arg(cell.cellSurfId);
        domainsREG += QString(" %1 %2").arg(cell.cellSurfId).arg(cell.nucSurfId);
        cellCounter++;
    }

    // optimized buffer cell (3000)
    if (is3DMode) {

        out << QString("3000  %1 -%2  -%3")
                   .arg(buffMatNo + 1)
                   .arg(Bufdensity)
                   .arg(containerId)
            << positiveCytoSurfaces << " $optimized buffer" << Qt::endl;

    } else {

        out << QString("3000  %1 -%2  (-%3 %4)")
                   .arg(buffMatNo + 1)
                   .arg(Bufdensity)
                   .arg(containerId)
                   .arg(pz_top)
            << positiveCytoSurfaces << " $optimized buffer" << Qt::endl;

    }

    /*
    out << QString("3000  %1 -%2  (-%3 %4)")
               .arg(buffMatNo + 1)
               .arg(Bufdensity)
               .arg(containerId)
               .arg(pz_top)
        << positiveCytoSurfaces << " $optimized buffer" << Qt::endl;


    if (!is3DMode) {

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

    }
*/

    // optimized outer world cell (4001)
    // complement operator (:) to define everything outside the RPP container

    /*
    out << QString("4001 %1 %2 -4000  ( %3 : -%4 ) #3001")
               .arg(num_materials + 1)
               .arg(Airdensity)
               .arg(containerId)
               .arg(pz_top)
        << " $optimized outer boundary" << Qt::endl;
*/


    if (is3DMode) {
        out << QString("4001 %1 %2 -4000  %3")
                   .arg(num_materials + 1)
                   .arg(Airdensity)
                   .arg(containerId)
            << " $optimized outer boundary" << Qt::endl;

    } else {

        out << QString("4001 %1 %2 -4000  ( %3 : -%4 )")
                   .arg(num_materials + 1)
                   .arg(Airdensity)
                   .arg(containerId)
                   .arg(pz_top)
            << " $optimized outer boundary" << Qt::endl;

    }

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

        if (is3DMode) {
            out << "c30 0.000" << Qt::endl;
        } else {
            out << "0.000 c30/40.0" << Qt::endl;
        }
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

        if (is3DMode) {
            out << "c30 0.000" << Qt::endl;
        } else {
            out << "0.000 c30" << Qt::endl;
        }
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

        if (is3DMode) {
            out << "c30 0.000" << Qt::endl;
        } else {
            out << "0.000 c30" << Qt::endl;
        }
        out << "axis = xy" << Qt::endl;
        out << "file = top_view_deposit.out" << Qt::endl;
        out << "part = all" << Qt::endl;
        out << "material = all" << Qt::endl;
        out << "output = dose" << Qt::endl;
        out << "unit = 2" << Qt::endl;
        out << "epsout = 1" << Qt::endl;

        out << Qt::scientific;
        out.setRealNumberPrecision(
            3); // precision

        out << "\n[ T - Deposit ]" << Qt::endl;
        out << "title = dose in cell constituents" << Qt::endl;
        out << "mesh = reg" << Qt::endl;
        out << "reg = " << domainsREG << Qt::endl;
        out << "volume" << Qt::endl;
        out << "reg      vol" << Qt::endl;

        for (const auto &cell : cells) {
            // nucleus is ALWAYS 3D ellipsoid
            double nucVol = M_PI * (4.0 / 3.0) * (cell.nrx * cell.nrx * cell.nrz);

            // cell volume depends on if it is a 3D full ellipsoid or a 2D dome
            double cellTotalVol = 0.0;
            if (is3DMode) {

                // full ellipsoid
                cellTotalVol = M_PI * (4.0 / 3.0) * (cell.rx * cell.rx * cell.rz);

            } else {

                // half ellipsoid cut at Z=0 plane
                cellTotalVol = M_PI * (2.0 / 3.0) * (cell.rx * cell.rx * cell.rz);

            }

            // to be precise, we remove the volume of the nucleus from the cytoplasm
            double cellVol = cellTotalVol - nucVol;

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


    int cellNo = ui->lineEdit_6->text().toInt();
    double cellSize = ui->lineEdit->text().toDouble();
    double nucleusSize = ui->lineEdit_10->text().toDouble();
    double CellPitch = ui->lineEdit_11->text().toDouble();

    // read ratios
    double cellZRatio = ui->lineEdit_3->text().toDouble();
    if (cellZRatio <= 0.0) cellZRatio = 0.5; // Safety fallback

    double nucZRatio = ui->lineEdit_4->text().toDouble();
    if (nucZRatio <= 0.0) nucZRatio = 0.8; // Safety fallback

    // calc. actual Z radii based on XY radii and user ratios
    double actualCellRz = (cellSize / 2.0) * cellZRatio;
    double actualNucRz = (nucleusSize / 2.0) * nucZRatio;

    double spacing = cellSize + CellPitch;
    double totalSpan = (cellNo - 1) * spacing;
    double halfSpan = totalSpan / 2.0;

    QString cytoplasmType = ui->comboBox->currentText();


    double userRadius = 0.0;
    if (cytoplasmType == "Random" || cytoplasmType == "3D in-vivo") {
        //userRadius = cellSize / 2.0;
        userRadius = ui->spinBoxRandomRadius->value();
    }

    const int maxAttempts = 1000;

    QVector<QPointF> placedCenters;
    bool placementOk = true;
    currentCellList.clear();


    // 3D in-vivo cell array definition
    if (cytoplasmType == "3D in-vivo") {
        ManualArrangeDialog dlg(this);

        dlg.setDefaultParams(cellSize, nucleusSize, cellNo, CellPitch, userRadius, cellZRatio, nucZRatio);

        if (dlg.exec() == QDialog::Accepted) {
            currentCellList = dlg.getFinalCells();

            if (currentCellList.isEmpty()) {
                ui->textBrowser->setText(tr("<font color='red'>Error: Generation aborted.</font>"));
                return;
            }

            is3DMode = true;
            renderManualCells();
        }
        return;
    }

    is3DMode = false;


    // 2D cell array data generation
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

            // calculated radii
            cc.rx = cellSize / 2.0;
            cc.rz = actualCellRz;

            cc.nrx = nucleusSize / 2.0;
            cc.nrz = actualNucRz;

            // set major axis X and Y to zero
            cc.majorX = 0.0;
            cc.majorY = 0.0;

            // valid Z space for 2D dome nucleus
            double minNz = cc.nrz + 0.000001;
            double maxNz = cc.rz * std::sqrt(1.0 - std::pow(cc.nrx / cc.rx, 2)) - cc.nrz - 0.000001;

            if (maxNz < minNz || cc.nrx >= cc.rx) {
                ui->textBrowser->setText(tr("<font color='red'>Error: Nucleus is too large or ratio causes protrusion.</font>"));
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
           "Version 1.0.0"));
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


    //double majorZ = 0.0;

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
                   z0, e0, currentCellList, bufH, cytoMatNo, nucMatNo,
                   buffMatNo, is3DMode);

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
        "<p>This utility automates the creation of complex 2D semi-ellipsoid and 3D full-ellipsoid multi cell geometries for PHITS Monte Carlo simulations. "
        "The tool ensures all components are mathematically contained within their respective boundaries, preventing geometry overlaps and errors.</p>"

        "<h3>1. Physical Geometry and Shape</h3>"
        "<ul>"
        "<li><b>Cell Shape:</b> Depending on the selected mode (under Cell Distribution dropdown menu), cells are modeled as either <b>2D semi-ellipsoids (domes)</b> resting on the Z=0 plane, or <b>3D full ellipsoids</b> submerged in tissue. "
        "The shape is controlled by a <b>'Z-to-XY Ratio'</b> (e.g., a ratio of 0.5 creates a realistic flattened cell based on the defined radius).</li>"
        "<li><b>Nucleus Shape & Placement:</b> The nucleus is modeled as a 3D <b>ellipsoid</b> using its own specific Z-to-XY Ratio. "
        "The software randomly positions it within the cytoplasm in full 3D space, mathematically ensuring it never protrudes through the cell membrane.</li>"
        "<li><b>Medium / Tissue Buffer:</b> For 2D <i>in-vitro</i> studies, the buffer extends upward from the floor. For 3D <i>in-vivo</i> modeling, "
        "it extends downwards from the tissue surface (Z=0) to completely submerge the cell array at user specified depth.</li>"
        "</ul>"

        "<h3>2. Array Distribution and Spacing</h3>"
        "<ul>"
        "<li><b>2D Array (uniform / random):</b> Generates a flat N x N grid of dome cells.</li>"
        "<li><b>3D In-Vivo Array:</b> Generates an N x N x Z multi layer grid of full ellipsoids. You can specify the exact <b>Depth In-Vivo</b> to submerge the highest cells below the surface.</li>"
        "<li><b>Uniform vs. Random:</b> 'Uniform' creates a perfect mathematical lattice. 'Random' safely shifts cells in 2D or 3D space within a defined 'Random Radius' "
        "to simulate realistic, non-idealized biological distributions. The algorithm rigorously checks to ensure random shifts never cause cell overlaps or breach the defined tissue depth.</li>"
        "</ul>"

        "<h3>3. Simulation Parameters & Scoring</h3>"
        "<ul>"
        "<li><b>Materials:</b> Assign materials to each region. The CellMaker automatically handles material IDs for the PHITS <code>[ Material ]</code> section.</li>"
        "<li><b>Source Configuration:</b> Define the particle beam as either a <b>Point</b> source or a <b>Disk</b> beam. The <code>maxcas</code> and <code>maxbch</code> parameters control statistical history generation.</li>"
        "<li><b>Tallies (Outputs):</b> Use checkboxes to automatically generate cards for geometry visualization (<code>[ T - Gshow ]</code>), particle tracking (<code>[ T - Track ]</code>), and highly precise microscopic dose deposition inside distinct cellular compartments (<code>[ T - Deposit ]</code>).</li>"
        "</ul>"

        "<h3>4. Workflow: Preview to Execution</h3>"
        "<ol>"
        "<li><b>Generate:</b> Process the mathematical model and visualize it in the 2D viewer.</li>"
        "<li><b>Verify:</b> Use the <b>XY (Top View)</b> to check lateral spacing and the <b>XZ (Side View)</b> for vertical profiles and buffer depths.</li>"
        "<li><b>Export:</b> Generate the <code>.inp</code> file with automatic <code>[ Surface ]</code> and <code>[ Cell ]</code> cards.</li>"
        "<li><b>Execute:</b> Run the local PHITS executable directly from the utility to begin the transport calculation.</li>"
        "</ol>"

        "<h3>Important Technical Notes</h3>"
        "<ul>"
        "<li><b>Geometry Validation:</b> If the user defined sizes or ratios result in a nucleus that is too large to fit safely inside the cell, or if random placement fails to resolve overlaps within 1000 attempts, generation will abort to protect the simulation physics.</li>"
        "<li><b>Units:</b> All inputs are in <b>micrometers (µm)</b>, but are automatically converted and exported to PHITS in <b>centimeters (cm)</b>.</li>"
        "</ul>"
    };

    // --- JAPANESE MANUAL ---
    ManualText japanese = {
                           "<h2 style='color: #0078D7;'>PHITS CellMaker の概要</h2>"
                           "<p>このユーティリティは、PHITSモンテカルロシミュレーション用の複雑な2次元半楕円体および3次元楕円体の多細胞ジオメトリの作成を自動化します。"
                           "このツールは、すべてのコンポーネントが数学的にそれぞれの境界内に収まることを保証し、ジオメトリの重複やエラーを防ぎます。</p>"

                           "<h3>1. 物理的ジオメトリと形状</h3>"
                           "<ul>"
                           "<li><b>細胞の形状:</b> （Cell Distribution ドロップダウンメニューで）選択したモードに応じて、細胞はZ=0平面上に配置された<b>2次元半楕円体（ドーム型）</b>、または組織内に沈められた<b>3次元完全楕円体</b>としてモデル化されます。"
                           "形状は<b>「Z-to-XY Ratio (Z/XY比)」</b>によって制御されます（例：比率を0.5にすると、指定した半径に基づいて現実的な扁平な細胞が作成されます）。</li>"
                           "<li><b>細胞核の形状と配置:</b> 細胞核は、独自のZ/XY比を用いた3次元<b>楕円体</b>としてモデル化されます。"
                           "ソフトウェアは、核を細胞質内の3次元空間にランダムに配置しますが、細胞膜から突き出ないよう数学的に保証しています。</li>"
                           "<li><b>媒質 / 組織バッファ:</b> 2次元の<i>in-vitro</i>（試験管内）研究では、バッファは底面から上方へ拡張します。3次元の<i>in-vivo</i>（生体内）モデリングでは、"
                           "組織表面（Z=0）から下方へ拡張し、ユーザーが指定した深さで細胞アレイを完全に沈めます。</li>"
                           "</ul>"

                           "<h3>2. アレイの配置と間隔</h3>"
                           "<ul>"
                           "<li><b>2D Array (均等 / ランダム):</b> ドーム型細胞の平坦な N x N グリッドを生成します。</li>"
                           "<li><b>3D In-Vivo Array:</b> 楕円体の N x N x Z 多層グリッドを生成します。最も高い位置にある細胞を表面から沈めるための正確な<b>Depth In-Vivo (生体内深度)</b>を指定できます。</li>"
                           "<li><b>均等 (Uniform) とランダム (Random):</b> 「均等」は完全な数学的格子を作成します。「ランダム」は、定義された「Random Radius (ランダム半径)」内で細胞を2次元または3次元空間で安全にシフトさせ、現実的で非理想化された生物学的分布をシミュレートします。アルゴリズムは、ランダムなシフトによって細胞の重複が発生したり、定義された組織の深さを超えたりしないよう厳密にチェックします。</li>"
                           "</ul>"

                           "<h3>3. シミュレーションパラメータとスコアリング</h3>"
                           "<ul>"
                           "<li><b>マテリアル (Materials):</b> 各領域にマテリアル（物質）を割り当てます。CellMakerは、PHITSの <code>[ Material ]</code> セクション用のマテリアルIDを自動的に処理します。</li>"
                           "<li><b>線源設定 (Source Configuration):</b> 粒子ビームを<b>点 (Point)</b>線源または<b>ディスク (Disk)</b>ビームとして定義します。<code>maxcas</code> および <code>maxbch</code> パラメータは、統計的ヒストリーの生成を制御します。</li>"
                           "<li><b>タリー (出力):</b> チェックボックスを使用すると、ジオメトリの可視化 (<code>[ T - Gshow ]</code>)、粒子のトラッキング (<code>[ T - Track ]</code>)、および個々の細胞コンパートメント内の高精度な微視的線量付与 (<code>[ T - Deposit ]</code>) のためのカードが自動的に生成されます。</li>"
                           "</ul>"

                           "<h3>4. ワークフロー: プレビューから実行まで</h3>"
                           "<ol>"
                           "<li><b>生成 (Generate):</b> 数学モデルを処理し、2Dビューアで可視化します。</li>"
                           "<li><b>確認 (Verify):</b> <b>XY (上面図)</b> で横方向の間隔を確認し、<b>XZ (側面図)</b> で垂直プロファイルとバッファの深さを確認します。</li>"
                           "<li><b>エクスポート (Export):</b> 自動生成された <code>[ Surface ]</code> および <code>[ Cell ]</code> カードを含む <code>.inp</code> ファイルを作成します。</li>"
                           "<li><b>実行 (Execute):</b> ユーティリティから直接ローカルのPHITS実行ファイルを起動し、輸送計算を開始します。</li>"
                           "</ol>"

                           "<h3>重要な技術的注意事項</h3>"
                           "<ul>"
                           "<li><b>ジオメトリの検証:</b> ユーザーが定義したサイズや比率によって、核が大きすぎて細胞内に安全に収まらない場合、またはランダム配置において1000回の試行で重複を解消できない場合、シミュレーションの物理的整合性を保護するために生成が中断されます。</li>"
                           "<li><b>単位:</b> すべての入力単位は<b>マイクロメートル (µm)</b>ですが、PHITSにエクスポートされる際に自動的に<b>センチメートル (cm)</b>に変換されます。</li>"
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

    // check view
    bool isCurrentlyXZ = false;
    if (ui->graphicsView->scene() != nullptr && ui->graphicsView->scene() == sceneXZ) {
        isCurrentlyXZ = true;
    }

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
    double minZRender = currentCellList[0].z - currentCellList[0].rz;

    for(const auto& cc : currentCellList) {
        minX = qMin(minX, cc.x - cc.rx);
        maxX = qMax(maxX, cc.x + cc.rx);
        minY = qMin(minY, cc.y - cc.rx);
        maxY = qMax(maxY, cc.y + cc.rx);
        maxZ = qMax(maxZ, cc.rz);
        minZRender = qMin(minZRender, cc.z - cc.rz);
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
    //sceneXZ->addLine(minX-40, 0, maxX+40, 0, QPen(QColor(150, 150, 150), 1.0))->setZValue(0);
    //sceneXZ->addRect(minX-20, -totalBufZ, (maxX-minX)+40, totalBufZ * 2, bufferOutline, bufferFill)->setZValue(0);
    sceneXZ->addLine(minX-40, 0, maxX+40, 0, QPen(QColor(150, 150, 150), 1.0))->setZValue(0);

    if (is3DMode) {

        double bottomZRender = minZRender - 20.0;
        sceneXZ->addRect(minX-20, 0, (maxX-minX)+40, std::abs(bottomZRender), bufferOutline, bufferFill)->setZValue(0);
    } else {

        sceneXZ->addRect(minX-20, -totalBufZ, (maxX-minX)+40, totalBufZ, bufferOutline, bufferFill)->setZValue(0);
    }

    for (const auto& cc : currentCellList) {
        if (is3DMode) {

            QLinearGradient cytoGradXZ(cc.x, -(cc.z + cc.rz), cc.x, -(cc.z - cc.rz));
            cytoGradXZ.setColorAt(0.0, QColor(0, 255, 255, 100));
            cytoGradXZ.setColorAt(1.0, QColor(0, 120, 220, 160));

            sceneXZ->addEllipse(cc.x - cc.rx, -(cc.z + cc.rz), cc.rx * 2, cc.rz * 2, cellOutline, QBrush(cytoGradXZ))->setZValue(1);

            double nucCenterZ = -(cc.nz);
            QRadialGradient nucGradXZ(cc.nx, nucCenterZ, cc.nrx);
            nucGradXZ.setFocalPoint(cc.nx - cc.nrx * 0.3, nucCenterZ - cc.nrz * 0.3);
            nucGradXZ.setColorAt(0.0, QColor(255, 120, 120));
            nucGradXZ.setColorAt(1.0, QColor(180, 0, 0));

            sceneXZ->addEllipse(cc.nx - cc.nrx, nucCenterZ - cc.nrz, cc.nrx * 2, cc.nrz * 2, nucleusOutline, QBrush(nucGradXZ))->setZValue(2);
        } else {

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
    }

    // user selected view on update
    if (isCurrentlyXZ) {
        ui->graphicsView->setScene(sceneXZ);
        ui->graphicsView->resetTransform();
        ui->graphicsView->fitInView(sceneXZ->itemsBoundingRect(), Qt::KeepAspectRatio);
        ui->graphicsView->scale(0.9, 0.9);
    } else {
        ui->graphicsView->setScene(sceneXY);
        ui->graphicsView->resetTransform();
        ui->graphicsView->fitInView(sceneXY->itemsBoundingRect(), Qt::KeepAspectRatio);
        //ui->graphicsView->scale(0.9, 0.9);
        ui->graphicsView->scale(0.9, -0.9);
    }
}



void cellmaker::on_actionSave_Model_triggered()
{
    on_pushButton_7_clicked();
}


void cellmaker::on_actionRun_PHITS_triggered()
{
    on_pushButton_2_clicked();
}


void cellmaker::on_actionManual_triggered()
{
    on_pushButton_8_clicked();
}
