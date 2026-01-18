#ifndef CELLMAKER_H
#define CELLMAKER_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui {
class cellmaker;
}
QT_END_NAMESPACE

class cellmaker : public QMainWindow
{
    Q_OBJECT
/*
    struct CellData {
        double x, y, z;
        double rx, ry, rz;
    };

    struct NucData {
        double nx,ny,nz;
        double nrx,nry,nrz;
    };
*/
    //cell data struct
    struct CompleteCell {
        //cytoplasm body
        double x, y, z;
        double rx, rz;
        double majorX, majorY;
        int cellSurfId;

        //nucleus (inside)
        double nx, ny, nz;
        double nrx, nrz;
        int nucSurfId;
    };


public:
    cellmaker(QWidget *parent = nullptr);
    ~cellmaker();

private slots:
    void on_pushButton_5_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_clicked();


    void on_pushButton_4_clicked();

    void on_pushButton_3_clicked();

    void on_actionQuit_triggered();

    void on_actionAbout_triggered();

    void on_actionGenerate_Cell_Array_triggered();

    void on_actionZoom_triggered();

    void on_actionZoom_2_triggered();

    void phitsScriptGen(const QString &path, const QString &maxcas, const QString &maxbch, const QString sourceType,
                        const QString proj, const QString r0, const QString z0, const QString e0,
                        QList<CompleteCell> cells, double bufH, double majorZ);




private:
    Ui::cellmaker *ui;
};
#endif // CELLMAKER_H
