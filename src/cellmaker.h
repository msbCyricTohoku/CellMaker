#ifndef CELLMAKER_H
#define CELLMAKER_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QWheelEvent>
#include "celldata.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class cellmaker;
}
QT_END_NAMESPACE

class cellmaker : public QMainWindow
{
    Q_OBJECT

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
                        QList<CompleteCell> cells, double bufH, double majorZ,int cytoMatNo,int nucMatNo,int buffMatNo, bool is3DMode);




    void on_pushButton_2_clicked();

    void on_pushButton_7_clicked();

    void on_checkBox_4_stateChanged(int arg1);

    void on_pushButton_8_clicked();

    void on_actionSave_Model_triggered();

    void on_actionRun_PHITS_triggered();

    void on_actionManual_triggered();

private:
    Ui::cellmaker *ui;
    QList<CompleteCell> currentCellList;
    bool is3DMode = false;

void renderManualCells();
    // pointers to hold our two 2D views
    QGraphicsScene *sceneXY = nullptr;
    QGraphicsScene *sceneXZ = nullptr;


protected:
    // mouse wheel to allow zooming
    bool eventFilter(QObject *obj, QEvent *event) override;
};
#endif // CELLMAKER_H
