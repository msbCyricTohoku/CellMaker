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

private:
    Ui::cellmaker *ui;
};
#endif // CELLMAKER_H
