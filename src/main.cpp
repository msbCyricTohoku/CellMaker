#include "cellmaker.h"
#include <QIcon>
#include <QSystemTrayIcon>
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QIcon myIcon(":/resource/icon.png");
    a.setWindowIcon(myIcon);


    QSystemTrayIcon trayIcon(myIcon);
    trayIcon.setToolTip("cellmaker");
    trayIcon.show();

    cellmaker w;
    w.show();
    return a.exec();
}
