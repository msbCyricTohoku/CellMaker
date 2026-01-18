#include "cellmaker.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    cellmaker w;
    w.show();
    return a.exec();
}
