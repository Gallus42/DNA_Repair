#include <QObject>
#include "mainwindow.h"

#include "DNARepairModel.hh"
#include "LogicHandler.hh"

#include <QApplication>


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow* window = new MainWindow();
    LogicHandler* logicHandler = new LogicHandler(window);
    logicHandler->setArgs(argc, argv);
    QObject::connect(window, &MainWindow::simulationRequested, logicHandler, &LogicHandler::onSimulationRequested);
    QObject::connect(window, &MainWindow::macSelectionRequested, logicHandler, &LogicHandler::onMacSelectionRequested);
    QObject::connect(window, &MainWindow::rootRequested, logicHandler, &LogicHandler::onROOTRequested);
    QObject::connect(logicHandler, &LogicHandler::rootGenerated, window, &MainWindow::onRootGenerated);
    
    QObject::connect(logicHandler, &LogicHandler::fileSelected, window, &MainWindow::onMacChanged);

    window->show();
    return a.exec();
}
