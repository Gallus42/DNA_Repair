#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QDialog>
#include <QLabel>
#include <QMainWindow>
#include <QPushButton>
#include <QLineEdit>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QRadioButton>
#include <QLineEdit>
#include <TabWidget.hh>


class RootSelectionDialog : public QDialog
{
    Q_OBJECT

public:
    RootSelectionDialog(QWidget* parent = nullptr);
    ~RootSelectionDialog() = default;
    
    void setSimPath();
    void setExpPath();

    QString simPath();
    QString expPath();
    QString legendTitle();
    int timeScale();

private:
    void createUi();

private:
    QString m_expDataPath;
    QString m_simDataPath;
    QLineEdit* m_ordinateLineEdit;
    QSpinBox* m_scaleSpinBox;
};

class MainWindow : public QDialog
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

signals:
    void simulationRequested(double a, double b, double Nirrep ,double Dz);
    void macSelectionRequested();
    void rootRequested(const QString& simDataPath, const QString& expDataPath, const QString& legendTitle, const int timeScale);

public slots:
    void onMacChanged(const QString& filePath);
    void onRootGenerated(const QString& plotPath);

private:
    void createUi();

    void createTab();

private:
    QPushButton*    m_simulationButton;
    QPushButton*    m_macSelectionButton;

    QLabel*         m_macLabel;
    QLabel*         m_ROOTGraphicWidget;

    QDoubleSpinBox*      m_aSpinBox;
    QDoubleSpinBox*      m_bSpinBox;
    QDoubleSpinBox*      m_nirrepSpinBox;
    QDoubleSpinBox*      m_DzSpinBox;

    TabWidget*      m_tabWidget;
};
#endif // MAINWINDOW_H
