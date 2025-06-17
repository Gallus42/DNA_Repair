#include <QGridLayout>
#include <QLabel>
#include <QButtonGroup>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QFileDialog>
#include <QDialogButtonBox>
#include <QLineEdit>

#include "mainwindow.h"


MainWindow::MainWindow(QWidget *parent)
    : QDialog(parent)
{
    setWindowTitle("DNAReparation");
    setFixedSize(870, 830);
    createUi();
}

MainWindow::~MainWindow()
{
}

void MainWindow::onMacChanged(const QString& filePath)
{
    m_macLabel->setText(tr("Mac file path:%1").arg(filePath));
    m_simulationButton->setDisabled(filePath.isEmpty());
}

void MainWindow::onRootGenerated(const QString& plotPath)
{
    QFileInfo filePath(plotPath);

    Tab* tab = dynamic_cast<Tab*>(m_tabWidget->currentWidget());
    m_tabWidget->setTabText(m_tabWidget->currentIndex(), filePath.fileName());
    tab->setFilePath(plotPath);
    auto layout = tab->layout();
    for (int index = 0; index < layout->count(); ++index)
    {
        QLabel* label = dynamic_cast<QLabel*>(layout->itemAt(index)->widget());
        if (nullptr == label)
            continue;

        label->setPixmap(QPixmap(plotPath));
        break;
    }
}

void MainWindow::createUi()
{
    auto mainLayout = new QGridLayout;

    m_simulationButton = new QPushButton("Start simulation");
    m_simulationButton->setDisabled(true);
    connect(m_simulationButton, &QPushButton::clicked, this, [this](){emit simulationRequested(m_aSpinBox->value(), m_bSpinBox->value(), m_nirrepSpinBox->value(), m_DzSpinBox->value());});

    m_macSelectionButton = new QPushButton("Select .in file");
    connect(m_macSelectionButton, &QPushButton::clicked, this, &MainWindow::macSelectionRequested);

    m_macLabel = new QLabel(tr("Mac file path:"));
    
    m_tabWidget = new TabWidget(this);
    m_tabWidget->setFixedSize(850, 700);
    createTab();
    connect(m_tabWidget, &TabWidget::createNewTabButtonClicked, this, &MainWindow::createTab);

    m_ROOTGraphicWidget = new QLabel;
    m_ROOTGraphicWidget->setFixedSize(800, 600);

    QLabel* aLbl = new QLabel("a Parameter:");
    QLabel* bLbl = new QLabel("b Parameter:");
    QLabel* nirrepLbl = new QLabel("Irreparable value:");
    QLabel* DzLbl = new QLabel("Dz value:");

    m_aSpinBox = new QDoubleSpinBox;
    m_aSpinBox->setDecimals(7);
    m_bSpinBox = new QDoubleSpinBox;
    m_bSpinBox->setDecimals(7);
    m_nirrepSpinBox = new QDoubleSpinBox;
    m_DzSpinBox = new QDoubleSpinBox;

    QWidget* parametrWidget = new QWidget;

    QHBoxLayout* parameterLayout = new QHBoxLayout;
    parameterLayout->addWidget(aLbl);
    parameterLayout->addWidget(m_aSpinBox);
    parameterLayout->addWidget(bLbl);
    parameterLayout->addWidget(m_bSpinBox);
    parameterLayout->addWidget(nirrepLbl);
    parameterLayout->addWidget(m_nirrepSpinBox);
    parameterLayout->addWidget(DzLbl);
    parameterLayout->addWidget(m_DzSpinBox);


    parametrWidget->setLayout(parameterLayout);

    mainLayout->addWidget(m_simulationButton, 0, 0);

    mainLayout->addWidget(parametrWidget, 1, 0);

    mainLayout->addWidget(m_macSelectionButton, 2, 0, Qt::AlignTop);
    mainLayout->addWidget(m_macLabel, 3, 0, Qt::AlignTop);
    mainLayout->addWidget(m_tabWidget, 4, 0, Qt::AlignTop);

    mainLayout->setRowStretch(4, 1);

    setLayout(mainLayout);
}

void MainWindow::createTab()
{
    Tab* newTab = new Tab;

    auto tabLayout = new QGridLayout;

    QPushButton* ROOTButton = new QPushButton("Generate ROOT plot");
    connect(ROOTButton, &QPushButton::clicked, this, [this]()
    {
        RootSelectionDialog* dialog = new RootSelectionDialog(this);
        
        connect(dialog, &RootSelectionDialog::accepted, this, [this, dialog] ()
        {
            emit rootRequested(dialog->simPath(), dialog->expPath(), dialog->legendTitle(), dialog->timeScale());
            dialog->deleteLater();
        });

        connect(dialog, &RootSelectionDialog::rejected, this, [dialog] ()
        {
            dialog->deleteLater();
        }); 

        dialog->show();
    });

    QLabel* m_ROOTGraphicWidget = new QLabel;
    m_ROOTGraphicWidget->setFixedSize(800, 600);

    tabLayout->addWidget(m_ROOTGraphicWidget, 0, 0);
    tabLayout->addWidget(ROOTButton, 1, 0);
    newTab->setLayout(tabLayout);

    m_tabWidget->addTab(newTab, "newTab");
}

RootSelectionDialog::RootSelectionDialog(QWidget* parent)
    : QDialog(parent)
{
    createUi();
}

void RootSelectionDialog::setSimPath()
{
    m_simDataPath = QFileDialog::getOpenFileName(dynamic_cast<QWidget*>(parent()), tr("Open Simulation File"), "C:", tr("Data Files (*.out *.csv *.dat)"));
}

void RootSelectionDialog::setExpPath()
{
    m_expDataPath = QFileDialog::getOpenFileName(dynamic_cast<QWidget*>(parent()), tr("Open Expiremental File"), "C:", tr("Data Files (*.out *.csv *.dat)"));
}

QString RootSelectionDialog::simPath()
{
    return m_simDataPath;
}

QString RootSelectionDialog::expPath()
{
    return m_expDataPath;
}

QString RootSelectionDialog::legendTitle()
{
    return m_ordinateLineEdit->text();
}

int RootSelectionDialog::timeScale()
{
    return m_scaleSpinBox->value();
}

void RootSelectionDialog::createUi()
{
    setModal(true);

    QPushButton* simulationButton = new QPushButton("Select simulation Data File");
    QPushButton* experimentalButton = new QPushButton("Select Experimental Data File");

    QLabel* simPathLabel = new QLabel("Simulation Path: ");
    QLabel* expPathLabel = new QLabel("ExperimentalPath: ");

    QLabel* ordinateLabel = new QLabel("Y-Axis legend:");
    m_ordinateLineEdit = new QLineEdit("Enzyme per cell (%)");

    QLabel* spinBoxLabel = new QLabel("Max Time(X):");
    m_scaleSpinBox = new QSpinBox;
    m_scaleSpinBox->setValue(25);
    m_scaleSpinBox->setMinimum(1);
    m_scaleSpinBox->setMaximum(25);

    QDialogButtonBox* btnBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QPushButton* okButton = btnBox->button(QDialogButtonBox::Ok);
    okButton->setDisabled(true);

    connect(btnBox, &QDialogButtonBox::rejected, this, &RootSelectionDialog::rejected);
    connect(okButton, &QPushButton::clicked, this, &RootSelectionDialog::accepted);

    connect(simulationButton, &QPushButton::clicked, this, [this, simPathLabel, okButton] ()
    {
        setSimPath();
        simPathLabel->setText("Simulation Path: " + m_simDataPath);
        okButton->setDisabled(m_simDataPath.isEmpty());
    });

    connect(experimentalButton, &QPushButton::clicked, this, [this, expPathLabel]()
    {
        setExpPath();
        expPathLabel->setText("ExperimentalPath: " + m_expDataPath);
    });

    auto layout = new QGridLayout;
    layout->addWidget(simulationButton, 0, 0, 1, 2);
    layout->addWidget(simPathLabel, 1, 0);
    layout->addWidget(experimentalButton, 2, 0, 1, 2);
    layout->addWidget(expPathLabel, 3, 0);
    layout->addWidget(ordinateLabel, 4, 0);
    layout->addWidget(m_ordinateLineEdit, 4, 1);
    layout->addWidget(spinBoxLabel, 5, 0);
    layout->addWidget(m_scaleSpinBox, 5, 1);

    layout->addWidget(btnBox, 6, 1);
    layout->setColumnStretch(1, 1);

    setLayout(layout);
}
