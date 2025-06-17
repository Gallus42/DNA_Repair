#include <QToolButton>
#include <QTabBar>
#include <QString>
#include <QHBoxLayout>
#include <QWidget>
#include <QVariant>
#include <QAction>
#include <QMenu>
#include <QLabel>
#include <QDialog>

#include <iostream>

#include "TabWidget.hh"


void Tab::setFilePath(const QString value)
{
    filePath = value;
}

QString Tab::getFilePath()
{
    return filePath;
}

TabWidget::TabWidget(QWidget *parent)
    : QTabWidget(parent)
{
    createUi();
}

TabWidget::~TabWidget()
{
}

void TabWidget::createUi()
{
    addCreateNewTabButton();
}

int TabWidget::addTab(QWidget* widget, const QString& title)
{
    const int ret = QTabWidget::addTab(widget, title);

    setCurrentIndex(ret);

    QTabBar* tabBar = this->tabBar();
    tabBar->setTabData(ret, QVariant::fromValue(QMap<int, QVariant>()));
    tabBar->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(tabBar, &QTabBar::customContextMenuRequested, this, [this, tabBar](const QPoint& point)
    {
        QMenu* menu = new QMenu;
        QAction* unpinAction = new QAction("Unpin Tab");
        connect(unpinAction, &QAction::triggered, this, &TabWidget::onUnpinTriggered);

        menu->addAction(unpinAction);
        menu->exec(mapToGlobal(point));
    });

    QToolButton* tabButton = new QToolButton(tabBar);
    tabButton->setToolButtonStyle(Qt::ToolButtonIconOnly);
    tabButton->setStyleSheet("QToolButton { background-color: transparent; border: none; }");
    tabButton->setIcon(QIcon("../9068678.png"));

    tabBar->setTabButton(ret, QTabBar::ButtonPosition::RightSide, tabButton);

    connect(tabButton, &QToolButton::clicked, this, [this, tabButton, tabBar]()
    {
        if (1 == tabBar->count())
            return;

        for (int tabIndex = 0; tabIndex < tabBar->count(); ++tabIndex)
        {
            QWidget* btn = tabBar->tabButton(tabIndex, QTabBar::ButtonPosition::RightSide);
            if (btn != tabButton)
                continue;

            removeTab(tabIndex);
            break;
        }
    });

    return ret;
}

void TabWidget::addCreateNewTabButton()
{
    QToolButton* addNewTabBtn = new QToolButton();
    addNewTabBtn->setIcon(QIcon("../add-button.png"));
    addNewTabBtn->setToolButtonStyle(Qt::ToolButtonIconOnly);
    addNewTabBtn->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_T));

    const QString styleSheet = "QToolButton { background-color: transparent;"
        "color: #bfbfbf;"
        "border-radius: 7px;"
        "border: 1px solid #E3E3FF;"
        "border-top: 1px solid #E3E3FF;"
        "border-left: 1px solid #E3E3FF;"
        "margin: 2px;"
        "width: 25px;"
        "height: 25px;"
        "selection-background-color: #18456d; }"

        "QToolButton:hover { background-color: #E3E3FF;"
        "color: #bfbfbf;"
        "border-radius: 7px;"
        "border: 1px solid #E3E3FF;"
        "border-top: 1px solid #E3E3FF;"
        "border-left: 1px solid #E3E3FF;"
        "margin: 2px;"
        "width: 25px;"
        "height: 25px;"
        "selection-background-color: #18456d; }";

    addNewTabBtn->setStyleSheet(styleSheet);

    QHBoxLayout* buttonsLayout = new QHBoxLayout();
    buttonsLayout->addWidget(addNewTabBtn);
    buttonsLayout->setContentsMargins(0, 0, 0, 0);
    buttonsLayout->setSpacing(0);

    QWidget* buttonsWidget = new QWidget();
    buttonsWidget->setLayout(buttonsLayout);

    setCornerWidget(buttonsWidget, Qt::TopRightCorner);

    connect(addNewTabBtn, &QToolButton::clicked, this, &TabWidget::createNewTabButtonClicked);
}

void TabWidget::onUnpinTriggered()
{
    if (count() < 2)
        return;

    QDialog* dialog = new QDialog;
    dialog->setWindowTitle(this->tabText(this->currentIndex()));

    Tab* tab = dynamic_cast<Tab*>(this->currentWidget());
    QGridLayout* layout = new QGridLayout;
    auto label = new QLabel();
    label->setFixedSize(800, 600);
    const QString filePath = tab->getFilePath();\
    if (!filePath.isEmpty())
        label->setPixmap(QPixmap(filePath));

    layout->addWidget(label);
    dialog->setLayout(layout);
    dialog->show();
    dialog->setContextMenuPolicy(Qt::CustomContextMenu);
    this->removeTab(this->currentIndex());
}