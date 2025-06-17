#ifndef TABWIDGET_H
#define TABWIDGET_H

#include <QTabWidget>

class Tab : public QWidget
{
public:
    Tab() = default;
    ~Tab() = default;

public:
    void setFilePath(const QString value);
    QString getFilePath();

private:
    QString filePath;
};

class TabWidget : public QTabWidget
{
    Q_OBJECT

public:
    TabWidget(QWidget *parent = nullptr);
    ~TabWidget();

    int addTab(QWidget* widget, const QString& title);
    void addCreateNewTabButton();

signals:
    void createNewTabButtonClicked();

private slots:
    void onUnpinTriggered();

private:
    void createUi();
};
#endif 
