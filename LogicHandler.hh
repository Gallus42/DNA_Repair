#include <QObject>


class LogicHandler : public QObject
{
    Q_OBJECT

public:
    
    enum EnfyzeType : int
    {

    };

    LogicHandler(QWidget* parent = nullptr);
    ~LogicHandler() = default;

    static int argc_main;
    static char** argv_main;

signals:
    void rootGenerated(const QString& plotPath);

public slots:
    void onSimulationRequested(double a, double b, double Nirrep ,double Dz);
    void onMacSelectionRequested();
    void onROOTRequested(const QString& simPath, const QString& expPath, const QString& legendTitle, const int timeScale);
    void setArgs(int argc, char** argv);

signals:
    void fileSelected(const QString& filePath);

private:
    QString m_macroFilePath;
};
