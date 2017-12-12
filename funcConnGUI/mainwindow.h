#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QDir>
#include<QJsonArray>
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();


private slots:
    void on_pushButton_chooseData_clicked();

    void on_pushButton_7_clicked();

    void on_pushButton_3_clicked();

//    void on_pushButton_6_clicked();

    void on_pushButton_12_clicked();

    void on_pushButton_18_clicked();

    void on_commandLinkButton_Reference_clicked();

    void on_commandLinkButton_OutDir_clicked();

    void on_commandLinkButton_ROIFile_clicked();

    void on_commandLinkButton_FSLDIR_clicked();

    void on_commandLinkButton_CorrFile_clicked();

    void on_help_ROIFile_clicked();

    void on_lineEdit_AnalysisName_textEdited(const QString &arg1);

    void on_radioButton_bw_Groups_clicked(bool checked);



    void on_radioButton_wt_Groups_clicked(bool checked);

    void on_chooseCombsButton_clicked();

    void on_spinBox_valueChanged(int arg1);

    void on_checkBox_AllCombs_clicked(bool checked);

    void writeStats(QJsonObject &json) const;

    void writeReferImgpath(QJsonObject &json) const;


    void writeAnalysisName(QJsonObject &json) const;

    void on_pushButton_Go_clicked();

    void on_radioButton_Unprocessed_clicked();

    void on_radioButton_PreprocwtFSL_clicked();

    void on_radioButton_PreprocwFSL_clicked();

    void on_checkBox_HighPass_clicked(bool checked);

    void on_checkBox_LowPass_clicked(bool checked);

    void on_commandLinkButton_Mask_clicked();

    void on_pushButton_Go3_clicked();

    void on_checkBox_BET_clicked(bool checked);

private:
    QJsonArray combinations;
    QString CommonDir = getenv("PWD");
    QString OutDir = getenv("PWD");
    QString FSLDIR = getenv("FSLDIR");
    QList<int> mirrored_checkStates;
    QList<QString> FeatFileNames;
    QList<QString> FunctionalFileNames;
    QList<QString> StructuralFileNames;
    QString OutChosenPath;
    QString WorkingDir = QDir::currentPath();
    int ProcessingWay =0;
    int spinBoxValue =1;
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
