#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "editlist.h"
#include "editfeat.h"
#include "choosecombinations.h"
#include <QSpinBox>
#include <QMessageBox>
#include<QFileDialog>
#include<QJsonDocument>
#include<QJsonObject>
#include<QJsonArray>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->radioButton_Unprocessed->setChecked(true);
    ui->checkBox->setChecked(false);
    ui->checkBox_2->setChecked(false);
    ui->checkBox_3->setChecked(false);
    ui->spinBox->setRange(1,5);
    ui->radioButton_wt_Groups->setChecked(true);
    ui->chooseCombsButton->hide();
    ui->label_Hypothesis->hide();
    ui->radioButton_G1_gr_G2->hide();
    ui->radioButton_G1_gr_G2->setChecked(true);
    ui->radioButton_G2_gr_G1->hide();
    ui->radioButton_NormFischer->setChecked(true);
    ui->radioButton_reg_y->setChecked(true);


    QString DefAnalysisName = "FConnectivityAnl";
    QString FSLDIRpos1 = "/usr/share/fsl/5.0/bin/";
    QString FSLDIRpos2 = "/usr/share/fsl/bin/";

    if (QDir(FSLDIRpos1).exists()){
        ui->lineEdit_FSLDIR->setText(FSLDIRpos1);
    }
    else if (QDir(FSLDIRpos2).exists()){
        ui->lineEdit_FSLDIR->setText(FSLDIRpos2);
    }

    ui->lineEdit_AnalysisName->setText(DefAnalysisName);
    ui->lineEdit_OutDir->setText(OutDir + DefAnalysisName);
}


MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_chooseData_clicked()
{
    int spin_1 = spinBoxValue;
    if (ui->radioButton_PreprocwFSL->isChecked()){
        editFeat editfeat;
        editfeat.setGroups(spin_1);
        editfeat.setModal(true);
        editfeat.exec();
        FeatFileNames= editfeat.get_FileNames();
    }
    else {
        EditList editlist;
        editlist.setGroups(spin_1);
        editlist.setModal(true);
        editlist.exec();
        StructuralFileNames = editlist.get_StructuralFileNames();
        FunctionalFileNames = editlist.get_FunctionalFileNames();
    }

}

void MainWindow::on_pushButton_7_clicked()
{
    QMessageBox::information(this,"Help", "Put the reference MNI structural file here.");
}

void MainWindow::on_pushButton_3_clicked()
{
    QApplication::quit();
}

void MainWindow::on_pushButton_6_clicked()
{
    QApplication::quit();
}

void MainWindow::on_pushButton_12_clicked()
{
    QApplication::quit();
}

void MainWindow::on_pushButton_18_clicked()
{
    QMessageBox::information(this,"Help", "If correlation files given, then there is no requirement of input.");

}

void MainWindow::on_commandLinkButton_Reference_clicked()
{
    QString new_input;
    QString already_input = ui->lineEdit_Reference->text();
    QFile Fout(already_input) ;

    if (Fout.exists())
    {
            new_input=QFileDialog::getOpenFileName(
                                                    this,
                                                    tr("Open File"),
                                                    already_input,
                                                    "All files (*.*);;Text File (*.txt)"
                                                       );
    }
    else
    {
        new_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open File"),
                                                CommonDir,
                                                "All files (*.*);;Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
    ui->lineEdit_Reference->setText(new_input);

    }
    else {
        ui->lineEdit_Reference->setText(already_input);
    }
}

void MainWindow::on_commandLinkButton_OutDir_clicked()
{
    QString new_Dir;

    if (QDir(OutDir).exists())
    {
        new_Dir=QFileDialog::getExistingDirectory(
                                                    this,
                                                    tr("Choose Directory"),
                                                    OutDir);
    }
    else {
        new_Dir =QFileDialog::getExistingDirectory(
                                                this,
                                                tr("Choose Directory"),
                                                CommonDir);
    }
    if (new_Dir!=NULL){
        ui->lineEdit_OutDir->setText(new_Dir +"/"+ ui->lineEdit_AnalysisName->text());
        OutDir = new_Dir+"/";
    }
    else{
        ui->lineEdit_OutDir->setText(OutDir + ui->lineEdit_AnalysisName->text());
    }
}

void MainWindow::on_commandLinkButton_ROIFile_clicked()
{
    QString new_input;
    QString already_input = ui->lineEdit_ROIFile->text();
    QFile Fout(already_input) ;

    if (Fout.exists())
    {
        already_input=QFileDialog::getOpenFileName(
                                                    this,
                                                    tr("Open File"),
                                                    already_input,
                                                    "All files (*.*);;Text File (*.txt)"
                                                       );
    }
    else
    {
        already_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open File"),
                                                CommonDir,
                                                "All files (*.*);;Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
        ui->lineEdit_ROIFile->setText(new_input);

    }
    else {
        ui->lineEdit_ROIFile->setText(already_input);
    }
}

void MainWindow::on_commandLinkButton_FSLDIR_clicked()
{

    QString filename=QFileDialog::getExistingDirectory(
                                                this,
                                                tr("Choose Directory"),
                                                "/home/deepak/");
    ui->lineEdit_FSLDIR->setText(filename);
}



void MainWindow::on_commandLinkButton_CorrFile_clicked()
{
    QString new_input;
    QString already_input = ui->lineEdit_CorrFile->text();
    QFile Fout(already_input) ;

    if (Fout.exists())
    {
      already_input=QFileDialog::getOpenFileName(
                                                    this,
                                                    tr("Open File"),
                                                    already_input,
                                                    "All files (*.*);;Text File (*.txt)"
                                                       );
    }
    else
    {
        already_input=QFileDialog::getOpenFileName(
                                                this,
                                                tr("Open File"),
                                                CommonDir,
                                                "All files (*.*);;Text File (*.txt)"
                                                   );
    }
    if (new_input!=NULL){
        ui->lineEdit_CorrFile->setText(new_input);

    }
    else {
        ui->lineEdit_CorrFile->setText(already_input);
    }
}

void MainWindow::on_help_ROIFile_clicked()
{
    QMessageBox::information(this,"Caution", "This file must be registered to reference.");
}


void MainWindow::on_lineEdit_AnalysisName_textEdited(const QString &arg1)
{
    ui->lineEdit_OutDir->setText(OutDir + arg1);

}

void MainWindow::on_radioButton_bw_Groups_clicked(bool checked)
{

    if (checked){
        if (spinBoxValue == 2){

            ui->label_Hypothesis->show();
            ui->radioButton_G1_gr_G2->show();
            ui->radioButton_G2_gr_G1->show();
        }
        else{
            ui->checkBox_AllCombs->show();
            ui->checkBox_AllCombs->setChecked(true);
        }
     }
}



void MainWindow::on_radioButton_wt_Groups_clicked(bool checked)
{
    if (checked){
        ui->chooseCombsButton->hide();
        ui->checkBox_AllCombs->hide();
        ui->label_Hypothesis->hide();
        ui->radioButton_G1_gr_G2->hide();
        ui->radioButton_G2_gr_G1->hide();
    }
}

void MainWindow::on_chooseCombsButton_clicked()
{
    /*
     * For testing purpose
    QList<int> list;
    list<<100<<500;
    */
    QString AB;
    chooseCombinations choosecombinations;
    choosecombinations.setCheckboxes(spinBoxValue);
    choosecombinations.setModal(true);
    choosecombinations.exec();
    //QMessageBox::information(this,"Help",QString::number(choosecombinations.checkedGroups[0]) );
    mirrored_checkStates = choosecombinations.checkedGroups;
    int GroupB = 2;
    int lastOne =1;
    for (int i = 0; i < mirrored_checkStates.size(); ++i) {
        if (mirrored_checkStates.at(i) == 2){

            AB+= "Group_" + QString::number(lastOne) + "_vs_Group_" + QString::number(GroupB) +",";
        }
        GroupB++;
        if (GroupB>spinBoxValue){
            lastOne++;
            GroupB =lastOne + 1;
        }
    }
    QMessageBox::information(this,"Help",AB);
}

void MainWindow::on_spinBox_valueChanged(int arg1)
{
    bool check_bw_groups = ui->radioButton_bw_Groups->isChecked();
    spinBoxValue = arg1;
    if (arg1==1){
        ui->radioButton_wt_Groups->setChecked(true);
        ui->radioButton_bw_Groups->hide();
        ui->checkBox_AllCombs->hide();
        ui->label_Hypothesis->hide();
        ui->radioButton_G1_gr_G2->hide();
        ui->radioButton_G2_gr_G1->hide();
    }
    else if (arg1==2){
        ui->radioButton_bw_Groups->show();
        ui->checkBox_AllCombs->hide();
        ui->chooseCombsButton->hide();
        if (check_bw_groups){
            ui->label_Hypothesis->show();
            ui->radioButton_G1_gr_G2->show();
            ui->radioButton_G2_gr_G1->show();
            ui->radioButton_G1_gr_G2->setChecked(true);
        }
    }
    else {
        ui->radioButton_bw_Groups->show();
        ui->label_Hypothesis->hide();
        ui->radioButton_G1_gr_G2->hide();
        ui->radioButton_G2_gr_G1->hide();
        if (check_bw_groups){
            ui->checkBox_AllCombs->show();
            ui->checkBox_AllCombs->setChecked(true);
        }

    }
}

void MainWindow::on_checkBox_AllCombs_clicked(bool checked)
{
    if (!checked){
        ui->chooseCombsButton->show();
    }
    else {
        ui->chooseCombsButton->hide();
    }
}

void MainWindow::writeReferImgpath(QJsonObject &json) const
{
    QString FilePath;
    QJsonObject ReferSummary;
    ReferSummary["Title"] = QString("It tells the path of the file used for reference.");
    ReferSummary["ReferImgPath"] = ui->lineEdit_Reference->text();
    json["ReferSummary"] = ReferSummary;

    QJsonObject ProcessingType;
    ProcessingType["Title"] = QString("The inputs can be of 3 types: "
                                      "#0 : Unprocessed "
                                      "#1 : Preprocessed without FSL "
                                      "#2 : Preprocessed with FSL");
    ProcessingType["ProcessingWay"] = ProcessingWay;
    json["ProcessingType"] = ProcessingType;

    if (ProcessingWay ==2){
        QJsonObject FeatList;
        FeatList["Info"] = QString("This is for case 2. We shall require the list of only the feat directories.");
        FeatList["No of Groups"] = spinBoxValue;
        QJsonArray FeatFilesArray;
        for (int i=1 ; i<=FeatFileNames.size() ; i++)
        {
            FilePath = OutChosenPath + "/"+ ui->lineEdit_AnalysisName->text() + "_Group_" + QString::number(i)+"_FeatDirs.txt";
            QFile::copy(FeatFileNames.at(i-1), FilePath );
            QFile file(FeatFileNames.at(i-1));
            file.remove();
            FeatFilesArray.append(FilePath);
        }
//        FeatFileNames.clear();
        FeatList["FeatFilePaths"] = FeatFilesArray;
        json["FeatFilesInfo"] = FeatList;
    }

    else{
        QJsonObject FilesInfo;
        FilesInfo["Info"] = QString("This is for case 0 and case 1. We shall require the list of both structural and functional files.");
        FilesInfo["No of Groups"] = spinBoxValue;
        QJsonArray StructuralArray;
        QJsonArray FunctionalArray;
        for (int i=1 ; i<=FunctionalFileNames.size() ; i++)
        {
            FilePath = OutChosenPath + "/"+ ui->lineEdit_AnalysisName->text() + "_Group_" + QString::number(i)+"_FunctionalFiles.txt";
            QFile::copy(FunctionalFileNames.at(i-1), FilePath);
            QFile file(FunctionalFileNames.at(i-1));
            file.remove();
            FunctionalArray.append(FilePath);
        }
        for (int i=1 ; i<=FunctionalFileNames.size() ; i++)
        {
            FilePath = OutChosenPath + "/"+ ui->lineEdit_AnalysisName->text() + "_Group_" + QString::number(i)+"_StructuralFiles.txt";
            QFile::copy(StructuralFileNames.at(i-1), FilePath);
            QFile file(StructuralFileNames.at(i-1));
            file.remove();
            StructuralArray.append(FilePath);
        }
//        FunctionalFileNames.clear();
//        StructuralFileNames.clear();
        FilesInfo["FunctionalFilePaths"] = FunctionalArray;
        FilesInfo["StructuralFilePaths"] = StructuralArray;
        json["FilesInfo"] = FilesInfo;
    }

}


void MainWindow::writeAnalysisName(QJsonObject &json) const
{
    json["Analysis Name"] = ui->lineEdit_AnalysisName->text();
    QJsonObject ReferenceImgpath;
    writeReferImgpath(ReferenceImgpath);
    json["AnalysisParams"] = ReferenceImgpath;


//    QJsonArray levelArray;
//    foreach (const Level level, mLevels) {
//        QJsonObject levelObject;
//        level.write(levelObject);
//        levelArray.append(levelObject);
//    }
//    json["levels"] = levelArray;
}
void MainWindow::on_pushButton_Go_clicked()
{
    OutChosenPath = ui->lineEdit_OutDir->text();
    if (!QDir(OutChosenPath).exists()){
        QDir dir;
        bool created = dir.mkdir(OutChosenPath);
        if (!created){
            QMessageBox::warning(this,"Help","Not able to create Directory: "+ OutChosenPath);
            return;
        }
    }

    QFile saveFile( OutChosenPath + QStringLiteral( "/Design.json"));
    if (!saveFile.open(QIODevice::WriteOnly)) {
        qWarning("Couldn't open save file.");
    }
    QJsonObject gameObject;
    writeAnalysisName(gameObject);
    QJsonDocument saveDoc(gameObject);
    saveFile.write(saveDoc.toJson());
}

void MainWindow::on_radioButton_Unprocessed_clicked()
{
    ProcessingWay = 0;
}

void MainWindow::on_radioButton_PreprocwtFSL_clicked()
{
    ProcessingWay = 1;
}

void MainWindow::on_radioButton_PreprocwFSL_clicked()
{
    ProcessingWay = 2;
}
